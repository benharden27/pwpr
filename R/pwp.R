#' Main PWP Program
#'
#' @param met_input_file
#' @param profile_input_file
#' @param dt time-step increment (seconds)
#' @param dz depth increment (meters)
#' @param days the number of days to run
#' @param depth the depth to run
#' @param dt_save time-step increment for saving to file (multiples of dt)
#' @param lat latitude (degrees)
#' @param rb critical bulk richardson number
#' @param rg critical gradient richardson number
#' @param rkz background diffusion
#' @param beta1 longwave extinction coefficient (m)
#' @param beta2 shortwave extinction coefficient (m)
#'
#' @return
#' @export
#'
#' @examples
pwp <- function(met_input_file, profile_input_file,
                dt = 900, dz = 10, days = 300, depth = 1000, dt_save = 4,
                lat = 60, rb = 0.64, rg = 0.25, rkz = 0,
                beta1 = 0.6, beta2 = 20) {

  # Setup additional variables
  g <- 9.81
  cpw	<- 4183.3
  f <- 2 * 7.29E-5 * sin(lat * pi / 180)

  # Format air-sea flux data --------------------------------------------------

  # Load air-sea forcing data from csv file
  met <- read.csv(met_input_file)

  # Interpolate the air/sea flux variables at dt resolution
  # Set nmet equal to the number of time increments using a resolution of dt

  nmet <- days * 8.64E4 / dt
  time <- met$time[1] + 0:(nmet-1) * dt / 8.64E4

  qi <- approx(met$time, met$sw, time)
  qo <- approx(met$time, met$lw + met$qlat + met$qsens, time)
  tx <- approx(met$time, met$tx, time)
  ty <- approx(met$time, met$ty, time)
  precip <- approx(met$time, met$precip, time)

  # Interpolate evaporation minus precipitation at dt resolution
  evap <- 0.03456 * approx(met$time, met$qlat, time) / (86400 * 1000)
  emp <- evap - precip


  # Format profile data -------------------------------------------------------

  # load profile data from csv file
  prof <- read.csv(profile_input_file)

  # Interpolate the profile variables at dz resolution
  # Set nz equal to the number of depth increments + 1 using a resolution of dz

  nz <- 1 + depth/dz
  z <- 0:(nz - 1) * dz

  t <- approx(prof$z, prof$t, z)
  s <- approx(prof$z, prof$s, z)
  d <- oce::swSigma(s,t)

  # Initialize additional profile variables at dz resolutions

  #u and v:   east and north current
  #absrb:   absorbtion factor

  u <- rep(0,nz)
  v <- rep(0,nz)
  absrb <- absorb(beta1,beta2)

  # Specify a simple "background" diffusion to be applied to the profiles

  dstab <- dt * rkz / dz^2

  if (dstab > 0.5){
    disp('Warning, this value of rkz will be unstable')
  }

  # Define the variables to be passes into pwp routine that will change at each time step
  pwp_in <- list(time = 0,
                 t = t,
                 s = s,
                 d = d,
                 u = u,
                 v = v)

  # save these as the first instance of pwp_output
  pwp_output <- pwp_in

  # set up list of params to pass to pwpgo
  # Get all variables in the envrionemtns
  defined_vars <- ls(envir = environment())
  # remove those that are already named in pwp_in
  defined_vars <- defined_vars[defined_vars %in% names(pwp_in)]
  # put all these variables in a list
  params <- mget(defined_vars)


  # Enter the PWP time loop -------------------------------------------------

  for (m in 1:nmet) {

    pwp_in <- pwpgo(pwp_in,params, m)

    if (m %% dt_save == 0) {
      pwp_output$time <- append(pwp_output$time, pwp_out$time)
      pwp_output$t <-  cbind(pwp_output$t, pwp_out$t)
      pwp_output$s <-  cbind(pwp_output$s, pwp_out$s)
      pwp_output$d <-  cbind(pwp_output$d, pwp_out$d)
      pwp_output$u <-  cbind(pwp_output$u, pwp_out$u)
      pwp_output$v <-  cbind(pwp_output$v, pwp_out$v)
    }


  }

  return(pwp_output)

  }


pwpgo <- function(pwp_in,params,m) {
  # ------------------------------------------------------------
  # Unpack variable lists into function environment
  # Does the same thing as the global calls in MATLAB scripts
  list2env(pwp_in, env = environment())
  list2env(params, env = environment())

  # Apply heat and fresh water fluxes to the top most grid cell
  t[1] <- t[1] + (qi * absrb[1] - qo) * dt / (dz * d[1] * cpw)
  s[1] <- s[1] / (1 - emp * dt/ dz)

  # Absorb solar radiation at depth

  t[2:nz] <- t[2:nz] + qi * absrb[2:nz] * dt / (dz * d[2:nz] *cpw)

  # Compute the density, and relieve static instability, if it occurs

  d <- swSigma(s,t)

  # remove_si ? function was here

  # At this point the density profile should be statically stable

  # Find the index of the surfaced mixed-layer right after the heat/salt fluxes

  ml_index <- min(which(diff(d) > 1E-4))

  # Get the depth of the surfaced mixed-layer

  ml_depth <- z[ml_index + 1]

  # Time step the momentum equation

  # Rotate the current throughout the water column through an
  # angle equal to inertial rotation for half of a time step

  ang <- -f * dt / 2

  rot(ang)

  # Apply the wind stress to the mixed layer as it now exists

  du <- (tx / (ml_depth * d[1])) * dt
  dv <- (ty / (ml_depth * d[1])) * dt
  u[1:ml_index] <- u[1:ml_index] + du
  v[1:ml_index] <- v[1:ml_index] + dv

  # Apply drag to the current (this is a horrible parameterization of
  # inertial-internal wave dispersion)

  if (ucon > 1E-10){
    u <- u * (1 - dt * ucon)
    v <- v * (1 - dt * ucon)
  }

  # Rotate another half time step

  rot(ang)

  # Finished with the momentum equation for this time step

  # Do the bulk Richardson number instability form of mixing (as in PWP)

  if (rb > 1E-5){
    bulk_mix(ml_index)
  }

  # Do the gradient Richardson number instability form of mixing

  if (rg > 0){
    grad_mix
  }
}

# --------------------------------------------------------------

bulk_mix <- function(ml_index) {
  rvc <- rb
  for (j in ml_index +1 :nz){
    h <- z[j]
    dd <- (d[j] - d[1]) / d[1]
    dv <- (u[j] - u[1])^2 + (v[j] - v[1])^2
    if (dv == 0){
      rv <- Inf
    } else {
      rv <- g * h * dd / dv
    }
    if (rv > rvc){
      break
    } else {
      mix5[j]
    }
  }
}

# ----------------------------------------------------------------

grad_mix <- function({

  # This function performs the gradient Richardson Number relaxation
  # by mixing adjacent cells just enough to bring them to a new
  # Richardson Number

  rc <- rg

  # Compute the gradient Richardson Number, taking care to avoid dividing by
  # zero in the mixed layer. The numerical values of the minimum allowable
  # density and velocity differences are entirely arbitrary, and should not
  # effect the calculations (except that on some occasions they evidently have)

  j1 <- 1
  j2 <- nz - 1

  while (1){
    for (j <- j1:j2){
      if (j <= 0){
        keyboard
      }
      dd <- (d[j + 1] - d[j]) / d[j]
      dv <- (u[j + 1] - u[j])^2 + (v[j + 1] - v[j])^2
      if (dv == 0){
        r[j] <- Inf
      } else {
        r[j] <- g * dz * dd / dv
      }
    }
    # Find the smallest value of r in profile

    rs <- min(r)
    js <- min(which(r == s))

    # Check to see whether the smallest r is critical or not

    if (rs > rc){
      return
    }

    # Mix the cells js and js+1 that had the smallest Richardson Number

    stir(rc, rs, js) #IDK WHAT STIR IS IN R?

    # Recompute the Richardson Number over the part of the profile that has changed

    j1 <- js - 2
    if (j1 < 1){
      j1 <- 1
    }
    j2 <- js + 2
    if (j2 > nz - 1){
      j2 <- nz - 1
    }
  }
}

# ------------------------------------------------------------------------

stir <- function(rc,r,j){

  # This subroutine mixes cells j and j+1 just enough so that
  # the Richardson Number after the mixing is brought up to
  # the value rnew. In order to have this mixing process converge,
  # rnew must exceed the critical value of the Richardson number where
  # mixing is presumed to start. If r critical <- rc <- 0.25 (the nominal value)
  # and r <- 0.20, then
  # rnew <- 0.3 would be reasonable. If r were smaller, then a
  # larger value of rnew - rc is used to hasten convergence.

  rcon <- 0.02 + (rc - r) / 2
  rnew <- rc + rcon / 5
  f <- 1 - r / rnew
  dt <- (t[j + 1] - t[j]) * f / 2
  t[j + 1] <- t[j + 1] - dt
  t[j] <- t[j] + dt
  ds <- (s[j + 1] - s[j]) * f / 2
  s[j + 1] <- s[j + 1] - ds
  s[j] <- s[j] + ds
  d[j:j + 1] <- swSigma(s[j:j + 1], t[j:j + 1])
  du <- (u[j +1] - u[j]) * f / 2
  u[j + 1] <- u[j + 1] - du
  u[j] <- u[j] + du
  dv <- (v[j + 1] - v[j]) * f / 2
  v[j +1] <- v[j + 1] - dv
  v[j] <- v[j] + dv
}

# ----------------------------------------------------------------------

mix5 <- function(j){

  # This subroutine mixes the arrays t, s, u, v down to level j

  t[1:j] <- mean(t[1:j])
  s[1:j] <- mean(s[1:j])
  d[1:j] <- swSigma(s[1:j], t[1:j])
  u[1:j] <- mean(u[1:j])
  v[1:j] <- mean(v[1:j])
}

# ------------------------------------------------------------------

rot <- function(ang){

  # This subroutine rotates the vector (u,v) through an angle, ang

  r <- (u + i * v) * exp (i * ang)
  u <- real(r)
  v <- imag(r)
}

# ---------------------------------------------------------------

remove_si <- function{

  # Find and relieve static instability that may occur in the
  # density array d. This simulates free convection.
  # ml_index is the index of the depth of the surface mixed layer
  # after adjustment

  while (1){
    ml_index <- min(which(diff(d) < 0))
    if (is.empty.model(ml_index)){
      break
    }
    mix5(ml_index + 1)
  }
}

# ---------------------------------------------------------------

absrb <- function(beta1, beta2){

  # Compute solar radiation absorbtion profile. This subroutine
  # assumes two wavelengths, and a double exponential depth for absorbtion

  # Subscript 1 is for red, non-penetrating light, and 2 is for blue,
  # penetrating light. rsl is the funcion assumed to be red

  rs1 <- 0.6
  rs2 <- 1.0 - rs1
  absrb <- rep(0,nz)
  z1 <- (0:nz-1) * dz
  z2 <- z1 + dz
  z1b1 <- z1 / beta1
  z2b1 <- z2 / beta1
  z1b2 <- z1 / beta2
  z2b2 <- z2 / beta2
  absrb <- (rs1 * (exp(-z1b1) - exp(-z2b1) + rs2 * (exp(-z1b2) - exp(-z2b2)))
}


