#' Main PWP Program
#'
#' @param met_input_file Input metereological data
#' @param profile_input_file Input profile data
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
  ucon <- 0.1*abs(f);

  # Format air-sea flux data --------------------------------------------------

  # Load air-sea forcing data from csv file
  met <- read.csv(met_input_file)

  # Interpolate the air/sea flux variables at dt resolution
  # Set nmet equal to the number of time increments using a resolution of dt

  nmet <- days * 8.64E4 / dt
  time <- met$time[1] + 0:(nmet-1) * dt / 8.64E4

  qi <- approx(met$time, met$sw, time)$y
  qo <- approx(met$time, met$lw + met$qlat + met$qsens, time)$y
  tx <- approx(met$time, met$tx, time)$y
  ty <- approx(met$time, met$ty, time)$y
  precip <- approx(met$time, met$precip, time)$y

  # Interpolate evaporation minus precipitation at dt resolution
  evap <- 0.03456 * approx(met$time, met$qlat, time)$y / (86400 * 1000)
  emp <- evap - precip

  # Format profile data -------------------------------------------------------

  # load profile data from csv file
  prof <- read.csv(profile_input_file)

  # Interpolate the profile variables at dz resolution
  # Set nz equal to the number of depth increments + 1 using a resolution of dz

  nz <- 1 + depth/dz
  z <- 0:(nz - 1) * dz

  t <- approx(prof$z, prof$t, z, rule = 2)$y
  s <- approx(prof$z, prof$s, z, rule = 2)$y
  d <- oce::swSigma(s, t, p = 0)

  # Initialize additional profile variables at dz resolutions

  #u and v:   east and north current
  #absrb:   absorbtion factor

  u <- rep(0,nz)
  v <- rep(0,nz)
  absrb <- absorb(beta1,beta2,nz,dz)

  # Specify a simple "background" diffusion to be applied to the profiles

  dstab <- dt * rkz / dz^2

  if (dstab > 0.5){
    disp('Warning, this value of rkz will be unstable')
  }

  # Define the variables to be passes into pwp routine that will change at each time step
  ts <- 0
  pwp_in <- mget(c("ts","t","s","d","u","v"))
  pwp_output <- pwp_in

  # set up list of params to pass to pwpgo
  # Get all variables in the envrionemtns
  defined_vars <- ls(envir = environment())
  # remove those that are already named in pwp_in
  defined_vars <- defined_vars[!(defined_vars %in% names(pwp_in))]
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


#' Run one loop of the PWP program
#'
#' @param pwp_in list of pwp input variables that will change through loop
#' @param params List of constants that will stay the same through loop
#' @param m particular time step index
#'
#' @return
#' @export
#'
#' @examples
pwpgo <- function(pwp_in, params, m) {
  # ------------------------------------------------------------
  # Unpack variable lists into function environment
  # Does the same thing as the global calls in MATLAB scripts
  pwp_names <- names(pwp_in)
  list2env(pwp_in, env = environment())
  param_names <- names(params)
  list2env(params, env = environment())

  # Apply heat and fresh water fluxes to the top most grid cell
  t[1] <- t[1] + (qi[m] * absrb[1] - qo[m]) * dt / (dz * d[1] * cpw)
  s[1] <- s[1] / (1 - emp[m] * dt/ dz)

  # Absorb solar radiation at depth
  t[2:nz] <- t[2:nz] + qi[m] * absrb[2:nz] * dt / (dz * d[2:nz] *cpw)

  # Compute the density, and relieve static instability, if it occurs

  d <- oce::swSigma(s,t, p=0)

  # remove_si ? function was here
  pwp_int <- mget(pwp_names)
  pwp_int <- remove_si(pwp_int)
  list2env(pwp_int, env = environment())

  # At this point the density profile should be statically stable

  # Find the index of the surfaced mixed-layer right after the heat/salt fluxes
  ml_index <- min(which(diff(d) > 1E-4))

  # Get the depth of the surfaced mixed-layer
  ml_depth <- z[ml_index + 1]

  # Time step the momentum equation

  # Rotate the current throughout the water column through an
  # angle equal to inertial rotation for half of a time step

  ang <- -f * dt / 2

  uv <- rot(u,v,ang)
  list2env(uv, env = environment())

  # Apply the wind stress to the mixed layer as it now exists

  du <- (tx[m] / (ml_depth * d[1])) * dt
  dv <- (ty[m] / (ml_depth * d[1])) * dt
  u[1:ml_index] <- u[1:ml_index] + du
  v[1:ml_index] <- v[1:ml_index] + dv

  # Apply drag to the current (this is a horrible parameterization of
  # inertial-internal wave dispersion)
  # ucon is the coefficient of inertial-internal wave dissipation (0) s^-1

  if (ucon > 1E-10){
    u <- u * (1 - dt * ucon)
    v <- v * (1 - dt * ucon)
  }

  # Rotate another half time step
  # NEED TO CREATE OUTPUT FROM THIS FUNCTION
  uv <- rot(u,v,ang)
  list2env(uv, env = environment())

  # Finished with the momentum equation for this time step

  pwp_int <- mget(pwp_names)
  params$ml_index <- ml_index

  # Do the bulk Richardson number instability form of mixing (as in PWP)

  # ADD ALL VARIABLES TO bulk_mix function definition
  if (rb > 1E-5){
    pwp_int <- bulk_mix(pwp_int, params)
  }

  # Do the gradient Richardson number instability form of mixing

  if (rg > 0){
    pwp_int <- grad_mix(pwp_int, params)
  }

  # MAKE SURE WE ADD A RETURN FUNCTION
  return(pwp_int)
}

# --------------------------------------------------------------

#' bulk_mix
#'
#' @param ml_index the index of the depth of the surface mixed layer after adjustment
#' @param rb
#' @param nz
#' @param z
#' @param d
#' @param u
#' @param v
#' @param g
#' @param t
#' @param s
#'
#' @return
#' @export
#'
#' @examples
bulk_mix <- function(pwp_int, params) {

  pwp_names <- names(pwp_int)
  list2env(pwp_int, env = environment())
  list2env(params, env = environment())

  rvc <- rb
  for (j in (ml_index+1):nz){
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
      mi <- mix5(j, t, s, d, u, v)
      list2env(mi)
    }
  }

  return(mget(pwp_names))
}

# ----------------------------------------------------------------

#' grad_mix
#'
#' @param rg
#' @param nz
#' @param d
#' @param u
#' @param v
#' @param g
#' @param dz
#' @param t
#' @param s
#'
#' @return
#' @export
#'
#' @examples
grad_mix <- function(pwp_int, params){

  pwp_names <- names(pwp_int)
  list2env(pwp_int, env = environment())
  list2env(params, env = environment())

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
  r <- rep(0,length(j1:j2)) #DONT KNOW IF TRANSLATED CORRECTLY FROM r = zeros(size(j1:j2))

  while (TRUE) {
    for (j in j1:j2){
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
    js <- min(which(r == rs))

    # Check to see whether the smallest r is critical or not

    if (rs > rc){
      break
    }

    # Mix the cells js and js+1 that had the smallest Richardson Number

    # MAKE SURE YOU CREATE AN OUTPUT FROM THE FUNCTION
    # LIKE WITH THE rot() FUNCTION
    st <- stir(rc, rs, js, s, t, d, u, v)
    list2env(st, env = environment())

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

  return(mget(pwp_names))

}

# ------------------------------------------------------------------------

#' stir
#'
#' @param rc
#' @param r
#' @param j
#' @param s
#' @param t
#' @param d
#' @param u
#' @param v
#'
#' @return
#' @export
#'
#' @examples
stir <- function(rc, rs, js, s, t, d, u, v){

  # This subroutine mixes cells j and j+1 just enough so that
  # the Richardson Number after the mixing is brought up to
  # the value rnew. In order to have this mixing process converge,
  # rnew must exceed the critical value of the Richardson number where
  # mixing is presumed to start. If r critical <- rc <- 0.25 (the nominal value)
  # and r <- 0.20, then
  # rnew <- 0.3 would be reasonable. If r were smaller, then a
  # larger value of rnew - rc is used to hasten convergence.

  rcon <- 0.02 + (rc - rs) / 2
  rnew <- rc + rcon / 5
  f <- 1 - rs / rnew
  dt <- (t[js + 1] - t[js]) * f / 2
  t[js + 1] <- t[js + 1] - dt
  t[js] <- t[js] + dt
  ds <- (s[js + 1] - s[js]) * f / 2
  s[js + 1] <- s[js + 1] - ds
  s[js] <- s[js] + ds
  d[js:(js + 1)] <- oce::swSigma(s[js:(js + 1)], t[js:(js + 1)], p = 0)
  du <- (u[js +1] - u[js]) * f / 2
  u[js + 1] <- u[js + 1] - du
  u[js] <- u[js] + du
  dv <- (v[js + 1] - v[js]) * f / 2
  v[js +1] <- v[js + 1] - dv
  v[js] <- v[js] + dv

  return(list(t = t, s = s, d = d, u = u, v = v))
}

# ----------------------------------------------------------------------

#' mix5
#'
#' @param j
#' @param t
#' @param s
#' @param d
#' @param u
#' @param v
#'
#' @return
#' @export
#'
#' @examples
mix5 <- function(j, t, s , d, u, v){

  # This subroutine mixes the arrays t, s, u, v down to level j

  t[1:j] <- mean(t[1:j])
  s[1:j] <- mean(s[1:j])
  d[1:j] <- oce::swSigma(s[1:j], t[1:j], p=0)
  u[1:j] <- mean(u[1:j])
  v[1:j] <- mean(v[1:j])

  return(list(t=t, s=s, d=d, u=u, v=v))
}

# ------------------------------------------------------------------

#' rot
#'
#' @param ang the angle the vector is rotated about
#' @param u
#' @param v
#'
#' @return
#' @export
#'
#' @examples
rot <- function(u,v,ang){

  # This subroutine rotates the vector (u,v) through an angle, ang

  r <-complex(1,u,v) * exp(complex(0,0,1) * ang)
  u <- Re(r)
  v <- Im(r)

  return(list(u = u, v = v))
}

# ---------------------------------------------------------------

#' remove_si
#'
#' @param t
#' @param s
#' @param d
#' @param u
#' @param v
#'
#' @return
#' @export
#'
#' @examples
remove_si <- function(pwp_int){

  # Find and relieve static instability that may occur in the
  # density array d. This simulates free convection.
  # ml_index is the index of the depth of the surface mixed layer
  # after adjustment
  pwp_names <- names(pwp_int)
  list2env(pwp_int, env = environment())

  while (1){
    mli <- which(diff(d) < 0)
    if(length(mli) == 0) {
      break
    }
    ml_index <- min(mli)
    a <- mix5(ml_index + 1, t, s , d, u, v)
    list2env(a,envir = environment())
  }
  return(mget(pwp_names))
}

# ---------------------------------------------------------------

#' absorb
#'
#' @param beta1 longwave extinction coefficient (m)
#' @param beta2 shortwave extinction coefficient (m)
#' @param nz
#' @param dz
#'
#' @return
#' @export
#'
#' @examples
absorb <- function(beta1, beta2, nz, dz){

  # Compute solar radiation absorbtion profile. This subroutine
  # assumes two wavelengths, and a double exponential depth for absorbtion

  # Subscript 1 is for red, non-penetrating light, and 2 is for blue,
  # penetrating light. rsl is the funcion assumed to be red

  rs1 <- 0.6
  rs2 <- 1.0 - rs1
  absrb <- rep(0,nz)
  z1 <- 0:(nz-1) * dz
  z2 <- z1 + dz
  z1b1 <- z1 / beta1
  z2b1 <- z2 / beta1
  z1b2 <- z1 / beta2
  z2b2 <- z2 / beta2
  absrb <- (rs1 * (exp(-z1b1) - exp(-z2b1) + rs2 * (exp(-z1b2) - exp(-z2b2))))
}


