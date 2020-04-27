met_input_file <- "inst/extdata/met_input_file_test_2.csv"
profile_input_file <- "inst/extdata/profile_input_file_test.csv"
dt = 900
dz = 10
days = 300
depth = 1000
dt_save = 4
lat = 60
rb = 0.64
rg = 0.25
rkz = 0
beta1 = 0.6
beta2 = 20

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

m<-1

vars <- ls()
vars <- vars[!(vars %in% c("pwp_in","params", "m"))]
rm(list = c(vars,"vars"))

pwp_output <- pwp_in

for (m in 1:params$nmet) {

  pwp_in <- pwpgo(pwp_in,params, m)

  if (m %% params$dt_save == 0) {
    pwp_output$time <- append(pwp_output$time, pwp_in$time)
    pwp_output$t <-  cbind(pwp_output$t, pwp_in$t)
    pwp_output$s <-  cbind(pwp_output$s, pwp_in$s)
    pwp_output$d <-  cbind(pwp_output$d, pwp_in$d)
    pwp_output$u <-  cbind(pwp_output$u, pwp_in$u)
    pwp_output$v <-  cbind(pwp_output$v, pwp_in$v)
  }


}


