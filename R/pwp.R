# Load libraries
library(tidyverse)
library(oce)

# Initialize the user-defined variables

# dt:   time-step increment (seconds)
# dz:   depth increment (meters)
# days:   the number of days to run
#depth:   the depth to run
#dt_save:   time-step increment for saving to file (multiples of dt)
#lat:   latitude (degrees)
#g:   gravity (9.8 m/s^2)
#cpw:   specific heat of water (4183.3 J/kgC)
#rb:  critical bulk richardson number (0.65)
#rg:  critical gradient richardson number (0.25)
#rkz:   background diffusion (0)
#beta1:   longwave extinction coefficient (0.6 m)
#beta2:   shortwave extinction coefficient (20 m)

pwp <- function(met_input_file, profile_input_file,
                dt = 900, dz = 10) {
days <- 300
depth <- 1000
dt_save <- 4
lat <- 60
g <- 9.8
cpw <- 4183.3
rb <- 0.65
rg <- 0.25
rkz <- 0
beta1 <- 0.6
beta2 <- 20

# Initialize additional variables

#f: coriolis parameter
#ucon: coefficient of inertial-internal wave dissipation

f <- 2 * 7.29E-5 * sin(lat * pi / 180)

# Format air-sea flux data --------------------------------------------------

# Load air-sea forcing data from csv file
met <- read_csv()

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
prof <- read.csv()

# Interpolate the profile variables at dz resolution
# Set nz equal to the number of depth increments + 1 using a resolution of dz

nz <- 1 + depth/dz
z <- 0:(nz - 1) * dz

t <- approx(prof$z, prof$t, z)
s <- approx(prof$z, prof$s, z)
d <- swSigma(s,t)

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

# Define the variables to be saved

pwp_output$dt <- dt
pwp_output$dz <- dz
pwp_output$lat <- lat
pwp_output$z <- z
pwp_output$time <- []
pwp_output$t <- []
pwp_output$s <- []
pwp_output$d <- []
pwp_output$u <- []
pwp_output$v <- []

return(pwp_output)
}

# Step through the PWP model

disp(['STATUS (out of' int2str(nmet)' steps):'])

for m = 1:nmet

# ------------------------------------------------------------

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

# --------------------------------------------------------------

# function rot(ang)
# This subroutine rotates the vector (u,v) through an angle, ang

r <- (u + i * v) * exp (i * ang)
u <- real(r)
v <- imag(r)

# ---------------------------------------------------------------

# function absrb = absorb(beta1,beta2)

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

#plz work
