library(R.matlab)
library(tidyverse)
library(oce)

# Met data

met <- lapply(readMat("inst/extdata/NCEP_timeseries_mooring.mat"), as.vector)
met <- as_tibble(met[2:length(met)])
NCEP_met_test <- mutate(met, time = as.POSIXct(time*60*60*24, tz = "UTC", origin = "0000-01-01"))

usethis::use_data(NCEP_met_test,overwrite = TRUE)


# Mooring data

mooring <- readMat("inst/extdata/IrmingerSeaMMPs/MMP_offshore_0203.mat")$MP
moorctd <- vector("list",dim(mooring[[1]])[1])
for (i in 1:dim(mooring[[1]])[1]) {
  mooradd <- tibble(time = as.POSIXct(mooring[[1]][i,]*60*60*24, origin = "0000-01-01"),
                    press = mooring[[2]][i,],
                    temp = mooring[[3]][i,],
                    sal = mooring[[4]][i,],
                    prof = i)
  moorctd[[i]] <- ctdDecimate(as.ctd(salinity = mooradd$sal,
                    temperature = mooradd$temp,
                    pressure = mooradd$press,
                    time=mooradd$time,
                    station = i),p=5)
  mooradd <- tibble(time = as.POSIXct(moorctd[[i]]@data$time, tz = "UTC", origin = "1970-01-01"),
                    press = moorctd[[i]]@data$pressure,
                    dep = swDepth(moorctd[[i]], lat = 60),
                    temp = moorctd[[i]]@data$temperature,
                    sal = moorctd[[i]]@data$salinity,
                    theta = swTheta(moorctd[[i]]),
                    sigma = swSigmaTheta(moorctd[[i]]))
  if (i == 1) {
    moor <- mooradd
  } else {
    moor <- bind_rows(moor, mooradd)
  }
}

mooring_test <- moor

usethis::use_data(mooring_test, overwrite = TRUE)

