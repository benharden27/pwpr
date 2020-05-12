# Create met and profile input files
library(tidyverse)
library(pwpr)

# profile data
ti <- 5
df <- filter(mooring_test, station == ti) %>%
  select(time, dep, temp, sal)

# fill to surface
ii <- which(!is.na(df$temp))[1] - 1
df$temp[1:ii] <- df$temp[ii+1]
df$sal[1:ii] <- df$sal[ii+1]
df$time[1:ii] <- df$time[ii+1]

timein <- df$time[1]


df <- transmute(df, z = dep, t = temp, s = sal) %>%
  mutate(d = oce::swSigma(s,t,p=0)) %>%
  arrange(d) %>%
  mutate(z = sort(z))

write_csv(df, "inst/extdata/irminger_initial_profile.csv")

# met data
# put all in a tibble
df <- tibble(dttm = NCEP_met_test$time,
             tx = NCEP_met_test$ums,
             ty = NCEP_met_test$vms,
             qlat = NCEP_met_test$lhf,
             qsens = NCEP_met_test$shf,
             lw = NCEP_met_test$nlwrs,
             sw = -NCEP_met_test$nswrs,
             precip = NCEP_met_test$prate/1000/3600)

#
ti <- sea::find_near(df$dttm, timein)
df <- mutate(df, time = ((1:nrow(df))-ti)/4) %>%
  filter(time >= 0) %>%
  mutate(precip = 0)

write_csv(df,"inst/extdata/ncep_forcing.csv")
