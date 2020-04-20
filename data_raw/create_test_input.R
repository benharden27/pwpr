# Create test csv files for using to troubleshoot pwp functions
library(tidyverse)
library(seapacific)

# Met CSV -----------------------------------------------------------------

# t is in days
nx = 100
df <- tibble(time = seq(0,nx-1,1),
                 sw = runif(100,0,1000),
                 lw = runif(100,0,1000),
                 qlat = runif(100,0,1000),
                 qsens = runif(100,0,1000),
                 tx = runif(100,0,1000),
                 ty = runif(100,0,1000),
                 precip = runif(100,0,1000)
                 )

write_csv(df, "inst/extdata/met_input_file_test.csv")



# Profile CSV -------------------------------------------------------------

df <- tibble(z = S285$ctd[[10]]@data$depth,
             t = S285$ctd[[10]]@data$temperature,
             s = S285$ctd[[10]]@data$salinity)
df <- df[6:nrow(df),]

write_csv(df, "inst/extdata/profile_input_file_test.csv")
