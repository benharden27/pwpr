# Create test csv files for using to troubleshoot pwp functions
library(tidyverse)
library(seapacific)

# Met CSV -----------------------------------------------------------------

# t is in days
nx = 100
df <- tibble(time = seq(0,nx-1,1),
                 sw = runif(100,0,50),
                 lw = runif(100,0,30),
                 qlat = runif(100,0,20),
                 qsens = runif(100,0,20),
                 tx = runif(100,0,5e-2),
                 ty = runif(100,0,5e-2),
                 precip = runif(100,0,1e-8)
                 )

write_csv(df, "inst/extdata/met_input_file_test.csv")

df <- tibble(time = seq(0,nx-1,1),
             sw = 0,
             lw = 50,
             qlat = 0,
             qsens = 0,
             tx = 0,
             ty = 0,
             precip = 0
)

write_csv(df, "inst/extdata/met_input_file_test_2.csv")

df <- tibble(time = seq(0,nx-1,1),
             sw = 0,
             lw = 50,
             qlat = 0,
             qsens = 0,
             tx = 0,
             ty = 0,
             precip = 0
)

write_csv(df, "inst/extdata/met_input_file_test_2.csv")



# Profile CSV -------------------------------------------------------------

df <- tibble(z = S285$ctd[[10]]@data$depth,
             t = S285$ctd[[10]]@data$temperature,
             s = S285$ctd[[10]]@data$salinity)
df <- df[6:nrow(df),]

write_csv(df, "inst/extdata/profile_input_file_test.csv")
