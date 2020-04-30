library(pwpr)

rm(list=ls())

met_input_file <- "inst/extdata/met_input_file_test.csv"
profile_input_file <- "inst/extdata/profile_input_file_test.csv"

pwp_output <- pwp(met_input_file,profile_input_file, days = 50, dz = 10)

oce::imagep(t(pwp_output$t), ylim = c(100,0))

ggplot(pwp_output$t)

