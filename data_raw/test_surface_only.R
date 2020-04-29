library(pwpr)

rm(list=ls())

met_input_file <- "inst/extdata/met_input_file_test_2.csv"
profile_input_file <- "inst/extdata/profile_input_file_test.csv"

pwp_output <- pwp(met_input_file,profile_input_file, days = 30)

oce::imagep(t(pwp_output$t), ylim = c(100,0))
