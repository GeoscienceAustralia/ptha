library(devtools)
document('.')
document('.')
build('.')

setwd('..')
system('R CMD check rptha_0.0.tar.gz')
