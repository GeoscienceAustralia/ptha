#
# Convert the 'se.dat' boundary condition file into a good format.
#
se_dat = read.table('se.dat')

new_dat = se_dat
new_dat[,1] = (se_dat[,1] - se_dat[1,1]) * 60.0
new_dat[,2] = new_dat[,2] - new_dat[1,2]
names(new_dat) = c("time", "stage")

write.csv(new_dat, file='se_dat_converted.csv', row.names=FALSE)
