#
# Basilisk wave problem
# 

x = read.table('db.txt')

# Indices at t=0
k = which(x[,1] == 0)
# Indices at t=tend
k1 = which(x[,1] > 0)

png('Basilisk_wave_after_100km_through_domain.png', width=10, height=4, units='in', res=300)
plot(x[k,2], x[k,3])
points(x[k1,2], x[k1,3], col='red')
title(paste0('Range of final wave as percentage of initial range = ', diff(range(x[k1,3]))/diff(range(x[k,3]))))
dev.off()

