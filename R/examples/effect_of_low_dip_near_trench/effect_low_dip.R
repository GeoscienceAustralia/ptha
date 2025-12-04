#
# Compare pure-thrust Okada unit sources at the trench with 5 degree and 10 degree dip
#

library(rptha)

xs = seq(-100.5,100.5,len=100)
ys = seq(-100.5,100.5,len=100)
mypts=expand.grid(xs, ys)*1000

width = 50
length = 50

dip = 5
centroid_depth = width/2*sin(dip/180*pi) # Rupture to trench
dip_5 = okada_tsunami(0,0,centroid_depth,0,dip,length,width,0,1, mypts[,1],mypts[,2])
dip_5_dz = dip_5[[3]]
dim(dip_5_dz) = c(length(xs), length(ys))

dip =10 
centroid_depth = width/2*sin(dip/180*pi) # Rupture to trench
dip_10 = okada_tsunami(0,0,centroid_depth,0,dip,length,width,0,1, mypts[,1],mypts[,2])
dip_10_dz = dip_10[[3]]
dim(dip_10_dz) = c(length(xs), length(ys))

dip =15 
centroid_depth = width/2*sin(dip/180*pi) # Rupture to trench
dip_15 = okada_tsunami(0,0,centroid_depth,0,dip,length,width,0,1, mypts[,1],mypts[,2])
dip_15_dz = dip_15[[3]]
dim(dip_15_dz) = c(length(xs), length(ys))


library(fields)
par(mfrow=c(1,3))
#par(oma=c(0,1,0,1))
ZLIM = c(-0.3, 0.3)
image.plot(xs, ys, dip_5_dz, zlim=ZLIM, main=paste0('Dip = 5 degrees, max = ', round(max(dip_5_dz), 2)), asp=1)
image.plot(xs, ys, dip_10_dz, zlim=ZLIM, main=paste0('Dip = 10 degrees, max = ', round(max(dip_10_dz), 2)), asp=1)
image.plot(xs, ys, dip_15_dz, zlim=ZLIM, main=paste0('Dip = 15 degrees, max = ', round(max(dip_15_dz), 2)), asp=1)

