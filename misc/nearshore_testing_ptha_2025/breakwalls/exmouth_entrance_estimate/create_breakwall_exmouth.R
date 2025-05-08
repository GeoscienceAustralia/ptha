library(rptha)
myshp = readOGR('breakwall_location/exmouth_entrance_estimate.shp')

mypts = approxSpatialLines(myshp, n=200)
mypts = coordinates(mypts)
mypts = data.frame(lon=mypts[,1], lat=mypts[,2], z=rep(3.1, nrow(mypts)))
write.csv(mypts, file='exmouth_entrance_breakwall_estimate.csv', row.names=FALSE)
