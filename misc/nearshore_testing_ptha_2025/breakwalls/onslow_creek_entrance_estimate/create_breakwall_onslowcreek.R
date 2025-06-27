library(rptha)
myshp = readOGR('creek_entrance_geo/creek_entrance_geo.shp')

mypts = approxSpatialLines(myshp, n=200)
mypts = coordinates(mypts)
mypts = data.frame(lon=mypts[,1], lat=mypts[,2], z=rep(2., nrow(mypts)))
write.csv(mypts, file='onslow_creek_entrance_breakwall_estimate.csv', row.names=FALSE)
