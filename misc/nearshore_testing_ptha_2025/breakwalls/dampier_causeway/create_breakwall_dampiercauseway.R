library(rptha)
myshp = readOGR('causeway_geo/causeway_geo.shp')

mypts = approxSpatialLines(myshp, n=200)
mypts = coordinates(mypts)
mypts = data.frame(lon=mypts[,1], lat=mypts[,2], z=rep(5., nrow(mypts)))
write.csv(mypts, file='dampier_causeway_breakwall_estimate.csv', row.names=FALSE)
