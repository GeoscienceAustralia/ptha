library(rgeos)
library(sf)
source('alternatives_rgeos.R')

# Test data
x = read_sf('hjort/hjort.shp')
y = as(x, 'Spatial')

# gBuffer byid=FALSE
t1 = rgeos::gBuffer(y, width=0.2)
t2 = gBuffer(y, width=0.2)
plot(t1, border='black')
plot(t2, border='red', lty='dashed', add=TRUE)

# gBuffer byid=TRUE
t1 = rgeos::gBuffer(y, width=0.2, byid=TRUE)
t2 = gBuffer(y, width=0.2, byid=TRUE)
plot(t1, border='black')
plot(t2, border='red', lty='dashed', add=TRUE)

# gCentroid, byid=FALSE
t1 = rgeos::gCentroid(y)
t2 = gCentroid(y)
plot(t1, col='black', pch=19)
plot(t2, col='red', pch=1, add=TRUE)

# gCentroid, byid=TRUE
t1 = rgeos::gCentroid(y, byid=TRUE)
t2 = gCentroid(y, byid=TRUE)
plot(t1, col='black', pch=19)
plot(t2, col='red', pch=1, add=TRUE)
