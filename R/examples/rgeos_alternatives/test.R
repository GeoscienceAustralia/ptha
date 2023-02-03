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

# gDistance
y2 = as(st_geometry(x) + 2, 'Spatial')
proj4string(y2) = proj4string(y)
t1 = rgeos::gDistance(y, y2)
t2 = gDistance(y, y2)
stopifnot(abs(t1 - t2) < 1e-08)

# gDistance byID
y2 = as(st_geometry(x) + 2, 'Spatial')
proj4string(y2) = proj4string(y)
t1 = rgeos::gDistance(y, y2, byid=TRUE)
t2 = gDistance(y, y2, byid=TRUE)
stopifnot(all(abs(t1 - t2) < 1e-08))

