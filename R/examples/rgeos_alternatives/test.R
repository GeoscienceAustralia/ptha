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

# gIntersection
y1 = as(st_geometry(x) + 1, 'Spatial')
proj4string(y1) = proj4string(y)
t1 = rgeos::gIntersection(y, y1)
t2 = gIntersection(y, y1)
plot(t1)
plot(t2, add=TRUE, border='red', lty='dotted')

# gIntersection by id
y1 = as(st_geometry(x) + 1, 'Spatial')
proj4string(y1) = proj4string(y)
t1 = rgeos::gIntersection(y, y1, byid=TRUE)
t2 = gIntersection(y, y1, byid=TRUE)
plot(t1)
plot(t2, add=TRUE, border='red', lty='dotted')

# gIntersection when there is no intersecting!!
y2 = as(st_geometry(x) + 2, 'Spatial')
proj4string(y2) = proj4string(y)
y2 = y2[1:2,]
t1 = rgeos::gIntersection(y, y2)
t2 = gIntersection(y, y2)
stopifnot(is.null(t1) & is.null(t2))
# As above with byid
t1 = rgeos::gIntersection(y, y2, byid=TRUE)
t2 = gIntersection(y, y2, byid=TRUE)


# gArea
t1 = rgeos::gArea(y)
t2 = gArea(y)
stopifnot(abs(t1 - t2) < 1e-08)

# gArea byID
t1 = rgeos::gArea(y, byid=TRUE)
t2 = gArea(y, byid=TRUE)
stopifnot(all(abs(t1 - t2) < 1e-08))

# gIntersects
y1 = as(st_geometry(x) + 1, 'Spatial')
proj4string(y1) = proj4string(y)
t1 = rgeos::gIntersects(y, y1)
t2 = gIntersects(y, y1)
stopifnot(all(t1 == t2) & t1)

# gIntersects by id
y1 = as(st_geometry(x) + 1, 'Spatial')
proj4string(y1) = proj4string(y)
y1 = y1[1:25,]
t1 = rgeos::gIntersects(y, y1, byid=TRUE)
t2 = gIntersects(y, y1, byid=TRUE)
stopifnot(all(t1 == t2))

# gContains
y1 = x[1:2,]
y1 = as(st_geometry(y1) + 1, 'Spatial')
proj4string(y1) = proj4string(y)
# Mixed example
y2 = x[1:25,] 
y2 = as(st_geometry(y2) + 1, 'Spatial') 
proj4string(y2) = proj4string(y)
t1 = rgeos::gContains(y, y1)
t2 = gContains(y, y1)
stopifnot((t1 == t2) & (!t1))
t1 = rgeos::gContains(y, y2)
t2 = gContains(y, y2)
stopifnot((t1 == t2) & (!t1))
# gContains, byid
t1 = rgeos::gContains(y, y1, byid=TRUE)
t2 = gContains(y, y1, byid=TRUE)
stopifnot(all(t1 == t2))
t1 = rgeos::gContains(y, y2, byid=TRUE)
t2 = gContains(y, y2, byid=TRUE)
stopifnot(all(t1 == t2))

# gCovers
y1 = x[1:2,]
y1 = as(st_geometry(y1) + 0.1, 'Spatial')
proj4string(y1) = proj4string(y)
y2 = gBuffer(gCentroid(y[1,]), width=0.01)

t1 = rgeos::gCovers(y, y1)
t2 = gCovers(y, y1)
stopifnot((t1 == t2) & (!t1))

t1 = rgeos::gCovers(y, y2)
t2 = gCovers(y, y2)
stopifnot((t1 == t2) & (t1))

t1 = rgeos::gCovers(y, y1, byid=TRUE)
t2 = gCovers(y, y1, byid=TRUE)
stopifnot(all((t1 == t2) & (!t1)))

t1 = rgeos::gCovers(y, y2, byid=TRUE)
t2 = gCovers(y, y2, byid=TRUE)
stopifnot(all((t1 == t2)) & t1[1] & !any(t1[-1]))

# gUnaryUnion
t1 = rgeos::gUnaryUnion(y)
t2 = gUnaryUnion(y)
plot(t1)
plot(t2, add=TRUE, border='red', lty='dotted')


t1 = rgeos::gUnaryUnion(y, id=y$dwndp_n)
t2 = gUnaryUnion(y, id=y$dwndp_n)
plot(t1)
plot(t2, add=TRUE, border='red', lty='dotted')

# gUnion
f1 = y[c(1,2),]
f2 = y[c(3,4),]
t1 = rgeos::gUnion(f1, f2)
t2 = gUnion(f1, f2)
plot(t1, border='black')
plot(t2, border='red', lty='dashed', add=TRUE)

# gUnion, byid=TRUE
t1 = rgeos::gUnion(f1, f2, byid=TRUE)
t2 = gUnion(f1, f2, byid=TRUE)
plot(t1, border='black')
plot(t2, border='red', lty='dashed', add=TRUE)

