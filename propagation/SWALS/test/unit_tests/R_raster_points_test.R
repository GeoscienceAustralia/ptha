library(raster)
rast = raster('data/PointData_elevation_c_INITIAL.tif')

x = 500000 + runif(10)*1000
y = 1610000 + runif(10)*1000

zL = extract(rast, cbind(x,y))
zB = extract(rast, cbind(x,y), method='bilinear')


## Make a new test raster

dx = 200
dy = 200
xmn = 499000
xmx = 509000
ymn = 1609000
ymx = 1619000

xs = seq(xmn, xmx - dx, by=dx) + dx/2
ys = seq(ymn, ymx - dy, by=dy) + dy/2

zs = matrix(NA, ncol = length(xs), nrow=length(ys))
for(i in 1:length(xs)){
    for(j in 1:length(ys)){
        zs[i,j] = xs[i] + 2*ys[j]
    }
}

test_rast = raster(zs, xmn, xmx, ymn, ymx, crs=CRS("+init=epsg:3123"))

writeRaster(test_rast, 'data/test_rast.tif', format='GTiff', overwrite=TRUE)

x1 = xmn + runif(10)*(xmx - xmn)
y1 = ymn + runif(10)*(ymx - ymn)

options(digits=12)

print(x1)
print(y1)

z1 = extract(test_rast, cbind(x1, y1))
print(z1)

z1 = extract(test_rast, cbind(x1, y1), method='bilinear')
print(z1)

