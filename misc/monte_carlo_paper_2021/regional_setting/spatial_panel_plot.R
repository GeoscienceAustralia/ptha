#
library(rptha) # for readOGR
library(raster)

# Background data 
kt_grid = readOGR('kermadectonga2.shp', layer='kermadectonga2')
## These datasets are not provided here -- the former is the DEM used in PTHA18, while the latter is a zero-contour
## derived from that.
dem   = raster('../../../../../AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')
zero_contour = '../../../../../AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/zero_contour/zero_contour.shp'
zc = readOGR(zero_contour, layer='zero_contour')

##
panel_plot<-function(){

    xrange = c(170, 195)
    yrange = c(-45, -5)

    par(mar=c(3,3,1,1))
    plot(c(0, 1), c(0, 1), col='white',
        cex=0.1, asp=1/cos(mean(yrange)/180*pi), 
        xlim=xrange, ylim=yrange, 
        main="", xlab="Lon", ylab="Lat", cex.lab=1.5, cex.axis=1.3)
    # Nice stuff
    image(dem, add=TRUE, zlim=c(-9000, 0), col=grey(seq(0,1,len=100), alpha=0.4), maxpixels=1e+07)
    plot(zc, add=TRUE)
    plot(kt_grid, col=rgb(1, 0, 0, alpha=0.5), border='black', add=TRUE)

    return(invisible(0))
}



#
# Make a plot comparing all the appproaches 
#
png('kermadectonga2_region.png', width=3.5, height=5, units='in', res=200)
panel_plot()
dev.off()



#
# Make a satellite plot with Tongatapu and sites of interest.
#
library(OpenStreetMap); library(sp)

upperLeft  = c(-20.8, 184.0)
lowerRight = c(-21.8, 185.4)
osm_backdrop = openmap(upperLeft = upperLeft, lowerRight= lowerRight, type='bing')
osm_backdrop_reproj = openproj(osm_backdrop, proj4string(CRS("+init=epsg:4326")))
osm_backdrop_reproj$bbox$p1[1] = osm_backdrop_reproj$bbox$p1[1] + 360
osm_backdrop_reproj$bbox$p2[1] = osm_backdrop_reproj$bbox$p2[1] + 360
osm_backdrop_reproj$tiles[[1]]$bbox$p1[1] = osm_backdrop_reproj$tiles[[1]]$bbox$p1[1] + 360
osm_backdrop_reproj$tiles[[1]]$bbox$p2[1] = osm_backdrop_reproj$tiles[[1]]$bbox$p2[1] + 360

png('Tongatapu_zoom.png', width=6, height=5.275, units='in', res=300)
par(mar=c(3,3,1.3,1))
plot(c(184.3, 185.3), c(-21.6, -20.8), col='white', asp=1/cos(-21.4/180*pi), xlab="", ylab="", xaxs='i', yaxs='i',
    cex.axis=1.4)
plot(osm_backdrop_reproj, add=TRUE)
plot(kt_grid, col=rgb(1, 0, 0, alpha=0.5), border=rgb(0,0,0,alpha=0.5), density=10, add=TRUE)

test_point_list = list(
    site_P = c(185.1239, -21.0888),
    test_point_1 = c(184.8292, -20.9286),
    test_point_2 = c(184.6758, -21.3586),
    test_point_3 = c(184.47748, -21.041364),
    test_point_4 = c(185.05833, -21.44512),
    test_point_NZ =  c(178.3945, -37.3940) 
    )

points(test_point_list$site_P[1], test_point_list$site_P[2], pch=19, col='white', cex=2)
text(test_point_list$site_P[1], test_point_list$site_P[2], 'Site P', col='white', cex=1.8, pos=3)
points(test_point_list$test_point_1[1], test_point_list$test_point_1[2], pch=19, col='white', cex=2)
text(test_point_list$test_point_1[1], test_point_list$test_point_1[2], 'Nearby site 1', col='white', cex=1.8, pos=3)
points(test_point_list$test_point_2[1], test_point_list$test_point_2[2], pch=19, col='white', cex=2)
text(test_point_list$test_point_2[1], test_point_list$test_point_2[2], 'Nearby site 2', col='white', cex=1.8, pos=3)
points(test_point_list$test_point_3[1], test_point_list$test_point_3[2], pch=19, col='white', cex=2)
text(test_point_list$test_point_3[1], test_point_list$test_point_3[2], 'Nearby site 3', col='white', cex=1.8, pos=3)
dev.off()
