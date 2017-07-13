# 
# Routines to extract data from e.g. the CMT catalogue, or the ISC-GEM catalogue, 
# within polygons that define our source-zones
#
library(rptha)
config = new.env()
source('config.R', local=config)

# Key inputs
cmt_catalogue_csv = config$cmt_catalogue_csv 
unit_source_grid_polygon_shapefiles = config$unit_source_grid_polygon_shapefiles 

# Properties of earthquake events we extract
mw_threshold = config$MW_MIN 
depth_threshold = config$depth_threshold #70 # depth < depth_threshold
buffer_width = config$buffer_width # Events are inside polygon, after polygon is buffered by 'buffer_width' degrees
# Events have rake_min <= rake1 <= rake_max; OR rake_min <= rake2 <= rake_max
rake_min = 90 - config$rake_deviation_thrust_events  
rake_max = 90 + config$rake_deviation_thrust_events  

#
# Read all unit-source grid shapefiles into a list, with names corresponding to
# source-zones
#
unit_source_grid_poly = lapply(
    as.list(unit_source_grid_polygon_shapefiles), 
    f<-function(x) readOGR(x, layer=gsub('.shp', '', basename(x))))

names(unit_source_grid_poly) = basename(dirname(dirname(dirname(
    unit_source_grid_polygon_shapefiles))))

#
# Read earthquake observational datasets
#
gcmt = read.csv(cmt_catalogue_csv)
days_in_year = 365.25
cmt_duration_years = diff(range(gcmt$julianDay1900))/days_in_year

#'
#' Determine whether lonlat coordinates are inside a given polygon
#'
#' Takes care of different longitude conventions (e.g. offsets by 360n).
#'
#' If buffer_width is provided, then poly is buffered by this amount prior
#' to testing inclusion. Buffer_width should be in the same units as poly's coordinates.
#'
#' @param lonlat 2 column matrix with longitude/latitude
#' @param poly SpatialPolygons object
#' @param buffer width thickness of buffer to apply to poly before testing the point inclusion.
#' @return logical vector with one entry for every row in lonlat.
#'
lonlat_in_poly<-function(lonlat, poly, buffer_width = 0){

    if(buffer_width != 0){
        poly = gBuffer(poly, width=buffer_width, byid=TRUE)
    }

    # Find point on poly, which is used to determine the longitude convention of 'poly'
    point_in_poly = as.numeric(coordinates(gCentroid(poly)))

    # Extend to same dimensions as lonlat
    point_in_poly = matrix(point_in_poly, nrow=length(lonlat[,1]), ncol=2, byrow=TRUE)

    # Make sure longitude convention is the same as for the polygon
    lonlat_near_poly = adjust_longitude_by_360_deg(lonlat, point_in_poly)

    lonlat_near_poly = SpatialPoints(lonlat_near_poly, proj4string=CRS(proj4string(poly)))

    inside_poly = over(lonlat_near_poly, as(poly, 'SpatialPolygons'))

    inside_poly = !is.na(inside_poly)

    return(inside_poly)

}

#
# Extract earthquake events > threshold for all sources
#

eq_events = vector(mode='list', length=length(unit_source_grid_poly))
for(i in 1:length(unit_source_grid_poly)){
   
    # Count as inside if EITHER the hypocentre, or the centroid, is inside 
    inside_events_hypo = lonlat_in_poly(
        gcmt[,c('hypo_lon', 'hypo_lat')], unit_source_grid_poly[[i]], buffer_width=buffer_width)
    inside_events_centroid = lonlat_in_poly(
        gcmt[,c('cent_lon', 'cent_lat')], unit_source_grid_poly[[i]], buffer_width=buffer_width)
   
    inside_events = (inside_events_hypo | inside_events_centroid)
 
    # Criterion for point selection - note we will keep the event if either rake1 or rake2 is within pi/4 of pure thrust
    inside_keep = which( inside_events & (gcmt$Mw >= mw_threshold) & (gcmt$depth <= depth_threshold) & 
        ((gcmt$rake1 >= rake_min & gcmt$rake1 <= rake_max) | (gcmt$rake2 >= rake_min & gcmt$rake2 <= rake_max) )
        )

    eq_events[[i]] = gcmt[inside_keep,]
}
names(eq_events) = names(unit_source_grid_poly)

#
# Get rate of events > threshold on source-zone
#

source_zone_exceedance_rate = sapply(eq_events, f<-function(x) length(x[,1])/cmt_duration_years)

# time between events [ignoring first/last points issues where the time is
# uncertain]
source_zone_time_differences = lapply(eq_events, f<-function(x) diff(x$julianDay1900)/days_in_year )

#
# EXPLORATORY ANALYSES BELOW HERE
#

# Gutenberg Richter cumulative distribution
GR<-function(x, b, mw_min) 1-10^(b*(mw_min - x))

# Gutenberg Richter pdf
gR<-function(x, b, mw_min) b*log(10)*10^(b*(mw_min - x))

negloglik_GR<-function(b, x, mw_min){
    -sum(log(gR(x, b, mw_min)))
}

## Find events within 0.2 degrees of a source
print_b_estimates<-function(poly){

    xx_c = lonlat_in_poly(gcmt[,c('cent_lon', 'cent_lat')], poly, buffer_width=0.2)
    xx_h = lonlat_in_poly(gcmt[,c('hypo_lon', 'hypo_lat')], poly, buffer_width=0.2)

    xx = (xx_c | xx_h)

    # Fit 'b' to the data for different possible Mw-min's
    for(mm in seq(5., 6.5, by=0.1)){

        kk = which(xx & gcmt$Mw > mm & gcmt$depth <= depth_threshold)

        out = optimize(negloglik_GR, lower=0.5, upper=1.5, x=gcmt$Mw[kk], mw_min=mm)
        out_exact = 1/(log(10)*mean(gcmt$Mw[kk] - mm))

        test_b_CI = seq(0.2, 2.5, by=0.01)

        out_multiple = sapply(test_b_CI, f<-function(b) negloglik_GR(b, x=gcmt$Mw[kk], mw_min=mm) )
        out_fun = approxfun(test_b_CI, out_multiple)
        out_fun_offset<-function(x) (out_fun(x) - out$objective - qchisq(0.95, df=1)/2 )
        out_lower = uniroot(out_fun_offset, lower=min(test_b_CI), upper=out$minimum)
        out_upper = uniroot(out_fun_offset, lower=out$minimum, upper=max(test_b_CI))

        # mag, N, solution, solution, lower profile likelihood CI, upper profile likelihood CI, analytical sd
        print(c(mm, length(kk), out_exact, out$minimum, out_lower$root, out_upper$root, out$minimum/sqrt(length(kk))))

    }
}

# Normalise time differences by rate of events [e.g. Geist et al, 2011, and
# much earthquake literature]
rate_normalised_sztd = source_zone_time_differences
for(i in 1:length(rate_normalised_sztd)){
    rate_normalised_sztd[[i]] = rate_normalised_sztd[[i]]*source_zone_exceedance_rate[i]
}

# Try exponential and gamma fit to the entire 'rate-normalised'
# time-between-events dataset
library(fitdistrplus)
mm = unlist(rate_normalised_sztd)
# One zero -- must be removed for fitting -- but check it
mm[which(mm==0)] = min(mm[which(mm>0)])
# exponential (poisson process)
exp_fit_mle = fitdist(data=mm, distr='exp', method='mle')
# gamma
gamma_fit_mme = fitdist(data=mm, distr='gamma', method='mme')
gamma_fit_mle = fitdist(data=mm, distr='gamma', start = as.list(gamma_fit_mme$estimate))
gamma_fit_mge = fitdist(data=mm, distr='gamma', method='mge')
# weibull
weibull_fit_mle = fitdist(data=mm, distr='weibull', method='mle', start=list(shape=1, scale=1))

#
# Goodness of fit tests
#
gofstat(gamma_fit_mle)
gofstat(weibull_fit_mle)
gofstat(exp_fit_mle)

pdf('test.pdf', width=10, height=8)
par(mfrow=c(2,2))
qqcomp(exp_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Exponential model QQ plot')
qqcomp(gamma_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Gamma model QQ plot')
qqcomp(weibull_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Weibull model QQ plot')

par(mfrow=c(3,3))
mw_store = c()
time_diff_store = c()
for(i in 1:length(eq_events)){
    try(plot(eq_events[[i]]$Mw, c(rate_normalised_sztd[[i]], NA), main=names(eq_events)[i]) )
    mw_store = c(mw_store, eq_events[[i]]$Mw)
    time_diff_store = c(time_diff_store, rate_normalised_sztd[[i]], NA)
}

#
# Is there a relation between the normalised time between events, and the magnitude of the first event?
#
par(mfrow=c(1,1))
plot(range(mw_store, na.rm=TRUE), range(time_diff_store[time_diff_store>0], na.rm=TRUE), col=0, 
    xlab='Mw', ylab='Rate normalised time between events', log='y')
for(i in 1:length(eq_events)){
    try(points(eq_events[[i]]$Mw, c(rate_normalised_sztd[[i]], NA), pch=i, col=rainbow(9)[i], cex=2))
}
legend('topright', names(eq_events), pch=1:8, col=rainbow(9)[1:8], pt.cex=2)

# Some relationship ?
mw_time_diff_test = cor.test(mw_store, time_diff_store, method='s')


#
# Spatial plots
#

library(RFOC)

for(i in 1:length(eq_events)){

    plot(unit_source_grid_poly[[i]], main=names(unit_source_grid_poly)[i], border='grey', axes=TRUE)
    points(eq_events[[i]]$cent_lon, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')
    # Make sure we don't miss points due to longitude convention
    points(eq_events[[i]]$cent_lon-360, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')
    points(eq_events[[i]]$cent_lon+360, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')

    #
    # Same plot as above, with beachballs
    #

    plot(unit_source_grid_poly[[i]], main=names(unit_source_grid_poly)[i], border='grey', axes=TRUE)

    pol_centroid = as.numeric(coordinates(gCentroid(unit_source_grid_poly[[i]])))

    for(j in 1:length(eq_events[[i]]$cent_lon)){

        lon_lat = c(eq_events[[i]]$cent_lon[j], eq_events[[i]]$cent_lat[j])
        lon_lat = adjust_longitude_by_360_deg(lon_lat, pol_centroid)
    
        # Beach-ball details
        mec = SDRfoc(eq_events[[i]]$strk1[j], eq_events[[i]]$dip1[j], eq_events[[i]]$rake1[j], PLOT=FALSE)

        fcol =  c( rgb(1,0,0, alpha=1), 'lightblue', 'green', 'orange', 'yellow', 'purple', 'black')[foc.icolor(mec$rake1)]

        par(lwd=0.2) # Try to avoid strong lines on beachballs
        try(
            justfocXY(mec, lon_lat[1], lon_lat[2],
                focsiz=4*((eq_events[[i]]$Mw[j]-6)/10)**2, xpd=TRUE, fcol=fcol,
                fcolback='white'),
            silent=TRUE
            )
        par(lwd=1)

    }
}

dev.off()

