#
# Plotting of BP09 results
#
image_flags = commandArgs(trailingOnly=TRUE)
if(length(image_flags) != 1){
    stop('Commandline arguments not recognized')
}

#
# Read field data. Note there are obvious spatial registration problems with
# these datasets relative to the DEMs used for modelling (discussed in the
# GEOCLAW validation report). Here we try to correct those, but it isn't perfect.
#

# Optionally treat "bad" observations, with a few different options
BAD_OBS_TREATMENT = 'Project'
stopifnot(any(BAD_OBS_TREATMENT == c('None', 'Drop', 'Project')))

# Optionally add source contours to the elevation plot, and make a zoomed elevation plot
PLOT_EXTRAS = FALSE

#
# Tohoku field data
#
field_data_tohoku = read.csv('../test_repository/BP09-FrankG-Okushiri_island/test_data/FieldData_Tohoku.csv')
# Approximate fix for systematic registration errors (see the GEOCLAW
# validation report)
tohoku_lonlat_adjust = c(0.001819111, -0.01124044) 

field = data.frame( lon = field_data_tohoku[,8] + field_data_tohoku[,9]/60 + field_data_tohoku[,10]/(60*60) + tohoku_lonlat_adjust[1],
                    lat = field_data_tohoku[,5] + field_data_tohoku[,6]/60 + field_data_tohoku[, 7]/(60*60) + tohoku_lonlat_adjust[2],
                    z_uncorrected = field_data_tohoku[,2], 
                    z_corrected = field_data_tohoku[,4],
                    tide = field_data_tohoku[,3],
                    quality = field_data_tohoku[,1] )
field[field==99999] = NA
# Convert cm to m
field$z_uncorrected = field$z_uncorrected * 1/100
field$z_corrected = field$z_corrected * 1/100

#
# UJNR field data
#
field_data_ujnr = read.csv('../test_repository/BP09-FrankG-Okushiri_island/test_data/FieldData_UJNR.csv', na.strings=c("", 'n/a'))
# Approximate fix for systematic registration errors
ujnr_lonlat_adjust = c(0.0002246666, -0.011496)
field_ujnr = data.frame(
        lon = field_data_ujnr[,6] + field_data_ujnr[,7]/60 + ujnr_lonlat_adjust[1],
        lat = field_data_ujnr[,4] + field_data_ujnr[,5]/60 + ujnr_lonlat_adjust[2],
        z = field_data_ujnr[,9])

#
# TSUJI field data -- note this contains some rows with high values, but no lon-lat. Perhaps they are 
# poorly constrained regional values? Unclear, but might be worth following up
#
field_data_tsuji = read.csv('../test_repository/BP09-FrankG-Okushiri_island/test_data/FieldData_Tsuji.csv', 
    na.strings=c("", 'n/a'))
tsuji_lonlat_adjust = c(0.005152444, -0.01318489 )
# Sometimes a range of elevations is provided - we use the max
field_tsuji = data.frame(
    lon = field_data_tsuji[,7] + field_data_tsuji[,8]/60 + field_data_tsuji[,9]/(60*60) + tsuji_lonlat_adjust[1],
    lat = field_data_tsuji[,4] + field_data_tsuji[,5]/60 + field_data_tsuji[,6]/(60*60) + tsuji_lonlat_adjust[2],
    # Use the max elevation value if a range is provided
    z_max = sapply(field_data_tsuji[,3], f<-function(x) max(as.numeric(strsplit(as.character(x), '-')[[1]]))))


#
# Get data
#

source('../../../plot.R')
eps=1.0e-03

MD_DIR = rev(Sys.glob(paste0('OUTPUTS/RUN_*')))[1]

get_domain_contour<-function(domain_number){

    most_recent_run = MD_DIR

    stage_elev_rasters = vector(mode='list', length=2)
    stage_elev_rasters[[1]] = merge_domains_nc_grids(multidomain_dir = most_recent_run, 
        domain_index=domain_number, desired_var = 'max_stage', return_raster=TRUE)
    stage_elev_rasters[[2]] = merge_domains_nc_grids(multidomain_dir = most_recent_run, 
        domain_index=domain_number, desired_var = 'elevation0', return_raster=TRUE)
    xs = xFromCol(stage_elev_rasters[[1]], 1:ncol(stage_elev_rasters[[1]])) 
    ys = yFromRow(stage_elev_rasters[[1]], nrow(stage_elev_rasters[[1]]):1) 

    wet_or_dry = as.matrix(stage_elev_rasters[[1]]) > (as.matrix(stage_elev_rasters[[2]]) + eps)
    # Flip around into the desired orientation
    wet_or_dry = t(wet_or_dry)
    wet_or_dry = wet_or_dry[,ncol(wet_or_dry):1]

    # Get wet-dry contours
    all_wet_dry_contours = contourLines(xs, ys, wet_or_dry, level=0.99999)
    
    # Look at the peak-stage along the contour
    peak_stage = lapply(all_wet_dry_contours, f<-function(x) extract(stage_elev_rasters[[1]], cbind(x$x, x$y), method='simple'))
    peak_elev = lapply(all_wet_dry_contours, f<-function(x) extract(stage_elev_rasters[[2]], cbind(x$x, x$y), method='simple'))

    return(environment())
}

#
# Plot of survey heights
#

centroid = c(139.47, 42.14) # Center for radial height plots
# Useful function to add radial arrows to plot
add_arrows<-function(lon, lat, len, centroid, coastal_points, ...){

    if(!is.null(coastal_points) & !(BAD_OBS_TREATMENT == 'None')){
        # Some data is clearly in an erronious position. There are 2 possible "Fixes"
        #    - Remove them by dropping points if the "distance" between the coast and the point is too big
        #    - Project their base point onto the nearest coastal point
        stopifnot(is(coastal_points, 'matrix') & (ncol(coastal_points) == 2))

        # Apply the fix
        tokeep = rep(TRUE, length(lon))
        for(i in 1:length(lon)){
            if(any(is.na(lon[i]) | is.na(lat[i]) | is.na(len[i]))) next

            if(BAD_OBS_TREATMENT == 'Drop'){
                # Throw away observed lon/lat that are too far away
                DROP_DATA_DIST = 300 # meters
                # Estimate distance, without depending on a spherical geometry package
                earth_rad_times_deg2rad = 6371000 * pi/180
                local_dist = earth_rad_times_deg2rad * 
                    sqrt( min( ((lon[i] - coastal_points[,1])*cos(42/180*pi))**2 + 
                                (lat[i] - coastal_points[,2])**2) )

                if(local_dist >= DROP_DATA_DIST) tokeep[i] = FALSE
            }else if(BAD_OBS_TREATMENT == 'Project'){
                # Replace observed lon/lat with nearest coastal point in model
                local_ind = which.min(
                    ((lon[i] - coastal_points[,1])*cos(42/180*pi))**2 + 
                    (lat[i]  - coastal_points[,2]                )**2)
                lon[i] = coastal_points[local_ind,1]
                lat[i] = coastal_points[local_ind,2]
            }
        }
        
        lon = lon[tokeep]
        lat = lat[tokeep]
        len = len[tokeep]
    }

    uv = cbind(lon - centroid[1], lat-centroid[2])
    uv = uv/sqrt(rowSums(uv**2))
    arrows(lon, lat, lon + uv[,1]*len, lat+uv[,2]*len, length=0, ...)
}

# Scale the runup bars (observed in meters, size plotted in degrees)
ARROW_LENGTH_SCALE = 1/1000

if(BAD_OBS_TREATMENT == 'None'){
    # This is not actually related to dropping bad points, but the options are
    # consistent with older behaviour of the script.
    PLOT_MODEL_ON_LAND_ONLY = TRUE
}else{
    PLOT_MODEL_ON_LAND_ONLY = FALSE
}

png(paste0('runup_heights_okushiri_', image_flags, '.png'), width=6.1, height=8, units='in', res=200)

# Setup plot frame with tohoku university group data
plot(field$lon, field$lat, pch=19, cex=0.2, asp=1/cos(42/180*pi), 
    ylim=c(42.02, 42.25), xlim=c(139.35, 139.6), col='white',
    xlab='Lon', ylab='Lat', cex.axis=1.5, cex.lab=1.7)
if(BAD_OBS_TREATMENT %in% c('Drop', 'None')){
    title('Okushiri runup (plotted radially from center point) \n Beware some data positioning errors', 
          cex.main=1.4, line=1)
}else if(BAD_OBS_TREATMENT == 'Project'){
    title('Okushiri runup (plotted radially from center point) \n Data projected to coast', 
          cex.main=1.4, line=1)
}

# High-ish res models
max_stage_store = list()
store_wd_points = matrix(NA, ncol=2, nrow=0)
for(j in 2:6){
    x5 = get_domain_contour(j)

    max_stage_store[[j]] = x5$stage_elev_rasters

    if(j == 2) next

    #plot(x5$stage_elev_rasters[[2]], zlim=c(0, 30), add=TRUE)
    for(i in 1:length(x5$all_wet_dry_contours)){
        lon = x5$all_wet_dry_contours[[i]]$x
        lat = x5$all_wet_dry_contours[[i]]$y
        stg = x5$peak_stage[[i]]
        elv = x5$peak_elev[[i]]
        iswet = (stg > (elv))

        len=stg*ARROW_LENGTH_SCALE # arrow length

        if(PLOT_MODEL_ON_LAND_ONLY){
            # Plot wet areas, with elevation > 0
            kk = which(elv > 0 & iswet)
        }else{
            # Or plot all wet areas at the coast
            kk = which(iswet)
        }
        add_arrows(lon[kk], lat[kk], len[kk], centroid=centroid, coastal_points = NULL, col='orange', lwd=3)
        points(lon[kk], lat[kk], pch=19, cex=0.2, col='purple')
        store_wd_points = rbind(store_wd_points, cbind(lon, lat))
    }
    plot(extent(x5$stage_elev_rasters[[2]]), add=TRUE, col='blue')
}

# Tohoku group
add_arrows(field$lon, field$lat, field$z_corrected*ARROW_LENGTH_SCALE, centroid=centroid, 
           coastal_points=store_wd_points)
#points(field$lon, field$lat, pch=19, cex=0.2)

# ujnr dataset
#points(field_ujnr$lon, field_ujnr$lat, cex=0.2, pch=19, col='green')
add_arrows(field_ujnr$lon, field_ujnr$lat, field_ujnr$z*ARROW_LENGTH_SCALE, centroid=centroid, 
           coastal_points=store_wd_points, col='green')

# Tsuji dataset (incomplete)
add_arrows(field_tsuji$lon, field_tsuji$lat, field_tsuji$z_max*ARROW_LENGTH_SCALE, centroid=centroid, 
           coastal_points=store_wd_points, col='red')

# Dense grid
abline(v=seq(139.3, 139.7,by=0.01), col='grey', lty='dotted')
abline(h=seq(42, 42.3, by=0.01), col='grey', lty='dotted')
points(centroid[1], centroid[2], pch=19)

legend('bottomright', c('Model', 'Data (Tohoku)', 'Data (UJNR)', 'Data (Tsuji)'), lty=rep(1, 4), cex=1.2,
       lwd=c(3,1,1,1), col=c('orange', 'black', 'green', 'red'), bty='n')

dev.off()

#
# Make similar max-stage plot
#
png(paste0('max_stage_okushiri_', image_flags, '.png'), width=6.1, height=8, units='in', res=800)
plot(field$lon, field$lat, pch=19, cex=0.2, asp=1/cos(42/180*pi), 
    ylim=c(42.02, 42.25), xlim=c(139.35, 139.6), col='white',
    xlab='Lon', ylab='Lat', cex.lab=1.7, cex.axis=1.5)
title('Okushiri Island modelled max-stage', cex.main=1.7)
for(j in 2:6){
    wet_stage = max_stage_store[[j]][[1]]
    elev = max_stage_store[[j]][[2]]
    wet_stage[wet_stage <= elev + eps] = NA
    plot(wet_stage, add=TRUE, zlim = c(0.1,31), 
         col=colorRampPalette(c('grey', 'lightblue', 'green', 'yellow', 'orange', 'red', 'purple'), bias=3)(255))
    plot(extent(wet_stage), add=TRUE, col='blue')
    contour(elev, levels=0, labels="", add=TRUE)

    #
    # Check max stage is reasonable
    # This will depend on the model resolution/setup
    #
    if(j == 5){
        # Monai domain -- high runup region
        desired_max_wet_stage = 31.7
        model_max_wet_stage = max(as.matrix(wet_stage), na.rm=TRUE)
        err_stat = abs(model_max_wet_stage - desired_max_wet_stage)/desired_max_wet_stage
        if((err_stat < 0.2) | (model_max_wet_stage > 25.0)){
            print(c('PASS', model_max_wet_stage))
        }else{
            print(c('FAIL', model_max_wet_stage))
        }

    }
}
dev.off()


#
# Check mass conservation in the domain
#
multidomain_log = readLines(Sys.glob(paste0(MD_DIR,'/multidomain*.log'))[1])
k = grep('unexplained ', multidomain_log)
final_mass_err = multidomain_log[k[length(k)]]
final_mass_err_val = as.numeric(strsplit(final_mass_err, ':')[[1]][2])
if(abs(final_mass_err_val) < 5){
    print(c('PASS', final_mass_err_val))
}else{
    print(c('FAIL', final_mass_err_val))
}

#
# Add an image of the elevation and domains
#
tmp = get_domain_interior_bbox_in_multidomain(MD_DIR)
COLPOW = 2

png(paste0('elevation_okushiri_', image_flags, '.png'), width=6.0, height=8, units='in', res=200)
multidomain_image(multidomain_dir=MD_DIR, variable='elevation0', 
    xlim=tmp$multidomain_xlim, ylim=tmp$multidomain_ylim, zlim=4000*c(-1,1), asp=1/cos(42/180*pi), 
    cols=c(rev(colorRampPalette(rev(c('darkblue', 'blue', 'steelblue', 'skyblue', 'lightblue')), bias=COLPOW)(255)), 
           colorRampPalette(c('lightgreen', 'green', 'darkgreen', 'orange', 'brown'), bias=COLPOW)(255) ),
    use_fields=TRUE)
mtext('Longitude', side=1, cex=1.5, line=2)
mtext('Latitude', side=2, cex=1.5, line=2)
title(main='Model elevation and domain nesting', cex.main=1.7)
# Add bounding boxes (without merging)
for(i in 1:length(tmp$domain_interior_bbox)){
  bb = tmp$domain_interior_bbox[[i]]
  points(rbind(bb, bb[1,]), t='l', col='red')
}

if(PLOT_EXTRAS){
    # Optionally add a contour with the initial deformation
    suppressMessages(library(raster))
    initial_condition = raster("../test_repository/BP09-FrankG-Okushiri_island/initial_condition_raster/HNO1993.tif")

    contour(initial_condition, level=setdiff(seq(-1, 5, by=0.5), 0), col='white', add=TRUE)    
    #contour(initial_condition, level=seq(-1, -2, -3, -4, -5), col='black', add=TRUE)    
}

dev.off()


if(PLOT_EXTRAS){
    # Zoom of the above plot

    png(paste0('elevation_okushiri_zoom_', image_flags, '.png'), width=6.1, height=8, units='in', res=200)
    multidomain_image(multidomain_dir=MD_DIR, variable='elevation0', asp=1/cos(42/180*pi),
        xlim=c(139.2, 139.8), ylim=c(41.7, 42.3), zlim=4000*c(-1,1), 
        cols=c(rev(colorRampPalette(rev(c('darkblue', 'blue', 'steelblue', 'skyblue', 'lightblue')), bias=COLPOW)(255)), 
               colorRampPalette(c('lightgreen', 'green', 'darkgreen', 'orange', 'brown'), bias=COLPOW)(255) ),
        use_fields=TRUE)
    mtext('Longitude', side=1, cex=1.5, line=2)
    mtext('Latitude', side=2, cex=1.5, line=2)
    title(main='Zoom around Okushiri Island', cex.main=2)
    # Add bounding boxes (without merging)
    for(i in 1:length(tmp$domain_interior_bbox)){
      bb = tmp$domain_interior_bbox[[i]]
      points(rbind(bb, bb[1,]), t='l', col='red')
    }

    suppressMessages(library(raster))
    initial_condition = raster("../test_repository/BP09-FrankG-Okushiri_island/initial_condition_raster/HNO1993.tif")

    contour(initial_condition, level=setdiff(seq(-1, 5, by=0.5), 0), col='white', add=TRUE)    

    dev.off()


}
