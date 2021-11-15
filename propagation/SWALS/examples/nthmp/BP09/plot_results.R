#
# Plotting of BP09 results
#
image_flags = commandArgs(trailingOnly=TRUE)
if(length(image_flags) != 1){
    stop('Commandline arguments not recognized')
}

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

get_domain_contour<-function(domain_number){

    #most_recent_run = rev(Sys.glob(paste0('OUTPUTS/RUN_*/RUN_*000', domain_number, '*')))[1]
    most_recent_run = rev(Sys.glob(paste0('OUTPUTS/RUN_*')))[1]

    #x5 = get_all_recent_results(most_recent_run)
    #stage_elev_rasters = make_max_stage_raster(x5, return_elevation=TRUE, na_outside_priority_domain=TRUE)
    #x5 = list()

    stage_elev_rasters = vector(mode='list', length=2)
    stage_elev_rasters[[1]] = merge_domains_nc_grids(multidomain_dir = most_recent_run, 
        domain_index=domain_number, desired_var = 'max_stage', return_raster=TRUE)
    stage_elev_rasters[[2]] = merge_domains_nc_grids(multidomain_dir = most_recent_run, 
        domain_index=domain_number, desired_var = 'elevation0', return_raster=TRUE)
    xs = xFromCol(stage_elev_rasters[[1]], 1:ncol(stage_elev_rasters[[1]])) 
    ys = yFromRow(stage_elev_rasters[[1]], nrow(stage_elev_rasters[[1]]):1) 

    # Get 'peak' elevation, which accounts for the fact that nesting + linear extrapolation
    # can change elevation in nesting regions. If not accounted for, we might mistakenly think
    # that nesting regions have wet peak stage in some areas
    #peak_elev = x5$elev[,,1]
    #for(i in 2:length(x5$elev[1,1,])){
    #    peak_elev = pmax(peak_elev, x5$elev[,,i])
    #}
    #setValues(stage_elev_rasters[[2]], c(peak_elev))

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
add_arrows<-function(lon, lat, len, centroid, ...){
    uv = cbind(lon - centroid[1], lat-centroid[2])
    uv = uv/sqrt(rowSums(uv**2))
    arrows(lon, lat, lon + uv[,1]*len, lat+uv[,2]*len, length=0, ...)
}

png(paste0('runup_heights_okushiri_', image_flags, '.png'), width=8, height=8, units='in', res=200)

# Setup plot frame with tohoku university group data
plot(field$lon, field$lat, asp=1, pch=19, cex=0.2,
    ylim=c(42.02, 42.25), xlim=c(139.35, 139.6),
    xlab='Lon', ylab='Lat', cex.axis=1.5, cex.lab=1.7)
title('Okushiri runup (plotted radially from center point) \n Beware some data positioning errors', 
      cex.main=1.7, line=1)

# High-ish res models
max_stage_store = list()
for(j in 3:6){
    x5 = get_domain_contour(j)
    #plot(x5$stage_elev_rasters[[2]], zlim=c(0, 30), add=TRUE)
    for(i in 1:length(x5$all_wet_dry_contours)){
        lon = x5$all_wet_dry_contours[[i]]$x
        lat = x5$all_wet_dry_contours[[i]]$y
        stg = x5$peak_stage[[i]]
        elv = x5$peak_elev[[i]]
        iswet = (stg > (elv))

        len=stg/1000 # arrow length
        # Plot wet areas, with elevation > 0 (to avoid getting most boundaries of interior nested grids)
        kk = which(elv > 0 & iswet)
        add_arrows(lon[kk], lat[kk], len[kk], centroid=centroid, col='orange', lwd=3)
        points(lon[kk], lat[kk], pch=19, cex=0.2, col='purple')
    }
    plot(extent(x5$stage_elev_rasters[[2]]), add=TRUE, col='blue')

    max_stage_store[[j]] = x5$stage_elev_rasters
}

# Tohoku group
add_arrows(field$lon, field$lat, field$z_corrected/(100*1000), centroid=centroid)
points(field$lon, field$lat, pch=19, cex=0.2)
# ujnr dataset
points(field_ujnr$lon, field_ujnr$lat, cex=0.2, pch=19, col='green')
add_arrows(field_ujnr$lon, field_ujnr$lat, field_ujnr$z/1000, centroid=centroid, col='green')
# Tsuji dataset (incomplete)
add_arrows(field_tsuji$lon, field_tsuji$lat, field_tsuji$z_max/1000, centroid=centroid, col='red')

# Dense grid
abline(v=seq(139.3, 139.7,by=0.01), col='grey', lty='dotted')
abline(h=seq(42, 42.3, by=0.01), col='grey', lty='dotted')
points(centroid[1], centroid[2], pch=19)

dev.off()

#
# Make similar max-stage plot
#
png(paste0('max_stage_okushiri_', image_flags, '.png'), width=7, height=8, units='in', res=800)
plot(field$lon, field$lat, asp=1, pch=19, cex=0.2,
    ylim=c(42.02, 42.25), xlim=c(139.35, 139.57), col='white',
    xlab='Lon', ylab='Lat', cex.lab=1.7, cex.axis=1.5)
title('Okushiri Island modelled max-stage', cex.main=2)
for(j in 3:6){
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
        # Monai domain
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
multidomain_log = readLines(rev(Sys.glob('OUTPUTS/RUN*/multidomain*.log'))[1])
k = grep('unexplained ', multidomain_log)
final_mass_err = multidomain_log[k[length(k)]]
final_mass_err_val = as.numeric(strsplit(final_mass_err, ':')[[1]][2])
if(abs(final_mass_err_val) < 5){
    print(c('PASS', final_mass_err_val))
}else{
    print(c('FAIL', final_mass_err_val))
}
