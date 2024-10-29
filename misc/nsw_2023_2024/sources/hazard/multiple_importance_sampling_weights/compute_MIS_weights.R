#
# Investigations of the "optimal weights" for multiple importance samples.
#
# We have N separate Monte Carlo samples, using different importance functions, which
# are in principle optimised differently for different sites. We can estimate the hazard
# by combining the samples IS1, IS2, IS3 ... with 
#    multiple_importance_sampling_exrate = sum_j{ w_j * exrate(ISj) } 
# where exrate(ISj) is an exceedance-rate of interest for the j'th importance sample
# for a scenario frequency model of interest, and
#    sum_j { w_j } = 1
#
# If the weights are determined without "peeking" at the samples, then this will be an unbiased
# estimate of the exceedance-rate of interest. The weights can vary spatially. 
# 
# Here we determine the weights that minimise the analytical variance of the multiple importance
# sampling solution at PTHA18 points, and extend them to rasters. 
#
# Note the calculations do not involve any 'peeking' at the actual samples.
# The weights are purely based on the analytical Monte Carlo variances for each
# importance sample at PTHA18 points, using PTHA18 wave heights. Conceptually it's
# analogous to the magnitude-bin optimization used in the stratified-importance-sampling.
#

# Get main function and dependencies (ptha18, isu, ...)
source('compute_MIS_weights_per_source_zone.R')

#
# INPUTS
#

# Find the importance samples (there are multiple)
all_sample_dirs = Sys.glob('../scenarios_ID*')
all_samples = basename(all_sample_dirs) # Used for file naming

# PTHA18 hazard point gauges at which we compute weights.
#model_point_gauges_file = 
all_ptha18_gauges = ptha18$get_all_gauges()
k = which(all_ptha18_gauges$lon > 140 & all_ptha18_gauges$lon < 180 &
          all_ptha18_gauges$lat > -57 & all_ptha18_gauges$lat < -22)
#model_ptha18_point_gauges = list(
    #'nsw' = read.csv('../../../point_gauges/point_gauges_2023_11_14.csv'),
    #'lordhowe' =     
model_ptha18_point_gauges = all_ptha18_gauges[k,]
    
# We will do calculations at stages matching these source-zone-specific
# exceedance-rates.
stage_exceedance_rate_thresholds = c(1/100, 1/500, 1/2500, 1/10000)

use_cache = TRUE # Read main calculations from pre-existing file if possible.
cache_file = 'compute_MIS_weights_per_source_zone_result.RDS'
fig_outdir = './FIG/' # Include trailing slash

# To make PNG's nice (focussed on NSW coast here)
plot_x_range = c(147.5, 155.5)
plot_y_range = c(-40, -26)

# Control pixels in output weights raster. Better to not have this too
# high-resolution (smooth variation in weights is good). Use interpolation
# when extracting raster values to avoid jumps at pixel boundaries.
weight_raster_xs = seq(140, 180, by=1)
weight_raster_ys = seq(-55, -25, by=1)
weight_raster_outdir = './importance_sample_weight_rasters/' # Include trailing slash
exceedance_rate_for_variance_reduction_plot = stage_exceedance_rate_thresholds[4]
    # In practice we'll probably use a single weight for each source zone.
    # The code will make a plot showing the expected variance reduction in this case (compared to equal weights).
    # This tag is used to find the correct raster to use for those plots.

#
# END INPUTS
#

# Read the data
if(!use_cache | !file.exists(cache_file)){
    mis_env = compute_MIS_weights_per_source_zone(ptha18, isu, 
        all_sample_dirs, all_samples, 
        model_ptha18_point_gauges, 
        stage_exceedance_rate_thresholds)

    if(use_cache) saveRDS(mis_env, cache_file)
}else{
    # Read from cache 
    mis_env = readRDS(cache_file)
}

# Unpack data
source_zones = mis_env$source_zones
szd = mis_env$szd

#
# Make lots of figures
#

# Restrict the range of points used for plotting, so we 
inside = (model_ptha18_point_gauges$lon > plot_x_range[1] & 
          model_ptha18_point_gauges$lon < plot_x_range[2] &   
          model_ptha18_point_gauges$lat > plot_y_range[1] &   
          model_ptha18_point_gauges$lat < plot_y_range[2] )
keepPlot = which(inside)
# An even more limited range for nearshore points
n0 = which(inside & 
    (model_ptha18_point_gauges$gaugeID - floor(model_ptha18_point_gauges$gaugeID) < 0.12) # Nearshore sites, detected by the non-integer part of their gauge ID
    )
# A limited range for offshore points (not nearshore)
not_n0 = which(inside &
    (model_ptha18_point_gauges$gaugeID - floor(model_ptha18_point_gauges$gaugeID) > 0.12) # Offshore sites, detected by the non-integer part of their gauge ID
    )
nexrate = length(stage_exceedance_rate_thresholds)
source('spatial_barplots.R')
dir.create(fig_outdir, showWarnings=FALSE)
for(sz in 1:length(source_zones)){

    source_zone = source_zones[sz]

    #
    # Spatial barplots with the optimal fraction for each site 
    #
    png(paste0(fig_outdir, '/Optimal_fraction_spatial_barplot_', source_zone, '.png'), width=10, height=4, units='in', res=300)
    barplot_col = c('red', 'green', 'blue')
    par(mfrow=c(1,4))
    par(mar=c(2,2,1,0.5))
    for(i in 1:nexrate){
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,i,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling weights
        #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[ ,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        spatial_barplots(model_ptha18_point_gauges$lon[keepPlot], model_ptha18_point_gauges$lat[keepPlot], bar_values=szd[[source_zone]]$best_weights[keepPlot,i,], 
            bar_width_scale=0.1, bar_height_scale=0.2, bar_col=barplot_col, asp=1, xlab='Lon', ylab='Lat')
        legend('topleft', all_samples, fill=barplot_col, title='Optimal weight', bty='n', cex=1.1)
        title(paste0(source_zone, ' Rate = ', stage_exceedance_rate_thresholds[i]), cex.main=1.4)
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

    #
    # Scatter plot of weights by latitude, zoom into the NSW coast
    #
    png(paste0(fig_outdir, '/Optimal_fraction_by_latitude_combined_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    par(mfrow=c(2,2))
    for(i in 1:nexrate){
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,i,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling best_weights
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[ ,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        pt_cex = szd[[source_zone]]$IS_effective_population_size[,i,]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

        plot(szd[[source_zone]]$best_weights[keepPlot,i,2], model_ptha18_point_gauges$lat[keepPlot], col=barplot_col[2], pch=19, xlim=c(0, 1), cex = pt_cex[keepPlot,2], 
            xlab='Importance sample weight', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
            main=paste0('Optimal sample weight by latitude, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
        points(szd[[source_zone]]$best_weights[keepPlot,i,3], model_ptha18_point_gauges$lat[keepPlot], col=barplot_col[3], pch=1, cex=pt_cex[keepPlot,3])
        points(szd[[source_zone]]$best_weights[keepPlot,i,1], model_ptha18_point_gauges$lat[keepPlot], col=barplot_col[1], pch=2, cex=pt_cex[keepPlot,1])
        legend('topright', all_samples, col=barplot_col, pch=c(2, 19, 1), title='Optimal weight', cex=1.4, bty='n')
        grid(col='orange')
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

    #
    # Scatter plot of weights by latitude, nearshore sites (inner most PTHA18 hazard point contours)
    #
    png(paste0(fig_outdir, '/Optimal_fraction_by_latitude_nearshore', '_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    par(mfrow=c(2,2))
    for(i in 1:nexrate){
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[n0,i,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling weights
        #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[n0,i,], na.rm=TRUE)
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[n0,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        pt_cex = szd[[source_zone]]$IS_effective_population_size[,i, ]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

        plot(szd[[source_zone]]$best_weights[n0,i,2], model_ptha18_point_gauges$lat[n0], col=barplot_col[2], pch=19, xlim=c(0, 1),
            cex = pt_cex[n0,2], 
            xlab='Importance sample weight', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
            main=paste0('Optimal sample weight by latitude nearshore, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
        points(szd[[source_zone]]$best_weights[n0,i,3], model_ptha18_point_gauges$lat[n0], col=barplot_col[3], pch=1, cex=pt_cex[n0,3])
        points(szd[[source_zone]]$best_weights[n0,i,1], model_ptha18_point_gauges$lat[n0], col=barplot_col[1], pch=2, cex=pt_cex[n0,1])
        legend('topright', all_samples, col=barplot_col, pch=c(2, 19, 1), title='Optimal weight', cex=1.4, bty='n')
        grid(col='orange')
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

    #
    # Scatter plot of weights by latitude, offshore sites (excluding nearshore sites as defined above)
    #
    png(paste0(fig_outdir, '/Optimal_fraction_by_latitude_offshore', '_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    par(mfrow=c(2,2))
    for(i in 1:nexrate){
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[not_n0,i,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling weights
        #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[not_n0,i,], na.rm=TRUE)
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[not_n0,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        pt_cex = szd[[source_zone]]$IS_effective_population_size[,i, ]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

        plot(szd[[source_zone]]$best_weights[not_n0,i,2], model_ptha18_point_gauges$lat[not_n0], col=barplot_col[2], pch=19, xlim=c(0, 1),
            cex=pt_cex[not_n0,2],
            xlab='Importance sample weight', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
            main=paste0('Optimal sample weight by latitude offshore, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
        points(szd[[source_zone]]$best_weights[not_n0,i,3], model_ptha18_point_gauges$lat[not_n0], col=barplot_col[3], pch=1, cex=pt_cex[not_n0,3])
        points(szd[[source_zone]]$best_weights[not_n0,i,1], model_ptha18_point_gauges$lat[not_n0], col=barplot_col[1], pch=2, cex=pt_cex[not_n0,1])
        legend('topright', all_samples, col=barplot_col, pch=c(2, 19, 1), title='Optimal weight', cex=1.4, bty='n')
        grid(col='orange')
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

    # How much does the exact variance reduce using the optimal weighting
    png(paste0(fig_outdir, '/Variance_reduction_by_latitude_rate_1in', '_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    par(mfrow=c(2,2))
    for(i in 1:4){
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,i,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling weights
        #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        pt_cex = szd[[source_zone]]$IS_effective_population_size[,i, ]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

        plot(szd[[source_zone]]$variance_reduction[keepPlot,i], model_ptha18_point_gauges$lat[keepPlot], col='black', pch=19, xlim=c(0, 1),
            cex=rowMeans(pt_cex[keepPlot,]),
            xlab='Variance reduction', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
            main=paste0('Optimal-variance / Equal-weight-variance, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
        grid(col='orange')
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

    #png(paste0(fig_outdir, '/Variance_reduction_with_modelled_wts_by_latitude_rate_1in', '_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    #par(mfrow=c(2,2))
    #for(i in 1:4){
    #    # Effective population (not sample) size in PTHA18
    #    efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,i,], na.rm=TRUE)
    #    # Effective population (not sample) size using the importance sampling weights
    #    #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)
    #    efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

    #    pt_cex = szd[[source_zone]]$IS_effective_population_size[,i, ]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

    #    plot(szd[[source_zone]]$modelled_weight_variance_reduction[,i], model_ptha18_point_gauges$lat, col='black', pch=19, xlim=c(0, 2),
    #        cex=rowMeans(pt_cex),
    #        xlab='Variance reduction', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
    #        main=paste0('Modelled-wt-variance / Equal-weight-variance, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
    #    grid(col='orange')
    #    title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    #}
    #dev.off()

    #png(paste0(fig_outdir, '/Variance_reduction_with_modelled_wts_1in100_by_latitude_rate_1in', '_', source_zone, '.png'), width=16, height=16, units='in', res=300)
    #par(mfrow=c(2,2))
    #for(i in 1:4){
    #    # Effective population (not sample) size in PTHA18
    #    efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,i,], na.rm=TRUE)
    #    # Effective population (not sample) size using the importance sampling weights
    #    #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)
    #    efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[,i,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

    #    pt_cex = szd[[source_zone]]$IS_effective_population_size[,i, ]/max(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)

    #    plot(szd[[source_zone]]$modelled_weight_100_variance_reduction[,i], model_ptha18_point_gauges$lat, col='black', pch=19, xlim=c(0, 2),
    #        cex=pt_cex,
    #        xlab='Variance reduction', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
    #        main=paste0('Modelled-wt-variance (rate=1/100) / Equal-weight-variance, rate= 1/', round(1/stage_exceedance_rate_thresholds[i]), '\n ', source_zone), cex.main=1.4)
    #    grid(col='orange')
    #    title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    #}
    #dev.off()

}

#
# Make rasters containing the weights
#

dir.create(weight_raster_outdir, showWarnings=FALSE)
source('spatial_weights_IDW.R')
weight_rasters = vector(mode='list', length=length(source_zones))
names(weight_rasters) = source_zones
for(sz in 1:length(source_zones)){

    source_zone = source_zones[sz]
    print(source_zone)

    weight_rasters[[source_zone]] = vector(mode='list', length=length(stage_exceedance_rate_thresholds))
    names(weight_rasters[[source_zone]]) = paste0('r_1in', round(1/stage_exceedance_rate_thresholds))

    # For each exceedance-rate, compute the weight rasters for each importance sample on the source zone.
    # The weight rasters sum to 1 everywhere (reflecting that our sample weights must sum to 1 everywhere).
    for(exr in 1:length(stage_exceedance_rate_thresholds)){

        # interp_weights_on_raster<-function(lonlat, weights, sample_names, raster_xs, raster_ys, radius=1){
        r_list =  interp_weights_on_raster(
            lonlat = as.matrix(model_ptha18_point_gauges[,1:2]),
            weights = szd[[source_zone]]$best_weights[,exr,],
            sample_names = all_samples,
            raster_xs = weight_raster_xs,
            raster_ys = weight_raster_ys,
            radius=2)

        # Store it
        weight_rasters[[source_zone]][[exr]] = r_list

        exrate_tag = names(weight_rasters[[source_zone]])[exr] # Useful name

        # Make a figure
        png(paste0(fig_outdir, '/Weight_rasters_', source_zone, '_', exrate_tag, '.png'),
            width=9, height=3, units='in', res=300)
        par(mfrow=c(1,length(r_list)))
        for(i in 1:length(r_list)){
            plot(r_list[[i]], range=c(0,1))
            points(model_ptha18_point_gauges[,1:2], pch='.')
            title(paste0(source_zone, '-', exrate_tag, '\n', 'weight of ', all_samples[i]))

        }
        dev.off()

        # Save the rasters as tifs
        for(i in 1:length(r_list)){
            output_file = paste0(weight_raster_outdir, 'weight_raster_', source_zone, '_', exrate_tag, '_', all_samples[i], '.tif')
            writeRaster(r_list[[i]], filename=output_file, overwrite=TRUE, gdal=c('COMPRESS=DEFLATE'))
        }

    }

    #
    # Make a variance reduction plot that will be more indicative of the real improvements
    #
    exr_ref = which(stage_exceedance_rate_thresholds == exceedance_rate_for_variance_reduction_plot)
    if(length(exr_ref) == 0){
        stop('exceedance_rate_for_variance_reduction_plot is not in stage_exceedance_rate_thresholds -- make sure they are equal at the floating point level')
    }
    wts = matrix(NA, nrow=nrow(model_ptha18_point_gauges), ncol=length(all_samples))
    # Compute the weights for each importance sample at model PTHA18 points
    for(i in 1:ncol(wts)){
        # Get weights from the raster.
        wts[,i] = extract(weight_rasters[[source_zone]][[exr_ref]][[i]], model_ptha18_point_gauges[,1:2],
                          method='bilinear', ID=FALSE, raw=TRUE)
        # Treat points outside the raster -- default to equal weights
        k = which(is.na(wts[,i]))
        if(length(k) > 0) wts[k,i] = 1/ncol(wts)
    }

    # Plot the variance reduction
    png(paste0(fig_outdir, '/Variance_reduction_with_unique_weights_per_source_by_latitude_',
                source_zone, '.png'),
                width=16, height=16, units='in', res=300)
    par(mfrow=c(2,2))
    for(exr in 1:length(stage_exceedance_rate_thresholds)){
        variance_reduction_with_selected_weight_raster =
            rowSums(wts**2 * szd[[source_zone]]$var[,exr,])/rowSums( (1/ncol(wts))**2 * szd[[source_zone]]$var[,exr,])
        # Effective population (not sample) size in PTHA18
        efps_ptha18 = median(szd[[source_zone]]$ptha18_effective_population_size[,exr,], na.rm=TRUE)
        # Effective population (not sample) size using the importance sampling weights
        #efps_IS = median(szd[[source_zone]]$IS_effective_population_size[,i,], na.rm=TRUE)
        efps_IS = apply(szd[[source_zone]]$IS_effective_population_size[,exr,], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE))

        pt_cex = szd[[source_zone]]$IS_effective_population_size[,exr, ]/max(szd[[source_zone]]$IS_effective_population_size[,exr,], na.rm=TRUE)

        plot(variance_reduction_with_selected_weight_raster[keepPlot], model_ptha18_point_gauges$lat[keepPlot],
            col='black', pch=19, xlim=c(0, 2), cex=pt_cex[keepPlot],
            xlab='Variance reduction', ylab='Latitude', cex.lab=1.4, cex.axis=1.4, 
            main=paste0('Modelled-wt-variance (rate=1/', round(1/exceedance_rate_for_variance_reduction_plot),
                        ') / Equal-weight-variance, rate= 1/', round(1/stage_exceedance_rate_thresholds[exr]),
                        '\n ', source_zone), 
            cex.main=1.4)
        grid(col='orange')
        title(sub=paste0('PTHA18 effective Pop N (median) = ', round(efps_ptha18), 
                         ', IS effective Pop N (median) = ', paste0(round(efps_IS), collapse=" ")))
    }
    dev.off()

}

# Check the variance reduction that we'll get if using the rasters for weights
