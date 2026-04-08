library(terra)

source('sites_at_sea_levels.R')

output_basedir = 'rasters_at_sites'
sites = dir(output_basedir)

# Rasters (varying elevation models) by site
sites_varying_elevation_rasts = lapply(sites, function(x) Sys.glob(paste0(output_basedir, '/', x, '/varying_elevation/*.tif')))
names(sites_varying_elevation_rasts) = sites

# Rasters (static sea level models) by site
sites_static_elevation_rasts = lapply(sites, function(x) Sys.glob(paste0(output_basedir, '/', x, '/static_sea_level_*/*.tif')))
names(sites_static_elevation_rasts) = sites

# We expect the depth to be very similar in each case. We expect the elevation
# and max-stage to differ by (almost) the chosen sea level, although the
# max-stage will also have effects due to the dynamics
events = unique(sapply(basename(sites_varying_elevation_rasts[[1]]), function(x) strsplit(x, split="-ambient_sea_level")[[1]][1]))

plot_site_and_event<-function(site, event){

    stopifnot(length(event)==1 & length(site)==1)

    # Get files matching site and event, varying elevation
    ii = grep(event, sites_varying_elevation_rasts[[site]])
    candidate_rasts_varying_elevation = sites_varying_elevation_rasts[[site]][ii]

    # Get files matching site and event, static elevation
    ii = grep(event, sites_static_elevation_rasts[[site]])
    candidate_rasts_static_elevation = sites_static_elevation_rasts[[site]][ii]

    # Get the desired extent of plots
    plot_bbox = sites_at_sea_levels[[site]]$plot_bbox
    plot_ext = ext(c(plot_bbox[1], plot_bbox[3], plot_bbox[2], plot_bbox[4]))

    # Read the rasters, cropped to the desired output extent
    tmplte = list(depth=NA, max_stage=NA, elevation0=NA)
    rasts = list(varying=tmplte, static=tmplte)   
    for(tag in c('depth', 'max_stage', 'elevation0')){
        filename = candidate_rasts_varying_elevation[grep(paste0('_', tag, '.tif'), candidate_rasts_varying_elevation)]
        rasts$varying[[tag]] = crop(rast(filename), plot_ext)
        filename = candidate_rasts_static_elevation[grep(paste0('_', tag, '.tif'), candidate_rasts_static_elevation)]
        rasts$static[[tag]] = crop(rast(filename), plot_ext)
    }

    # Plots
    par(mfrow=c(2,2))
    par(mar=c(4,4,5,2))
    # Difference in max depth
    depth_diff = rasts$varying$depth - rasts$static$depth
    rng_depth = range(as.matrix(depth_diff), na.rm=TRUE)
    percentiles_depth_diff = quantile(na.omit(as.matrix(depth_diff)), probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
    plot(depth_diff, main=paste0('Difference in depths (range ', round(rng_depth[1], 3), ' ', round(rng_depth[2], 3), ')'))
    title(main=event, line=1, cex.main=0.5)

    # Max stage (varying) plus a wet/dry contour derived from the static model
    plot(rasts$varying$max_stage, main='Max stage (varying) with dry limit (static)', col=rainbow(255)[1:200])
    contour(is.na(rasts$static$max_stage), levels=0.5, add=TRUE, labels="")
    title(main=event, line=1, cex.main=0.5)

    # Difference in elevation
    elev_diff = rasts$varying$elevation0 - rasts$static$elevation0
    rng_elev = range(as.matrix(elev_diff), na.rm=TRUE)
    plot(elev_diff, main=paste0('Difference in elevation (range ', round(rng_elev[1], 3), ' ', round(rng_elev[2], 3), ')'))
    title(sub=event)

    # Histogram of differences in depth
    hist(na.omit(as.matrix(depth_diff)), n=101); abline(v=0, col='red')

    max_stage_diff = rasts$varying$max_stage - rasts$static$max_stage - mean(elev_diff, na.rm=TRUE)
    rng_max_stage = range(as.matrix(max_stage_diff), na.rm=TRUE)
    #plot(max_stage_diff, main='Difference in max-stage (minus mean(elev_diff))')

    return(invisible(list(site=site, event=event, rng_depth=rng_depth, rng_elev=rng_elev, rng_max_stage=rng_max_stage, percentiles_depth=percentiles_depth_diff)))
}

# Store output stats
site_results = vector(mode='list', length=length(sites))
names(site_results) = sites
for(site in sites){
    print(paste0('Plotting ', site))

    # Store output stats
    site_results[[site]] = vector(mode='list', length=length(events))
    names(site_results[[site]]) = events

    output_pdf = paste0(output_basedir, '/', site, '/Compare_maxima_', site, '.pdf')
    pdf(output_pdf, width=10, height=10)
    for(event in events){
        site_results[[site]][[event]] = plot_site_and_event(site, event)
    }
    dev.off()
}

saveRDS(site_results, 'site_results.RDS')

#
# Plot the difference in depth at different sites
#
pc_depth_diff = lapply(site_results, function(x) do.call(rbind, lapply(x, function(y) y$percentiles_depth)))
pdf('Differences_in_max_depth_by_site.pdf', width=9, height=6)
for(site in names(site_results)){
    ylimit = max(abs((pc_depth_diff[[site]])))*c(-1,1)
    par(mar=c(5,6,4,3))
    plot(pc_depth_diff[[site]][,1], t='o', xlab='Scenario index', ylab='Depth difference (m) \n (varying elevation - static)', 
        main=paste0(site, ': Differences in depth with "varying elevation" vs "static" tides.'), ylim=ylimit,
        cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
    for(i in 2:5){
        points(pc_depth_diff[[site]][,i], t='o', col=i)
    }
    legend('top', colnames(pc_depth_diff[[site]]), col=1:5, lty='solid', pch=1, horiz=TRUE, bty='n')
    grid(col='orange')
}
dev.off()


