
#
# Inputs
#

# All the data is now available here
source('plot_parameters.R', local=TRUE)

#
# End inputs
#


# Read shapefiles
usp = readOGR(unit_source_polygon_shapefile, 
    layer=gsub('.shp', '', basename(unit_source_polygon_shapefile)))

zc = readOGR(zero_contour, gsub('.shp', '', basename(zero_contour)))

# Read tables with earthquake metadata
stoch_meta = list()
unif_meta = list()
vari_unif_meta = list()

stoch_meta$events_with_Mw = read_table_from_netcdf(stochastic_slip_events)
unif_meta$events_with_Mw = read_table_from_netcdf(uniform_slip_events)
vari_unif_meta$events_with_Mw = read_table_from_netcdf(variable_uniform_slip_events)

stoch_meta$unit_source_statistics = read_table_from_netcdf(unit_source_statistics)
unif_meta$unit_source_statistics = read_table_from_netcdf(unit_source_statistics)
vari_unif_meta$unit_source_statistics = read_table_from_netcdf(unit_source_statistics)

# Unit source statistics
uss = stoch_meta$unit_source_statistics

# Crop the DEM for fast plotting
if(is.null(plot_region)){
    plot_region = extent(usp) #extent(138, 145, 33, 44)
}
plot_poly = as(plot_region, 'SpatialPolygons')
# Extend the raster a bit so it goes to the plot boundaries
tmp = plot_region
tmp@xmin = tmp@xmin-10
tmp@ymin = tmp@ymin-10
tmp@xmax = tmp@xmax+10
tmp@ymax = tmp@ymax+10
smalldem = crop(dem, tmp)

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)


#
# unit sources "close enough" to the target event,
# where the definition of "close enough" varies according to the magnitude of the
# observed earthquake
#
nearby_uss = find_unit_sources_near_hypocentre(event_hypocentre, usp, uss, 
    event_mw, scaling_relation_type='Strasser')

uniform_slip_eis = strsplit(unif_meta$events_with_Mw$event_index_string, split='-')
keep_uniform = mapply(f<-function(mw, eis, rate){
    out = (any(eis %in% nearby_uss) & (abs(mw - event_mw) < magnitude_difference) & (rate>0))
    return(out)
    },
    unif_meta$events_with_Mw$Mw, uniform_slip_eis, unif_meta$events_with_Mw$rate_annual)
keep_stochastic = (stoch_meta$events_with_Mw$uniform_event_row %in% which(keep_uniform))
keep_vari_unif = (vari_unif_meta$events_with_Mw$uniform_event_row %in% which(keep_uniform))

#
# Stochastic slip plots
#
# @param i row_index in the event table to plot
# @param sim_meta  -- one of stoch_meta, unif_meta, vari_unif_meta
# @param base start of filename for each figure
# @param dem_zlim zlim for DEM in plot
# @param max_slip_val before plotting, clip all slip values to below this
# @param surface_def_max before plotting, clip all surface deformation values to below this
# @param surface_def_min before plotting, clip all surface deformation values to above this
# @param uniform_slip flag denoting uniform slip values
# @return 0 (but makes a plot as a side-effect)
single_plot<-function(i, sim_meta, base, dem_zlim, max_slip_val, surface_def_max, 
    surface_def_min, uniform_slip){ 

    # Skip 'zero rate' events
    if(sim_meta$events_with_Mw$rate_annual[i] == 0) return(0)
    
    png_name = paste0(out_dir, '/', base, '_', 100000 + i, '.png')
    print(png_name)

    # Get included unit-sources (so we know which ones to colour)
    included_uss = get_unit_source_indices_in_event(sim_meta$events_with_Mw[i,])
    included_uss_dd = uss$downdip_number[included_uss]
    included_uss_as = uss$alongstrike_number[included_uss]
    included_uss_usp = (usp$dwndp_n %in% included_uss_dd & 
        usp$alngst_ %in% included_uss_as)

    # Get slip
    if(!uniform_slip){ 
        included_uss_slip = as.numeric(strsplit(
            sim_meta$events_with_Mw$event_slip_string[i], '_')[[1]]) 
    }else{ 
        included_uss_slip = included_uss * 0 + sim_meta$events_with_Mw$slip[i]
    }

    # Make deformation raster
    local_raster_files = uss$initial_condition_file[included_uss]
    r1 = ptha$get_initial_condition_for_event(sim_meta, event_ID=i)

    # Map slip onto polygons
    included_slip_vals = rep(0, length(usp))
    for(j in 1:length(included_uss)){
        matching_ind = which((usp$dwndp_n %in% included_uss_dd[j] & 
            usp$alngst_ %in% included_uss_as[j]))
        included_slip_vals[matching_ind] = included_uss_slip[j]
    }

    # Make an informative title
    if(!uniform_slip){
        title_extra = paste0('Mw = ', round(sim_meta$events_with_Mw$Mw[i], 1), 
            '; Peak slip = ', 
            round(max(as.numeric(strsplit(sim_meta$events_with_Mw$event_slip_string[i], '_')[[1]])), 1))
    }else{
        title_extra = paste0('Mw = ', round(sim_meta$events_with_Mw$Mw[i], 1), '; Peak slip = ', 
            round(sim_meta$events_with_Mw$slip[i], 1))
    }

    png(png_name, width=8, height=16, res=75, units='in')

    par(mfrow=c(2,1))

    #
    # First panel of plot
    #
    plot(plot_region, col='white', asp=1)
    plot(usp, add=TRUE)
    image(smalldem, add=TRUE, col=grey(seq(0,1,len=100), alpha=0.3), 
        zlim=dem_zlim, maxpixels=1e+08)
    plot(zc, add=TRUE)

    ncol = 100
    colz = c('green', rev(heat.colors(ncol)))
    slip2col = ceiling(pmin(included_slip_vals/max_slip_val, 1)*(ncol)) + 1
    plot(usp, add=TRUE, col=colz[slip2col], 
        density=50)
    #title(paste0('Kurils-Japan stochastic slip events'), cex.main=2)
    title(title_extra, cex.main=2)
    # Add a legend to the plot
    plot(smalldem, add=TRUE, legend.only=TRUE, zlim=c(0, max_slip_val),
        col = colz)
    if(!all(is.null(wave_sites))){
        points(wave_sites, pch=19, cex=2)
        text(wave_sites[,1], wave_sites[,2], paste0('S', 1:2), pos=2, col='red', cex=1.5)
    }
        
    #
    # Second panel of plot
    #

    mycolz = rev(rainbow(255))[50:255]
    # Make values near zero transparent
    myind = approx(c(surface_def_min, surface_def_max), c(1, length(mycolz)), 
        xout=c(-surface_def_max/50, surface_def_max/50))$y
    mycolz[floor(min(myind)):ceiling(max(myind))] = rgb(1, 1, 1, alpha=0)

    plot(plot_region, col='white', asp=1)
    image(smalldem, add=TRUE, col=grey(seq(0,1,len=100), alpha=0.3), 
        zlim=dem_zlim, maxpixels=1e+08)
    plot(usp, add=TRUE)
    plot(max(min(r1, surface_def_max), surface_def_min), 
        zlim=c(surface_def_min, surface_def_max), 
        col=mycolz, add=TRUE)
    title(main='Ocean Surface Deformation (m)', cex.main=2)
    plot(usp, add=TRUE)
    plot(zc, add=TRUE)
    if(!all(is.null(wave_sites))){
        points(wave_sites, pch=19, cex=2)
        text(wave_sites[,1], wave_sites[,2], paste0('S', 1:2), pos=4, col='red', cex=1.5)
    }

    dev.off()

    return(0)
}

#
# Do lots of plots in parallel
#
# @param sim_meta  -- one of stoch_meta, unif_meta, vari_unif_meta
# @param base start of filename for each figure
# @param dem_zlim zlim for DEM in plot
# @param max_slip_val before plotting, clip all slip values to below this
# @param surface_def_max before plotting, clip all surface deformation values to below this
# @param surface_def_min before plotting, clip all surface deformation values to above this
# @param uniform_slip flag denoting uniform slip values
# @param desired_indices event indices to plot (row_indices in sim_meta$events_with_Mw)
# @return 0 (but makes a plot as a side-effect)
many_plots<-function(sim_meta, base='stochastic_slip',
    dem_zlim = c(-6500, 0), max_slip_val=50, 
    surface_def_max=12, surface_def_min=-8,
    uniform_slip=FALSE, desired_indices=NULL){

    if(is.null(desired_indices)){
        desired_indices = length(sim_meta$events_with_Mw[,1])
    }

    try_single_plot<-function(...){
        try(single_plot(...))
    }
   
    library(parallel)
    error_catch = mclapply(as.list(desired_indices), try_single_plot, 
        sim_meta=sim_meta, base=base,
        dem_zlim = dem_zlim, max_slip_val=max_slip_val, 
        surface_def_max=surface_def_max, surface_def_min=surface_def_min,
        uniform_slip=uniform_slip,
        mc.cores=MC_CORES)

    return(error_catch)
}

##
## Uncomment these to remake the plots
##
#many_plots(unif_meta, base='uniform_slip', uniform_slip=TRUE, desired_indices = which(keep_uniform))
errs = many_plots(stoch_meta, base='stochastic_slip', uniform_slip=FALSE, desired_indices = which(keep_stochastic),
    surface_def_min=SURFACE_DEF_MIN, surface_def_max=SURFACE_DEF_MAX, max_slip_val=MAX_SLIP_VAL,
    dem_zlim = DEM_ZLIM)
#many_plots(vari_unif_meta, base='variable_uniform_slip', uniform_slip=FALSE, desired_indices = which(keep_vari_unif))


# Stand-alone plot of earthquake slip
slip_panel_plot<-function(i, sim_meta, dem_zlim, max_slip_val, uniform_slip){

    # Get included unit-sources (so we know which ones to colour)
    included_uss = get_unit_source_indices_in_event(sim_meta$events_with_Mw[i,])
    included_uss_dd = uss$downdip_number[included_uss]
    included_uss_as = uss$alongstrike_number[included_uss]
    included_uss_usp = (usp$dwndp_n %in% included_uss_dd & 
        usp$alngst_ %in% included_uss_as)

    if(!uniform_slip){ 
        included_uss_slip = as.numeric(strsplit(
            sim_meta$events_with_Mw$event_slip_string[i], '_')[[1]]) 
    }else{ 
        included_uss_slip = included_uss * 0 + sim_meta$events_with_Mw$slip[i]
    }

    # Make deformation raster
    #local_raster_files = uss$initial_condition_file[included_uss]
    #r1 = ptha$get_initial_condition_for_event(sim_meta, event_ID=i)

    # Map slip onto polygons
    included_slip_vals = rep(0, length(usp))
    for(j in 1:length(included_uss)){
        matching_ind = which((usp$dwndp_n %in% included_uss_dd[j] & 
            usp$alngst_ %in% included_uss_as[j]))
        included_slip_vals[matching_ind] = included_uss_slip[j]
    }

    # Make the panel
    aspect_ratio = 1/cos(0.5*(plot_region@ymax + plot_region@ymin)/180 * pi)
    plot(usp, asp=aspect_ratio, axes=TRUE, cex.axis=1.5, border='grey')
    image(smalldem, add=TRUE, col=grey(seq(0,0.8,len=100), alpha=0.3), 
        zlim=dem_zlim, maxpixels=1e+08)
    plot(zc, add=TRUE)

    ncol = 100
    colz = c('green', rev(heat.colors(ncol)))
    slip2col = ceiling(pmin(included_slip_vals/max_slip_val, 1)*(ncol)) + 1
    plot(usp, add=TRUE, col=colz[slip2col], 
        density=50, border='white')
    for(i in 1:length(usp)){
        if(colz[slip2col][i] != 'green') plot(usp[i,], border='black', add=TRUE)
    }
    #title(paste0('Kurils-Japan stochastic slip events'), cex.main=2)
    # Add a legend to the plot
    plot(smalldem, add=TRUE, legend.only=TRUE,zlim=c(0, max_slip_val),
        col = colz, smallplot = c(0.8, 0.83, 0.3, 0.7), 
        legend.args=list(text='Slip (m)', side=2, line=0.5, cex=1.5))

}

##  #
##  #
##  #
##  png('slip_types_plot.png', width=12, height=9, units='in', res=300)
##  par(mfrow=c(2,3))
##  par(mar=c(3,3,1,2))
##  # FAUS
##  slip_panel_plot(3168, unif_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=TRUE)
##  title('FAUS', cex=2, line=-2, cex.main=2.5)
##  text(140, 56, 'a)', cex=3)
##  # VAUS
##  slip_panel_plot(47525, vari_unif_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=FALSE)
##  title('VAUS', line=-2, cex.main=2.5)
##  text(140, 56, 'c)', cex=3)
##  # HS
##  slip_panel_plot(47525, stoch_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=FALSE)
##  title('HS', line=-2, cex.main=2.5)
##  text(140, 56, 'e)', cex=3)
##  # FAUS
##  slip_panel_plot(3184, unif_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=TRUE)
##  title('FAUS', line=-2, cex.main=2.5)
##  text(140, 56, 'b)', cex=3)
##  # VAUS
##  slip_panel_plot(47741, vari_unif_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=FALSE)
##  title('VAUS', line=-2, cex.main=2.5)
##  text(140, 56, 'd)', cex=3)
##  # HS
##  slip_panel_plot(47634, stoch_meta, dem_zlim=c(-9500, 0), max_slip_val = 50, uniform_slip=FALSE)
##  title('HS', line=-2, cex.main=2.5)
##  text(140, 56, 'f)', cex=3)
##  dev.off()

mytime = as.character(julian(Sys.time(), units='secs'))
save.image(paste0('earthquake_types_plot_image_', mytime, '.RData'))
