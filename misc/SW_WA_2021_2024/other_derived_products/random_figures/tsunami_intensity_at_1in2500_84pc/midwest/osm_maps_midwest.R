library(rptha)
library(cptcity)
library(OpenStreetMap)
library(magick)
library(sf)
library(cartography)

# For running in multple batches
batch_number = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[1])
batch_count = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[2])

#
# INPUT PARAMETERS
#

# Regions for zoom plots. These were manually created
ZOOM_PLOTS_REGIONS = readOGR('plot_zoom_regions', verbose=FALSE) 
# Order from south to north to make identification easier
k = order(coordinates(ZOOM_PLOTS_REGIONS)[,2])
ZOOM_PLOTS_REGIONS = ZOOM_PLOTS_REGIONS[k,]

# Zones for maps
NO_THREAT_ZONE = readOGR('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/jatwc_inundation_zones/Geraldton-Coast/Geraldton-Coast_no_threat_with-PTHA-exrate-limit_84pc_4e-04', verbose=FALSE)
MARINE_WARNING_ZONE = readOGR('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/jatwc_inundation_zones/Geraldton-Coast/Geraldton-Coast_marine_warning_with-PTHA-exrate-limit_84pc_4e-04', verbose=FALSE)
LAND_WARNING_ZONE = readOGR('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/jatwc_inundation_zones/Geraldton-Coast/Geraldton-Coast_land_warning_with-PTHA-exrate-limit_84pc_4e-04', verbose=FALSE)

# Max pixels in image plot
MAXPIXELS=1e+07

# Do not plot tifs with resolution coarser than this. Can speed things up
# depending on multidomain design
MIN_TIF_RES = 1.1/(60*5)

##
## Various tifs derived from the model
##
TIF_FILES = list()
TIF_FILES$max_stage_1in2500_84pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/max_stage_1in2500_84pc/*.tif')
TIF_FILES$max_depth_1in2500_84pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/max_depth_1in2500_84pc/*.tif') 
TIF_FILES$max_flux_1in2500_84pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/max_flux_1in2500_84pc/*.tif') 
TIF_FILES$max_speed_1in2500_84pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/max_speed_1in2500_84pc/*.tif') 
TIF_FILES$elevation = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/elevation_in_model/*.tif')
TIF_FILES$inundation_rate_ltm = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/inundation_rate_logic_tree_mean/*.tif')
TIF_FILES$inundation_rate_84pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/inundation_rate_84pc/*.tif')
TIF_FILES$inundation_rate_16pc = Sys.glob('../../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/inundation_rate_16pc/*.tif')

# We use different kinds of plot for different variables. Group like this
FLOW_VARIABLE_NAMES = c(paste0('max_', c('stage', 'depth', 'flux', 'speed'), '_1in2500_84pc'), 'elevation')
INUNDATION_RATE_NAMES = names(TIF_FILES)[grep('inundation', names(TIF_FILES))]

# Store the tif extents
TIF_EXTENTS = vector(mode='list', length=length(TIF_FILES))
names(TIF_EXTENTS) = names(TIF_FILES)
for(nm in names(TIF_FILES)){
    TIF_EXTENTS[[nm]] = lapply(TIF_FILES[[nm]], function(tif_file) extent(raster(tif_file)))
}

# Store the tif resolution
TIF_RES = vector(mode='list', length=length(TIF_FILES))
names(TIF_RES) = names(TIF_FILES)
for(nm in names(TIF_FILES)){
    TIF_RES[[nm]] = unlist(lapply(TIF_FILES[[nm]], function(tif_file) res(raster(tif_file))[1]))
}

##
## Color palettes for flow variables
##
pal0 = cpt("ukmo_wow_humidity", n=255)[1:220]
pal1 = colorRampPalette(pal0, bias=1.8)(1024)
tmp = col2rgb(pal1)
pal1_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=160, maxColorValue=255)
# Another alternative
palDS = c('steelblue1', 'blue1', 'blue4')
palDS1 = colorRampPalette(palDS, bias=3)(1024)
tmp = col2rgb(palDS1)
palDS1_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=150, maxColorValue=255)
# Yet another alternative
palDS2 = cpt("arendal_temperature", n=400)
tmp = col2rgb(palDS2)
palDepthStage_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# Yet another alternative
palDS3 = cpt("cb_seq_YlOrRd_05", n=400)
tmp = col2rgb(palDS3)
palDS3_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# And another
palSpeed = cpt("h5_jet", n=400)
tmp = col2rgb(palSpeed)
palSpeed_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255) 
# And another
palFlux = colorRampPalette(cpt("oc_sst", n=400), bias=1.8)(400) # Concentrate variation around low values
tmp = col2rgb(palFlux)
palFlux_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255) 

#
# Classes to colour the inundation rates
#
inundation_classmat = matrix(
    # Lower, upper, integerID
    c(1/100  , 9999999, 1, 
      1/500  , 1/100  , 2, 
      1/2500 , 1/500  , 3, 
      1/10000, 1/2500 , 4, 
      -1     , 1/10000, 5), 
    ncol=3, byrow=TRUE)
inundation_classmat_labels = c('>1/100', '1/500 - 1/100', '1/2500 - 1/500', '1/10000 - 1/2500', '< 1/10000')
# Colours with transparency
inundation_classmat_col = c('red', 'orange', 'yellow', 'grey', NA)
icct = col2rgb(inundation_classmat_col, alpha=TRUE)
inundation_classmat_col = rgb(icct[1,], icct[2,], icct[3,], alpha=icct[4,]*0.5, maxColorValue=255)


# For each tif associated with the plot var, return TRUE if it intersects with the
# bounding box defined by lower_left and upper_right, and FALSE otherwise.
is_tif_touching_bbox<-function(lower_left, upper_right, plotvar='max_stage'){

    rast_extents = TIF_EXTENTS[[plotvar]]
    local_extent = extent(c(lower_left[1], upper_right[1], lower_left[2], upper_right[2]))
    dbsize = max(upper_right - lower_left)
    local_extent = extend(local_extent, dbsize) # To avoid missing things near the edges

    is_touching = lapply(rast_extents, function(x){
        intersect(x, local_extent)
    })
    is_touching = !unlist(lapply(is_touching, is.null))
    return(is_touching)
}

# Get the OSM data, with caching to speed it up
get_osm_backdrop_reproj<-function(lower_left, upper_right, dbox, backdrop){
    osm_cachefile = paste0('.osm_tile_cache_', upper_right[1], upper_right[2], lower_left[1], lower_left[2], '.RDS')
    dmax = max(dbox)
    if(!file.exists(osm_cachefile)){
        # Get a piece that is bigger than we need, to avoid gaps in plots when
        # the plot limits are extended by R (which depends on the figure size
        # too)
        osm_backdrop = try(openmap(upperLeft = c(upper_right[2] + dmax, lower_left[1]  - dmax), 
                               lowerRight= c(lower_left[2]  - dmax, upper_right[1] + dmax), 
                               type=backdrop, zoom=16))
        if(is(osm_backdrop, 'try-error')){
            stop('Download failed')
        }
        # Reproject it
        osm_backdrop_reproj = openproj(osm_backdrop, proj4string(CRS("+init=EPSG:4326")))
        saveRDS(osm_backdrop_reproj, file=osm_cachefile)
    }else{
        # Read from cache
        osm_backdrop_reproj = readRDS(osm_cachefile)
    }
    return(osm_backdrop_reproj)
}

#
# Make a 'zoomed-in' plot with a raster, and a backdrop
#
regional_plot<-function(lower_left, upper_right, backdrop='osm', plotvar='max_stage', colpal = pal1_trans, 
    zlim=c(0, 10)){

    # Slightly reduce the LHS white-space margin. This was to reduce a gap
    # (whitespace) between combined figures
    orig_mar = par('mar')
    newmar = c(5.1, 2.5, 4.1, 2.1)
    par('mar' = newmar)
    on.exit(par('mar' = orig_mar))

    plot_ext = extent(c(lower_left[1], upper_right[1], lower_left[2], upper_right[2]))
    plot_asp = 1/cos(0.5*(plot_ext@ymin + plot_ext@ymax)/180 * pi)
    r1 = raster(plot_ext, nrows=50, ncols=50) # to consistently initialise plot

    if(plotvar %in% names(TIF_FILES)){

        if(plotvar %in% FLOW_VARIABLE_NAMES){
            # Plot flow variables. These have a legend on the panel
            r1 = setValues(r1, rep(zlim, 50*50/2))
            plot(max(min(r1, zlim[2]), zlim[1]), ext=plot_ext, 
                zlim=zlim, asp=plot_asp, maxpixels=MAXPIXELS, col=colpal,
                xlab="", ylab="", cex.axis=1.2, 
                xlim=c(lower_left[1], upper_right[1]), ylim=c(lower_left[2], upper_right[2]), 
                legend.shrink=1, axis.args=list(cex.axis=1.6))

        }else if(plotvar %in% INUNDATION_RATE_NAMES){
            # Plot inundation exceedance-rates. These do not have a legend on the panel
            r1 = setValues(r1, rep(c(0, 1), 50*50/2)) # Fake data
            image(r1, 
                asp=plot_asp, maxpixels=MAXPIXELS, col='white',
                xlab="", ylab="", cex.axis=1.2, 
                xlim=c(lower_left[1], upper_right[1]), ylim=c(lower_left[2], upper_right[2]), 
                cex.axis=1.2)
        }

        # Find rasters to plot, searching for those overlapping with the figure limits
        tmp_ll = par('usr')[c(1,3)]
        tmp_ur = par('usr')[c(2,4)]
        #tmp_dx = tmp_ur - tmp_ll # Extra buffer to get more tifs. Having problems with gaps.
        #tmp_ll = tmp_ll - tmp_dx/3
        #tmp_ur = tmp_ur + tmp_dx/3
        raster_inds = which(is_tif_touching_bbox(tmp_ll, tmp_ur, plotvar) & (TIF_RES[[plotvar]] < MIN_TIF_RES))
        raster_files_to_plot = TIF_FILES[[plotvar]][ raster_inds ]

    }else if(plotvar == 'JATWC_2_inundation_zones'){
        # Plotting shapefiles
        r1 = setValues(r1, rep(as.factor(c(1,2,3)), length.out=50*50))
        plot(r1, legend=FALSE, col=c('green', 'blue', 'red'))
    }else{
        stop(paste0('unknown plotvar ', plotvar))
    }

    tmp_ll = par('usr')[c(1,3)]
    tmp_ur = par('usr')[c(2,4)]

    # Get the OSM data
    dbox = tmp_ur - tmp_ll #upper_right - lower_left
    osm_backdrop_reproj = get_osm_backdrop_reproj(tmp_ll, tmp_ur, dbox, backdrop)
    #osm_backdrop_reproj = get_osm_backdrop_reproj(lower_left, upper_right, dbox, backdrop)
    # Add to plot
    plot(osm_backdrop_reproj, add=TRUE)

    # Make transparent a bit
    #rect(lower_left[1]-dbox[1]*1.5, lower_left[2]-dbox[2]*1.5, upper_right[1]+dbox[1]*1.5, upper_right[2]+dbox[2]*1.5, 
    rect(tmp_ll[1], tmp_ll[2], tmp_ur[1], tmp_ur[2], 
        border=NA, col=rgb(1,1,1,alpha=0.5,maxColorValue=1))

    if(plotvar %in% FLOW_VARIABLE_NAMES){
        # Overplot with the flow variable rasters
        for(i in 1:length(raster_files_to_plot)){
            r1 = raster(raster_files_to_plot[i])
            image(max(min(r1, zlim[2]), zlim[1]), 
                  zlim=zlim, add=TRUE, maxpixels=MAXPIXELS, col=colpal)
        }
    }else if(plotvar %in% INUNDATION_RATE_NAMES){
        # Overplot with the inudation rasters, after classifing into discrete bands
        for(i in 1:length(raster_files_to_plot)){
            r1 = raster(raster_files_to_plot[i])
            r1b = reclassify(r1, inundation_classmat)
            image(r1b, add=TRUE, col=inundation_classmat_col, 
                # Need to specify zlim here to prevent unfortunate defaults for single-class plot
                zlim=c(1, length(inundation_classmat_col)))
        }
    }else if(plotvar == 'JATWC_2_inundation_zones'){

        # Overplot with the shapefiles. BEWARE THIS CAN GET THE "INSIDE" WRONG ON OCCASION
        #plot(LAND_WARNING_ZONE, col='red', border='red', density=20, add=TRUE)
        #plot(MARINE_WARNING_ZONE, col='blue', border='blue', density=20, add=TRUE)
        #plot(NO_THREAT_ZONE, col='green', border='green', density=20, add=TRUE)

        # Work around occasional problems with hatching in the routines above
        plot(LAND_WARNING_ZONE, col=NA, border='red', add=TRUE)
        plot(MARINE_WARNING_ZONE, col=NA, border='blue', add=TRUE)
        plot(NO_THREAT_ZONE, col=NA, border='green', add=TRUE)
        hatchedLayer(st_as_sf(LAND_WARNING_ZONE), col='red', pattern='right2left', density=4, add=TRUE)
        hatchedLayer(st_as_sf(MARINE_WARNING_ZONE), col='blue', pattern='right2left', density=4, add=TRUE)
        hatchedLayer(st_as_sf(NO_THREAT_ZONE), col='green', pattern='right2left', density=4, add=TRUE)

        par(xpd=TRUE)
        legend('right', legend=c('No Threat', 'Marine W.', 'Land W.'), 
            fill = c('green', 'blue', 'red'), inset=-0.20)
        par(xpd=FALSE)
    }

    # Scale bar (outside the plot)
    par(xpd=TRUE)
    d_km = pointDistance(matrix(lower_left, ncol=2), cbind(upper_right[1], lower_left[2]), 
        lonlat=TRUE)/2 * 1/1000
    d_km = signif(d_km, 1)
    xy_lower = c(plot_ext@xmin + 0.6*(plot_ext@xmax-plot_ext@xmin), 
                plot_ext@ymin  - 0.26*(plot_ext@ymax-plot_ext@ymin))
    scalebar(d=d_km, xy=xy_lower, type='bar', divs=2, lonlat=TRUE, below='km')
    par(xpd=FALSE)

}

# Make a few plots for the ith polygonal region defined in ZOOM_PLOTS_REGIONS
make_region_plot<-function(i, out_dir='.'){

    plot_region = extent(ZOOM_PLOTS_REGIONS[i,])

    dx = c(plot_region@xmax - plot_region@xmin, plot_region@ymax - plot_region@ymin)

    plot_centre = c(plot_region@xmin, plot_region@ymin) + dx/2

    plot_ll = plot_centre - max(dx)/2
    plot_ur = plot_centre + max(dx)/2

    siteName = gsub(' ', '-', ZOOM_PLOTS_REGIONS$Site[i])

    # Plot dimensions  # 80% of A4 in inches
    P_W = 7 # 8.27 * 0.8
    P_H = 5 # 11.69 * 0.8
    # resolution
    P_RES = 200 #300


    plot_index = as.character(100 + i)
    plot_index = substring(plot_index, 2, nchar(plot_index))

    # Plot depth
    depth_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_depth_1in250084pc.png')
    png(depth_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_depth_1in2500_84pc', 
        zlim=c(0,4), backdrop='osm', colpal=palDepthStage_trans)
    title(paste0(siteName, ': max-depth (m) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot max-stage
    stage_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_stage_1in250084pc.png')
    png(stage_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_stage_1in2500_84pc', 
        zlim=c(0,6), backdrop='osm', colpal=palDepthStage_trans)
    title(paste0(siteName, ': max-stage (m) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot max-speed
    speed_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_speed_1in250084pc.png')
    png(speed_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_speed_1in2500_84pc', 
        zlim=c(0,5), backdrop='osm', colpal=palSpeed_trans)
    title(paste0(siteName, ': max-speed (m/s) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot max-flux
    flux_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_flux_1in250084pc.png')
    png(flux_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_flux_1in2500_84pc', 
        zlim=c(0,50), backdrop='osm', colpal=palFlux_trans)
    title(paste0(siteName, ': max-flux (m.m/s) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    ## Combine into a single plot
    four_images = image_read(c(depth_file, stage_file, speed_file, flux_file))
    iminfo = image_info(four_images[1])
    combined_figs = image_montage(four_images, tile="2x2", geometry=paste0(iminfo$width, 'x', iminfo$height))
    image_write(combined_figs, path=paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_combined_fig.png'), format='png')
    

    # Plot model based warning zones
    png(paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_JATWC_zones.png'), 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='JATWC_2_inundation_zones', 
        zlim=c(1,5), backdrop='osm', colpal=palDepthStage_trans)
    title(paste0(siteName, ': JATWC inundation zones'), cex.main=1.6)
    dev.off()

    # Plot inundation exceedance-rates, logic tree mean
    inundation_ltm_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_inundation_rate_ltm.png')
    png(inundation_ltm_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, plotvar='inundation_rate_ltm', 
        zlim=c(NA, NA), backdrop='osm', colpal=NA)
    title(paste0(siteName, ': Inundation rate \n logic tree mean'), cex.main=1.9)
    dev.off()

    # Plot inundation exceedance-rates, 84th percentile
    inundation_84pc_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_inundation_rate_84pc.png')
    png(inundation_84pc_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, plotvar='inundation_rate_84pc', 
        zlim=c(NA, NA), backdrop='osm', colpal=NA)
    title(paste0(siteName, ': Inundation rate \n 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot inundation exceedance-rates, 16th percentile
    inundation_16pc_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_inundation_rate_16pc.png')
    png(inundation_16pc_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, plotvar='inundation_rate_16pc', 
        zlim=c(NA, NA), backdrop='osm', colpal=NA)
    title(paste0(siteName, ': Inundation rate \n 16th percentile'), cex.main=1.9)
    dev.off()

    # Legend for inundation rate
    inundation_rate_legend_file = paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_inundation_rate_legend.png')
    png(inundation_rate_legend_file, width=P_W, height=P_H, units='in', res=P_RES)
    plot(c(0,1), c(0, 1), col='white', axes=FALSE, frame.plot=FALSE, ann=FALSE) 
    k = 1:4 # Colours that are visible
    legend('center', legend=inundation_classmat_labels[k], fill=inundation_classmat_col[k], title='Inundation rate (events/year)', 
        cex=2, bg='white', border='black')
    dev.off()

    ## Combine into a single plot
    four_images = image_read(c(inundation_ltm_file, inundation_84pc_file, inundation_16pc_file, inundation_rate_legend_file))
    iminfo = image_info(four_images[1])
    combined_figs = image_montage(four_images, tile="2x2", geometry=paste0(iminfo$width, 'x', iminfo$height))
    image_write(combined_figs, path=paste0(out_dir, '/plot_region_', plot_index, '_osm_', siteName, '_inundation_combined_fig.png'), format='png')
}
library(parallel)
my_inds = splitIndices(length(ZOOM_PLOTS_REGIONS), batch_count)[[batch_number]]
out_dir = './test_figs'
dir.create(out_dir, showWarnings=FALSE)
for(i in my_inds) make_region_plot(i, out_dir=out_dir)

