library(rptha)
library(cptcity)
library(OpenStreetMap)
library(magick)
library(sf)
library(cartography)
library(httr)
library(jsonlite)
library(ggspatial)
library(ggpattern)
library(ggplot2)

# For running in multple batches
batch_number = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[1])
batch_count = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[2])

#
# INPUT PARAMETERS
#

# Regions for zoom plots. These were manually created
ZOOM_PLOTS_REGIONS = readOGR('../probabilistic_inundation/plots/plot_regions', verbose=FALSE) 
# Order from south to north to make identification easier
k = order(coordinates(ZOOM_PLOTS_REGIONS)[,2])
ZOOM_PLOTS_REGIONS = ZOOM_PLOTS_REGIONS[k,]

# Max pixels in image plot
MAXPIXELS=1e+07

# Do not plot tifs with resolution coarser than this. Can speed things up
# depending on multidomain design
MIN_TIF_RES = 1.1/(60*5)

##
## Various tifs derived from the model, by sea level (names in millimetres)
##
TIF = list(sl_000=list(), sl_800=list())

TIF$sl_000$max_stage$files = Sys.glob('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.0/RUN_20241219_164414473/max_stage*.tif')
TIF$sl_000$max_speed$files = Sys.glob('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.0/RUN_20241219_164414473/max_speed*.tif')

TIF$sl_800$max_stage$files = Sys.glob('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.8/RUN_20241219_164342781/max_stage*.tif')
TIF$sl_800$max_speed$files = Sys.glob('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.8/RUN_20241219_164342781/max_speed*.tif')

TIF$sl_800$max_stage_diff$files = Sys.glob('./diff/sl_800mm_max_stage_domain*_diff.tif')
TIF$sl_800$max_speed_diff$files = Sys.glob('./diff/sl_800mm_max_speed_domain*_diff.tif')

# check that none of them are empty
for (sl in names(TIF)){
    for (var in names(TIF[[sl]])){
        if(length(TIF[[sl]][[var]]$files) == 0){
            stop(paste0('No files found for ', sl, ' ', var))
        }
    }
}

# This has a strong influence on the quality of the background image, and the
# time it takes to download. Probably needs experimentation in different cases.
# Consider using library(basemaps) instead. Note hacking of this inside the main function.
BACKDROP_ZOOM_DEFAULT = 14

# Store the tif extents
for (sl in names(TIF)){
    for (var in names(TIF[[sl]])){
        # store the extents as a list of of raster::extent objects 
        TIF[[sl]][[var]][['extents']] = lapply(TIF[[sl]][[var]]$files, function(tif_file) extent(raster(tif_file)))

    }
}

# Store the tif resolution
for (sl in names(TIF)){
    for (var in names(TIF[[sl]])){
        # store the resolution as a list of of raster::extent objects 
        TIF[[sl]][[var]][['res']] = unlist(lapply(TIF[[sl]][[var]]$files, function(tif_file) res(raster(tif_file))[1]))
    }
}

# Color palettes for flow variables
source('colour_palette.R')

# For each tif associated with the plot var, return TRUE if it intersects with the
# bounding box defined by lower_left and upper_right, and FALSE otherwise.
is_tif_touching_bbox<-function(lower_left, upper_right, plotvar='max_stage', sl='sl_800'){

    rast_extents = TIF[[sl]][[plotvar]]$extents
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
get_osm_backdrop_reproj<-function(lower_left, upper_right, dbox, backdrop, backdrop_zoom){
    osm_dir = paste0('.osm_tile_cache_', backdrop)
    if(!dir.exists(osm_dir)){
        dir.create(osm_dir)
    }
    osm_cachefile = paste0('.osm_tile_cache_', backdrop, "/", upper_right[1], upper_right[2], lower_left[1], lower_left[2], '.RDS')
    dmax = max(dbox)
    if(!file.exists(osm_cachefile)){
        # Get a piece that is bigger than we need, to avoid gaps in plots when
        # the plot limits are extended by R (which depends on the figure size
        # too)
        print(paste0('Downloading OSM data for ', backdrop))
        osm_backdrop = try(openmap(upperLeft = c(upper_right[2] + dmax, lower_left[1]  - dmax), 
                               lowerRight= c(lower_left[2]  - dmax, upper_right[1] + dmax), 
                               type=backdrop, zoom=backdrop_zoom))
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

# Fetch the tsunami evacuation zones
fetch_tsunami_zones<-function(){
        tsunami_evac_url = "https://services1.arcgis.com/vkTwD8kHw2woKBqV/arcgis/rest/services/Tsunami_LGA_AEIP_Data_V2/FeatureServer/0/query"
        query_params <- list(
            where = "LGA_NAME IN ('LIVINGSTONE SHIRE', 'GLADSTONE REGIONAL')",
            outFields = "Evac_zone",
            f = "geojson"
        )
        response <- httr::GET(url = tsunami_evac_url, query = query_params)
        if (status_code(response) == 200) {
            geojson_data <- content(response, as = "text")
            polygon_data <- st_read(geojson_data, quiet = TRUE)
        } else {
            stop("Failed to retrieve data from the ArcGIS server")
        }
    return(polygon_data)
}

# Make a 'zoomed-in' plot with a raster, and a backdrop
regional_plot<-function(lower_left, upper_right, backdrop='osm', plotvar='max_stage', colpal = pal1_trans, 
    zlim=c(0, 10), sl){

    # Scale the zoom for the size of the area
    if(prod(upper_right - lower_left) > 0.02){
        backdrop_zoom = BACKDROP_ZOOM_DEFAULT - 2
    }else{
        backdrop_zoom = BACKDROP_ZOOM_DEFAULT
    }

    # Slightly reduce the LHS white-space margin. This was to reduce a gap
    # (whitespace) between combined figures
    orig_mar = par('mar')
    newmar = c(4.9, 2.3, 3.9, 1.9)
    par('mar' = newmar)
    on.exit(par('mar' = orig_mar))

    plot_ext = extent(c(lower_left[1], upper_right[1], lower_left[2], upper_right[2]))
    plot_asp = 1/cos(0.5*(plot_ext@ymin + plot_ext@ymax)/180 * pi)
    r1 = raster(plot_ext, nrows=50, ncols=50) # to consistently initialise plot

    # Plot flow variables. These have a legend on the panel
    r1 = setValues(r1, rep(zlim, 50*50/2))
    plot(max(min(r1, zlim[2]), zlim[1]), ext=plot_ext, 
        zlim=zlim, asp=plot_asp, maxpixels=MAXPIXELS, col=colpal,
        xlab="", ylab="", cex.axis=1.2, 
        xlim=c(lower_left[1], upper_right[1]), ylim=c(lower_left[2], upper_right[2]), 
        legend.shrink=1, axis.args=list(cex.axis=1.6))

    # Find rasters to plot, searching for those overlapping with the figure limits
    tmp_ll = par('usr')[c(1,3)]
    tmp_ur = par('usr')[c(2,4)]
    #tmp_dx = tmp_ur - tmp_ll # Extra buffer to get more tifs. Having problems with gaps.
    #tmp_ll = tmp_ll - tmp_dx/3
    #tmp_ur = tmp_ur + tmp_dx/3
    raster_inds = which(is_tif_touching_bbox(tmp_ll, tmp_ur, plotvar, sl) & (TIF[[sl]][[plotvar]]$res < MIN_TIF_RES))
    raster_files_to_plot = TIF[[sl]][[plotvar]][['files']][ raster_inds ]

    tmp_ll = par('usr')[c(1,3)]
    tmp_ur = par('usr')[c(2,4)]

    # Get the OSM data
    dbox = tmp_ur - tmp_ll
    osm_backdrop_reproj = get_osm_backdrop_reproj(tmp_ll, tmp_ur, dbox, backdrop, backdrop_zoom)

    # Add to plot
    plot(osm_backdrop_reproj, add=TRUE)

    # Make transparent a bit
    rect(tmp_ll[1], tmp_ll[2], tmp_ur[1], tmp_ur[2], 
        border=NA, col=rgb(1,1,1,alpha=0.25,maxColorValue=1))

    # Overplot with the flow variable rasters
    for(i in 1:length(raster_files_to_plot)){
        r1 = raster(raster_files_to_plot[i])
        image(max(min(r1, zlim[2]), zlim[1]), 
                zlim=zlim, add=TRUE, maxpixels=MAXPIXELS, col=colpal)
        }
    # plot(ZERO_CONTOUR, add=TRUE, col='black')

    # Scale bar (inside the plot)
    d_km = pointDistance(matrix(lower_left, ncol=2), cbind(upper_right[1], lower_left[2]), 
        lonlat=TRUE)/2 * 1/1000
    d_km = signif(d_km, 1)
    xy_lower = c(plot_ext@xmin + 0.05*(plot_ext@xmax-plot_ext@xmin), 
                plot_ext@ymin  + 0.05*(plot_ext@ymax-plot_ext@ymin))
    scalebar(d=d_km, xy=xy_lower, type='bar', divs=2, lonlat=TRUE, below='km')
    par(xpd=FALSE)

}

# Make a few plots for the ith polygonal region defined in ZOOM_PLOTS_REGIONS
make_region_plot<-function(i, out_dir='.', prep_dir='./prep', sl=0){
    plot_region = extent(ZOOM_PLOTS_REGIONS[i,])

    SCALE_FACTOR = 1.0 # Increase/decrease to enlarge/shrink the plot domain
    dx = c(plot_region@xmax - plot_region@xmin, plot_region@ymax - plot_region@ymin) * SCALE_FACTOR

    plot_centre = c(plot_region@xmin, plot_region@ymin) + dx/2

    plot_ll = plot_centre - max(dx)/2
    plot_ur = plot_centre + max(dx)/2

    siteName = gsub(' ', '-', ZOOM_PLOTS_REGIONS$Site[i])

    # Plot dimensions  # 80% of A4 in inches
    P_W = 7 # 8.27 * 0.8
    P_H = 7 # 11.69 * 0.8
    # resolution
    P_RES = 200 #300

    plot_index = as.character(100 + i)
    plot_index = substring(plot_index, 2, nchar(plot_index))

    # 1. Plot stage 0m
    stage_000_file = paste0(prep_dir, siteName, '_sl000_stage.png')
    png(
        stage_000_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_stage', 
        zlim=c(0,5),
        colpal=palDepthStage_trans,
        sl='sl_000'
    )
    title(sub="Max stage (m) SL[0.0m]", cex.main=1.9, cex.sub=1.5)
    dev.off()

    # 2. Plot max-speed 0m
    speed_000_file = paste0(prep_dir, siteName, '_sl000_speed.png')
    png(
        speed_000_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_speed', 
        zlim=c(0,5),
        colpal=palSpeed_trans,
        sl='sl_000'
    )
    title(sub="Max Speed (m/s) SL[0.0m]", cex.main=1.9, cex.sub=1.5)
    dev.off()

    # # 3. Plot max-stage 800mm
    stage_800_file = paste0(prep_dir, siteName, '_sl800_stage.png')
    png(
        stage_800_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_stage', 
        zlim=c(0,5),
        colpal=palDepthStage_trans,
        sl='sl_800'
    )
    title(main=siteName, sub="Max stage (m) SL[0.8m]", cex.main=1.9, cex.sub=1.5)
    dev.off()

    # 4. Plot max-speed 800mm
    speed_800_file = paste0(prep_dir, siteName, '_sl800_speed.png')
    png(
        speed_800_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_speed', 
        zlim=c(0,5),
        colpal=palSpeed_trans,
        sl='sl_800'
    )
    title(main=siteName, sub="Max Speed (m/s) SL[0.8m]", cex.main=1.9, cex.sub=1.5)
    dev.off()

    # 5. Plot difference in stage from sl_000
    stage_diff_file = paste0(prep_dir, siteName, '_stage_diff.png')
    png(
        stage_diff_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_stage_diff', 
        zlim=c(-1,1)+0.8,
        colpal=palDiff,
        sl='sl_800'
    )
    title(sub = expression('Stage Difference (m) SL[0.8m] - SL[0m]'), cex.main=1.9, cex.sub=1.5)
        dev.off()

    # 6. Plot difference in speed from sl_000
    speed_diff_file = paste0(prep_dir, siteName, '_speed_diff.png')
    png(
        speed_diff_file, 
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(
        lower_left = plot_ll,
        upper_right=plot_ur,
        plotvar='max_speed_diff', 
        zlim=c(-2, 2),
        colpal=palDiff,
        sl='sl_800'
    )
    title(sub = expression('Speed Difference (m/s) SL[0.8m] - SL[0m]'), cex.main=1.9, cex.sub=1.5)
    dev.off()

    # 7. Combine stage 000, stage 800 and stage diff
    images_list = image_read(c(stage_000_file, stage_800_file, stage_diff_file))
    iminfo = image_info(images_list[1])
    combined_figs = image_montage(images_list, tile="3x1", geometry=paste0(iminfo$width, 'x', iminfo$height))
    image_write(combined_figs, path=paste0(out_dir, siteName, '_stage.png'), format='png')

    # 8. combine speed 000, speed 800 and speed diff
    images_list = image_read(c(speed_000_file, speed_800_file, speed_diff_file))
    iminfo = image_info(images_list[1])
    combined_figs = image_montage(images_list, tile="3x1", geometry=paste0(iminfo$width, 'x', iminfo$height))
    image_write(combined_figs, path=paste0(out_dir, siteName, '_speed.png'), format='png')
}


out_dir = './figures_unadjusted/'
dir.create(out_dir, showWarnings=FALSE)
prep_dir = './figures_unadjusted/prep/'
dir.create(prep_dir, showWarnings=FALSE)
# create figures for all regions
for(i in 1:nrow(ZOOM_PLOTS_REGIONS)){
    print(ZOOM_PLOTS_REGIONS$Site[i])
    make_region_plot(i, out_dir=out_dir, prep_dir=prep_dir)
}

library(terra)
mins<-sapply(Sys.glob('diff/sl_800mm_max_stage_domain*_diff.tif'), function(x) min(as.matrix(rast(x))))
all_min = min(mins)
