library(rptha)
library(cptcity)
library(OpenStreetMap)
library(magick)
library(sf)
library(cartography)
library(httr)
library(jsonlite)
library(parallel)


# For running in multple batches
batch_number = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[1])
batch_count = 1 #as.numeric(commandArgs(trailingOnly=TRUE)[2])

#
# INPUT PARAMETERS
#

# Regions for zoom plots. These were manually created
ZOOM_PLOTS_REGIONS = readOGR('plot_regions', verbose=FALSE) 
# Order from south to north to make identification easier
k = order(coordinates(ZOOM_PLOTS_REGIONS)[,2])
ZOOM_PLOTS_REGIONS = ZOOM_PLOTS_REGIONS[k,]

ATWS_COASTAL_ZONES = readOGR('../../jatwc_to_inundation/ATWS_ZONES/ATWS_Zones_V2_2_4/ATWS_Zones_V2_2_4.shp')
ATWS_COASTAL_ZONES$ATWS_Zones = gsub(" ", "-", ATWS_COASTAL_ZONES$ATWS_Zones) # Replace spaces with '-' to help list naming
# ATWS_COASTAL_ZONES is missing projection info
crs(ATWS_COASTAL_ZONES) = crs(ZOOM_PLOTS_REGIONS)

#
# Mapped zones matching ATWS coastal zones
#
jatwc_zones_mapped = basename(Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/JATWC_inundation_zones/Inundation_zones/*'))
NO_THREAT_ZONE = vector(mode='list', length=length(jatwc_zones_mapped))
names(NO_THREAT_ZONE) = jatwc_zones_mapped
MARINE_WARNING_ZONE = NO_THREAT_ZONE
LAND_WARNING_ZONE = NO_THREAT_ZONE
for(i in 1:length(jatwc_zones_mapped)){
    nm = jatwc_zones_mapped[i]
    NO_THREAT_ZONE[[nm]] = readOGR(paste0('../../../data_package/gladstone_tsunami_modelling_project_data_package/JATWC_inundation_zones/Inundation_zones/', 
        nm, '/', nm, '_no_threat_with-PTHA-exrate-limit_84pc_4e-04'), verbose=FALSE)
    MARINE_WARNING_ZONE[[nm]] = readOGR(paste0('../../../data_package/gladstone_tsunami_modelling_project_data_package/JATWC_inundation_zones/Inundation_zones/', 
        nm, '/', nm, '_marine_warning_with-PTHA-exrate-limit_84pc_4e-04'), verbose=FALSE)
    LAND_WARNING_ZONE[[nm]] = readOGR(paste0('../../../data_package/gladstone_tsunami_modelling_project_data_package/JATWC_inundation_zones/Inundation_zones/', 
        nm, '/', nm, '_land_warning_with-PTHA-exrate-limit_84pc_4e-04'), verbose=FALSE)
}

# Max pixels in image plot
MAXPIXELS=1e+07

# Do not plot tifs with resolution coarser than this. Can speed things up
# depending on multidomain design
MIN_TIF_RES = 1.1/(60*5)

# ZERO_CONTOUR = readOGR('../../../elevation/countour_0/Australian_contour.shp')

##
## Various tifs derived from the model
##
TIF_FILES = list()

#
# 1/2500, 84th percentile
#
TIF_FILES$max_stage_1in2500_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/tsunami_flow_variables/flow_variables_1in2500_at_84th_percentile/max_stage_with_exceedance_rate_1in2500_at_84th_percentile/*.tif')
# Use depth because initial condition is 0 m. Depth above  above IC reduces to stage.
TIF_FILES$max_depth_1in2500_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/tsunami_flow_variables/flow_variables_1in2500_at_84th_percentile/max_depth_with_exceedance_rate_1in2500_at_84th_percentile/*.tif') 
TIF_FILES$max_flux_1in2500_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/tsunami_flow_variables/flow_variables_1in2500_at_84th_percentile/max_flux_with_exceedance_rate_1in2500_at_84th_percentile/*.tif')
TIF_FILES$max_speed_1in2500_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/tsunami_flow_variables/flow_variables_1in2500_at_84th_percentile/max_speed_with_exceedance_rate_1in2500_at_84th_percentile/*.tif')
TIF_FILES$min_stage_1in2500_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/tsunami_flow_variables/flow_variables_1in2500_at_84th_percentile/min_stage_with_exceedance_rate_1in2500_at_84th_percentile/*.tif')
TIF_FILES$elevation = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/elevation_in_model/*.tif')

#
# Rate of inundation
#
TIF_FILES$inundation_rate_ltm = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/probabilistic_inundation_hazard/inundation_rate/logic_tree_mean/*.tif')
TIF_FILES$inundation_rate_84pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/probabilistic_inundation_hazard/inundation_rate/84th_percentile/*.tif')
TIF_FILES$inundation_rate_16pc = Sys.glob('../../../data_package/gladstone_tsunami_modelling_project_data_package/probabilistic_inundation_hazard/inundation_rate/16th_percentile/*.tif')

# This has a strong influence on the quality of the background image, and the
# time it takes to download. Probably needs experimentation in different cases.
# Consider using library(basemaps) instead. Note hacking of this inside the main function.
BACKDROP_ZOOM_DEFAULT = 14

# We use different kinds of plot for different variables. Group like this
FLOW_VARIABLE_NAMES = c(
    paste0('max_', c('stage', 'depth', 'flux', 'speed'), '_1in2500_84pc'), 
    'elevation', 
    'min_stage_1in2500_84pc'
)

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
colorless = rgb(255, 255, 255, alpha=0, maxColorValue=255)
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
palDepthStage_trans[1] = colorless
# reversed arendal temperature
palDS2_rev = rev(palDS2)
tmp = col2rgb(palDS2_rev)
palDepthStage_trans_rev = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palDepthStage_trans_rev[length(palDepthStage_trans_rev)] = colorless
# Yet another alternative
palDS3 = cpt("cb_seq_YlOrRd_05", n=400)
tmp = col2rgb(palDS3)
palDS3_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# And another
palSpeed = cpt("h5_jet", n=400)
tmp = col2rgb(palSpeed)
palSpeed_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palSpeed_trans[1] = colorless
# And another
palFlux = colorRampPalette(cpt("oc_sst", n=400), bias=1.8)(400) # Concentrate variation around low values
tmp = col2rgb(palFlux)
palFlux_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palFlux_trans[1] = colorless

# Colours for speed that will split it into 4 hazard categories using zlim=c(0, 6)
palSpeedHazard = rep(c('beige', 'orange', 'red', 'black'), each=100)
tmp = col2rgb(palSpeedHazard)
palSpeedHazard_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# make the first colour transparent
palSpeedHazard_trans[1] = colorless

#
# Classes to colour the inundation rates
#
inundation_classmat = matrix(
    # Lower, upper, integerID
    c(1    ,   999999999999999,     1,
      1/100  , 1,       2, 
      1/500  , 1/100  , 3, 
      1/2500 , 1/500  , 4, 
      1/10000, 1/2500 , 5, 
      -1     , 1/10000, 6), 
    ncol=3, byrow=TRUE)
inundation_classmat_labels = c('>1', '1 - 1/100', '1/500 - 1/100', '1/2500 - 1/500', '1/10000 - 1/2500', '< 1/10000')
# Colours with transparency
inundation_classmat_col = c('cornflowerblue','red', 'orange', 'yellow', 'grey', NA)
icct = col2rgb(inundation_classmat_col, alpha=TRUE)
# the blue should be slightly transparent, and the grey
icct[4,1] = 0.7*icct[4,1]
icct[4,5] = 0.7*icct[4,5]
inundation_classmat_col = rgb(icct[1,], icct[2,], icct[3,], icct[4,], maxColorValue=255)

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
            where = "LGA_NAME IN ('LIVINGSTONE SHIRE', 'GLADSTONE REGIONAL', 'ROCKHAMPTON REGIONAL')",
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


#
# Make a 'zoomed-in' plot with a raster, and a backdrop
#
regional_plot<-function(lower_left, upper_right, backdrop='esri-imagery', plotvar='max_stage', colpal = pal1_trans, 
    zlim=c(0, 10), atws_zone = ""){

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

    }else if(plotvar == 'JATWC_inundation_zones'){
        # Plotting shapefiles
        r1 = setValues(r1, rep(as.factor(c(1,2,3)), length.out=50*50))
        # colours taken from background colours on JATWC website
        cols = c(
            rgb(51,211,57,255,maxColorValue=255),
            rgb(58,59,248,255,maxColorValue=255),
            rgb(218,83,83,255,maxColorValue=255)
        )
        plot(r1, legend=FALSE, col=cols)
    }else{
        stop(paste0('unknown plotvar ', plotvar))
    }

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
    }else if(plotvar == 'JATWC_inundation_zones'){
        bbox_poly <- st_as_sf(
            sf::st_sfc(
            sf::st_polygon(list(
                matrix(c(
                tmp_ll[1], tmp_ll[2],
                tmp_ur[1], tmp_ll[2],
                tmp_ur[1], tmp_ur[2],
                tmp_ll[1], tmp_ur[2],
                tmp_ll[1], tmp_ll[2]
                ), ncol=2, byrow=TRUE)
            )),
            crs = 4326
            )
        )

        # get the  current evacuation zones
        current_evac_polygon <- fetch_tsunami_zones()
        current_evac_polygon <- st_make_valid(current_evac_polygon)
        # for some reason, sometimes holes in the polygon only display correctly by taking the compliment twice
        try({
            current_evac_polygon_compliment <- st_difference(bbox_poly, current_evac_polygon)
            current_evac_polygon <- st_difference(bbox_poly, current_evac_polygon_compliment)
        })

        # plot the current evacuation
        grey_transp <- rgb(0.5, 0.5, 0.5, 0.25, maxColorValue=1)
        colourless <- rgb(255, 255, 255, alpha=0, maxColorValue=255)
        plot(current_evac_polygon, col=grey_transp, add=TRUE, border='darkgrey', lwd=2)

        # Work around for occasional problems 
        plot(LAND_WARNING_ZONE[[atws_zone]], col=NA, border=cols[3], add=TRUE)
        plot(MARINE_WARNING_ZONE[[atws_zone]], col=NA, border=cols[2], add=TRUE)
        plot(NO_THREAT_ZONE[[atws_zone]], col=NA, border=cols[1], add=TRUE)
        # Overplot with the shapefiles.
        hatchedLayer(st_as_sf(LAND_WARNING_ZONE[[atws_zone]]), col=cols[3], pattern='right2left', density=4, add=TRUE)
        hatchedLayer(st_as_sf(MARINE_WARNING_ZONE[[atws_zone]]), col=cols[2], pattern='right2left', density=4, add=TRUE)
        hatchedLayer(st_as_sf(NO_THREAT_ZONE[[atws_zone]]), col=cols[1], pattern='right2left', density=4, add=TRUE)

        par(xpd=TRUE)
        legend(
            'topright',
            legend=c('No Threat', 'Marine Warning', 'Land Warning', 'Current Evacuation Zone'), 
            fill = c(NA, NA, NA, grey_transp),
            border=c(cols[1], cols[2], cols[3], 'darkgrey'),
            inset=0.005,
        )
        par(xpd=FALSE)
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
make_region_plot<-function(i, out_dir='.', prep_dir='./prep'){
    plot_region = extent(ZOOM_PLOTS_REGIONS[i,])

    SCALE_FACTOR = 1.0 # Increase/decrease to enlarge/shrink the plot domain
    dx = c(plot_region@xmax - plot_region@xmin, plot_region@ymax - plot_region@ymin) * SCALE_FACTOR

    plot_centre = c(plot_region@xmin, plot_region@ymin) + dx/2

    plot_ll = plot_centre - max(dx)/2
    plot_ur = plot_centre + max(dx)/2

    siteName = gsub(' ', '-', ZOOM_PLOTS_REGIONS$Site[i])

    # Find the ATWS zone  
    intersect_zoom_with_zone = st_intersection(st_as_sf(ZOOM_PLOTS_REGIONS[i,]), st_make_valid(st_as_sf(ATWS_COASTAL_ZONES)))
    if(nrow(intersect_zoom_with_zone) == 0){
         stop('error: Could not find intersecting ATWS zone')
    }else if(nrow(intersect_zoom_with_zone) == 1){
        # A single zone matches
        atws_zone = intersect_zoom_with_zone$ATWS_Zones[1]
    }else{
        if(all(intersect_zoom_with_zone$Site ==  intersect_zoom_with_zone$ATWS_Zones[1])){
            # A single zone in multiple pieces
            atws_zone = intersect_zoom_with_zone$ATWS_Zones[1]
        }else{
            # Several different zones match, need to choose one or throw an error
            stop('error: Intersection with multiple ATWS zones, need to implement a method to select them and get the right model shapefiles (noting they can overlap).')
        }
    }


    # Plot dimensions  # 80% of A4 in inches
    P_W = 7 # 8.27 * 0.8
    P_H = 7 # 11.69 * 0.8
    # resolution
    P_RES = 200 #300


    plot_index = as.character(100 + i)
    plot_index = substring(plot_index, 2, nchar(plot_index))

    #
    # 1/2500 84th percentile
    #

    # Plot depth
    depth_file = paste0(prep_dir, siteName, '_depth_1in2500_84pc.png')
    png(depth_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_depth_1in2500_84pc', 
        zlim=c(0,5), colpal=palDepthStage_trans, atws_zone = atws_zone)
    title(paste0(siteName, ': Max Depth (m) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # # Plot max-stage
    stage_file = paste0(prep_dir, siteName, '_stage_1in2500_84pc.png')
    png(stage_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_stage_1in2500_84pc', 
        zlim=c(0,3), colpal=palDepthStage_trans, atws_zone = atws_zone)
    title(paste0(siteName, ': Max Stage (m) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot max-speed
    speed_file = paste0(prep_dir, siteName, '_speed_1in2500_84pc.png')
    png(speed_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_speed_1in2500_84pc', 
        zlim=c(0,5), colpal=palSpeed_trans, atws_zone = atws_zone)
    title(paste0(siteName, ': Max Speed (m/s) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot max-flux
    flux_file = paste0(prep_dir, siteName, '_flux_1in2500_84pc.png')
    png(flux_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_flux_1in2500_84pc', 
        zlim=c(0,35), colpal=palFlux_trans, atws_zone = atws_zone)
    title(paste0(siteName, ': Max Flux (m.m/s) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Plot min-stage
    min_stage_file = paste0(prep_dir, siteName, '_min_stage_1in2500_84pc.png')
    png(min_stage_file, 
        width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='min_stage_1in2500_84pc',
        zlim=c(-3,0), colpal=palDepthStage_trans_rev, atws_zone = atws_zone)
    title(paste0(siteName, ': Min Stage (m) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # # Speed hazard categories
    speed_file_cat = paste0(prep_dir, siteName, '_speedcategories_1in2500_84pc.png')
    png(speed_file_cat, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='max_speed_1in2500_84pc', 
        zlim=c(0,6), colpal=palSpeedHazard_trans, atws_zone = atws_zone)
    title(paste0(siteName, ': Current hazard (m/s) \n 1/2500, 84th percentile'), cex.main=1.9)
    dev.off()

    # Combine marine hazard: speed categories, stage and min-stage
    images_list = image_read(c(speed_file_cat, stage_file, min_stage_file))
    iminfo = image_info(images_list[1])
    combined_figs = image_montage(images_list, tile="3x1", geometry=paste0(iminfo$width, 'x', iminfo$height))
    marine_dir = paste0(out_dir, 'marine_hazard/')
    dir.create(marine_dir, showWarnings=FALSE)
    image_write(combined_figs, path=paste0(marine_dir, siteName, '_marine_hazard_2500_84pc.png'), format='png')

    # Combine onshore hazard: depth, speed, flux
    images_list = image_read(c(depth_file, speed_file, flux_file))
    iminfo = image_info(images_list[1])
    combined_figs = image_montage(images_list, tile="3x1", geometry=paste0(iminfo$width, 'x', iminfo$height))
    onshore_dir = paste0(out_dir, 'onshore_hazard/')
    dir.create(onshore_dir, showWarnings=FALSE)
    image_write(combined_figs, path=paste0(onshore_dir, siteName, '_onshore_hazard_2500_84pc.png'), format='png')

    #
    # Plot model based warning zones
    #
    jatwc_dir = paste0(out_dir, 'jatwc/')
    dir.create(jatwc_dir, showWarnings=FALSE)
    png(
        paste0(jatwc_dir, siteName, '_JATWC_zones.png'),
        width=P_W,
        height=P_H,
        units='in',
        res=P_RES
    )
    regional_plot(lower_left = plot_ll, upper_right=plot_ur, plotvar='JATWC_inundation_zones', 
        zlim=c(1,5), colpal=palDepthStage_trans, atws_zone = atws_zone, backdrop='osm')
    title(paste0(siteName, ': JATWC Inundation Zones'), cex.main=1.6)
    dev.off()
    #
    # Plot inundation exceedance-rates, logic tree mean
    #
    inundation_ltm_file = paste0(prep_dir, siteName, '_inundation_rate_ltm.png')
    png(inundation_ltm_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, backdrop='osm', plotvar='inundation_rate_ltm', 
        zlim=c(NA, NA), colpal=NA, atws_zone = atws_zone)
    title(paste0(siteName, ': Inundation Rate \n Logic Tree Mean'), cex.main=1.9)
    dev.off()

    #
    # Plot inundation exceedance-rates, 84th percentile
    #
    inundation_84pc_file = paste0(prep_dir, siteName, '_inundation_rate_84pc.png')
    png(inundation_84pc_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, backdrop='osm', plotvar='inundation_rate_84pc', 
        zlim=c(NA, NA), colpal=NA, atws_zone = atws_zone)
    title(paste0(siteName, ': Inundation Rate \n 84th Percentile'), cex.main=1.9)
    dev.off()

    #
    # Plot inundation exceedance-rates, 16th percentile
    #
    inundation_16pc_file = paste0(prep_dir, siteName, '_inundation_rate_16pc.png')
    png(inundation_16pc_file, width=P_W, height=P_H, units='in', res=P_RES)
    regional_plot(lower_left=plot_ll, upper_right=plot_ur, backdrop='osm', plotvar='inundation_rate_16pc', 
        zlim=c(NA, NA), colpal=NA, atws_zone = atws_zone)
    title(paste0(siteName, ': Inundation Rate \n 16th Percentile'), cex.main=1.9)
    dev.off()

    # Legend for inundation rate
    inundation_rate_legend_file = paste0(prep_dir, siteName, '_inundation_rate_legend.png')
    png(inundation_rate_legend_file, width=P_W, height=P_H, units='in', res=P_RES)
    plot(c(0,1), c(0, 1), col='white', axes=FALSE, frame.plot=FALSE, ann=FALSE) 
    k = 1:5 # Visible colour indices #which(!is.na(inundation_classmat_col))
    legend(
        'topright',
        legend=inundation_classmat_labels[k],
        fill=inundation_classmat_col[k],
        title='Inundation Rate (events/year)', 
        cex=1,
        bg='white',
        border='black',
        inset=0.005
        )
    dev.off()

    # Combine into a single plot
    three_images = image_read(c(inundation_16pc_file, inundation_ltm_file, inundation_84pc_file))
    # legend_image = image_read(inundation_rate_legend_file)
    iminfo = image_info(three_images[1])
    combined_figs = image_montage(three_images, tile="3x1", geometry=paste0(iminfo$width, 'x', iminfo$height))
    # overlay the legend on top of the 16th percentile in the bottom left. Make it smaller as an inset
    # legend_image = image_scale(legend_image, '40%')
    # combined_figs = image_composite(combined_figs, legend_image, gravity = 'southwest', offset = '+100+100')
    inundation_dir = paste0(out_dir, 'inundation/')
    dir.create(inundation_dir, showWarnings=FALSE)
    image_write(combined_figs, path=paste0(inundation_dir, siteName, '_inundation_16_LTM_84.png'), format='png')

}
my_inds = splitIndices(length(ZOOM_PLOTS_REGIONS), batch_count)[[batch_number]]
out_dir = './figures/'
dir.create(out_dir, showWarnings=FALSE)
prep_dir = './figures/prep/'
dir.create(prep_dir, showWarnings=FALSE)
for(i in my_inds) {
    print(ZOOM_PLOTS_REGIONS$Site[i])
    make_region_plot(i, out_dir=out_dir, prep_dir=prep_dir)
}
