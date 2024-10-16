#
# Create rasters containing "Merged models" with the 84th percentile inundation rate
#

library(sf)
library(stars)
library(terra)

# Polygons showing where we should make rasters, and their resolution
raster_regions = st_read('merged_model_raster_output_zones/merged_model_raster_output_zones.shp')

# Rasters we extract results from
model_rasters = list(
    'GreaterPerth' = rast(Sys.glob('../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/inundation_rate_84pc/*.vrt')),
    'BunburyBusselton' = rast(Sys.glob('../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/bunbury_busselton/inundation_rate_84pc/*.vrt')),
    'Midwest' = rast(Sys.glob('../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/inundation_rate_84pc/*.vrt')),
    'BunburyBusseltonAlternate' = rast(Sys.glob('../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/bunbury_busselton_alternative_with_bunbury_floodgate_shut/inundation_rate_84pc/*.vrt'))
    )

# Where we write the rasters
output_dir = 'merged_inundation_rate_84pc'
# Where we write the contours
output_dir_contour = 'merged_inundation_rate_contour_1in2500_1in500_1in100_84pc'


# Determine the preferred raster from which to extract outputs
make_rasters_to_use<-function(){
    rasters_to_use = rep("", nrow(raster_regions))

    # - Contained in the bunker bay polygon --> Greater Perth
    bunkerBay_poly = st_read('priority_model_zones/bunkerBay/bunkerBay.shp')
    k = st_contains_properly(bunkerBay_poly, raster_regions)[[1]]
    rasters_to_use[k] = 'GreaterPerth'

    # - Contained in a region toward the north of perth that will otherwise be in the midwest model -- > Greater Perth
    north_perth_poly = st_read('priority_model_zones/Perth/Perth.shp')
    k = st_contains_properly(north_perth_poly, raster_regions)[[1]]
    rasters_to_use[k] = "GreaterPerth"

    # - Contained in a region near bunbury --> Max of the 2 Bunbury models
    bunbury_poly = st_read('priority_model_zones/bunbury_2models/bunbury_2models.shp')
    k = st_contains_properly(bunbury_poly, raster_regions)[[1]]
    rasters_to_use[k] = "max_Bunbury_or_BunburyAlternate"

    # - Contained in a region near Geraldton --> Max of Greater Perth or Midwest. This was requested by Adrian Brannigan because
    # the midwest model contained buildings in the topography, which no-where else did, so he preferred to use a conservative maximum.
    Geraldton_poly = st_read('priority_model_zones/midwestOrPerth/midwestOrPerth.shp')
    k = st_contains_properly(Geraldton_poly, raster_regions)[[1]]
    rasters_to_use[k] = "max_Midwest_or_GreaterPerth"

    # Deal with remaining rasters
    raster_centroids = st_coordinates(st_centroid(raster_regions))

    # - Lat < -33 --> BunburyBusselton
    k = which(rasters_to_use == "" & raster_centroids[,2] < -33)
    rasters_to_use[k] = 'BunburyBusselton'

    # - Lat > -33 & Lat < -31.5 --> GreaterPerth
    k = which(rasters_to_use == "" & raster_centroids[,2] < -31.5)
    rasters_to_use[k] = "GreaterPerth"

    # - Lat > -31.5 --> Midwest
    k = which(rasters_to_use == "")
    rasters_to_use[k] = "Midwest"

    return(rasters_to_use)
}
rasters_to_use = make_rasters_to_use()

# Populate the new raster with cell values from one or more original rasters.
resample_to_new_raster<-function(raster_poly, raster_to_use, output_raster_name){

    # Make a raster
    dx = raster_poly$dx
    rast_ext = ext(raster_poly)
    ncols = round((rast_ext$xmax - rast_ext$xmin)/dx)
    nrows = round((rast_ext$ymax - rast_ext$ymin)/dx) 
    r1 = rast(rast_ext, nrows=nrows, ncols=ncols, crs='EPSG:4326')

    # Check the resolution matches intention (prevent mistake with nrows/ncols etc)
    stopifnot(all( abs(res(r1)-dx) < (1e-04*dx)))

    # Populate the raster cells
    if(raster_to_use == 'GreaterPerth'){
        r1 = resample(model_rasters$GreaterPerth, r1, method='near')
    }else if(raster_to_use == 'BunburyBusselton'){
        r1 = resample(model_rasters$BunburyBusselton, r1, method='near')
    }else if(raster_to_use == 'Midwest'){
        r1 = resample(model_rasters$Midwest, r1, method='near')
    }else if(raster_to_use == 'max_Midwest_or_GreaterPerth'){
        t0 = resample(model_rasters$Midwest, r1, method='near')
        t1 = resample(model_rasters$GreaterPerth, r1, method='near')
        t0[is.na(t0)] = 0
        t1[is.na(t1)] = 0
        r1 = max(t0, t1) 
        r1[r1 == 0] = NA
        rm(t0, t1); gc()
    }else if(raster_to_use == 'max_Bunbury_or_BunburyAlternate'){
        t0 = resample(model_rasters$BunburyBusselton, r1, method='near')
        t1 = resample(model_rasters$BunburyBusseltonAlternate, r1, method='near')
        t0[is.na(t0)] = 0
        t1[is.na(t1)] = 0
        r1 = max(t0, t1) 
        r1[r1 == 0] = NA
        rm(t0, t1); gc()
    }else{
        stop(paste0('Unknown value of raster_to_use: ', raster_to_use))
    }

    # Mask the output at domains finer than the current domain
    k = which(raster_regions$dx < dx*0.99)
    if(length(k) > 0){
        r1 = mask(r1, raster_regions[k,], touches=FALSE, inverse=TRUE)
    }

    writeRaster(r1, file=output_raster_name, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    return(invisible(0))
}

# To create a single raster tile 
parallel_fun<-function(i){
   try(resample_to_new_raster(raster_regions[i,], rasters_to_use[i], paste0(output_dir, '/', raster_regions$outfile[i])))
}

dir.create(output_dir, showWarnings=FALSE)
# Make rasters in parallel
library(parallel)
parallel_run = mclapply(1:nrow(raster_regions), parallel_fun, mc.cores=16, mc.preschedule=FALSE)

# Create a VRT to view them all at once
mydir = getwd()
setwd(output_dir)
system('gdalbuildvrt -resolution highest merged_models_inundation_rate_84th_percentile.vrt *.tif')
setwd(mydir)

#
# Now make a contour
#

make_contour<-function(raster_file, rate=c(1/2500,1/500, 1/100)){

    r1 = read_stars(raster_file)
    # For some reason I can't do all the contour levels at once. Loop instead
    for(i in 1:length(rate)){
        mycont = st_contour(r1, contour_lines=TRUE, breaks=rate[i])
        names(mycont)[1] = 'rate'
        if(i == 1){ 
            contours = mycont
        }else{
            contours = rbind(contours, mycont) 
        }
    }         

    return(contours)
}
# Make contours in parallel. Contours will not be joined at raster boundaries,
# but it probably doesn't matter
parallel_run = mclapply(Sys.glob(paste0(output_dir, '/*.tif')), make_contour, 
    mc.cores=16, mc.preschedule=FALSE)
all_contours = do.call(rbind, parallel_run)

st_write(all_contours, dsn=output_dir_contour, 
    layer=paste0(basename(output_dir_contour), '.shp'),
    driver='ESRI Shapefile', 
    append=FALSE)
