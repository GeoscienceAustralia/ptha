#
# Get the max-stage rasts from a Greater Perth model run
# - clip them in areas that were never wet (having zero max-flux). Nowadays we can use the SWALS option
#

# Store the adjusted tifs here
output_dir = './Sumatra2004_max_stage_greater_perth_model/'

# Source of the tifs
max_stage_rasts = Sys.glob('/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/greater_perth_revised2023/swals/OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20230831_172007573/max_stage*.tif')
max_flux_rasts = gsub('max_stage', 'max_flux', max_stage_rasts)


library(terra)

for(i in 1:length(max_stage_rasts)){
    print(i)
    x = rast(max_stage_rasts[i])
    y = rast(max_flux_rasts[i])

    # NA where max-flux is zero (i.e. areas that are always dry)
    x[y == 0] = NA

    new_file = paste0(output_dir, 
        gsub('max_stage', 'DRY_CLIPPED_max_stage', basename(max_stage_rasts[i]))
    writeRaster(x, new_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}

