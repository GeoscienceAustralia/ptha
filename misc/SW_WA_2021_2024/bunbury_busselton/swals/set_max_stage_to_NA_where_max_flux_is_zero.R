#
# In areas where the elevation is changing over time, the current max-stage calculations 
# can lead to ground areas appearing like they are wet (if the ground subsides).
#
# This is not right, and can affect our large scale plots (near the earthquake source).
#
# A cheap fix is to set the value to NA where max_flux == 0
#

#max_stage_rasts = Sys.glob('OUTPUTS/Fuji_andaman2004_24hrs_domain301122-full-ambient_sea_level_0.0/RUN_20230309_174201139/max_stage*.tif') 
max_stage_rasts = Sys.glob('OUTPUTS/Fuji_sumatra2005_24hrs_domain301122-full-ambient_sea_level_0.0/RUN_20230308_144745735/max_stage*.tif')

max_flux_rasts = gsub('max_stage', 'max_flux', max_stage_rasts)

library(terra)

for(i in 1:length(max_stage_rasts)){
    print(i)
    x = rast(max_stage_rasts[i])
    y = rast(max_flux_rasts[i])

    # NA where max-flux is zero (i.e. areas that are always dry)
    x[y == 0] = NA

    new_file = gsub('max_stage', 'NA_filtered_max_stage', max_stage_rasts[i])
    writeRaster(x, new_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}
