# The computed rasters with a flow quantity at an exceedance-rate/percentile were
# calculated using root-finding between specified lower and upper limits. They can
# be confusing to interpret near the lower-bound of the search space. Depending
# on the flow variable and root-finding tolerance, a value near the lower bound
# should be missing (i.e. reflecting no flow at the given exceedance-rate/percentile). 
#
# This script creates new rasters that are easier to interpret near the lower
# bounds.
library(terra)
library(parallel)

raster_data_list = commandArgs(trailing=TRUE)[1]
if(raster_data_list == "") raster_data_list = '84pc'

if(raster_data_list == '50pc'){
    #
    # INPUTS, 50th percentile case
    #
    rasters_to_process = list(

        # Max-depth rasters
        max_depth=list(
            rasters = Sys.glob('ptha_highres_domains_depth_at_epistemic_uncertainty_50pc/*.tif'),
            NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_depth_1in250_50pc',
            output_vrt = 'ptha_max_depth_1in250_50pc.vrt'),

        # Max-stage rasters
        max_stage = list(
            rasters = Sys.glob('ptha_highres_domains_max_stage_at_epistemic_uncertainty_50pc/*.tif'),
            NA_below = 1.1015, # No lower limit required #0.6 + 0.0015, # Lower limit was 0.6. Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_stage_1in250_50pc',
            output_vrt = 'ptha_max_stage_1in250_50pc.vrt'),

        # Max-flux rasters
        max_flux = list(
            rasters = Sys.glob('ptha_highres_domains_max_flux_at_epistemic_uncertainty_50pc/*.tif'),
            NA_below = 0.015, # Compare to root-finding tolerance of 1e-02
            output_dir = './ptha_max_flux_1in250_50pc',
            output_vrt = 'ptha_max_flux_1in250_50pc.vrt'),

        # Max-speed rasters
        max_speed = list(
            rasters = Sys.glob('ptha_highres_domains_max_speed_at_epistemic_uncertainty_50pc/*.tif'),
            NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_speed_1in250_50pc',
            output_vrt = 'ptha_max_speed_1in250_50pc.vrt')

        )
}else if(raster_data_list == '84pc'){
    #
    # INPUTS, 84th percentile case
    #
    rasters_to_process = list(

        # Max-depth rasters
        max_depth=list(
            rasters = Sys.glob('ptha_highres_domains_depth_at_epistemic_uncertainty_84pc/*.tif'),
            NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_depth_1in2500_84pc',
            output_vrt = 'ptha_max_depth_1in2500_84pc.vrt'),

        # Max-stage rasters
        max_stage = list(
            rasters = Sys.glob('ptha_highres_domains_max_stage_at_epistemic_uncertainty_84pc/*.tif'),
            NA_below = 1.1015, # No lower limit required #0.6 + 0.0015, # Lower limit was 0.6. Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_stage_1in2500_84pc',
            output_vrt = 'ptha_max_stage_1in2500_84pc.vrt'),

        # Max-flux rasters
        max_flux = list(
            rasters = Sys.glob('ptha_highres_domains_max_flux_at_epistemic_uncertainty_84pc/*.tif'),
            NA_below = 0.015, # Compare to root-finding tolerance of 1e-02
            output_dir = './ptha_max_flux_1in2500_84pc',
            output_vrt = 'ptha_max_flux_1in2500_84pc.vrt'),

        # Max-speed rasters
        max_speed = list(
            rasters = Sys.glob('ptha_highres_domains_max_speed_at_epistemic_uncertainty_84pc/*.tif'),
            NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
            output_dir = './ptha_max_speed_1in2500_84pc',
            output_vrt = 'ptha_max_speed_1in2500_84pc.vrt')

        )

}else{
    stop(paste0('unknown input ', raster_data_list))
}

MC_CORES = 48

#
# END INPUTS
#
for(varname in names(rasters_to_process)){

    rasters = rasters_to_process[[varname]]$rasters
    NA_below = rasters_to_process[[varname]]$NA_below
    output_dir = rasters_to_process[[varname]]$output_dir
    output_vrt = rasters_to_process[[varname]]$output_vrt

    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)        

    threshold_raster<-function(filename, output_dir, NA_below){
        r1 = rast(filename)
        output_file = paste0(output_dir, '/', basename(filename))        
        r1[r1 < NA_below] = NA
        writeRaster(r1, file=output_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(r1); gc()
        return(output_file)
    }

    mclapply(rasters, threshold_raster, output_dir=output_dir, NA_below=NA_below, mc.cores=MC_CORES)

    # Make a vrt file
    mydir = getwd()
    setwd(output_dir)
    vrt_command = paste0('gdalbuildvrt -resolution highest ', output_vrt, ' *.tif')
    system(vrt_command)
    setwd(mydir)
}