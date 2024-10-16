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

#
# INPUTS
#
rasters_to_process = list(

    # Max-depth rasters
    max_depth=list(
        rasters = Sys.glob('midwest_revised_highres_domains_depth_at_epistemic_uncertainty_84pc/*.tif'),
        NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
        output_dir = './midwest_max_depth_1in2500_84pc',
        output_vrt = 'midwest_max_depth_1in2500_84pc.vrt'),

    # Max-stage rasters
    max_stage = list(
        rasters = Sys.glob('midwest_revised_highres_domains_max_stage_at_epistemic_uncertainty_84pc/*.tif'),
        NA_below = -9999, # No lower bound needed #0.6 + 0.0015, # Lower limit was 0.6. Compare to root-finding tolerance of 1e-03
        output_dir = './midwest_max_stage_1in2500_84pc',
        output_vrt = 'midwest_max_stage_1in2500_84pc.vrt'),

    # Max-flux rasters
    max_flux = list(
        rasters = Sys.glob('midwest_revised_highres_domains_max_flux_at_epistemic_uncertainty_84pc/*.tif'),
        NA_below = 0.015, # Compare to root-finding tolerance of 1e-02
        output_dir = './midwest_max_flux_1in2500_84pc',
        output_vrt = 'midwest_max_flux_1in2500_84pc.vrt'),

    # Max-speed rasters
    max_speed = list(
        rasters = Sys.glob('midwest_revised_highres_domains_max_speed_at_epistemic_uncertainty_84pc/*.tif'),
        NA_below = 0.0015, # Compare to root-finding tolerance of 1e-03
        output_dir = './midwest_max_speed_1in2500_84pc',
        output_vrt = 'midwest_max_speed_1in2500_84pc.vrt')

    )

MC_CORES = 48

print('WARNING: This code includes a hack to work around a problem with initial condition of the WA Midwest model at Sandy Cape, which leads to flow in a small patch of land far from the coast (that would definitely have no flow if the initial condition were set correctly). You do not want to include this in any other models!!!!')
# This model has a minor problem where the initial stage at Sandy Cape is
# set to a low value.  By mistake there were cells on land with low
# elevation (far from wet areas) on the boundary of the polygon defining
# the area with low max-stage. This led to a "dam-break" initial condition
# at a site in domain 454! Hence there is non-zero flow that seems weird.  
#
# Because the area is on land, and completely disconnected from the real
# model solutions, we can just remove these cells.
raster_needs_fixing<-function(raster_name){
    # Only domain 454. The first grep is just to reduce the risk of this hack being
    # applied by accident in future models
    (grepl('midwest_revised_highres', raster_name, fixed=TRUE) & 
    grepl("domain_index_454.tif", raster_name, fixed=TRUE))
}
raster_cells_to_NA_in_domain_454_polygon = vect(
    # Polygon containing cells in domain 454 that should be set to NA
    matrix(c(
            115.055702, -30.2864 ,
            115.05709 , -30.2864 ,
            115.05709 , -30.284936,
            115.055702, -30.284936), 
        ncol=2, byrow=TRUE),
    type='polygons',
    crs='EPSG:4326')


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

        # Workaround for problem in domain 454
        if(raster_needs_fixing(filename)){
            r1 = mask(r1, raster_cells_to_NA_in_domain_454_polygon, inverse=TRUE)
        }

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
