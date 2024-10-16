#
# Make rasters in the multidomain directory. 
#
# Run as
#    Rscript make_rasters.R multidomain_directory_glob NUMBER_OF_CORES_TO_USE variables_to_make_rasters_for
# where
#    1. Sys.glob(multidomain_directory_glob) gives the multidomain directories to work on.
#    2. NUMBER_OF_CORES_TO_USE should be an integer. Use smaller counts if running out of memory.
#    3. variables_to_make_rasters_for is one or more variables from SWALS's nontemporal_grids_to_store, such as
#           elevation0 
#       or 
#           max_stage max_depth elevation0 arrival_time max_speed max_flux time_of_max_stage
#       and we also support these variables if UH/VH time grids were saved
#           last_timestep_UH last_timestep_VH
#

# Get the SWALS plot scripts
file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

parallel_fun<-function(i, multidomain_dir, raster_variable){

    print(paste0('Making ', raster_variable, ' raster for domain ', i))

    if(raster_variable %in% c('max_stage', 'max_depth')){
        # These also require elevation (to set max_stage to NA in dry areas, or
        # to compute depth)

        # Max-stage
        ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage', 
            domain_index=i, return_raster=TRUE)
        # Elevation
        elev = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0', 
            domain_index=i, return_raster=TRUE)
        output_file = paste0(multidomain_dir, '/elevation0_domain_', i, '.tif')

        # Masked dry areas
        ms[ms < elev + 1.0e-03] = NA

        if(raster_variable == 'max_stage'){
            output_file = paste0(multidomain_dir, '/max_stage_domain_', i, '.tif')
            writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    rm(ms, elev); gc()
        }else if(raster_variable == 'max_depth'){
            depth = ms - elev
            output_file = paste0(multidomain_dir, '/depth_as_max_stage_minus_elevation0_domain_', i, '.tif')
            writeRaster(depth, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
            rm(ms, elev, depth); gc()
        }

    }else if(raster_variable %in% c('last_timestep_UH', 'last_timestep_VH')){
        # Export UH or VH at the last time.
        # Idea is that we might notice nesting artefacts
        all_times = get_multidomain_output_times(multidomain_dir)
        N = length(all_times)
        tm = round(all_times[N])

        if(raster_variable == 'last_timestep_UH'){
            output_file = paste0(multidomain_dir, '/UH_time_', round(tm), 's_domain_', i, '.tif')
            ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='uh',
                desired_time_index=N, domain_index=i, return_raster=TRUE)
        }else if(raster_variable == 'last_timestep_VH'){
            output_file = paste0(multidomain_dir, '/VH_time_', round(tm), 's_domain_', i, '.tif')
            ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='vh',
                desired_time_index=N, domain_index=i, return_raster=TRUE)
        }

        writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(ms)
        gc()

    }else{
        # Standard variables

        desired_var = raster_variable
        gridded_var = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var=desired_var, 
            domain_index=i, return_raster=TRUE)
        output_file = paste0(multidomain_dir, '/', desired_var, '_domain_', i, '.tif')
        writeRaster(gridded_var, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(gridded_var); gc()

    }

    return(0)
}

input_args = commandArgs(trailingOnly=TRUE)
stopifnot(length(input_args) >= 3)
# First argument is the multidomain directory, or something that globs a few of them.
MULTIDOMAIN_DIRS = Sys.glob(input_args[1])
# Second argument is the integer number of cores
MC_CORES = as.numeric(input_args[2]); stopifnot(is.finite(MC_CORES)); stopifnot(MC_CORES >= 1)
# Third and later arguments give the rasters to create
RASTER_VARIABLES = input_args[-c(1,2)]

try_parallel_fun <- function(i, multidomain_dir, raster_variable){
    try(parallel_fun(i, multidomain_dir, raster_variable))
}

all_domain_inds = get_domain_indices_in_multidomain(MULTIDOMAIN_DIRS[1])

all_inputs = expand.grid(all_domain_inds, MULTIDOMAIN_DIRS, RASTER_VARIABLES, stringsAsFactors=FALSE)

library(parallel)
cl = makeForkCluster(MC_CORES)
clusterExport(cl, varlist=ls(all=TRUE))
clusterMap(cl, try_parallel_fun, i=all_inputs[,1], multidomain_dir=all_inputs[,2], raster_variable=all_inputs[,3], 
    .scheduling='dynamic')
stopCluster(cl)

