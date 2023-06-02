#
# Make rasters in the multidomain directory. Run as
#    Rscript make_rasters.R path_to_multidomain_directory
#

# Get the SWALS plot scripts [works on NCI or GD home machine]
file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

parallel_fun<-function(i, multidomain_dir){
    print(c('Making raster for domain ', i))
    # Max-stage
    ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage', 
        domain_index=i, return_raster=TRUE)
    # Elevation
    elev = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0', 
        domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'elevation0_domain_', i, '.tif')
    writeRaster(elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    ## Difference between max-stage and elevation
    #max_stg_less_elev = ms - elev
    #max_stg_less_elev[max_stg_less_elev < 1.0e-03] = NA
    #output_file = paste0(multidomain_dir, 'depth_as_max_stage_minus_elevation0_domain_', i, '.tif')
    #writeRaster(max_stg_less_elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    #rm(max_stg_less_elev)
    #gc()

    ## Max stage, masked in dry areas
    ms[ms < elev + 1.0e-03] = NA
    output_file = paste0(multidomain_dir, 'max_stage_domain_', i, '.tif')
    writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(ms, elev)
    gc()

    for(desired_var in c('arrival_time', 'max_speed', 'max_flux')){
        ## Arrival time
        gridded_var = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var=desired_var, 
            domain_index=i, return_raster=TRUE)
        output_file = paste0(multidomain_dir, desired_var, '_domain_', i, '.tif')
        writeRaster(gridded_var, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        rm(gridded_var); gc()
    }

    ## Export UH at the last time.
    ## Idea is that we might notice nesting artefacts
    #all_times = get_multidomain_output_times(multidomain_dir)
    #N = length(all_times)
    #tm = round(all_times[N])
    #output_file = paste0(multidomain_dir, 'UH_time_', round(tm), 's_domain_', i, '.tif')
    #ms = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='uh',
    #    desired_time_index=N, domain_index=i, return_raster=TRUE)
    #writeRaster(ms, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    #rm(ms)
    #gc()
}

try_parallel_fun<-function(i, multidomain_dir){
    try(parallel_fun(i, multidomain_dir))
}

#multidomain_dirs = commandArgs(trailingOnly=TRUE)
#multidomain_dirs = paste0(Sys.glob('OUTPUTS/VAUS*/RUN*'), '/')

multidomain_dirs = c(
    "OUTPUTS/Fuji_andaman2004_24hrs_domain301122-test_load_balance-ambient_sea_level_0.0/RUN_20221216_145827463/"
)

all_domain_inds = get_domain_indices_in_multidomain(multidomain_dirs[1])

all_inputs = expand.grid(all_domain_inds, multidomain_dirs, stringsAsFactors=FALSE)

library(parallel)
mcmapply(try_parallel_fun, 
    i=all_inputs[,1], 
    multidomain_dir=all_inputs[,2],
    SIMPLIFY=FALSE,
    mc.cores=16,
    mc.preschedule=TRUE)

