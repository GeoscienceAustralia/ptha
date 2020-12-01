# Get the SWALS plot scripts [works on NCI or GD home machine]
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
file_home = '~/Code_Experiments/fortran/Structured_shallow_water/plot.R'
if(file.exists(file_nci)){
    source(file_nci)
}else{
    source(file_home)
}


# Create the max-depth
make_max_depth_rasters_on_multidomain<-function(multidomain_dir){

    all_domain_inds = get_domain_indices_in_multidomain(multidomain_dir)

    for(i in all_domain_inds){

        # Max-stage
        r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage',
                    domain_index=i, return_raster=TRUE)
        #output_file = paste0(multidomain_dir, '/max_stage_domain_', i, '.tif')
        #writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
        ms = r1
        rm(r1)

        # Elevation
        r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0',
                    domain_index=i, return_raster=TRUE)
        #output_file = paste0(multidomain_dir, '/elevation0_domain_', i, '.tif')
        #writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

        # Difference between max-stage and elevation
        max_stg_less_elev = ms - r1
        max_stg_less_elev[max_stg_less_elev < 1.0e-03] = NA
        output_file = paste0(multidomain_dir, '/depth_as_max_stage_minus_elevation0_domain_', i, '.tif')
        writeRaster(max_stg_less_elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }
    rm(r1, ms, max_stg_less_elev)
    gc()

    return(multidomain_dir)
}


# Main program
#all_md_dirs = Sys.glob('OUTPUTS/ptha18_tonga_MSL0/ptha*/RUN*')
#all_md_dirs = Sys.glob('OUTPUTS/ptha18_tonga_MSL0.8/ptha*/RUN*')
#all_md_dirs = Sys.glob('OUTPUTS/ptha18_tonga_MSL0_meshrefine2/ptha*/RUN*')
all_md_dirs = Sys.glob('OUTPUTS/ptha18_tonga_MSL*/ptha*/RUN*')
MC_CORES=48
library(parallel)

# Protect parallel run against occasional failures
try_make_max_depth_rasters_on_multidomain<-function(x) try(make_max_depth_rasters_on_multidomain(x))

mclapply(all_md_dirs, try_make_max_depth_rasters_on_multidomain, mc.cores=MC_CORES)
