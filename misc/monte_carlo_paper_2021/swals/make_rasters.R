#
# Make rasters in the multidomain directory. Run as
#    Rscript make_rasters.R path_to_multidomain_directory
#

# Get the SWALS plot scripts [works on NCI or GD home machine]
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
file_home = '~/Code_Experiments/fortran/Structured_shallow_water/plot.R'
if(file.exists(file_nci)){
    source(file_nci)
}else{
    source(file_home)
}

multidomain_dir = commandArgs(trailingOnly=TRUE)[1]
all_domain_inds = get_domain_indices_in_multidomain(multidomain_dir)
all_grid_times = get_multidomain_output_times(multidomain_dir)

for(i in all_domain_inds){
    print(c('Making raster for domain ', i))
    # Max-stage
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage', 
				domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'max_stage_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    ms = r1
    rm(r1)

    # Elevation
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0', 
				domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'elevation0_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    # Difference between max-stage and elevation
    max_stg_less_elev = ms - r1
    max_stg_less_elev[max_stg_less_elev < 1.0e-03] = NA
    output_file = paste0(multidomain_dir, 'depth_as_max_stage_minus_elevation0_domain_', i, '.tif')
    writeRaster(max_stg_less_elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    ## As above subtracting 0.7
    #max_stg_less_elev = ms - r1 - 0.7
    #max_stg_less_elev[max_stg_less_elev < 1.0e-03] = NA
    #output_file = paste0(multidomain_dir, 'depth_as_max_stage_minus_elevation0_minus70cm_domain_', i, '.tif')
    #writeRaster(max_stg_less_elev, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(r1, ms, max_stg_less_elev)

    # Initial stage
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='stage', 
                desired_time_index=1, domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'initial_stage_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    rm(r1)

    # Final stage
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='stage', 
                desired_time_index=length(all_grid_times), domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'final_stage_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    rm(r1)
}


