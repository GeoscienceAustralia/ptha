source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')

multidomain_dir = commandArgs(trailingOnly=TRUE)[1]
all_domain_inds = get_domain_indices_in_multidomain(multidomain_dir)

for(i in all_domain_inds){
    print(c('Making raster for domain ', i))
    # Max-stage
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='max_stage', 
				domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'max_stage_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    rm(r1)
    # Elevation
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0', 
				domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'elevation0_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    rm(r1)
}


