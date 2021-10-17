#
# Make a single set of elevation rasters for domains 3-7 that are
# not perturbed by any earthquake
#

# Get the SWALS plot scripts [works on NCI or GD home machine]
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
file_home = '~/Code_Experiments/fortran/Structured_shallow_water/plot.R'
if(file.exists(file_nci)){
    source(file_nci)
}else{
    source(file_home)
}

# Use the Chile2015 model to make grids for domains 3-7.
# The earthquake is far-field so this elevation is unperturbed
multidomain_dir = Sys.glob('../../swals/OUTPUTS/Chile2015*/RUN*')[1]
all_domain_inds = 3:7

dir.create('initial_elevation_grids/', showWarnings=FALSE)

for(i in all_domain_inds){
    print(c('Making raster for domain ', i))

    # Elevation
    r1 = merge_domains_nc_grids(multidomain_dir = multidomain_dir, desired_var='elevation0',
                domain_index=i, return_raster=TRUE)
    output_file = paste0(multidomain_dir, 'elevation0_domain_', i, '.tif')
    writeRaster(r1, file=output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

}
