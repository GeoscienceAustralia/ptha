#
# Make a shapefile depicting the multidomain structure. Run as
#    Rscript make_domains_shapefile.R path_to_multidomain_directory
#

file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

md_dir = commandArgs(trailingOnly=TRUE)[1]
if(!file.exists(md_dir)) stop(paste0('Could not find file ', md_dir))

all_domains = get_domain_interior_bbox_in_multidomain(md_dir, include_SpatialPolygonsDataFrame=TRUE)
# Write out a SpatialPolygonsDataFrame containing the directory names
ad2 = all_domains$all_domain_spdf
ad2@data = cbind(ad2@data, data.frame(folder=basename(all_domains$domain_folders)))
# Sort so that lower-res domains are on top
k = rev(order(ad2@data$dx))
ad2 = ad2[k,]          

library(sf)
outdir = paste0(dirname(all_domains$domain_folders[1]), '/domains_shapefile/')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
st_write(st_as_sf(ad2), dsn=outdir, layer='domains_shapefile', driver='ESRI Shapefile', overwrite=TRUE)
