#
# We have computed rasters with exceedance-rates for max-depth in
# c(0.001, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10)
#
# In this script, for each domain, we count the rasters where that exceedance-rate
# is greater than some threshold. This lets us approximate the depth with the given exceedance-rate
#
library(raster)

# This matches each unique domain in the files_to_search
domains_char = paste0('_domain_', seq(1, 559), '_')
files_to_search = Sys.glob('summed_results*/*.tif')


# For each domain, we have exceedance-rates for the following max-depth
max_depth_values = c(0.001, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10)

# We will estimate the max-depth with (exceedance-rate > threshold_exrate)
threshold_exrate = 1/10000 #1/2500

output_dir = paste0('binned_max_depth_exceeding_', threshold_exrate)
dir.create(output_dir, showWarnings=FALSE)

count_files_exceeding_threshold_exrate<-function(i){

    domain_char = domains_char[i]

    matching_files = files_to_search[ grep(domain_char, files_to_search, fixed=TRUE) ]

    stopifnot(length(matching_files) == length(max_depth_values))

    # Count the number of files having exceedance-rate above the threshold
    output = 1*( raster(matching_files[1]) > threshold_exrate )
    for(i in 2:length(matching_files)){
        output = output + 1*( raster(matching_files[i]) > threshold_exrate )
    }

    # Sites that score zero in all cases are likely to be halos or regions that are 
    # never wet, and we want them to be NA so that the plotting works OK.
    output[output == 0] = NA

    # Convert the 'number of exceedences' to the depth thresholds
    output_vals = getValues(output)
    output = setValues(output, max_depth_values[output_vals])

    output_file = paste0(output_dir, '/binned_max_depth', domain_char, '_exceeding_', threshold_exrate, '.tif')

    writeRaster(output, output_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(output)
    gc()
    return(1)
}

# Protect parallel runs against failure of a single case by wrapping in 'try'
try_count_files_exceeding_threshold_exrate<-function(i){
    try(count_files_exceeding_threshold_exrate(i))
}

library(parallel)
mclapply(1:length(domains_char), try_count_files_exceeding_threshold_exrate, mc.cores=12)
