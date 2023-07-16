#
# We have computed rasters with exceedance-rates for max-stage in
# 0.6 + c(0.001, 1:10)
#
# In this script, for each domain, we count the rasters where that exceedance-rate
# is greater than some threshold. This lets us approximate the stage with the given exceedance-rate
#
library(raster)

# This matches each unique domain in the files_to_search
domains_char = paste0('_domain_', seq(1, 526), '_')
files_to_search = Sys.glob('summed_results*/*.tif')


# For each domain, we have exceedance-rates for the following max-stages
# Note they ALMOST differ by 1 -- FIXME: the code logic below assumes this
max_stage_values = 0.6 + c(0.001, 1:10)

# We will estimate the max-stage with (exceedance-rate > threshold_exrate)
threshold_exrate = 1/100 #1/500 #1/10000 #1/2500

output_dir = paste0('binned_max_stage_exceeding_', threshold_exrate)
dir.create(output_dir, showWarnings=FALSE)

count_files_exceeding_threshold_exrate<-function(i){

    domain_char = domains_char[i]

    matching_files = files_to_search[ grep(domain_char, files_to_search, fixed=TRUE) ]

    stopifnot(length(matching_files) == 11)

    # Count the number of files having exceedance-rate above the threshold
    output = 1*( raster(matching_files[1]) > threshold_exrate )
    for(i in 2:length(matching_files)){
        output = output + 1*( raster(matching_files[i]) > threshold_exrate )
    }

    # Sites that score zero in all cases are likely to be halos or regions that are 
    # never wet, and we want them to be NA so that the plotting works OK.
    output[output == 0] = NA

    # To change the 'raster counts' to the 'tsunami maxima'.
    # In this case we just convert 1, 2, 3, 11 --> 0, 1, 2, 3, ... 10
    # FIXME: This logic is fragile if other max-stage values are used.
    # For now just force a failure in that case.
    stopifnot(all(abs((max_stage_values - 0.6) - seq(0, 10)) < 2.0e-03))
    output = output - 1.0

    output_file = paste0(output_dir, '/binned_max_stage', domain_char, '_exceeding_', threshold_exrate, '.tif')

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
