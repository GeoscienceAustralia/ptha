# Suppose we have computed rasters with exceedance-rates for max-stage in a set of thresholds, such as
# 1.101, 2.1, 3.1, 4.1, ... 10.1
# For any specified exceedance-rate ER, we can compute an "approximate" corresponding max-stage by searching
# the above rasters for the largest max-stage with exceedance-rate exceeding ER. 
#
# This script does that (and similarly for other variables and sets of thresholds)
# Run with
#   Rscript compute_binned_thresholds_matching_exrate_from_set_of_exrate_rasters.R directory/with/rasters variable exceedance_rate
# e.g. (for 1/500 max-stage)
#   Rscript compute_binned_thresholds_matching_exrate_from_set_of_exrate_rasters.R ptha18-GreaterPerth2023-sealevel60cm/highres_max_stage_with_variance/ptha18-GreaterPerth2023-sealevel60cm-max_stage-LogicTreeMean-sum_of_source_zones max_stage 0.002
#
# The result is a raster with the highest threshold (in the input set) with
# exceedance_rate less than the input exceedance-rate (0.002 in the example above).
#
library(terra)
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)

# INPUTS
folder_with_rasters_for_multiple_thresholds = commandArgs(trailingOnly=TRUE)[1]
variable_of_interest = commandArgs(trailingOnly=TRUE)[2]
threshold_exrate = as.numeric(commandArgs(trailingOnly=TRUE)[3])

# Find files containing the exceedance rates (not the Monte Carlo variances or the upper 95% CI)
files_to_search = Sys.glob(paste0(folder_with_rasters_for_multiple_thresholds, '/*',
    variable_of_interest, '_domain__exceedance_rate_with_threshold_*.tif'))
stopifnot(length(files_to_search) > 1)

# Extract the 'threshold' from each file
file_thresholds = unlist(lapply(
    strsplit(basename(files_to_search), split="_threshold_"), 
    function(x) gsub('.tif', '', x[2], fixed=TRUE)))
unique_thresholds_sorted = sort(as.numeric(unique(file_thresholds)))

# Extract the domains_char from each file
file_domains_char = unlist(lapply(
    strsplit(gsub('summed_HS', '', basename(files_to_search)), 
             split=paste0(variable_of_interest, '_')),
    function(x) x[1]))
domains_char = unique(file_domains_char)

# The files-to-search should have every domain at every threshold, and nothing else.
stopifnot(length(domains_char)*length(unique_thresholds_sorted) == length(files_to_search))
stopifnot(all(table(file_domains_char) == length(unique_thresholds_sorted)))

output_dir = paste0(dirname(folder_with_rasters_for_multiple_thresholds), 
    '/binned_', variable_of_interest, '_that_has_exrate_greater_than_', threshold_exrate)
dir.create(output_dir, showWarnings=FALSE)

count_files_exceeding_threshold_exrate<-function(domain_ind){

    domain_char = domains_char[domain_ind]

    matching_files = files_to_search[ grep(domain_char, files_to_search, fixed=TRUE) ]

    stopifnot(length(matching_files) == length(unique_thresholds_sorted))

    # Count the number of files having exceedance-rate above the threshold
    count_above_exrate = 1*(rast(matching_files[1]) > threshold_exrate)
    count_above_exrate[is.na(count_above_exrate)] = 0 # NA values mean "zero exrate" (or halos)
    for(i in 2:length(matching_files)){
        tmp = 1*( rast(matching_files[i]) > threshold_exrate )
        tmp[is.na(tmp)] = 0 # NA values mean "zero exrate" (or halos)
        count_above_exrate = count_above_exrate + tmp
    }

    # Sites that score zero in all cases are likely to be halos or regions that are 
    # never wet, and we want them to be NA so that the plotting works OK.
    count_above_exrate[count_above_exrate == 0] = NA
    threshold_below_exrate = count_above_exrate * 0
    for(i in 1:length(unique_thresholds_sorted)){
        ti = unique_thresholds_sorted[i]
        threshold_below_exrate[count_above_exrate == i] = ti
    }

    output_file = paste0(output_dir, '/binned_', variable_of_interest, domain_char, 
        'that_has_exrate_greater_than_', threshold_exrate, '.tif')

    writeRaster(threshold_below_exrate, output_file, gdal=c('COMPRESS=DEFLATE'), overwrite=TRUE)

    rm(output)
    gc()
    return(1)
}

# Protect parallel runs against failure of a single case by wrapping in 'try'
try_count_files_exceeding_threshold_exrate<-function(i){
    try(count_files_exceeding_threshold_exrate(i))
}

library(parallel)
result = mclapply(1:length(domains_char), try_count_files_exceeding_threshold_exrate, mc.cores=asfm$DEFAULT_MC_CORES)
