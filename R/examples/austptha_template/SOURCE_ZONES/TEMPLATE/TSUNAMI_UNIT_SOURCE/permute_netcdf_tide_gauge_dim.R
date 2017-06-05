#
# Code to permute the dimensions of all tide-gauge netcdf unit-source files
# which have been run for a source-zone
#
# This changes the files to support fast access to time-series at single
# stations -- so can greatly speed up later analyses.
library(parallel)

#
# INPUT DATA
#

# Name of source. In this template, that should be the name of the directory
# above the working directory
source_name = basename(dirname(getwd()))

# Vector of all netcdf files that we wish to permute
tide_gauge_files = Sys.glob(
    paste0('/g/data/w85/tsunami/AustPTHA/version1/unit_sources/', 
    source_name, '/unit_source_tsunami/*/*/Gauge*.nc')) 

# Number of cores to use in parallel (shared memory only)
mc_cores = 16

#
# END INPUT DATA
#


#' Permute dimensions of netcdf tidal gauge file
#'
#' R interface to ncpdq for permuting tide gauge file dimensions.
#' This needs the nco module to be loaded. Permuting dimensions is
#' useful to enable high-speed access to single station time-series.
#'
#' @param infile Full path to the tide gauge netcdf file
#' @param outfile If not NULL, the full path where the permuted tide gauge file
#' will be written. If NULL, then the infile will be overwritten.
#' @return Nothing
permute_netcdf_tide_gauge_dimensions<-function(infile, outfile=NULL, verbose=TRUE){

    # Only overwrite infile if outfile is NULL
    overwrite_file = FALSE
    if(is.null(outfile)){
        outfile = gsub('.nc', '_permuted_dim.nc', infile)
        overwrite_file = TRUE
    }

    permute_command = paste0('ncpdq -a station,time ', infile, ' ', outfile)
    print(permute_command)
    # The system command will return 0 if it works, and 1 if it fails
    call_failed = system(permute_command)
 
    if(call_failed) stop()
    
    if(overwrite_file){ 
        if(verbose) print('... Overwriting infile')
        file.rename(outfile, infile)
    }

    return(invisible())
}

# Main computation here
permute_all_files = mclapply(as.list(tide_gauge_files), 
    permute_netcdf_tide_gauge_dimensions, 
    mc.cores = mc_cores, mc.preschedule=TRUE)

