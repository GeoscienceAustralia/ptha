check_for_failed_rasters<-function(expected_files){

    library(terra)
    read_files = lapply(expected_files, function(x) try(rast(x)))
    is_try_error = unlist(lapply(read_files, function(x) is(x, 'try-error')))
    k = which(is_try_error)
    return(expected_files[k])
}

check_threshold_epistemic_uncertainty_runs<-function(MIN_DOMAIN, MAX_DOMAIN){
    dirs = Sys.glob('nsw_full_coast_MIS_highres_domains_*epistemic_uncertainty*')
    all_failures = c()
    for(checkdir in dirs){
        all_tifs = Sys.glob(paste0(checkdir, '/*.tif'))
        expected_files = sapply(seq(MIN_DOMAIN, MAX_DOMAIN), function(x) gsub(as.character(MIN_DOMAIN), as.character(x), all_tifs[1]))
        files_that_failed = check_for_failed_rasters(expected_files)
        all_failures = c(all_failures, files_that_failed)
        print(files_that_failed)
    }
    return(all_failures)
}
