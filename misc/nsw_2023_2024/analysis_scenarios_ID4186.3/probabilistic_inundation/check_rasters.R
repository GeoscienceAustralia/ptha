check_for_failed_rasters<-function(expected_files){

    library(terra)
    read_files = lapply(expected_files, function(x) try(rast(x)))
    is_try_error = unlist(lapply(read_files, function(x) is(x, 'try-error')))
    k = which(is_try_error)
    return(expected_files[k])
}

