make_vrt<-function(all_files){
    temp_file = paste0(tempfile(), '.vrt')
    gdal_command = paste0('gdalbuildvrt -resolution highest ', temp_file, ' ', paste(all_files, collapse=" "))
    system(gdal_command, intern=TRUE)
    return(temp_file)
}

