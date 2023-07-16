# Workaround control vrt resolution -- does not seem to work in 'terra':https://github.com/rspatial/terra/issues/747 
# Actually after writing this code I noticed it can work with stars:
#     stars_vrt = read_stars(st_mosaic(all_files, options=c('-resolution', 'highest')))
# It can also work in terra
#     terra_vrt = vrt(all_files, options=c('-resolution', 'highest'))
make_vrt<-function(all_files){
    temp_file = paste0(tempfile(), '.vrt')
    gdal_command = paste0('gdalbuildvrt -resolution highest ', temp_file, ' ', paste(all_files, collapse=" "))
    system(gdal_command, intern=TRUE)
    return(temp_file)
}

