fj6_base_dir = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/'

source_zones = basename(list.dirs('.', recursive=FALSE))
dirs_to_copy = paste0(source_zones, '/EQ_SOURCE_bathyslope2025/Unit_source_data/', source_zones)

for(i in 1:length(source_zones)){

    sz = source_zones[i]
    sz_dir_local = dirs_to_copy[i]
    stopifnot(file.exists(sz_dir_local))

    sz_dir_fj6 = paste0(fj6_base_dir, sz_dir_local)
    print(c(sz_dir_local, sz_dir_fj6))
    dir.create(sz_dir_fj6, recursive=TRUE, showWarnings=FALSE)

    tiff_files_to_copy = Sys.glob(paste0(sz_dir_local, '/*_vertical_inc_bathyslope_kajiura_*.tif'))
    related_RDS_files = Sys.glob(paste0(sz_dir_local, '/*.RDS'))

    if(length(tiff_files_to_copy) != length(related_RDS_files)){
        stop(paste0('Not every RDS file for ', sz, ' has a corresponding tif file'))
    }

    copy_worked = file.copy(tiff_files_to_copy, sz_dir_fj6)
    if(!all(copy_worked)) print(paste0('Not all copies worked for ', sz ))
}
