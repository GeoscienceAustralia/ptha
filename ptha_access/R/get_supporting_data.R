#
# Codes to copy unit_source_grid geometries, from NCI to this repository
#
if(!exists('config_env')){
    config_env = new.env()
    source('R/config.R', local=config_env)
}

#nci_thredds_source_zones_base = 'https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/'
nci_thredds_source_zones_base = paste0(config_env$.GDATA_HTTP_BASE_LOCATION, 'SOURCE_ZONES/')

install_base = getwd()

# Useful shortcut
source_names_all = config_env$source_names_all

#'
#' Download one or more unit source grid shapefiles to output_base_dir
#'    
#'
download_unit_source_grid_shapefile<-function(source_names_all, 
    output_base_dir=paste0(install_base, '/SOURCE_ZONES/')){

    # Expand the source names into all shapefile components
    all_shp_files_no_extension = paste0(nci_thredds_source_zones_base, 
        source_names_all, '/EQ_SOURCE/unit_source_grid/', source_names_all)
    all_shp_files = paste0(all_shp_files_no_extension, '.shp')
    all_shx_files = paste0(all_shp_files_no_extension, '.shx')
    all_prj_files = paste0(all_shp_files_no_extension, '.prj')
    all_dbf_files = paste0(all_shp_files_no_extension, '.dbf')

    # Download
    for(i in 1:length(source_names_all)){
        source_name = source_names_all[i]
        # Directory where we save the shapefile
        shp_dir = paste0(output_base_dir, source_name, '/EQ_SOURCE/unit_source_grid/')
        dir.create(shp_dir, recursive=TRUE, showWarnings=FALSE)
        download.file(all_shp_files[i], destfile = paste0(shp_dir, '/', basename(all_shp_files[i])), mode='wb')
        download.file(all_shx_files[i], destfile = paste0(shp_dir, '/', basename(all_shx_files[i])), mode='wb')
        download.file(all_prj_files[i], destfile = paste0(shp_dir, '/', basename(all_prj_files[i])))
        download.file(all_dbf_files[i], destfile = paste0(shp_dir, '/', basename(all_dbf_files[i])), mode='wb')
    }
}

#
# This will download shapefiles for all sources in config$source_names_all
# Useful for doing the interactive map
#
download_all_unit_source_grid_shapefiles<-function(output_base_dir=paste0(install_base, '/SOURCE_ZONES/')){

    download_unit_source_grid_shapefile(source_names_all, output_base_dir=output_base_dir)
}


#'
#' Download the unit_source_statistics_SOURCEZONE.nc file for one or more
#' sourcezones. 
#'
download_unit_source_statistics_netcdf<-function(source_names_all,
    output_base_dir=paste0(install_base, '/SOURCE_ZONES/')){

    all_unit_source_files = paste0(nci_thredds_source_zones_base, 
        source_names_all, '/TSUNAMI_EVENTS/unit_source_statistics_', 
        source_names_all, '.nc')

    for(i in 1:length(all_unit_source_files)){

        unit_source_file = all_unit_source_files[i]
        uss_dir = paste0(output_base_dir, source_name, '/TSUNAMI_EVENTS/')
        dir.create(uss_dir, recursive=TRUE, showWarnings=FALSE)
        download.file(unit_source_file, destfile = paste0(uss_dir, '/', basename(unit_source_file)), mode='wb')

    }

}

