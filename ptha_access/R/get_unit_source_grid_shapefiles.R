#
# Codes to copy unit_source_grid geometries, from NCI to this repository
#
config_env = new.env()
source('R/config.R', local=config_env)

nci_thredds_source_zones_base = 'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/'

install_base = getwd()

# Useful shortcut
source_names_all = config$source_names_all

# Expand the source names into all shapefile components
all_shp_files_no_extension = paste0(nci_thredds_source_zones_base, 
    source_names_all, '/EQ_SOURCE/unit_source_grid/', source_names_all)
all_shp_files = paste0(all_shp_files_no_extension, '.shp')
all_shx_files = paste0(all_shp_files_no_extension, '.shx')
all_prj_files = paste0(all_shp_files_no_extension, '.prj')
all_dbf_files = paste0(all_shp_files_no_extension, '.dbf')

dir.create('SOURCE_ZONES', showWarnings=FALSE)

# Download
for(i in 1:length(source_names_all)){
    source_name = source_names_all[i]
    # Directory where we save the shapefile
    shp_dir = paste0('SOURCE_ZONES/', source_name, '/EQ_SOURCE/unit_source_grid/')
    dir.create(shp_dir, recursive=TRUE)
    download.file(all_shp_files[i], destfile = paste0(shp_dir, '/', basename(all_shp_files[i])))
    download.file(all_shx_files[i], destfile = paste0(shp_dir, '/', basename(all_shx_files[i])))
    download.file(all_prj_files[i], destfile = paste0(shp_dir, '/', basename(all_prj_files[i])))
    download.file(all_dbf_files[i], destfile = paste0(shp_dir, '/', basename(all_dbf_files[i])))
}

