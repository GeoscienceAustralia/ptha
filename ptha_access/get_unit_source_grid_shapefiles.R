#
# Codes to copy unit_source_grid geometries, from NCI to this repository
#

nci_thredds_source_zones_base = 'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/'

install_base = getwd()

source_names = c(
   'alaskaaleutians',
   'arutrough',
   'banda_detachment',
   'cascadia',
   'flores',
   'hjort',
   'izumariana',
   'kermadectonga',
   'kurilsjapan',
   'macquarienorth',
   'makran',
   'manokwari',
   'manus',
   'mexico',
   'moresby_trough',
   'mussau',
   'newguinea',
   'newhebrides',
   'north_sulawesi',
   'outer_rise_timor',
   'outerrise_kermadectonga',
   'outerrise_puysegur',
   'outerrisenewhebrides',
   'outerrisesolomon',
   'outerrisesunda',
   'philippine',
   'puysegur',
   'ryuku',
   'sandwich',
   'sangihe',
   'sangihe_backthrust',
   'se_sulawesi',
   'seram_thrust',
   'seramsouth',
   'solomon',
   'southamerica',
   'sunda',
   'tanimbar',
   'timor',
   'tolo_thrust',
   'trobriand')


all_shp_files_no_extension = paste0(nci_thredds_source_zones_base, 
    source_names, '/EQ_SOURCE/unit_source_grid/', source_names)
all_shp_files = paste0(all_shp_files_no_extension, '.shp')
all_shx_files = paste0(all_shp_files_no_extension, '.shx')
all_prj_files = paste0(all_shp_files_no_extension, '.prj')
all_dbf_files = paste0(all_shp_files_no_extension, '.dbf')

dir.create('SOURCE_ZONES', showWarnings=FALSE)

for(i in 1:length(source_names)){
    source_name = source_names[i]
    # Directory where we save the shapefile
    shp_dir = paste0('SOURCE_ZONES/', source_name, '/EQ_SOURCE/unit_source_grid/')
    dir.create(shp_dir, recursive=TRUE)
    download.file(all_shp_files[i], destfile = paste0(shp_dir, '/', basename(all_shp_files[i])))
    download.file(all_shx_files[i], destfile = paste0(shp_dir, '/', basename(all_shx_files[i])))
    download.file(all_prj_files[i], destfile = paste0(shp_dir, '/', basename(all_prj_files[i])))
    download.file(all_dbf_files[i], destfile = paste0(shp_dir, '/', basename(all_dbf_files[i])))
}

## # Copy over the unit-source grids
## setwd(project_base)
## all_unit_source_grids = dirname(Sys.glob('SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp'))
## new_unit_source_grids = paste0(install_base, '/', all_unit_source_grids)
## for(i in 1:length(new_unit_source_grids)){
##     dir.create(new_unit_source_grids[i], recursive=TRUE, showWarnings=FALSE)
##     file.copy(Sys.glob(paste0(all_unit_source_grids[i], '/*')), new_unit_source_grids[i])
## }
## setwd(install_base)
## 
## dir.create('DATA', showWarnings=FALSE)
## dir.create('DATA/HAZARD_POINTS', showWarnings=FALSE)
## file.copy(paste0(project_base, 'DATA/HAZARD_POINTS/merged_hazard_points.csv'), 'DATA/HAZARD_POINTS/')
library(ncdf4)

