suppressPackageStartupMessages(library(rgdal))

#
# Should we re-download all map data, even if the data is locally cached?
#
REFRESH_MAP = FALSE

#
# All source names
#
source_names_all = c(
   'alaskaaleutians',
   'arutrough',
   'banda_detachment',
   'cascadia',
   'floreswetar',
   'hjort',
   'izumariana',
   'kermadectonga',
   'kurilsjapan',
   'macquarieislandnorth',
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
   'timortrough',
   'tolo_thrust',
   'trobriand')


#
# Define the initial part of key web-addresses where our data can be accessed
#
if(file.exists('/g/data/fj6/PTHA/AustPTHA_1')){
    # We are on the NCI filesystem, read the data locally
    .GDATA_OPENDAP_BASE_LOCATION = '/g/data/fj6/PTHA/AustPTHA_1/'
    .GDATA_HTTP_BASE_LOCATION = '/g/data/fj6/PTHA/AustPTHA_1/'

}else{
    # We are not on the nci filesystem
    # These locations allow remote query/subsetting
    #
    # We include a flag that ensures strings are long enough
    # FIXME: Can reduce stringlength
    max_stringlength = '[stringlength=4096]'
    #
    .GDATA_OPENDAP_BASE_LOCATION = paste0(max_stringlength, 
        'http://dapds00.nci.org.au/thredds/dodsC/fj6/PTHA/AustPTHA_1/')
    #
    # Non-netcdf data -- this location only allows download
    .GDATA_HTTP_BASE_LOCATION = 'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/'
}

#'
#' Some files on gdata contain paths like /g/data/fj6/PTHA/AustPTHA_1/....
#' We want to replace the start of this path with the .GDATA_OPENDAP_BASE_LOCATION
#' (or perhaps the .GDATA_HTTP_BASE_LOCATION).
#' This function does that
#'
#' @param filepath The path of the file, starting with /g/data/fj6/....
#' @param prefix_type either 'opendap' or 'http'
#' @return The same filepath, with prefix suitable for remote access
#'
adjust_path_to_gdata_base_location<-function(filepath, prefix_type='opendap'){

    if(length(filepath) != 1) stop('adjust_path_to_gdata_base_location is not vectorized')

    split_path = strsplit(filepath,'/')[[1]]
    k = grep('AustPTHA_1', split_path)
    ls = length(split_path)
    if(length(k) == 0) stop('did not find AustPTHA_1 in path')

    if(k < ls){
        new_base_path = paste(split_path[(k+1):ls], collapse='/')
    }else{
        new_base_path = ''
    }

    if(prefix_type == 'opendap'){
        new_base_path = paste0(.GDATA_OPENDAP_BASE_LOCATION, new_base_path)
    }else if(prefix_type == 'http'){
        new_base_path = paste0(.GDATA_HTTP_BASE_LOCATION, new_base_path)
    }else{
        stop(paste0('Unknown value of prefix_type = ', prefix_type))
    }

    return(new_base_path)

}

# Use this variable to determine if we have already sourced config.R
.HAVE_SOURCED_CONFIG=TRUE
