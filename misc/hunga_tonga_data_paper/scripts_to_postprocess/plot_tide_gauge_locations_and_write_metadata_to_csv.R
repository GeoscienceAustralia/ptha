#
# Make files containing the pressure gauge names, locations, and data file-paths
# Also make some plots.
#
library(sp)
source('get_simple_world_map_data.R')
source('global_variables.R')

#
# Tide gauges
#
# We need to combined a few data sources:
#     BOM
#     MHL
#     TAS
#     DES
#     IOC
#     NSW Port Authority [in this case we cannot release the original data, just de-tided data]
#     AAD
#

# Open a file connection to report on skipped gauges
SKIPPED_GAUGES_FILECON = file(IGNORED_TIDEGAUGE_FILE, open='w')
writeLines("List of tide gauges that were not post-processed:", SKIPPED_GAUGES_FILECON)

extract_tas_revised_tide_gauge_locations<-function(){
    all_TAS_files = Sys.glob('../original/01_tide_gauges/2022_Tsunami_TAS/*.csv')
    all_TAS_headers = lapply(all_TAS_files, function(x) readLines(x, n=27))
    names(all_TAS_headers) = gsub('.csv', '', basename(all_TAS_files), fixed=TRUE)
    lon = unlist(lapply(all_TAS_headers, function(x) as.numeric(strsplit(x[4], split=",")[[1]][2]) ))
    lat = unlist(lapply(all_TAS_headers, function(x) as.numeric(strsplit(x[3], split=",")[[1]][2]) ))
    # Lats are mostly missing the minus sign, correct this
    k = which(lat > 0)
    if(length(k) > 0) lat[k] = -lat[k]

    tas_site_name = gsub('_aest_ahd.csv', '', basename(all_TAS_files), fixed=TRUE)

    output = data.frame(name=tas_site_name, 
                        lon=lon, lat=lat, file=all_TAS_files, Dataset = rep('TAS', length(lon)),
                        data_file = all_TAS_files)
    rownames(output) = NULL
    return(output)
}

extract_macquarie_island_tide_gauge_locations<-function(){ 
    # Only one file in this case
    mac_file = Sys.glob('../original/01_tide_gauges/Macquarie_Island/Macca*.csv')
    metadata = readLines('../original/01_tide_gauges/Macquarie_Island/MaccaUTC_WGS84_JAN2022_metadata.txt')
    lon = as.numeric(strsplit(metadata[grep('Longitude:', metadata, fixed=TRUE)], ':')[[1]][2])
    lat = as.numeric(strsplit(metadata[grep('Latitude:', metadata, fixed=TRUE)], ':')[[1]][2])
    output = data.frame(name='Macquarie_Island', lon=lon, lat=lat, file=mac_file, Dataset='AAD',
            data_file=mac_file)
    rownames(output) = NULL
    return(output)
}

extract_NSWPortAuth_tide_gauge_locations<-function(){

    all_NSW_files = Sys.glob('../original/01_tide_gauges/NSW_Port_Authority_DO_NOT_DISTRIBUTE/*.csv')

    all_NSW_headers = lapply(all_NSW_files, function(x) readLines(x, n=12))
    names(all_NSW_headers) = gsub('.csv', '', basename(all_NSW_files), fixed=TRUE)

    lon = unlist(lapply(all_NSW_headers, function(x) as.numeric(strsplit(gsub(":", ",", x[7]), split=",")[[1]][3]) ))
    lat = unlist(lapply(all_NSW_headers, function(x) as.numeric(strsplit(gsub(":", ",", x[7]), split=",")[[1]][4]) ))

    # Make a name, dealing with special characters
    tmp = unlist(lapply(names(all_NSW_headers), function(x) strsplit(x, '_')[[1]][2]))
    tmp = gsub('(', '', tmp, fixed=TRUE)
    tmp = gsub(')', '', tmp, fixed=TRUE)
    tmp = gsub('-', '_', tmp, fixed=TRUE)
    tmp = gsub(' ', '_', tmp, fixed=TRUE)
    tmp = gsub('___', '_', tmp, fixed=TRUE)
    gauge_name = tmp

    output = data.frame(name=gauge_name, lon=lon, lat=lat, 
        file=all_NSW_files, Dataset = rep('PANSW', length(lon)),
        data_file = all_NSW_files)
    rownames(output) = NULL
    return(output)
}

extract_DES_tide_gauge_locations<-function(){
    file_match = '../original/01_tide_gauges/DES_QGHL_data/*.csv'
    all_DES_files = Sys.glob(file_match)
    all_DES_headers = lapply(all_DES_files, function(x) readLines(x, n=36))
    names(all_DES_headers) = gsub('.csv', '', basename(all_DES_files), fixed=TRUE)
    dms_to_dec<-function(x){
        N = gregexpr("'", x, fixed=TRUE)[[1]][1]
        if(N < 1) N = nchar(x) # Workaround needed for one file that lacked "'"
        as.numeric(substring(x, 1, 3)) + as.numeric(substring(x, 4, N-1))/60
    }
    # Here I had trouble with special characters in a few files -- the iconv(..) call takes care of that
    lon = unlist(lapply(all_DES_headers, function(x) dms_to_dec(strsplit(iconv(x[8], sub=""), split=",")[[1]][2]) ))
    lat = -1 * unlist(lapply(all_DES_headers, function(x) dms_to_dec(strsplit(iconv(x[7], sub=""), split=",")[[1]][2]) ))
    output = data.frame(name=unlist(lapply(names(all_DES_headers), function(x) strsplit(x, '_')[[1]][1])), 
                        lon=lon, lat=lat, file=all_DES_files, Dataset = rep('DES', length(lon)),
                        data_file = all_DES_files)
    rownames(output) = NULL
    return(output)
}

extract_ioc_tide_gauge_locations<-function(){
    file = '../original/01_tide_gauges/ioc_sealevelmonitoring/ioc_australian_station_list.csv'
    site_data = read.csv(file)
    lon = site_data$lon
    lat = site_data$lat
    name = site_data$Code
    longer_name = gsub(' ', '_', site_data$Location)

    # Location of the data
    data_files = paste0('../original/01_tide_gauges/ioc_sealevelmonitoring/australian_data/',
        gsub(" ", "", name), '.csv')
    k = which(!file.exists(data_files))
    data_files[k] = ""

    output = data.frame(name=longer_name, shortname = name, lon=lon, lat=lat, 
        Dataset=rep('IOC', length(name)), data_file=data_files)
    rownames(output) = NULL
    return(output)
}

extract_BOM_and_MHL_tide_gauge_locations<-function(){

    # This file initially contained an incorrect location for Thursday Island (original file fixed during QC)
    bom_mhl_tide_gauge_locations = read.csv('../original/01_tide_gauges/BOM_and_MHL_tide_gauge_metadata_Table.csv')
    n = ncol(bom_mhl_tide_gauge_locations)
    names(bom_mhl_tide_gauge_locations)[2:n] = tolower(names(bom_mhl_tide_gauge_locations)[2:n])
    rownames(bom_mhl_tide_gauge_locations) = NULL

    # Find the data files matching each site.
    target_files = c(
        Sys.glob('../original/01_tide_gauges/BOM_tidegaugedata/BOM*.csv'),
        Sys.glob('../original/01_tide_gauges/MHL_data/*.csv') )

    # Idea -- split of the text at the end of the csv filename, and find it (in
    # lower-case) in bom_mhl_tide_gauge_locations$name
    # Note this doesn't work for the modified MHL data files -- we work around this below
    target_files_tag = strsplit(basename(target_files), "_")
    target_files_startname = unlist(lapply(target_files_tag, function(x) x[1]))
    target_files_endname = gsub('.csv', '', unlist(lapply(target_files_tag, function(x) x[2])) )
    target_files_endname = gsub('-', ' ', target_files_endname) # Empty for corrected MHL files

    # NOTE: This function is potentially fragile -- it need to be checked if any input files change!
    # This partly reflects that the data organisation is not so clean.
    match_file_to_name = function(endname, startname){

            if(grepl('Level1', startname)){
                # MHL data
                possible_match = grep(tolower(substring(startname, 1, 4)), 
                    tolower(substring(bom_mhl_tide_gauge_locations$name, 1, 4)))
                if(length(possible_match) != 1){
                    # This works because of our ordering of the files, but is fragile
                    possible_match = max(possible_match)
                }
                output = possible_match

            }else{

                possible_match = grep(endname, tolower(bom_mhl_tide_gauge_locations$name), fixed=TRUE)

                # Workarounds for some special cases
                if(endname == 'sydney'){
                    # Sydney workaround
                    possible_match = which("Sydney" == bom_mhl_tide_gauge_locations$name)                
                }
                if(endname == "cocos island"){
                    possible_match = which("Cocos Is.  " == bom_mhl_tide_gauge_locations$name)
                }
                if(endname == "groote eylandt (milner bay)"){
                    possible_match = which("Milner Bay (Groote Eylandt) " == bom_mhl_tide_gauge_locations$name)
                }
                if(endname == 'thursday island'){
                    possible_match = which("Thursday Is. " == bom_mhl_tide_gauge_locations$name)
                }
                output = possible_match
            }
            return(output)
        }

    # Find the indices of bom_mhl_tide_gauge_locations$name that match each file.
    # There can be no match (e.g. we prevent matches to ioc files, as they will not be used)
    matching_site_inds = mapply(match_file_to_name,
        target_files_endname, 
        target_files_startname)
    matching_site_names = unlist(lapply(matching_site_inds, 
        function(x){
            if(length(x) == 0){
                return("")
            }else{
                return(bom_mhl_tide_gauge_locations$name[x])
            }
        }))
    ii = match(bom_mhl_tide_gauge_locations$name, matching_site_names)
    bom_target_files = target_files[ii]
    bom_target_files[is.na(bom_target_files)] = ''

    bom_mhl_tide_gauge_locations = cbind(bom_mhl_tide_gauge_locations, 
        data.frame(data_file = bom_target_files))

    return(bom_mhl_tide_gauge_locations)
}

plot_tide_gauges<-function(){
    
    bom_mhl_tide_gauge_locations = extract_BOM_and_MHL_tide_gauge_locations()
    #tas_tide_gauge_locations = extract_tas_tide_gauge_locations()
    tas_tide_gauge_locations = extract_tas_revised_tide_gauge_locations()
    des_tide_gauge_locations = extract_DES_tide_gauge_locations()
    ioc_tide_gauge_locations = extract_ioc_tide_gauge_locations()
    nswPorts_tide_gauge_locations = extract_NSWPortAuth_tide_gauge_locations()
    macquarieisland_tide_gauge_locations = extract_macquarie_island_tide_gauge_locations()

    # Combine the gauge lcoations
    all_gauge_locations = rbind(
        bom_mhl_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')],
        tas_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')],
        macquarieisland_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')],
        des_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')],
        ioc_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')],
        nswPorts_tide_gauge_locations[c('name', 'lon', 'lat', 'Dataset', 'data_file')]
        )

    # Rename 'data_file' to 'original_data_file' for clarity in the output
    k = which(names(all_gauge_locations) == 'data_file')
    names(all_gauge_locations)[k] = 'original_data_file'

    # Improve the gauge names
    all_gauge_locations$name = gsub(' ', '_', all_gauge_locations$name) # Remove spaces
    all_gauge_locations$name = gsub('(', '', all_gauge_locations$name, fixed=TRUE) # Remove parenthesis
    all_gauge_locations$name = gsub(')', '', all_gauge_locations$name, fixed=TRUE) # Remove parenthesis
    all_gauge_locations$name = gsub('.', '', all_gauge_locations$name, fixed=TRUE) # Remove fullstop
    all_gauge_locations$name = gsub('_$', '', all_gauge_locations$name) # Remove trailing underscores

    # The above changes could lead to non-unique names -- make them unique by adding in the Dataset
    all_gauge_locations$name = paste0(all_gauge_locations$name, '_', all_gauge_locations$Dataset)

    all_gauge_locations$name = gsub('/', '_', all_gauge_locations$name, fixed=TRUE) # Remove '/'

    # Remove gauges that are missing data files. This happens due to a mismatch
    # between stations recorded in some metadata files, and the data we
    # actually receive.
    to_remove = which(all_gauge_locations$original_data_file == "")
    print(c('Removing rows without matching files: ', 
        paste0("    ", all_gauge_locations$name[to_remove])))
    writeLines(paste0(all_gauge_locations$name[to_remove], ' in Dataset ', 
            all_gauge_locations$Dataset[to_remove], 
            ' (Station name included in metadata, but not matching time series file was found)'), 
        SKIPPED_GAUGES_FILECON)
    if(length(to_remove) > 0)  all_gauge_locations = all_gauge_locations[-to_remove,]

    # Remove gauges with "clearly problematic" data in terms of our ability to see the tsunami
    k = which(all_gauge_locations$name == "Port_Giles_BOMPorts")
    if(length(k) > 0){
        # Record that we skipped the gauge
        writeLines(paste0(all_gauge_locations$name[k], ' with time-series file ', 
            all_gauge_locations$original_data_file[k], ' (QC identified problematic time-series)'), 
            SKIPPED_GAUGES_FILECON)
        print(paste0('Removing noisy gauge:', all_gauge_locations$name[k]))
        all_gauge_locations = all_gauge_locations[-k,]
    }

    # Make a name that will hold the output file
    all_gauge_locations = cbind(all_gauge_locations, 
        data.frame(postprocessed_file=paste0(OUTPUT_TIDE_DIR, all_gauge_locations$name, '.csv')))

    png(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_locations.png'), width=9, height=5.2, units='in', res=300)
    par('mar' = c(3.1, 4.1, 3.1, 2.1))
    plot(wrld_simpl, 
         xlim=c(min(all_gauge_locations$lon, na.rm=TRUE), 185), 
         ylim=range(all_gauge_locations$lat, na.rm=TRUE), 
         asp=1/cos(-30/180*pi), col='gray', bg='white', border='darkgray',
         axes=TRUE, xlab='', ylab='', cex.axis=1.6, cex.lab=1.5, las=1,
         main='B) Tide gauges', cex.main = 2.5)
    points(all_gauge_locations$lon, all_gauge_locations$lat, col='darkred', pch=19, cex=0.8)
    points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)
    dev.off()
    return(all_gauge_locations)

}

tide_gauge_locations = plot_tide_gauges()
write.csv(tide_gauge_locations, file=TIDEGAUGE_METADATA_TABLE_FILE, row.names=FALSE)

close(SKIPPED_GAUGES_FILECON)
