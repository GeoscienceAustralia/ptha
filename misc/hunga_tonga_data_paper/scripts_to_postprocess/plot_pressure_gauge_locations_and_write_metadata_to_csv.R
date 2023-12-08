#
# Make files containing the pressure gauge names, locations, and data file-paths
# Also make some plots.
#
library(sp)
source('get_simple_world_map_data.R')
source('global_variables.R')
source('create_README_in_postprocessed_folder.R')
source('create_README_in_postprocessed_graphical_checks_folder.R')

# Make the output folders
dir.create(OUTPUT_DIR, showWarnings=FALSE)
dir.create(OUTPUT_GRAPHICS_DIR, showWarnings=FALSE)

# Open a file connection to report on skipped gauges
SKIPPED_GAUGES_FILECON = file(IGNORED_MSLP_FILE, open='w')
writeLines("List of MSLP sensor time-series files that were skipped in post-processing:", SKIPPED_GAUGES_FILECON)

#
# BOM Pressure gauges
#
get_BOM_gauge_locations<-function(){

    gauge_locations = read.csv('../original/02_mslp_sensors/HD01D_StnDet.csv', header=FALSE)
    names(gauge_locations)[4] = 'station'
    names(gauge_locations)[7] = 'lat'
    names(gauge_locations)[8] = 'lon'
    names(gauge_locations)[2] = 'id'
    names(gauge_locations)[11] = 'height_m'
    names(gauge_locations)[12] = 'height_instrument_m'
    names(gauge_locations)[10] = 'state'
    names(gauge_locations)[16] = 'percent_complete'

    # Find the data files for each gauge
    string_matches = paste0('-', substring(1e+07 + as.numeric(gauge_locations$id), 3, 8), '_')
    data_files = lapply(string_matches, function(x){
        mtch = Sys.glob(paste0('../original/02_mslp_sensors/BOM_mslpdata/mslp', x, '*.csv'))
        if(length(mtch) != 1){
            return("")
        }else{
            return(mtch)
        }

    })
    gauge_locations$data_file = unlist(data_files)

    # Extract start/end times
    start_and_end_time = lapply(data_files, function(filename){
        if(filename == ""){
            start_time = NA
            end_time = NA
        }else{
            record = read.csv(filename)
            start_time = record[1,1]
            end_time = record[nrow(record), 1]
        }
        return(list(start_time=start_time, end_time=end_time))
    })
    gauge_locations$start_time = unlist(lapply(start_and_end_time, function(x) x$start_time))
    gauge_locations$end_time = unlist(lapply(start_and_end_time, function(x) x$end_time))

    # Return the named data columns as metadata
    gauge_colnames_to_keep = c('station', 'id', 'lat', 'lon', 'state',
        'height_m', 'height_instrument_m', 'start_time', 'end_time', 'percent_complete', 'data_file')
    gauge_locations = gauge_locations[, gauge_colnames_to_keep]

    # Remove locations that don't have a data file
    to_remove = which(gauge_locations$data_file == "")
    if(length(to_remove) > 0) gauge_locations = gauge_locations[-to_remove,]

    # We will put the "post-processed" data file here
    gauge_locations$postprocessed_file = paste0(OUTPUT_MSLP_DIR, 'BOM_', basename(gauge_locations$data_file))
    # 
    gauge_locations$Dataset = rep('MSLPBOM', length(gauge_locations$postprocessed_file))

    # Change the name 'data_file' to 'original_data_file' for greater clarity
    k = which(names(gauge_locations) == 'data_file')
    names(gauge_locations)[k] = 'original_data_file'

    # Remove a couple of stations where missing data prevents them from recording the HTHH waves,
    # or (for Birdsville) where a major discontinuity produces large artefacts in our the short-period wave filtering, 
    # or (for Douglas River and Thredbow) where there is lots of missing data that obscures part of the HTHH record.
    # This was determined from an initial version of the time-series plots
    IDs_to_remove = c('077094','083084', '038026', '014901', '071032')
    to_remove = sapply(IDs_to_remove, function(x) grep(x, gauge_locations$original_data_file))
    # Record that we skipped the gauges
    writeLines(gauge_locations$original_data_file[to_remove], SKIPPED_GAUGES_FILECON)
    # Print a message
    print('Deliberately removing stations that miss the volcanic pressure wave (or have lots of missing data in that time), and the Birdsville station with a major discontinuity: ')
    print(paste0('    ', IDs_to_remove, ', ', gauge_locations$station[to_remove]))
    gauge_locations = gauge_locations[-to_remove,]

    # Make ID include the organisation, now we also have DES pressure
    gauge_locations$id = paste0('BOM_', gauge_locations$id)

    return(gauge_locations)
}

#
# DES pressure gauge
#
get_DES_gauge_locations<-function(){

    gauge_xy = read.csv('../original/02_mslp_sensors/DES_station_locations_all.csv')

    # The gauge_xy file has some rows beneath the actual station location data,
    # which we remove
    k = which(gauge_xy[,1] == "") 
    if(length(k) != 1 | k <= 1){
        stop('ERROR in reading DES pressure gauge locations')
    }
    gauge_xy = gauge_xy[1:(k-1),]

    # Convert the coordinates into decimal degrees
    parse_coord<-function(coord_string){
        coord_parts = strsplit(coord_string, '° ')[[1]]
        n2 = nchar(coord_parts[2])
        coord = as.numeric(coord_parts[1]) + as.numeric(substring(coord_parts[2], 1, n2-2))/60
        # If it is 'degrees south', then make negative
        if(substring(coord_parts[2], n2) == 'S') coord = -coord
        return(coord)
    }
    lon = sapply(gauge_xy[,4], parse_coord, USE.NAMES=FALSE)
    lat = sapply(gauge_xy[,3], parse_coord, USE.NAMES=FALSE)
    gauge_name = gauge_xy[,1]

    # Some name changes are required to match filenames with station names
    # Below we ensure that gauge_name is the same as the start of the associated file
    k = which(gauge_name == "Gladstone Southtrees"); stopifnot(length(k) == 1)
    gauge_name[k] = "Gladstone South Trees"
    k = which(gauge_name == "GC Seaway") ; stopifnot(length(k) == 1)
    gauge_name[k] = "Gold Coast Seaway"
    k = which(gauge_name == "Laguna Quays") ; stopifnot(length(k) == 1)
    gauge_name[k] = "Laguna Quay"
    k = which(gauge_name == "Lucinda inshore") ; stopifnot(length(k) == 1)
    gauge_name[k] = 'Lucinda Inshore'
    k = which(gauge_name == 'Noosa River Sand Jetty'); stopifnot(length(k) == 1)
    gauge_name[k] = 'Noosa River Sand Jetty STG'
    k = which(gauge_name == 'Shorncliffe'); stopifnot(length(k) == 1)
    gauge_name[k] = 'Shorncliffe Pier'
    k = which(gauge_name == 'Lucinda'); stopifnot(length(k) == 1)
    gauge_name[k] = 'Lucinda Offshore'
    k = which(gauge_name == 'Tweed sand bypass'); stopifnot(length(k) == 1)
    gauge_name[k] = "Tweed Sand Bypass Jetty"
    
    # Find a match between the locations, and the files 
    pressure_files = Sys.glob('../original/02_mslp_sensors/DES_QGHL_pressure_data/*.csv')
    file_metadata = lapply(strsplit(basename(pressure_files), split='_'), function(x) list(station=x[1]))
    for(i in 1:length(file_metadata)){

        # Populate fields, with names matching the BOM data.
        file_metadata[[i]]$original_data_file = pressure_files[i]
        file_metadata[[i]]$postprocessed_file = paste0(OUTPUT_MSLP_DIR, 'DES_mslp-', 
            gsub(' ', '-', basename(file_metadata[[i]]$original_data_file)))
        file_metadata[[i]]$state = 'QLD'
        file_metadata[[i]]$id = paste0('DES_', i) # This is a nominal ID, because the BOM data has relevant values

        # Matching file, lon, lat
        k = match(file_metadata[[i]]$station, gauge_name)
        stopifnot(length(k) == 1)
        file_metadata[[i]]$lon = lon[k]
        file_metadata[[i]]$lat = lat[k]
        file_metadata[[i]]$height_m = 0
        file_metadata[[i]]$height_instrument_m = 0

        # Get start/end times in UTC
        pd = read.csv(pressure_files[i])
        start_end_times = as.difftime(-10, units='hours') + # 10h offset of AEST from UTC
            strptime(c(pd[1,1], pd[length(pd[,1]), 1]), format='%Y-%m-%dT%H:%M', tz='Etc/UTC')
        file_metadata[[i]]$start_time = format(start_end_times[1], '%Y-%m-%d %H:%M:%S')
        file_metadata[[i]]$end_time = format(start_end_times[2], '%Y-%m-%d %H:%M:%S')

        # Percent complete
        file_metadata[[i]]$percent_complete = mean(is.finite(pd[,2]))*100
    
    }
    

    # Make into a table, with row ordering matching the BOM data
    for(i in 1:length(file_metadata)) file_metadata[[i]] = as.data.frame(file_metadata[[i]])
    gauge_locations = do.call(rbind, file_metadata)
    gauge_colnames_to_keep = c('station', 'id', 'lat', 'lon', 'state',
        'height_m', 'height_instrument_m', 'start_time', 'end_time', 'percent_complete', 
        'original_data_file', 'postprocessed_file')
    gauge_locations = gauge_locations[, gauge_colnames_to_keep]
    # Add in a 'Dataset' column 
    gauge_locations$Dataset = rep('MSLPDES', length(gauge_locations$postprocessed_file))

    # Remove gauges that were judged to be too noisy (based on graphical checks)
    # - Geraldton
    # - Cardwell
    k = which(gauge_locations$station %in% c('Bundaberg', 'Cardwell'))
    # Record that we skipped the gauges
    writeLines(gauge_locations$original_data_file[k], SKIPPED_GAUGES_FILECON)
    # Print a related message
    print(paste0('Removing DES gauge at (too noisy) @ ', gauge_locations$station[k]))
    gauge_locations = gauge_locations[-k,]
    

    return(gauge_locations)
}

plot_pressure_gauges<-function(){

    gauge_locations_BOM = get_BOM_gauge_locations()
    gauge_locations_DES = get_DES_gauge_locations()
    gauge_locations = rbind(gauge_locations_BOM, gauge_locations_DES)

    # Make a plot of locations
    png(paste0(OUTPUT_GRAPHICS_DIR, '/BOM_MSLP_gauge_locations.png'), width=9, height=5.2, units='in', res=300)
    #par('mar') = c(5.1, 4.1, 4.1, 2.1)
    par('mar' = c(3.1, 4.1, 3.1, 2.1))

    plot(wrld_simpl, 
         xlim=c(min(gauge_locations$lon), 185),  ylim=range(gauge_locations$lat), 
         asp=1/cos(-30/180*pi), col='grey', bg='white', border='darkgray',
         axes=TRUE, xlab='', ylab='', cex.axis=1.5, cex.lab=1.5, las=1, 
         main='A) MSLP sensors', cex.main = 2.5)
    points(gauge_locations$lon, gauge_locations$lat, col='darkred', pch=19, cex=0.8)
    points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)

    # Add location of Broome Airport
    ba_i = grep('BROOME AIRPORT', gauge_locations$station)
    lonlat = c(gauge_locations$lon[ba_i], gauge_locations$lat[ba_i])
    arrows(lonlat[1], lonlat[2], 110, -17, len=0)
    text(110, -17, 'Broome Airport', adj=1, cex=1.8)

    # Add location of gauge with outlier pressure maxima time (due to oscillations)
    dw_i = grep("DWELLINGUP", gauge_locations$station)     
    lonlat = c(gauge_locations$lon[dw_i], gauge_locations$lat[dw_i])
    arrows(lonlat[1], lonlat[2], 108, -35, len=0)
    text(108, -35, 'Dwellingup', adj=1, cex=1.8)


    dev.off()


    return(gauge_locations)
}

# Make a plot, and return all the MSLP gauge metadata
pressure_gauge_locations = plot_pressure_gauges()
# Save the metadata to a CSV
write.csv(pressure_gauge_locations, MSLP_METADATA_TABLE_FILE, row.names=FALSE)

# Make the README
make_OUTPUT_DIR_README()
make_OUTPUT_GRAPHICS_DIR_README()

# Finish writing the file with skipped pressure gauges
close(SKIPPED_GAUGES_FILECON)
