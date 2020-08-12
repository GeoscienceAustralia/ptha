#
# Record gauges data files on a per-site / per-event basis, and provide a
# simple interface to read these files (which are diverse).
#
# Each site has a name (e.g. 'Eden_1min_DPIE'), an associated coordinate, and
# an associated list of de-tided tsunami data files, where the 'list name'
# corresponds to the GMT date on which the earthquake happened.
#
# The data files are in a range of formats -- so we'll need a suite of
# functions to read them all.
#


#
# Read the 1 min data that was processed from NSW tidal data (Department of
# Planning Infrastructure & Environment).
#
read_tsunami_DPIE<-function(filename){
    # Here the data gives a julian time. The 'raw' data from DPIE was in AEST, 
    # but the code which produced the file 'filename' converted to GMT time.
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    output = data.frame(time=tide_obs$time, juliant = tide_obs$juliant, 
                        height = tide_obs$height, resid=tide_obs$resid)
    return(output)
}

#
# Read the 1960 Chile tsunami record @ Fort Denison
#
read_tsunami_Chile1960_FortDenison<-function(filename){
    # This data was digitized by Kaya Wilson from a scan of the archives (by Dave Hanslow).
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    # The start-time was written as 24/05/1960 17:00. However the earthquake occurred
    # on 1960-06-22 19:11:20 GMT. From the data & expected arrival time, I think the
    # date must actually be 23/05/1960 17:00 in AEST.
    # Convert to GMT here (10 hours back)
    start_time = strptime('1960-05-23 07:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
    all_times = start_time + as.difftime(tide_obs[,2], units='secs')

    all_julian_times = julian(all_times)

    output = data.frame(time=all_times, juliant = all_julian_times, height=tide_obs[,5],
                        resid = tide_obs[,6])
    return(output)
}

#
# Read the 1960 Chile tsunami record @ Cronulla CSIRO Fisheries lab.
#
read_tsunami_Chile1960_Cronulla<-function(filename){
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    # The time-zone is Sydney local time (10 hours ahead of GMT). Here we
    # convert to GMT time
    times = strptime(tide_obs[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT') -
        as.difftime(10, units='hours')
    juliant = julian(times)

    output = data.frame(time=times, juliant=juliant, height=tide_obs$observed_stage_m,
                        resid=tide_obs$resid_m)
    return(output)
}


#
# Read some data for the Andaman 2004 event, provided to GA by BOM 
#
read_tsunami_BOM2004<-function(filename){

    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    # The time is in UTC
    all_times = strptime(tide_obs$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
    all_julian_times = julian(all_times)

    output = data.frame(time=all_times, juliant = all_julian_times, height=tide_obs$height,
                        resid = tide_obs$resid)
    return(output)
}


# Post-processed versions of files received by BOM -- I extract month-long
# subsets around the tsunami events and de-tide.
read_tsunami_BOM_1min<-function(filename){
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    all_times = strptime(tide_obs$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
    juliant = tide_obs$juliant
    height = tide_obs$stage
    resid = tide_obs$tsunami
    output = data.frame(time = all_times, juliant = juliant, height=height, resid=resid)
    return(output)
}

# Work on post-processed files from WA DOT.
read_tsunami_WADoT_5min<-function(filename){
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    all_times = strptime(tide_obs$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
    juliant = tide_obs$juliant
    height = tide_obs$stage
    resid = tide_obs$tsunami
    output = data.frame(time = all_times, juliant = juliant, height=height, resid=resid)
    return(output)
}


#
# Read some post-processed (de-tided) data from the Port Authority. We cannot distribute
# the raw data but can freely distribute this de-tided data. Hence 'height' is all NA values
#
read_detided_portauthority_1min<-function(filename){
    detided_obs = read.csv(filename, stringsAsFactors=FALSE)
    all_times = strptime(detided_obs[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
    all_julian_times = julian(all_times)
    output = data.frame(time=all_times, juliant=all_julian_times, height=detided_obs$observed_stage, 
        resid=detided_obs$residual)
    return(output)
}

GAUGE_DATA = list(

    #$
    # Twofold Bay (Eden) 1min BOM, after 2007$
    #$
    'TwofoldBay_1min_BOM' = list(
        # Coordinate from BOM metadata
        coord = c(149 + 54/60 + 27.86128/(60*60), -(37 + 4/60 + 25.17661/(60*60))),
        events = list(
            # Again there is a fair bit of activity going on around this
            '2010-02-27' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2010-02-27.csv',
            # Tohoku is quite clear, but there is something else big around March 22 -- but super high-frequency 
            # content [artefacts or just very short motions?]
            '2011-03-11' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2011-03-11.csv',
            # Pretty hard to detect.
            '2015-09-16' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2015-09-16.csv'
            ),
        readfun = read_tsunami_BOM_1min
        ),


    #
    # Batemans bay 1 min data
    #
    'BatemansBay_PrincessJetty_1min_DPIE' = list(
        coord = c(150.177831, -35.703811),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Princess Jetty 1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Ulladulla 1 min data
    #                   
    'Ulladulla_1min_DPIE' = list(
        coord = c(150.476525, -35.357672),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Ulladulla  1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Jervis Bay 1 min data
    #
    'JervisBay_1min_DPIE' = list(
        coord = c(150.707439, -35.121953),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Jervis Bay 1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Bundeena 1 min data
    #
    'Bundeena_1min_DPIE' = list(
        coord = c(151.1509, -34.082683),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Bundeena 1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Sydney Middle Harbour 1 min data
    #
    'Sydney_MiddleHarbour_1min_DPIE' = list(
        coord = c(151.258533, -33.825461),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Sydney 1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Sydney, Fort Dension, 1960 tidal gauge copied from the archives (by David Hanslow?) and digitized by Kaya Wilson
    #
    'Sydney_FortDenison1960' = list(
        coord = c(151.223, -33.851),
        events = list(
            '1960-05-22' = './DATA/TIDES/NSW_TIDAL_GAUGES/Fort_Denison_Chile1960/Chile1960_eventdata_4Gareth/Chile1960_tide_tsunami_residual.csv'
            ),
        readfun = read_tsunami_Chile1960_FortDenison
        ),

    #
    # Cronulla (CSIRO Fisheries) data digitized from an old Fisheries publication
    #
    'Cronulla_CSIRO_Fisheries_1960' = list(
        coord = c(151.1474, -34.0728), # From Becarri (2009) report on Chile 1960 tsunami in NSW.
        events = list(
            '1960-05-22' = './DATA/TIDES/NSW_TIDAL_GAUGES/Cronulla_Chile1960/Cronulla_Fisheries_1960_05_23_tidalgauge.csv'
            ),
        readfun = read_tsunami_Chile1960_Cronulla
        ),
    

    # Sydney data from Port Authority
    'Sydney_FortDenison_1min_PA' = list(
        coord = c(151 + 13.55/60, -(33 + 51.28/60)),
        events = list(
            '2004-12-26' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20041226_-_20050102_Fort_Denison_1_minute.rpt.csv',
            '2010-02-27' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20100227_-_20100306_Fort_Denison_1_minute.rpt.csv',
            '2011-03-11' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20110311_-_20110318_Fort_Denison_1_minute.rpt.csv',
            '2015-09-16' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20150916_-_20150923_Fort_Denison_1_minute.rpt.csv'
            ),
        readfun = read_detided_portauthority_1min
        ),

    # Botany bay data from Port Authority		  
    'Sydney_BotanyBay_1min_PA' = list(
        coord = c(151 + 12/60 + 41.97/(60*60), -(33 + 58/60 + 26.33/(60*60))),
        events = list(
            '2004-12-26' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20041226_-_20050102_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2010-02-27' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20100227_-_20100306_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2011-03-11' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20110311_-_20110318_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2015-09-16' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20150916_-_20150923_Bulk_Liquids_Berth_1_minute.rpt.csv'
            ),
        readfun = read_detided_portauthority_1min
        ),

    #
    # Hawkesbury (Patonga) 1 min data
    #
    'Hawkesbury_Patonga_1min_DPIE' = list(
        coord = c(151.274619 ,-33.550983),
        events = list(
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Patonga 1min.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    #
    # Port Kembla -- Andaman 2004
    #                   
    'PortKembla_BOM_1min_2004' = list(
        # Coordinate from Port Authority web-page (location in 2020)
        coord = c(150 + 54.71/60, -(34 + 28.43/60) ),
        events = list(
            # Sumatra 2004, from a record provided by BOM
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004/detided_obs_2004/port_kembla_2004_tsunami.csv"
            ),
        readfun = read_tsunami_BOM2004
        ),


   #
    # Port Kembla 1min BOM, after 2007
    #
    'PortKembla_1min_BOM' = list(
        # Coordinate from Port Authority web-page (location in 2020)
        coord = c(150.9118, -34.4734),
        events = list(
            '2010-02-27' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2010-02-27.csv',
            '2011-03-11' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2011-03-11.csv',
            '2015-09-16' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2015-09-16.csv'
            ),
        readfun = read_tsunami_BOM_1min
        ),

    #
    # Portland -- Andaman 2004
    #
    'Portland_BOM_1min_2004' = list(
        coord = c(141 + 36.8/60, -(38 + 20.6/60) ),
        events = list(
            # Sumatra 2004, from a record provided by BOM
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004/detided_obs_2004/portland_2004_tsunami.csv"
            ),
        readfun = read_tsunami_BOM2004
        ),

    #
    # Hillarys -- Andaman 2004
    #
    'Hillarys_BOM_1min_2004' = list(
        coord = c(115 + 44.3/60, -(31 + 49.5/60) ),
        events = list(
            # Sumatra 2004, from a record provided by BOM
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004/detided_obs_2004/hillarys_2004_tsunami.csv"
            ),
        readfun = read_tsunami_BOM2004
        ),

    # Freemantle Harbour
    'Freemantle_WADoT_5min_2004' = list(
        coord = c( 115.748056, -32.065556 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Freemantle.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    # Mangles Bay
    'ManglesBay_WADoT_5min_2004' = list(
        coord = c( 115.703333, -32.274444 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/ManglesBay.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    # Barrak Street
    'BarrackStreet_WADoT_5min_2004' = list(
        coord = c(115.85725, -31.959694),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/BarrackSt.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),
        
    'NOTHING' = list(
        coord = c(NA, NA),
        events = list()
        )

    )

if(FALSE){
    # Move the files to the current directory so it can be distributed for the paper
    files_to_create = unlist(lapply(GAUGE_DATA, f<-function(x) unlist(x$events)))
    folders_to_create = dirname(files_to_create)
    for(i in 1:length(folders_to_create)){
        dir.create(folders_to_create[i], showWarnings=FALSE, recursive=TRUE)
        target_file = paste0('/g/data/w85/tsunami/', 
    			 substring(folders_to_create[i], 3, nchar(folders_to_create[i])), 
    		 '/', basename(files_to_create[i]))
        file.copy(target_file, folders_to_create[i])
    }
    stop()
}


all_files_exist = lapply(GAUGE_DATA, f<-function(x) all(unlist(lapply(x$events, file.exists))))
if(!all(unlist(all_files_exist))){
    print('WARNING: Some files do not exist')
    print(all_files_exist)
}



# Store the directory, so if we source this from another location we can still find the gauges
BASE_DIR = getwd()

# Print a warning if not all the files are found -- this helps port the code elsewhere.
all_files_exist = lapply(GAUGE_DATA, f<-function(x) all(unlist(lapply(x$events, file.exists))))
if(!all(unlist(all_files_exist))){
    print('WARNING: Some files do not exist')
    print(all_files_exist)
}


#' Read tidal data corresponding to an event
#'
#' @param date a character string of such as '2011-03-11', matching that used
#' to encode events in the GAUGE_DATA list.
#' @return a data structure containing the data
get_data_for_event<-function(date){

    # Find gauges that include an event with name matching "date"
    matching_events = lapply(GAUGE_DATA, f<-function(x) which(names(x$events) == date))
    gauges_with_data = which(unlist(lapply(matching_events, f<-function(x) length(x) > 0)))

    if(length(gauges_with_data) == 0) stop(paste0('No gauges have data matching date ', date))

    # Copy the GAUGE_DATA structure to be modified for output
    output = GAUGE_DATA
    # Store the matching event index
    for(i in 1:length(output)){
        output[[i]]$matching_event_ind = matching_events[[i]]
    }
    # Remove events without data
    output = output[gauges_with_data]
    # Remove non-matching events
    for(i in 1:length(output)){
        output[[i]]$events = output[[i]]$events[output[[i]]$matching_event_ind]
        output[[i]]$events = paste0(BASE_DIR, '/', output[[i]]$events)
    }
    # Append the data
    for(i in 1:length(output)){
        output[[i]]$obs = output[[i]]$readfun(output[[i]]$events)
    }

    return(output)

}

