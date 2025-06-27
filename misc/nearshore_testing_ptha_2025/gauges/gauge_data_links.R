#
# Record gauges data files on a per-site / per-event basis, and provide a
# simple interface to read these files (which are diverse).
#
# Each site has a name (e.g. 'Eden_1min_DPIE'), an associated coordinate, and
# an associated list of de-tided tsunami data files, where the 'list name'
# corresponds to the GMT date on which the earthquake happened.
#
# The data files are in a range of formats -- so we need a suite of functions
# to read them all. But all functions return a data.frame with the entries:
#     time -- the time as a string ('year-month-day hour-minute-seconds') in GMT time
#     juliant -- the time in 'days since 1970-01-01 00:00:00 GMT' 
#     height -- the measured series
#     resid -- an "empirically de-tided" time-series representing the tsunami
#


#
# Read the 1 min and 15 min data that was processed from NSW tidal data
# (DPIE = Department of Planning Infrastructure & Environment). This also works
# for the 6 minute BOM data. All of these were subject to variants of my
# dart-buoy processing scripts, and so the output format is the same.
#
read_tsunami_DPIE<-function(filename){
    # Here the data gives a julian time. The 'raw' data from DPIE was in AEST, 
    # but the code which produced the file 'filename' converted to GMT time.
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)

    # Work-around for the fact that we have a few different file formats
    if('resid' %in% names(tide_obs) & 'height' %in% names(tide_obs)){

        # Much of the post-processed DPIE data is like this
        output = data.frame(time=tide_obs$time, juliant = tide_obs$juliant, 
                            height = tide_obs$height, resid=tide_obs$resid)

    }else if('tsunami' %in% names(tide_obs) & 'stage' %in% names(tide_obs)){

        # The 2021 events have a different post-processed file structure
        output = data.frame(time=tide_obs$time, juliant = tide_obs$juliant, 
                            height = tide_obs$stage, resid=tide_obs$tsunami)

    }else{
        stop('Unknown DPIE post-processed file')
    }
    return(output)
}

#
# Read the 1960 Chile tsunami record @ Fort Denison
#
read_tsunami_Chile1960_FortDenison<-function(filename){
    # This data was digitized by Kaya Wilson from a scan of the archives (that
    # I think Dave Hanslow did?).
    tide_obs = read.csv(filename, stringsAsFactors=FALSE)
    # The start-time was written as 24/05/1960 17:00. However the earthquake occurred
    # on 1960-05-22 19:11:20 GMT. From the data & expected arrival time, the
    # date must actually be 23/05/1960 17:00 in AEST.
    # Convert to GMT here (10 hours back).
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
# Read some data for the Andaman 2004 event, provided to GA by BOM, and subsequently
# further post-processed to remove the long-period component. 
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

    output = data.frame(time=all_times, juliant=all_julian_times, 
        height=detided_obs$observed_stage, resid=detided_obs$residual)

    return(output)
}

#
#
#
read_WADOT_revised_2024_detiding<-function(filename){
    detided_obs = read.csv(filename, comment.char='#', header=TRUE)     
    # These datasets are in WA standard time, convert to GMT
    all_times = strptime(detided_obs$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT') - as.difftime(8, units='hours')
   
    output = data.frame(time=all_times, juliant=julian(all_times),
        height=detided_obs$height, resid=detided_obs$highfreq)
    return(output)
}



# Given a gauge observation with numeric vectors 'time, stage', compute a time-smoothed stage
# with the given smooth_time
smooth_stage_by_time<-function(time, stage, smooth_time){
    
    f = approxfun(time, stage, rule=2)
    # Make a bunch of lag-times, and average over them
    ts = seq(-smooth_time/2, smooth_time/2, len=15)
    smooth_stage = 0
    for(i in 1:length(ts)){
        smooth_stage = smooth_stage + f(time + ts[i])
    }
    smooth_stage = smooth_stage/length(ts)

    return(smooth_stage)
}

smooth_model_1min<-function(time, stage){
    return(smooth_stage_by_time(time, stage, 60))
}

smooth_model_6min<-function(time, stage){
    return(smooth_stage_by_time(time, stage, 60*6))
}

GAUGE_DATA = list(

    #
    # The first few gauges below are from a particular MHL delivery. 
    # For some of these gauges we also have other datasets -- those cases are dealt with separately, deeper in the script.
    # 

    "BallinaBreakwall_1min_DPIE" = list(
        coord = c(153.584429,-28.875377),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_BallinaBreakwall.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_BallinaBreakwall.Level1.csv", 
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_BallinaBreakwall.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_BallinaBreakwall.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_BallinaBreakwall.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_BallinaBreakwall.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_BallinaBreakwall.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_BallinaBreakwall.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "Bermagui_1min_DPIE" = list(
        coord = c(150.071478,-36.426333),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Bermagui.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Bermagui.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Bermagui.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Bermagui.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Bermagui.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Bermagui.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Bermagui.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Bermagui.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "BrunswickHeads_1min_DPIE" = list(
        coord = c(153.552769,-28.537025),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_BrunswickHeads.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_BrunswickHeads.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_BrunswickHeads.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_BrunswickHeads.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_BrunswickHeads.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_BrunswickHeads.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_BrunswickHeads.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_BrunswickHeads.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_BrunswickHeads.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

   "CoffsHarbourInnerPumpoutJetty_1min_DPIE" = list(
        coord = c(153.146144,-30.302869),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_CoffsHarbourInnerPumpoutJetty.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_CoffsHarbourInnerPumpoutJetty.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "CrookhavenHeads_1min_DPIE" = list(
        coord = c(150.759397,-34.90534),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_CrookhavenHeads.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_CrookhavenHeads.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_CrookhavenHeads.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_CrookhavenHeads.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_CrookhavenHeads.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_CrookhavenHeads.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_CrookhavenHeads.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_CrookhavenHeads.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_CrookhavenHeads.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "CrowdyHeadFishermansWharf_1min_DPIE" = list(
        coord = c(152.750014,-31.838706),
        events = list(
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_CrowdyHeadFishermansWharf.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_CrowdyHeadFishermansWharf.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_CrowdyHeadFishermansWharf.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_CrowdyHeadFishermansWharf.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_CrowdyHeadFishermansWharf.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_CrowdyHeadFishermansWharf.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_CrowdyHeadFishermansWharf.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_CrowdyHeadFishermansWharf.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "Forster_1min_DPIE" = list(
        coord = c(152.508206,-32.173986),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Forster.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Forster.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Forster.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Forster.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Forster.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Forster.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Forster.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Forster.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Forster.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "LordHoweIsland_1min_DPIE" = list(
        coord = c(159.05815,-31.523638),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_LordHoweIsland.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_LordHoweIsland.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_LordHoweIsland.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_LordHoweIsland.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_LordHoweIsland.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_LordHoweIsland.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_LordHoweIsland.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_LordHoweIsland.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_LordHoweIsland.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    "PortMacquarie_1min_DPIE" = list(
        coord = c(152.911125,-31.426825),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_PortMacquarie.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_PortMacquarie.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_PortMacquarie.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_PortMacquarie.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_PortMacquarie.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_PortMacquarie.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_PortMacquarie.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_PortMacquarie.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_PortMacquarie.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


   "ShoalBay_1min_DPIE" = list(
        coord = c(152.17565,-32.719672),
        events = list(
            # The 2014-04-01 event occurred just before the gauge started
            #"2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_ShoalBay.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_ShoalBay.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_ShoalBay.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_ShoalBay.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_ShoalBay.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_ShoalBay.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_ShoalBay.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_ShoalBay.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    "TweedEntranceSouth_1min_DPIE" = list(
        coord = c(153.551186,-28.170639),
        events = list(
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_TweedEntranceSouth.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_TweedEntranceSouth.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_TweedEntranceSouth.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_TweedEntranceSouth.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_TweedEntranceSouth.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_TweedEntranceSouth.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_TweedEntranceSouth.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    "Yamba_1min_DPIE" = list(
        coord = c(153.362061,-29.428958),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Yamba.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Yamba.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Yamba.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Yamba.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Yamba.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Yamba.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Yamba.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Yamba.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Yamba.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Eden 1 min data
    #
    'Eden_1min_DPIE' = list(
        coord = c(149.908288, -37.071241),
        events = list(
            # Sunda 2004
            '2004-12-26' = "./DATA/TIDES/NSW_TIDAL_GAUGES/Eden_MHL/de_tided_data/Eden_1min_Dec_2004.csv",
            # Solomon 2007
            '2007-04-01' = "./DATA/TIDES/NSW_TIDAL_GAUGES/Eden_MHL/de_tided_data/Eden_1min_Apr_2007.csv",
            # Puysegur 2009
            '2009-07-15' = "./DATA/TIDES/NSW_TIDAL_GAUGES/Eden_MHL/de_tided_data/Eden_1min_Jul_2009.csv",
            # Chile 2010
            '2010-02-27' = "./DATA/TIDES/NSW_TIDAL_GAUGES/Eden_MHL/de_tided_data/Eden_1min_Feb_2010.csv",
            # Tohoku 2011
            '2011-03-11' = "./DATA/TIDES/NSW_TIDAL_GAUGES/Eden_MHL/de_tided_data/Eden_1min_Mar_2011.csv",
            # Santa Cruz 2013
            '2013-02-06' = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Eden.Level1.csv",
            # South America, 2014 Mw 8.2
            '2014-04-01' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Eden 1min.csv",
            # South America, 2015 Mw 8.3
            '2015-09-16' = "./DATA/TIDES/NSW_TIDAL_GAUGES/tsunami_extract_revised/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Eden 1min.csv",
            # Solomon event
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Eden.Level1.csv",
            # Mexico event
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Eden.Level1.csv",
            # New Hebrides, 2021 Mw 7.7
            '2021-02-10' = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW_2021_Feb_March/detided_data_near_events/Eden.csv",
            # Kermadec, 2021 Mw 8.1 -- same file as above
            '2021-03-04' = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW_2021_Feb_March/detided_data_near_events/Eden.csv",
            # Hunga tonga
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Eden.Level1.csv",
            # New Hebrides
            '2023-05-19' = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Eden.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Twofold Bay (Eden) 1min BOM, after 2007$
    #
    'TwofoldBay_1min_BOM' = list(
        # This coordinate was is from the ioc sea-level website: https://www.ioc-sealevelmonitoring.org/station.php?code=tbwc
        # and is consistent with information Yuelong Miao sent me. It was previously wrong but corrected in 2021 or 2022.
        coord = c(149.9266, -37.1003), 
        events = list(
            '2009-07-15' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2009-07-15.csv',
            '2009-09-29' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2009-09-29.csv',
            '2010-02-27' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2010-02-27.csv',
            '2011-03-11' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2011-03-11.csv',
            '2014-04-01' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2014-04-01.csv',
            '2015-09-16' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/twofold_near_2015-09-16.csv',
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/TwofoldBay.csv',
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/TwofoldBay.csv'
            ),
        readfun = read_tsunami_BOM_1min
        ),

    #
    # Gold Coast Sand Bypass Jetty (from IOC website)
    #
    'GoldCoastSandBypass' = list(
        coord = c(153.4326, -27.9387), # From IOC website
        events = list(
            # New Hebrides, 2021 Mw 7.7
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/GoldCoastSandBypass.csv',
            # Kermadec-Tonga, 2021 Mw 8.1, same dataset as above
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/GoldCoastSandBypass.csv'
            ),
        readfun = read_tsunami_DPIE
        ),


    #
    # Batemans bay 1 min data
    #

    "BatemansBay_PrincessJetty_1min_DPIE" = list(
        coord = c(150.177831,-35.703811),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_PrincessJetty.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_PrincessJetty.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_PrincessJetty.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_PrincessJetty.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_PrincessJetty.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_PrincessJetty.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_PrincessJetty.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_PrincessJetty.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_PrincessJetty.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    #
    # Ulladulla 1 min data
    #                   
    "Ulladulla_1min_DPIE" = list(
        coord = c(150.476525,-35.357672),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Ulladulla.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Ulladulla.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Ulladulla.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Ulladulla.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Ulladulla.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Ulladulla.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Ulladulla.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Ulladulla.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Ulladulla.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    #
    # Jervis Bay 1 min data
    #
    "JervisBay_1min_DPIE" = list(
        coord = c(150.707439,-35.121953),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_JervisBay.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_JervisBay.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_JervisBay.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_JervisBay.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_JervisBay.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_JervisBay.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_JervisBay.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_JervisBay.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_JervisBay.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Bundeena 1 min data
    #
    "Bundeena_1min_DPIE" = list(
        coord = c(151.1509,-34.082683),
        events = list(
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Bundeena.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Bundeena.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Bundeena.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Bundeena.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Bundeena.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Bundeena.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Bundeena.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    #
    # Sydney Middle Harbour 1 min data
    #
    "Sydney_MiddleHarbour_1min_DPIE" = list(
        coord = c(151.258533,-33.825461),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Sydney.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Sydney.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Sydney.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Sydney.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Sydney.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Sydney.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Sydney.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Sydney.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Sydney.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),


    #
    # Sydney, Fort Dension, 1960 tidal gauge copied from the archives (by David Hanslow?) and digitized by Kaya Wilson
    #
    'Sydney_FortDenison1960' = list(
        # This coordinate was originally incorrect (off by a few hundred meters), fixed 2023/06/06
        coord = c(151 + 13.55/60, -(33 + 51.28/60)),
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
        ## Minor update based on Port Authority email during 2023 data purchase
        #coord = c(151.22577, -33.854651),
        events = list(
            '2004-12-26' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20041226_-_20050102_Fort_Denison_1_minute.rpt.csv',
            '2007-04-01' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20070401_-_20070408_Fort_Denison_1_minute.rpt.csv',
            '2009-07-15' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20090715_-_20090722_Fort_Denison_1_minute.rpt.csv',
            '2009-09-29' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20090929_-_2009106_Fort_Denison_1_minute.rpt.csv',
            '2010-02-27' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20100227_-_20100306_Fort_Denison_1_minute.rpt.csv',
            '2011-03-11' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20110311_-_20110318_Fort_Denison_1_minute.rpt.csv',
            '2014-04-01' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20140401_-_20140408_Fort_Denison_1_minute.rpt.csv',
            '2015-09-16' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20150916_-_20150923_Fort_Denison_1_minute.rpt.csv'
            ),
        readfun = read_detided_portauthority_1min
        ),

    # More Port Authority data at Fort Denison -- but it uses a different file format to the batch above.
    'Sydney_FortDenison_1min_PA_b' = list(
        # Coordinate is slightly different to site above, based on Port Authority email during 2023 data purchase.
        # But I think they should be considered "the same". 
        # Also the Port Authority confirmed that their "FDT01" and "FDT02" sites are the same.
        coord = c(151.22577, -33.854651),
        events = list(
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/fort_denison_newHebrides_2021_02_10_detided_WITHOUT_RAW_HEIGHT.csv',
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/fort_denison_kermadec_2021_03_04_detided_WITHOUT_RAW_HEIGHT.csv',
            '2023-05-19' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/fort_denison_newHebrides_2023_05_19_detided_WITHOUT_RAW_HEIGHT.csv'
            ),
        readfun = read_tsunami_DPIE
        ),

    # Botany bay data from Port Authority -- Bulk Liquids Berth
    'Sydney_BotanyBay_1min_PA' = list(
        coord = c(151 + 12/60 + 41.97/(60*60), -(33 + 58/60 + 26.33/(60*60))),
        events = list(
            '2004-12-26' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20041226_-_20050102_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2007-04-01' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20070401_-_20070408_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2009-07-15' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20090715_-_20090722_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2009-09-29' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20090929_-_2009106_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2010-02-27' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20100227_-_20100306_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2011-03-11' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20110311_-_20110318_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2014-04-01' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20140401_-_20140408_Bulk_Liquids_Berth_1_minute.rpt.csv',
            '2015-09-16' = './DATA/TIDES/NSW_TIDAL_GAUGES/port_authority/data/DETIDED_DATA/TSUNAMI_ONLY_20150916_-_20150923_Bulk_Liquids_Berth_1_minute.rpt.csv'
            ),
        readfun = read_detided_portauthority_1min
        ),

    # Botany Bay (Pilot Jetty)
    'Sydney_Botany_Bay_Pilot_Jetty_1min_PA' = list(
        coord = c(151.220267, -33.967903),
        events = list(
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/botany_pilot_jetty_newHebrides_2021_02_10_detided_WITHOUT_RAW_HEIGHT.csv',
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/botany_pilot_jetty_kermadec_2021_03_04_detided_WITHOUT_RAW_HEIGHT.csv',
            '2023-05-19' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/botany_pilot_jetty_newHebrides_2023_05_19_detided_WITHOUT_RAW_HEIGHT.csv'
        ),
        readfun = read_tsunami_DPIE
        ),

    # Botany Bay (Kurnell Wharf)
    'Sydney_Botany_Bay_Kurnell_Wharf_1min_PA' = list(
        coord = c(151.211088, -34.004186),
        events = list(
            '2023-05-19' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/botany_kurnell_wharf_newHebrides_2023_05_19_detided_WITHOUT_RAW_HEIGHT.csv'
        ),
        readfun = read_tsunami_DPIE
    ),

    #
    # Hawkesbury (Patonga) 1 min data
    #
    "Hawkesbury_Patonga_1min_DPIE" = list(
        coord = c(151.274619,-33.550983),
        events = list(
            "2013-02-06" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/santaCruz_2013_02_06_Mw7.9/santaCruz_2013_02_06_Mw7.9_Patonga.Level1.csv",
            "2014-04-01" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2014_04_01_Mw8.2/southamerica_2014_04_01_Mw8.2_Patonga.Level1.csv",
            "2015-09-16" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/southamerica_2015_09_16_Mw8.3/southamerica_2015_09_16_Mw8.3_Patonga.Level1.csv",
            "2016-12-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/solomon_2016_12_08_Mw7.8/solomon_2016_12_08_Mw7.8_Patonga.Level1.csv",
            "2017-09-08" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/mexico_2017_09_08_Mw8.2/mexico_2017_09_08_Mw8.2_Patonga.Level1.csv",
            "2021-02-10" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2021_02_10_Mw7.7/newhebrides_2021_02_10_Mw7.7_Patonga.Level1.csv",
            "2021-03-04" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/kermadec_2021_03_04_Mw8.0/kermadec_2021_03_04_Mw8.0_Patonga.Level1.csv",
            "2022-01-15" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/hungatonga_2022_01_15_volcano/hungatonga_2022_01_15_volcano_Patonga.Level1.csv",
            "2023-05-19" = "./DATA/TIDES/NSW_TIDAL_GAUGES/NSW-MHL-1min-2011-2023/detiding/detided_MHL_1min_data/newhebrides_2023_05_19_Mw7.7/newhebrides_2023_05_19_Mw7.7_Patonga.Level1.csv"
            ),
        readfun = read_tsunami_DPIE
        ),

    'Newcastle_east_1min_PA' = list(
        coord = c(151.788944, -32.923806),
        events = list(
            '2011-03-11' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_tohoku_2011_03_11_detided_WITHOUT_RAW_HEIGHT.csv',
            '2014-04-01' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_chile_2014_04_01_detided_WITHOUT_RAW_HEIGHT.csv',
            '2015-09-16' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_chile_2015_09_16_detided_WITHOUT_RAW_HEIGHT.csv',
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_newHebrides_2021_02_10_detided_WITHOUT_RAW_HEIGHT.csv',
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_kermadec_2021_03_04_detided_WITHOUT_RAW_HEIGHT.csv',
            '2023-05-19' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/newcastle_east_newHebrides_2023_05_19_detided_WITHOUT_RAW_HEIGHT.csv'
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
            '2007-04-01' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2007-04-01.csv',
            '2009-07-15' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2009-07-15.csv',
            '2009-09-29' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2009-09-29.csv',
            '2010-02-27' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2010-02-27.csv',
            '2011-03-11' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2011-03-11.csv',
            '2014-04-01' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2014-04-01.csv',
            '2015-09-16' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/BOM_1minute_gauges/detided_data_near_events/port_kembla_near_2015-09-16.csv',
            # New Hebrides, 2021 Mw 7.7
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/PortKembla.csv',
            # Kermadec-Tonga, 2021 Mw 8.1, same dataset as above
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/IOC_2021_Feb_March/detided_data_near_events/PortKembla.csv'
            ),
        readfun = read_tsunami_BOM_1min
        ),

    #
    # Grain terminal at Port Kembla. This is after they started collecting 1min
    # data. Earlier data was 6min only -- it was used in Allen and Greenslade
    # (2016) but visually doesn't look to sample enough to capture the tsunami.
    #
    'PortKembla_GrainTerminal_1min_PA' = list(
        # From Port Authority email in which they supplied data
        coord = c(150.89358, -34.45558),
        events = list(
            '2021-02-10' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/port_kembla_grain_terminal_newHebrides_2021_02_10_detided_WITHOUT_RAW_HEIGHT.csv',
            '2021-03-04' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/port_kembla_grain_terminal_kermadec_2021_03_04_detided_WITHOUT_RAW_HEIGHT.csv',
            '2023-05-19' = './DATA/TIDES/NSW_TIDAL_GAUGES/NSW_Port_Authority_Purchase_2023_11/detiding/detided_obs_without_original_data/port_kembla_grain_terminal_newHebrides_2023_05_19_detided_WITHOUT_RAW_HEIGHT.csv'
            ),
        readfun = read_tsunami_DPIE
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

    #
    # Hillarys -- Sumatra 2005 
    #
    'Hillarys_BOM_1min_2005' = list(
        coord = c(115 + 44.3/60, -(31 + 49.5/60) ),
        events = list(
            # Sumatra 2005, from a record provided by BOM
            '2005-03-28' = './DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_Hillarys/detided_obs_2005/hillarys_2005_tsunami.csv'
            ),
        readfun = read_tsunami_DPIE # Turns out this one will work.
        ),

    # Hillarys again -- but with a different detided file format
    'Hillarys_BOM_1min' = list(
        coord = c(115 + 44.3/60, -(31 + 49.5/60) ),
        events = list(
            '2021-08-12' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/2021_08_Sandwich/detided_data_near_events/Hillarys.csv"
            ),
        readfun = read_tsunami_DPIE # Turns out this one will work.
        ),

    # Fremantle Harbour
    'Freemantle_WADoT_5min_2004' = list(
        coord = c( 115.748056, -32.065556 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Freemantle.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Freemantle.csv",
            '2021-08-12' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/2021_08_Sandwich/detided_data_near_events/Freemantle.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    # Fremantle Harbour -- there is a second gauge.
    'FreemantleAlternate_WADoT_5min_2004' = list(
        coord = c( 115.748056, -32.065556 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/FreemantleAlternate.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/FreemantleAlternate.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    # Fremantle Inner Harbour -- this gauge has a very long history
    'FremantleInnerHarbour_5min' = list(
        coord=c(115.73950, -32.0541666666),
        events = list(
            '2004-12-26' = "./DATA/TIDES/WA_TIDAL_GAUGES/FremantlePorts/detiding/detided_data/FremantleInner_Sumatra2004.csv",
            '2005-03-28' = "./DATA/TIDES/WA_TIDAL_GAUGES/FremantlePorts/detiding/detided_data/FremantleInner_Sumatra2005.csv"
            ),
        readfun = read_tsunami_DPIE
    ),

    # Mangles Bay
    'ManglesBay_WADoT_5min_2004' = list(
        coord = c( 115.703333, -32.274444 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/ManglesBay.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/ManglesBay.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    # Barrak Street
    'BarrackStreet_WADoT_5min_2004' = list(
        coord = c(115.85725, -31.959694),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/BarrackSt.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/BarrackSt.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

   'MandurahFishermansJetty_WADoT_5min_2004' = list(
        coord = c(115.715278, -32.528611),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/MandurahFJ.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/MandurahFJ.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'PeelInlet_WADoT_5min_2004' = list(
        coord = c(115.713889 , -32.592222),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/PeelInlet.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/PeelInlet.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'CapeBouvard_WADoT_5min_2004' = list(
        coord = c( 115.629722 , -32.601389),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/CapeBouvard.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/CapeBouvard.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Canarvon_WADoT_5min_2004' = list(
        coord = c( 113.65103,  -24.89875),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Canarvon.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Canarvon.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Caddup_WADoT_5min_2004' = list(
        coord = c( 115.643056,  -32.609444),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Caddup.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Caddup.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Harvey_WADoT_5min_2004' = list(
        coord = c( 115.676141, -32.68353),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Harvey.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Harvey.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Bunbury_WADoT_5min_2004' = list(
        coord = c(115.641065, -33.309671),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Bunbury.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Bunbury.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'BunburyInner_WADoT_5min_2004' = list(
        # coord = c(115.659988, -33.323476),
        # Deliberately shift the coord slightly off of "land"
        coord =  c(115.65967834,-33.32290673), 
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/BunburyInner.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/BunburyInner.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'BusseltonPortGeographe_WADoT_5min_2004' = list(
        #coord = c(115.394398, -33.630437),
        # Deliberately shift location to catch gauge in water
        coord = c(115.39472867, -33.63057753),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/BusseltonPortGeographe.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/BusseltonPortGeographe.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Lancelin_WADoT_5min_2004' = list(
        coord = c(115.327222, -31.014444),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Lancelin.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Lancelin.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'JurianBay_WADoT_5min_2004' = list(
        coord = c(115.042778, -30.287222 ),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/JurianBay.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/JurianBay.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'Geraldton_WADoT_5min_2004' = list(
        coord = c( 114.602147, -28.775919),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/Geraldton.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/Geraldton.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),

    'GeraldtonAlternate_WADoT_5min_2004' = list(
        coord = c( 114.602147, -28.775919),
        events = list(
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004_WADOT/detided_data_near_events/GeraldtonAlternate.csv",
            '2005-03-28' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Sumatra_2005_WADOT/detided_data_near_events/GeraldtonAlternate.csv"
            ),
        readfun = read_tsunami_WADoT_5min
        ),
    #
    # Cocos Island -- Andaman 2004
    #
    'CocosIsland_BOM_1min_2004' = list(
        coord = c(96 + 53.5/60, -(12 + 7.0/60) ),
        events = list(
            # Sumatra 2004, from a record provided by BOM
            '2004-12-26' = "./DATA/TIDES/MISC_TIDAL_GAUGES/tide_gauge/Andaman_2004/detided_obs_2004/cocos_island_2004_tsunami.csv"
            ),
        readfun = read_tsunami_BOM2004
        ),

    #
    # Lots of auto-generated data for WA
    #
    "CNCAR02_01_2004_WADOT_5min" = list(
        coord = c(113.65103 ,-24.89875),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/CNCAR02_01/CNCAR02_01_2004.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2004a_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2004a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2004b_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2004b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "EXEXM02_01_2004_WADOT_5min" = list(
        coord = c(114.140889 ,-21.954861),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/EXEXM02_01/EXEXM02_01_2004.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "NWONS01_02_2004_WADOT_5min" = list(
        coord = c(115.13153 ,-21.64967),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/NWONS01_02/NWONS01_02_2004.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "PWLAM01_03_2004_WADOT_5min" = list(
        coord = c(117.185833 ,-20.587778),
        events = list( 
            "2004-12-26" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/PWLAM01_03/PWLAM01_03_2004.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),


    "CNCAR02_01_2005_WADOT_5min" = list(
        coord = c(113.65103 ,-24.89875),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/CNCAR02_01/CNCAR02_01_2005.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2005a_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2005a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2005b_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2005b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "EXEXM02_01_2005_WADOT_5min" = list(
        coord = c(114.140889 ,-21.954861),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/EXEXM02_01/EXEXM02_01_2005.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "NWONS01_02_2005_WADOT_5min" = list(
        coord = c(115.13153 ,-21.64967),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/NWONS01_02/NWONS01_02_2005.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "PWLAM01_03_2005_WADOT_5min" = list(
        coord = c(117.185833 ,-20.587778),
        events = list( 
            "2005-03-28" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/PWLAM01_03/PWLAM01_03_2005.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),

    "BUBNY01_01_2006_WADOT_5min" = list(
        coord = c(115.659988 ,-33.323476),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/BUBNY01_01/BUBNY01_01_2006.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "CNCAR02_01_2006_WADOT_5min" = list(
        coord = c(113.65103 ,-24.89875),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/CNCAR02_01/CNCAR02_01_2006.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2006a_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2006a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2006b_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2006b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "EXEXM02_01_2006_WADOT_5min" = list(
        coord = c(114.140889 ,-21.954861),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/EXEXM02_01/EXEXM02_01_2006.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "GNGER02_02_2006a_WADOT_5min" = list(
        coord = c(114.602147 ,-28.775919),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/GNGER02_02/GNGER02_02_2006a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "GNGER02_02_2006b_WADOT_5min" = list(
        coord = c(114.602147 ,-28.775919),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/GNGER02_02/GNGER02_02_2006b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "NWONS01_02_2006_WADOT_5min" = list(
        coord = c(115.13153 ,-21.64967),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/NWONS01_02/NWONS01_02_2006.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "PWLAM01_03_2006_WADOT_5min" = list(
        coord = c(117.185833 ,-20.587778),
        events = list( 
            "2006-07-17" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/PWLAM01_03/PWLAM01_03_2006.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),

    "BSGEO01_01_2007_WADOT_5min" = list(
        coord = c(115.394398 ,-33.630437),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/BSGEO01_01/BSGEO01_01_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "CNCAR02_01_2007_WADOT_5min" = list(
        coord = c(113.65103 ,-24.89875),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/CNCAR02_01/CNCAR02_01_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2007a_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2007a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "DPKBY01_03_2007b_WADOT_5min" = list(
        coord = c(116.748889 ,-20.623611),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/DPKBY01_03/DPKBY01_03_2007b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "ESESP03_01_2007_WADOT_5min" = list(
        coord = c(121.897162 ,-33.870692),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/ESESP03_01/ESESP03_01_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "EXEXM02_01_2007_WADOT_5min" = list(
        coord = c(114.140889 ,-21.954861),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/EXEXM02_01/EXEXM02_01_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "GNGER02_02_2007a_WADOT_5min" = list(
        coord = c(114.602147 ,-28.775919),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/GNGER02_02/GNGER02_02_2007a.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "GNGER02_02_2007b_WADOT_5min" = list(
        coord = c(114.602147 ,-28.775919),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/GNGER02_02/GNGER02_02_2007b.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "NWONS01_02_2007_WADOT_5min" = list(
        coord = c(115.13153 ,-21.64967),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/NWONS01_02/NWONS01_02_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),
    "PWLAM01_03_2007_WADOT_5min" = list(
        coord = c(117.185833 ,-20.587778),
        events = list( 
            "2007-09-12" = "./DATA/TIDES/WA_TIDAL_GAUGES/revised_detiding_2024/processed_data/PWLAM01_03/PWLAM01_03_2007.txt"
            ), 
        readfun = read_WADOT_revised_2024_detiding 
        ),


    'NOTHING' = list(
        coord = c(NA, NA),
        events = list(),
        readfun = NULL
        )

    )

# Store the directory, so if we source this from another location we can still find the gauges
BASE_DIR = getwd()

# Print a warning if not all the files are found -- this helps catch errors when adding files
# or porting the code elsewhere.
all_files_exist = lapply(GAUGE_DATA, function(x) all(unlist(lapply(x$events, file.exists))))
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
    matching_events = lapply(GAUGE_DATA, function(x) which(names(x$events) == date))
    gauges_with_data = which(unlist(lapply(matching_events, function(x) length(x) > 0)))

    # Sanity check: we should never have more than 1 gauge at a site (unless
    # I've entered the dates incorrectly!).
    if(any(unlist(lapply(matching_events, function(x) length(x) > 1)))){
        stop('Error: Multiple gauges detected at a single site. This should not happen')
    }

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

# Dump the gauge names and coordinates to the csv file gauge_coords.csv
export_coordinates<-function(){
    all_coords = do.call(rbind, lapply(GAUGE_DATA, function(x) x$coord))
    colnames(all_coords) = c('lon', 'lat')

    site_name = rep(NA, length(GAUGE_DATA))
    for(i in 1:length(site_name)) site_name[i] = match_site_to_station_name[[names(GAUGE_DATA)[i]]]

    all_events = unlist(lapply(GAUGE_DATA, function(x) paste0(names(x$events), collapse=" & ")))
    stopifnot(length(all_events) == nrow(all_coords))

    all_coords = cbind(all_coords, data.frame(site_name = site_name, site_data_list_id=rownames(all_coords), events=all_events))

    data_provided_by = sapply(all_coords$site_data_list_id, match_data_source_to_name, USE.NAMES=FALSE)
    all_coords = cbind(all_coords, data.frame(data_provided_by=data_provided_by))

    write.csv(all_coords, file='gauge_coords.csv', row.names=FALSE, quote=TRUE)
}

match_data_source_to_name<-function(site_name){
    # Here site_name is an entry of names(GAUGE_DATA)
    if(grepl('WADOT', site_name) | grepl('WADoT', site_name)){
        return('WA Department of Transport')
    }else if(grepl('DPIE', site_name)){
        return('NSW Department of Planning, Industry and Environment')
    }else if(grepl('BOM', site_name)){
        return('Bureau of Meteorology')
    }else if(grepl('_PA', site_name)){
        return('Port Authority of NSW')
    }else if(grepl('FremantleInnerHarbour_5min', site_name)){
        return("Fremantle Ports")
    }else if(grepl("GoldCoastSandBypass", site_name)){
        return('Bureau of Meteorology')
    }else if(grepl("Sydney_FortDenison1960", site_name)){
        return('Wilson et al. (2018) https://doi.org/10.1038/s41598-018-33156-w')
    }else if(grepl("Cronulla_CSIRO_Fisheries_1960", site_name)){
        return('Davies et al. (2020) https://www.frontiersin.org/article/10.3389/feart.2020.598235')
    }else{
        return(NA)
    }
}

match_site_to_station_name = list(
    "BallinaBreakwall_1min_DPIE" = 'Ballina Breakwall',
    "Bermagui_1min_DPIE" = 'Bermagui',                    
    "BrunswickHeads_1min_DPIE" = 'Brunswick Heads',               
    "CoffsHarbourInnerPumpoutJetty_1min_DPIE" = 'Coffs Harbour Inner Pumpout Jetty',
    "CrookhavenHeads_1min_DPIE" = 'Crookhaven Heads',              
    "CrowdyHeadFishermansWharf_1min_DPIE" = 'Crowdy Head Fishermans Wharf',    
    "Forster_1min_DPIE" = 'Forster',                      
    "LordHoweIsland_1min_DPIE" = 'Lord Howe Island',               
    "PortMacquarie_1min_DPIE" = 'Port Macquarie',                 
    "ShoalBay_1min_DPIE" = 'Shoal Bay',                     
    "TweedEntranceSouth_1min_DPIE" = 'Tweed Entrance South',           
    "Yamba_1min_DPIE" = 'Yamba',                        
    "Eden_1min_DPIE" = 'Eden',                          
    "TwofoldBay_1min_BOM" = 'Twofold Bay',                    
    "GoldCoastSandBypass" = 'Gold Coast Sand Bypass',                    
    "BatemansBay_PrincessJetty_1min_DPIE" = 'Batemans Bay Princess Jetty',    
    "Ulladulla_1min_DPIE" = 'Ulladulla',                    
    "JervisBay_1min_DPIE" = 'Jervis Bay',                    
    "Bundeena_1min_DPIE" = 'Bundeena',                     
    "Sydney_MiddleHarbour_1min_DPIE" = 'Sydney Middle Harbour',         
    "Sydney_FortDenison1960" = 'Sydney Fort Denison',                 
    "Cronulla_CSIRO_Fisheries_1960" = 'Cronulla',          
    "Sydney_FortDenison_1min_PA" = 'Sydney Fort Denison',             
    "Sydney_FortDenison_1min_PA_b" = 'Sydney Fort Denison',           
    "Sydney_BotanyBay_1min_PA" = 'Botany Bay Bulk Liquids Berth',                
    "Sydney_Botany_Bay_Pilot_Jetty_1min_PA" = 'Botany Bay Pilot Jetty', 
    "Sydney_Botany_Bay_Kurnell_Wharf_1min_PA" = 'Botany Bay Kurnell Wharf',
    "Hawkesbury_Patonga_1min_DPIE" = 'Patonga',          
    "Newcastle_east_1min_PA"  = 'Newcastle East',                
    "PortKembla_BOM_1min_2004" = 'Port Kembla',                
    "PortKembla_1min_BOM" = 'Port Kembla',                   
    "PortKembla_GrainTerminal_1min_PA"  = 'Port Kembla Grain Terminal',      
    "Portland_BOM_1min_2004"  = 'Portland',                
    "Hillarys_BOM_1min_2004" = 'Hillarys',                
    "Hillarys_BOM_1min_2005" = 'Hillarys',                
    "Hillarys_BOM_1min" = 'Hillarys',                     
    "Freemantle_WADoT_5min_2004" = 'Fremantle Boat Harbour',            
    "FreemantleAlternate_WADoT_5min_2004" = 'Fremantle Boat Harbour',   
    "FremantleInnerHarbour_5min" = 'Fremantle Inner Harbour',            
    "ManglesBay_WADoT_5min_2004" = 'Mangles Bay',            
    "BarrackStreet_WADoT_5min_2004" = 'Barrack Street',         
    "MandurahFishermansJetty_WADoT_5min_2004" = "Mandurah Fishermans Jetty",
    "PeelInlet_WADoT_5min_2004" = 'Peel Inlet',             
    "CapeBouvard_WADoT_5min_2004" = 'Cape Bouvard',           
    "Canarvon_WADoT_5min_2004" = 'Carnarvon',              
    "Caddup_WADoT_5min_2004" = 'Caddup',                
    "Harvey_WADoT_5min_2004" = 'Harvey',                
    "Bunbury_WADoT_5min_2004" = 'Bunbury',                
    "BunburyInner_WADoT_5min_2004" = 'Bunbury Inner Harbour',           
    "BusseltonPortGeographe_WADoT_5min_2004" = 'Port Geographe', 
    "Lancelin_WADoT_5min_2004" = 'Lancelin',               
    "JurianBay_WADoT_5min_2004" = 'Jurien Bay',             
    "Geraldton_WADoT_5min_2004" = 'Geraldton',             
    "GeraldtonAlternate_WADoT_5min_2004" = 'Geraldton',     
    "CocosIsland_BOM_1min_2004" = 'Cocos Island',             
    "CNCAR02_01_2004_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2004a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2004b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2004_WADOT_5min" = 'Exmouth',            
    "NWONS01_02_2004_WADOT_5min" = 'Onslow Beadon Creek',
    "PWLAM01_03_2004_WADOT_5min" = 'Cape Lambert',            
    "CNCAR02_01_2005_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2005a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2005b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2005_WADOT_5min" = 'Exmouth',            
    "NWONS01_02_2005_WADOT_5min" = 'Onslow Beadon Creek',            
    "PWLAM01_03_2005_WADOT_5min" = 'Cape Lambert',            
    "BUBNY01_01_2006_WADOT_5min" = 'Bunbury Inner Harbour',            
    "CNCAR02_01_2006_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2006a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2006b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2006_WADOT_5min" = 'Exmouth',            
    "GNGER02_02_2006a_WADOT_5min" = 'Geraldton',           
    "GNGER02_02_2006b_WADOT_5min" = 'Geraldton',           
    "NWONS01_02_2006_WADOT_5min"  = 'Onslow Beadon Creek',           
    "PWLAM01_03_2006_WADOT_5min" = 'Cape Lambert',            
    "BSGEO01_01_2007_WADOT_5min" = 'Port Geographe',            
    "CNCAR02_01_2007_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2007a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2007b_WADOT_5min" = 'King Bay',           
    "ESESP03_01_2007_WADOT_5min" = 'Esperance',            
    "EXEXM02_01_2007_WADOT_5min" = 'Exmouth',            
    "GNGER02_02_2007a_WADOT_5min" = 'Geraldton',           
    "GNGER02_02_2007b_WADOT_5min" = 'Geraldton',           
    "NWONS01_02_2007_WADOT_5min" = 'Onslow Beadon Creek',            
    "PWLAM01_03_2007_WADOT_5min" = 'Cape Lambert',
    "NOTHING" = 'Deliberate blank entry'
)
