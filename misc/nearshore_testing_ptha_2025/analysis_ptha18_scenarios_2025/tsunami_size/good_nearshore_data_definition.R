# Determine whether a tide gauge has "good nearshore data".
# For most analyses we only want to use good nearshore data, among other
# restrictions (e.g. model gauge close to data gauge, both in high-resolution
# model zone, etc).
is_good_nearshore_data<-function(event_stats){
    #
    # We won't want to use all of the data. We should skip;
    # - "double-ups (2 gauges @ 1 site)"
    # - "clearly erronious data"
    # - "very small data" (where the signal is more likely noise than tsunami)
    # - data with infrequent sampling (e.g. 15 minute in NSW)
    #
    # For the nearshore data, use 1-minute data -- or 5 minute data (for WA). Except
    # for Chile 1960 where we take whatever we have (digitized scans of old tidal
    # gauges). 
    #
    # See comments below to explain why gauges were removed

    # Skip gauges where both "data max" and "abs(data min)" are below this threshold
    # This choice will remove Peel Inlet and Harvey from Sumatra 2005, and Peel Inlet from
    # Sumatra 2004. While the models work fine at these sites too, the datasets have clear artefacts
    # due to the number of decimal places in the data. Since the waves are tiny anyway it seems
    # reasonable to drop these sites.
    DATA_SIZE_THRESHOLD = 0.025 #1e-02

    good_nearshore_data = (
        # One minute, or 5 minute, or anything for Chile 1960
        ( grepl('_1min_', event_stats$sites) | endsWith(event_stats$sites, '_1min') | 
          grepl('_5min_', event_stats$sites) | endsWith(event_stats$sites, '_5min') | 
          (event_stats$event_name == 'chile1960') ) & 

        # Twofold 1min data for Puysegur 2009 is clearly noisy, although it also captures the tsunami somewhat.
        # Anyway now we have good data nearby at Eden from MHL, so remove this.
        ( !grepl("puysegur2009_TwofoldBay_1min_BOM", event_stats$site_and_event) ) &

        # For Sumatra 2004, there are 2 gauges at Fremantle Boat harbour which give quite different results.
        # They are in the same location (according to the metadata) and the gauges give similar results for the
        # subsequent 2005 Sumatra tsunami. (? Maybe Sumatra 2004 made operators realise the gauges needed maintenance ?).
        # It isn't clear that either record is 'better'. One has a more extreme maxima, the other a more extreme minima.
        # I suspect both have problems. 
        ( !grepl("sumatra2004_Freemantle_WADoT_5min_2004", event_stats$site_and_event) ) &
        ( !grepl("sumatra2004_FreemantleAlternate_WADoT_5min_2004", event_stats$site_and_event) ) &
        # Two gauges at Geraldton, one has obvious problems, the other looks good.
        ( !grepl("sumatra2004_Geraldton_WADoT_5min_2004", event_stats$site_and_event)) &
        # Two gauges at King Bay, just keep one
        ( !grepl("sumatra2004_DPKBY01_03_2004a_WADOT_5min", event_stats$site_and_event)) &

        # For Sumatra 2005, we have two basically identical records for Fremantle boat harbour and Geraldton.
        # Change it so we only have one at each site. 
        ( !grepl("sumatra2005_FreemantleAlternate_WADoT_5min_2004", event_stats$site_and_event) ) &
        ( !grepl("sumatra2005_Geraldton_WADoT_5min_2004"          , event_stats$site_and_event) ) &
        # Two gauges at King Bay, just keep one
        ( !grepl("sumatra2005_DPKBY01_03_2005b_WADOT_5min", event_stats$site_and_event)) &

        # Java 2006, we have a double-up at King Bay (Dampier), only keep one
        ( !grepl("java2006_DPKBY01_03_2006b_WADOT_5min", event_stats$site_and_event)) &
        # Double up at Geraldton
        ( !grepl("java2006_GNGER02_02_2006b_WADOT_5min", event_stats$site_and_event)) &

        # Sumatra 2007, we have a double up at King Bay (Dampier) 
        ( !grepl("sumatra2007_DPKBY01_03_2007b_WADOT_5min", event_stats$site_and_event)) &
        # Double up at Geraldton
        ( !grepl("sumatra2007_GNGER02_02_2007b_WADOT_5min", event_stats$site_and_event)) &

        # Sumatra 2010 wasn't modelled as it's hard to see at WA tide-gauges

        # Some gauges are problematic for southamerica2014
        ( !grepl("southamerica2014_TwofoldBay_1min_BOM", event_stats$site_and_event) ) & # Noisy
        # The following gauge could probably be cleaned and then included, but not so easy to clean.
        ( !grepl("southamerica2014_Sydney_FortDenison_1min_PA", event_stats$site_and_event) ) & 

        # This one has 'staircasing' which might be fixed with downsampling, although it's not too bad.
        #( !grepl("southamerica2014_PortKembla_1min_BOM", event_stats$site_and_event) ) & 

        # For southamerica 2015 
        ( !grepl("southamerica2015_Eden_1min_DPIE", event_stats$site_and_event) ) &  # Noisy and missing

        # For newhebrides 2021
        ( !grepl("newhebrides2021_TwofoldBay_1min_BOM", event_stats$site_and_event)) & # Noisy
        # In Sydney and Hawkesbury, we also remove late times (seiche) in manually_despike_data.R

        # For Kermadec 2021
        # In this case we remove later-time waves at several gauges in manually_despike_data.R 

        # For south sandwich 2021 the waves in WA are quite short period -- only the 1min data at Hillarys seems to capture it.
        ( !(grepl('sandwich2021', event_stats$site_and_event) & !grepl('Hillarys', event_stats$site_and_event) )) &

        # Skip records where the data is very small
        (! (pmax(event_stats$data_max, -event_stats$data_min) < DATA_SIZE_THRESHOLD) )

        )

}
