#
# Store the DART files we use for testing
#

if(!exists('variable_mu')){
    stop('variable_mu must be defined (either TRUE or FALSE) before sourcing config_DART_test_files.R')
}

#
# Get a vector with the Rdata files produced by running
# ../SOURCE_ZONES/sourcezone/TSUNAMI_EVENTS/plots/gauge_summary_statistics.R
# for each of the historical comparison events
#
if(variable_mu){
    #all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*varyMu.Rdata')

    #
    # Safest to hard-code the paths we need to use.
    # Note that events are skipped if they are A) Not thrust; B) Mw < 7.7 in GCMT, C) On a source-zone that was updated to use SLAB2
    #

    all_Rdata = c(
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1_varyMu.Rdata",
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0_varyMu.Rdata",
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7_varyMu.Rdata",
        "../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1_varyMu.Rdata",
        "../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0_varyMu.Rdata",
        #"../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7_varyMu.Rdata", 
        "../SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3_varyMu.Rdata", 
        "../SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1_varyMu.Rdata", 
        #"../SOURCE_ZONES/newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
        #"../SOURCE_ZONES/newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
        "../SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
        "../SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
        #"../SOURCE_ZONES/outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76_varyMu.Rdata",
        #"../SOURCE_ZONES/puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
        "../SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
        #"../SOURCE_ZONES/solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81_varyMu.Rdata", 
        #"../SOURCE_ZONES/solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78_varyMu.Rdata",
        "../SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81_varyMu.Rdata", 
        "../SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78_varyMu.Rdata",
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8_varyMu.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80_varyMu.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8_varyMu.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2_varyMu.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3_varyMu.Rdata",
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8_varyMu.Rdata",
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9_varyMu.Rdata",
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5_varyMu.Rdata", 
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8_varyMu.Rdata",
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9_varyMu.Rdata", 
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5_varyMu.Rdata",
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8_varyMu.Rdata"
        ) 

}else{

    #all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*[0-9].Rdata')

    #
    # Safest to hard-code the paths we need to use.
    # Note that events are skipped if they are A) Not thrust; B) Mw < 7.7 in GCMT, C) On a source-zone that was updated to use SLAB2
    #

    all_Rdata = c(
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1.Rdata",
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0.Rdata",
        #"../SOURCE_ZONES/kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7.Rdata",
        "../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata",
        "../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0.Rdata",
        #"../SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7.Rdata", 
        "../SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3.Rdata", 
        "../SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata", 
        #"../SOURCE_ZONES/newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
        #"../SOURCE_ZONES/newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
        "../SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
        "../SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
        #"../SOURCE_ZONES/outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76.Rdata",
        #"../SOURCE_ZONES/puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8.Rdata", 
        "../SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata", 
        #"../SOURCE_ZONES/solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81.Rdata", 
        #"../SOURCE_ZONES/solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78.Rdata",
        "../SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata", 
        "../SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78.Rdata",
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata", 
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata",
        "../SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8.Rdata",
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9.Rdata",
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5.Rdata", 
        #"../SOURCE_ZONES/sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8.Rdata",
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9.Rdata", 
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5.Rdata",
        "../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8.Rdata"
        ) 

}

