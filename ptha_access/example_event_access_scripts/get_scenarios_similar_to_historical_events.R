#
# Code to get 'good fitting' scenario for a bunch of historical events from the
# PTHA18 results {see GA Record for explanation}
#

# See ptha/ptha_access/README.md for a tutorial on using these functions

# Get the ptha_access functions (from the ptha repository)
source('../get_PTHA_results.R', chdir=TRUE)

#
# The following data structure gives the scenario row indices that look most
# like events of interest. It was made using the 'find_desired_event_rows.R'
# code. 
# Where the source-zone has been updated, the old versions is commented out.
#

scenarios_like_events = list(
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(24481, 24555, 24559)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1_varyMu.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(24559, 26445, 26446)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(20848, 20874, 22838)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0_varyMu.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(22838, 22862, 24812)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(16983, 16987, 17032)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7_varyMu.Rdata", 
    #   source_zone = "kermadectonga", desired_event_rows = c(10106, 16983, 16987)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(24563, 28606, 28678)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1_varyMu.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(30715, 32759, 32783)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(22713, 22745, 26925)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0_varyMu.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(29004, 30983, 31008)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(13653, 13663, 16946)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7_varyMu.Rdata", 
        source_zone = "kermadectonga2", desired_event_rows = c(20033, 22852, 22941)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3.Rdata", 
        source_zone = "kurilsjapan", desired_event_rows = c(35025, 37356, 37490)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3_varyMu.Rdata", 
        source_zone = "kurilsjapan", desired_event_rows = c(35025, 37356, 37490)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata", 
        source_zone = "kurilsjapan", desired_event_rows = c(46994, 47004, 47634)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1_varyMu.Rdata", 
        source_zone = "kurilsjapan", desired_event_rows = c(46994, 47004, 47634)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
    #   source_zone = "newhebrides", desired_event_rows = c(5504, 5976, 6446)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
    #   source_zone = "newhebrides", desired_event_rows = c(5976, 6906, 6913)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
    #   source_zone = "newhebrides", desired_event_rows = c(5933, 5964, 5967)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
    #   source_zone = "newhebrides", desired_event_rows = c(5933, 6420, 6429)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
        source_zone = "newhebrides2", desired_event_rows = c(7216, 7223, 7768)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
        source_zone = "newhebrides2", desired_event_rows = c(8821, 8831, 9335)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
        source_zone = "newhebrides2", desired_event_rows = c(6066, 7063, 7084)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
        source_zone = "newhebrides2", desired_event_rows = c(6066, 7063, 7084)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76.Rdata", 
        source_zone = "outerrise_kermadectonga", desired_event_rows = c(3723, 4511, 4527)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    ## This one is removed because there is no 'variable-shear-modulus' for
    ## normal fault events. The computations used the thrust-event profile
    ## (since they were automated), but we do not want that
    #structure(list(stats_file = "outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76_varyMu.Rdata", 
    #    source_zone = "outerrise_kermadectonga", desired_event_rows = c(3723, 
    #4511, 4527)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8.Rdata", 
    #    source_zone = "puysegur", desired_event_rows = c(2351, 2689, 
    #2695)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
    #    source_zone = "puysegur", desired_event_rows = c(2102, 2108, 
    #2351)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata", 
        source_zone = "puysegur2", desired_event_rows = c(1567, 1579, 1782)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
        source_zone = "puysegur2", desired_event_rows = c(1567, 1579, 1782)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81.Rdata", 
    #    source_zone = "solomon", desired_event_rows = c(11671, 12275, 
    #12285)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81_varyMu.Rdata", 
    #    source_zone = "solomon", desired_event_rows = c(10999, 11020, 
    #12291)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78.Rdata", 
    #    source_zone = "solomon", desired_event_rows = c(9498, 9543, 
    #9589)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78_varyMu.Rdata", 
    #    source_zone = "solomon", desired_event_rows = c(8195, 9543, 
    #9589)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata", 
        source_zone = "solomon2", desired_event_rows = c(10634, 10637, 11216)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81_varyMu.Rdata", 
        source_zone = "solomon2", desired_event_rows = c(10607, 10637, 11216)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78.Rdata", 
        source_zone = "solomon2", desired_event_rows = c(7549, 8739, 8757)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78_varyMu.Rdata", 
        source_zone = "solomon2", desired_event_rows = c(8739, 9822, 9831)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(128450, 128607, 128608)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(128450, 128607, 128608)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(81177, 81215, 81269)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(81177, 81215, 81269)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(61720, 71520, 71579)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(53117, 53123, 53131)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(87595, 87689, 87719)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(80289, 80343, 80427)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(108269, 108275, 108314)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(108269, 108275, 108314)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(65317, 75224, 75333)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8_varyMu.Rdata", 
        source_zone = "southamerica", desired_event_rows = c(55822, 55830, 65383)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(57154, 57278, 57310)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9_varyMu.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(57278, 57310, 62677)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(78610, 78650, 78731)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5_varyMu.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(73358, 73474, 78610)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(39777, 39805, 39813)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8_varyMu.Rdata", 
    #    source_zone = "sunda", desired_event_rows = c(39777, 39805, 39813)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(53116, 59402, 59403)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9_varyMu.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(64616, 64665, 69873)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(80013, 80152, 80242)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5_varyMu.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(80152, 80242, 89672)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(46769, 46844, 53856)), .Names = c("stats_file", "source_zone", "desired_event_rows")), 
    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8_varyMu.Rdata", 
        source_zone = "sunda2", desired_event_rows = c(46769, 46844, 53856)), .Names = c("stats_file", "source_zone", "desired_event_rows"))
    )



#
# In this commented out code, I accidently made a list with the 'least similar' events! Kept because it might be quasi-informative.
#
##scenarios_like_events = list(
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(20625, 
##    #    22629, 22637)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1_varyMu.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(22637, 
##    #    24607, 26491)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(18905, 
##    #    20946, 22965)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2006_05_03_Mw8.0_varyMu.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(22965, 
##    #    24929, 27648)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(12845, 
##    #    14894, 16989)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7_varyMu.Rdata", 
##    #    source_zone = "kermadectonga", desired_event_rows = c(15016, 
##    #    15058, 16935)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(26562, 
##        26624, 26648)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_samoa_2009_09_29_Mw8.1_varyMu.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(22391, 
##        26562, 26648)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(22823, 
##        26945, 26973)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2006_05_03_Mw8.0_varyMu.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(26945, 
##        26973, 31068)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(13639, 
##        13646, 16845)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kermadectonga2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kermadectonga2_tonga_2009_03_19_Mw7.7_varyMu.Rdata", 
##        source_zone = "kermadectonga2", desired_event_rows = c(13836, 
##        20164, 22848)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3.Rdata", 
##        source_zone = "kurilsjapan", desired_event_rows = c(37366, 
##        37639, 37823)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_kuril_2006_11_15_Mw8.3_varyMu.Rdata", 
##        source_zone = "kurilsjapan", desired_event_rows = c(27299, 
##        37316, 37366)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata", 
##        source_zone = "kurilsjapan", desired_event_rows = c(46859, 
##        46938, 48374)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "kurilsjapan/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1_varyMu.Rdata", 
##        source_zone = "kurilsjapan", desired_event_rows = c(46859, 
##        46938, 47532)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
##    #    source_zone = "newhebrides", desired_event_rows = c(6393, 
##    #    6403, 6416)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
##    #    source_zone = "newhebrides", desired_event_rows = c(6393, 
##    #    6403, 6863)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
##    #    source_zone = "newhebrides", desired_event_rows = c(5884, 
##    #    5886, 5895)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "newhebrides/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
##    #    source_zone = "newhebrides", desired_event_rows = c(5884, 
##    #    5886, 6906)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9.Rdata", 
##        source_zone = "newhebrides2", desired_event_rows = c(7760, 
##        8274, 8292)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_northNewHebrides_2013_02_06_Mw7.9_varyMu.Rdata", 
##        source_zone = "newhebrides2", desired_event_rows = c(7760, 
##        8292, 8797)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8.Rdata", 
##        source_zone = "newhebrides2", desired_event_rows = c(7024, 
##        7028, 7031)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "newhebrides2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_newhebrides2_vanuatu_north_2009_10_07_Mw7.8_varyMu.Rdata", 
##        source_zone = "newhebrides2", desired_event_rows = c(7028, 
##        8195)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76.Rdata", 
##        source_zone = "outerrise_kermadectonga", desired_event_rows = c(4492, 
##        4496, 4497)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    ## NOTE: Variable shear modulus is not relevant for these outer-rise
##    ## events. This code would have applied subduction like variability, which
##    ## is not what we want. So do not include this event.
##    #structure(list(stats_file = "outerrise_kermadectonga/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_outerrise_kermadectonga_outerrise_kermadectonga_2011_07_06_Mw76_varyMu.Rdata", 
##    #    source_zone = "outerrise_kermadectonga", desired_event_rows = c(4492, 
##    #    4496, 4497)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8.Rdata", 
##    #    source_zone = "puysegur", desired_event_rows = c(2677, 2681, 
##    #    2710)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "puysegur/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
##    #    source_zone = "puysegur", desired_event_rows = c(3039, 3045, 
##    #    3454)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    structure(list(stats_file = "puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8.Rdata", 
##        source_zone = "puysegur2", desired_event_rows = c(1578, 1783, 
##        1976)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "puysegur2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_puysegur2_puysegur_2009_07_15_Mw7.8_varyMu.Rdata", 
##        source_zone = "puysegur2", desired_event_rows = c(1973, 2203, 
##        2229)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81.Rdata", 
##    #    source_zone = "solomon", desired_event_rows = c(11681, 12200, 
##    #    12208)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2007_04_01_Mw81_varyMu.Rdata", 
##    #    source_zone = "solomon", desired_event_rows = c(11681, 14069, 
##    #    14075)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78.Rdata", 
##    #    source_zone = "solomon", desired_event_rows = c(9612, 9634, 
##    #    9660)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "solomon/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon_solomon_2016_12_08_Mw78_varyMu.Rdata", 
##    #    source_zone = "solomon", desired_event_rows = c(8316, 9612, 
##    #    11502)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81.Rdata", 
##    #    source_zone = "solomon2", desired_event_rows = c(10531, 11175, 
##    #    11228)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2007_04_01_Mw81_varyMu.Rdata", 
##        source_zone = "solomon2", desired_event_rows = c(10531, 12256, 
##        12810)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78.Rdata", 
##        source_zone = "solomon2", desired_event_rows = c(8785, 8801, 
##        8806)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "solomon2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_solomon2_solomon_2016_12_08_Mw78_varyMu.Rdata", 
##        source_zone = "solomon2", desired_event_rows = c(8761, 8797, 
##        8801)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(123569, 
##        128274, 128277)), .Names = c("stats_file", "source_zone", 
##    "desired_event_rows")), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(114143, 
##        118877, 118940)), .Names = c("stats_file", "source_zone", 
##    "desired_event_rows")), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(81303, 
##        88702, 88709)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_08_15_Mw80_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(81303, 
##        88687, 95997)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(53306, 
##        53348, 53353)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2007_11_14_Mw7.8_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(53302, 
##        53348, 53353)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(87613, 
##        87647, 87654)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(87490, 
##        87492, 109257)), .Names = c("stats_file", "source_zone", 
##    "desired_event_rows")), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(93675, 
##        100958, 101357)), .Names = c("stats_file", "source_zone", 
##    "desired_event_rows")), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(79447, 
##        86376, 86779)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(55853, 
##        55855, 55856)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "southamerica/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_southamerica_southamerica_2016_04_16_Mw7.8_varyMu.Rdata", 
##        source_zone = "southamerica", desired_event_rows = c(82742, 
##        82795, 90044)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(57089, 57121, 
##    #    57122)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9_varyMu.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(57089, 57121, 
##    #    57122)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(78938, 83405, 
##    #    86949)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5_varyMu.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(78938, 83405, 
##    #    86949)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(39890, 52643, 
##    #    52650)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    #structure(list(stats_file = "sunda/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8_varyMu.Rdata", 
##    #    source_zone = "sunda", desired_event_rows = c(52643, 58001, 
##    #    63443)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    #)), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(52816, 52888, 
##        59188)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_mentawai_2010_10_25_Mw7.9_varyMu.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(64438, 64439, 
##        64474)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(85514, 85520, 
##        85531)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2007_09_12_Mw8.5_varyMu.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(85520, 85531, 
##        93265)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(39922, 39929, 
##        54025)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )), 
##    structure(list(stats_file = "sunda2/TSUNAMI_EVENTS/plots/gauge_summary_stats_session_sunda2_sumatra_2010_04_06_Mw7.8_varyMu.Rdata", 
##        source_zone = "sunda2", desired_event_rows = c(54018, 54025, 
##        60131)), .Names = c("stats_file", "source_zone", "desired_event_rows"
##    )))


source_zone_data = vector(mode='list', length=length(scenarios_like_events))
for(i in 1:length(scenarios_like_events)){

    # Get the scenario metadata for the desired events
    source_zone = scenarios_like_events[[i]]$source_zone
    desired_events = scenarios_like_events[[i]]$desired_event_rows
    source_zone_data[[i]] = get_source_zone_events_data(source_zone, desired_event_rows = desired_events)
    deformation_rasters = vector(mode='list', length=length(desired_events))

    # Get the deformation rasters for the desired_events
    for(j in 1:length(desired_events)){
        deformation_rasters[[j]] = get_initial_condition_for_event(
            source_zone_data[[i]], j)
    }

    # For convenience, store the rasters inside the source_zone_data
    source_zone_data[[i]]$deformation_rasters = deformation_rasters
}

# Save the image so we can easily investigate again later
save.image('scenarios_like_historic_events.Rdata')

# Make plots
pdf('scenario_initial_conditions.pdf', width=10, height=10)
for(i in 1:length(scenarios_like_events)){

    event_title = basename(scenarios_like_events[[i]]$stats_file)

    for(j in 1:length(source_zone_data[[i]]$deformation_rasters)){
        par(oma=c(0,0, 2, 0))
        plot(source_zone_data[[i]]$deformation_rasters[[j]])
        title(event_title, outer=TRUE, line=-3)
    }
}
dev.off()

