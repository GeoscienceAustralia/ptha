event_stats = readRDS('event_stats.RDS')

# Figure out durations related to a few tide gauges that were truncated due to a seiche 

#
# New Hebrides 2021
#
k = which(event_stats$event_name == 'newhebrides2021')[1]
# Match definition in manually_despike_data
site = c('Sydney_FortDenison_1min_PA_b', 'Sydney_MiddleHarbour_1min_DPIE', "Hawkesbury_Patonga_1min_DPIE", "Sydney_Botany_Bay_Pilot_Jetty_1min_PA")
truncated_gauges_time = as.numeric( 24*(c(18670.1, 18670.4, 18670.5, 18670.25) - event_stats$model_start[k]) )
# truncated_gauges_time
# [1] 37.06806 44.26806 46.66806 40.66806

#
# Kermadec 2021 -- we use the same truncation everywhere
#
k = which(event_stats$event_name == 'kermadec2021' & event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
#> summary(as.numeric(event_stats$data_end[k] - event_stats$model_start[k]))*24
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  33.31   33.31   33.31   33.31   33.31   33.31 

