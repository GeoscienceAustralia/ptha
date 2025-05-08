#
# For plotting it's nice to show the allowed along-strike range of the
# centroids, but I didn't store it earlier, so this script extracts it.
#

library(rptha)
source('set_range_of_mw_and_centroid/compare_with_data_environment.R', chdir=TRUE)
all_events = read.csv('target_scenarios_data_frame.csv')

alongstrike_lower_ind = rep(NA, nrow(all_events))
alongstrike_upper_ind = rep(NA, nrow(all_events))
for(i in 1:nrow(all_events)){

    tsunami_event_dir = all_events$tsunami_event_dir[i]
    source_zone_env = make_scenario_selection_environment(tsunami_event_dir)

    unit_source_geometry = source_zone_env$unit_source_geometry
    unit_source_statistics = source_zone_env$unit_source_statistics

    event_point1 = c(all_events$event_point_1_lon[i], all_events$event_point_1_lat[i])
    event_point2 = c(all_events$event_point_2_lon[i], all_events$event_point_2_lat[i])
    
    i1 = find_unit_source_index_containing_point(event_point1, unit_source_geometry, unit_source_statistics)
    i2 = find_unit_source_index_containing_point(event_point2, unit_source_geometry, unit_source_statistics)
    alongstrike_lower_ind[i] = min(unit_source_statistics$alongstrike_number[c(i1, i2)])
    alongstrike_upper_ind[i] = max(unit_source_statistics$alongstrike_number[c(i1, i2)])

}

all_events = cbind(all_events, 
    data.frame(alongstrike_lower_ind = alongstrike_lower_ind, 
               alongstrike_upper_ind=alongstrike_upper_ind))

write.csv(all_events, 'target_scenarios_data_frame_with_alongstrike_index.csv', row.names=FALSE)
