#
# Extract spatial coordinates of "good nearshore" gauges in the highres domains
# that are close enough to our points of interest to be included. 
#

event_stats = readRDS('event_stats.RDS')

# Find sites that are actually used in our statistics
k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain & event_stats$distance_to_gauge < 200)
event_stats_k = event_stats[k,]

#
# Get each "site name" only once. This does not completely avoid double-ups
# because often we use different names for the same site at different times.
#
unique_gauges = unique(event_stats_k$sites)
N = match(unique_gauges, event_stats_k$sites)
output_df = data.frame(sites = event_stats_k$sites[N], lon=event_stats_k$obs_lon[N], lat=event_stats_k$obs_lat[N])
write.csv(output_df, file='gauge_location_info.csv', row.names=FALSE)

#
# Count gauges observations for each event
#
unique_sites_events = unique(event_stats_k$site_and_event)
sites_only = unlist(lapply(unique_sites_events, function(x) strsplit(x, split="_")[[1]][1]))

print(unique_sites_events)
print(table(sites_only))
