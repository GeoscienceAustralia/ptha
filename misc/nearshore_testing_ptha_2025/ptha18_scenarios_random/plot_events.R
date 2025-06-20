#
# Plot the locations of the historical earthquakes that we use.
#

# Get data to make the plots
source('data_for_plots.R')

png('target_scenario_locations.png', width=8.5, height=5, units='in', res=200)
# Plot them
plot(c(-40, 320), c(-75, 75), asp=1, col='white',
    xlab="", ylab="", xaxs='i', yaxs='i')
#plot(zero_contour, add=TRUE, col='grey')
plot(tm_wb, col='grey', border='grey', add=TRUE)
plot(shift(tm_wb, dx=-360), col='grey', border='grey', add=TRUE)
plot(shift(tm_wb, dx=360), col='grey', border='grey', add=TRUE)
for(i in 1:length(uss)){
    plot(uss[[i]], add=TRUE, border=NA, col='red')
}

# Plot the focal mechanisms
library(RFOC)
plot_order = rev(order(all_events_focal_mech$Mw)) # Plot high Mw first so low Mw display beneath.
for(i in plot_order){
    # To avoid double-plotting some focal mechanisms we skip 'batch 2'
    # (since the same scenarios were run in the first batch, for 30 samples,
    # followed by another 30 in batch 2)
    if(grepl('batch2', all_events_focal_mech$event_name[i])) next

    # Get the mechanism
    s = all_events_focal_mech$strk1[i]
    d = all_events_focal_mech$dip1[i]
    r = all_events_focal_mech$rake1[i]
    lon = all_events_focal_mech$hypo_lon[i]
    lon = lon + 360*(lon < -40)
    lat = all_events_focal_mech$hypo_lat[i]
    Mw = all_events_focal_mech$Mw[i]
    size = (Mw - 6.5)*0.15
    ball_col = c('orange', 'green', 'blue')[floor(Mw) - 6]
    MECH = SDRfoc(s, d, r, PLOT=FALSE)
    justfocXY(MECH, x=lon, y=lat, focsiz=size, fcol=ball_col)
}
dev.off()

#
# Alongstrike ranges
#
target_events_with_alongstrike = read.csv('target_scenarios_data_frame_with_alongstrike_index.csv')

source_zone = basename(dirname(target_events_with_alongstrike$tsunami_event_dir))

png('target_scenario_locations_with_centroid_ranges.png', width=8.5, height=5, units='in', res=300)
# Plot them
plot(c(-40, 320), c(-75, 75), asp=1, col='white',
    xlab="", ylab="", xaxs='i', yaxs='i')
#plot(zero_contour, add=TRUE, col='grey')
plot(tm_wb, col='grey', border='grey', add=TRUE)
plot(shift(tm_wb, dx=-360), col='grey', border='grey', add=TRUE)
plot(shift(tm_wb, dx=360), col='grey', border='grey', add=TRUE)
for(i in 1:length(uss)){
    plot(uss[[i]], add=TRUE, border=NA, col='red')
}

store_borders = list()
counter=0
for(i in 1:nrow(target_events_with_alongstrike)){
    # Avoid double up of batch 2 events (since that should be part of the first batch to have 60 scenarios total)
    if(grepl('batch2', target_events_with_alongstrike$event_name[i])) next
    counter=counter+1

    sz = source_zone[i]
    unit_source = uss[[sz]]
    alongstrike_lower = target_events_with_alongstrike$alongstrike_lower_ind[i]
    alongstrike_upper = target_events_with_alongstrike$alongstrike_upper_ind[i]
    k = which(unit_source$alngst_ >= alongstrike_lower  & 
              unit_source$alngst_ <= alongstrike_upper )

    orange_trans = rgb(1,1,0, alpha=0.7)
    plot(unit_source[k,], col=orange_trans, add=TRUE, border=NA)
    
    region_border = gUnaryUnion(unit_source[k,])
    store_borders[[counter]] = region_border
}
for(i in 1:length(store_borders)){
    plot(store_borders[[i]], add=TRUE, border='black', col=NA, lwd=0.5)
}
dev.off()
