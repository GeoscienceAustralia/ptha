#
# Make some plots used in the paper.
#
source('global_variables.R')
source('parse_gauge_data.R')
library(geosphere)

# Time corresponding to 'Julian day zero' in R's julian() function
R_JULIAN_ORIGIN = as.POSIXlt("1970-01-01", tz = "Etc/UTC")

pressure_metadata = read.csv(MSLP_METADATA_TABLE_FILE)

#
# Pressure gauge time-series plot
#

# Broome record
pg_ind = grep('BROOME AIRPORT', pressure_metadata$station)

pg1 = read.csv(pressure_metadata$postprocessed_file[pg_ind], comment.char='#')
site_lon = pressure_metadata$lon[pg_ind]
site_lat = pressure_metadata$lat[pg_ind]
dist_tonga = distHaversine(matrix(c(site_lon, site_lat), ncol=2), matrix(hunga_tonga_volcano, ncol=2))

XLIM = c(19006, 19012)
XLIM_ZOOM = c(19007., 19008.)
MAINSIZE = 1.7
AXISSIZE = 1.3
LABSIZE = 1.5
# We plot with Julian days, but label axes with a date/time format
XLAB_JULIAN = seq(XLIM[1], XLIM[2])
XLAB_JULIAN_label = as.difftime(XLAB_JULIAN, units='days') +  R_JULIAN_ORIGIN
XLAB_JULIAN_ZOOM = seq(XLIM_ZOOM[1], XLIM_ZOOM[2], by=1/6)
XLAB_JULIAN_ZOOM_label = as.difftime(XLAB_JULIAN_ZOOM, units='days') + R_JULIAN_ORIGIN

png(paste0(OUTPUT_GRAPHICS_DIR, '/pressure_gauge_processing_example.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(4.5,6,3,2))
# Panel 1: Main pressure data
plot(pg1$juliant, pg1$pressure, t='l', 
    xlab='Time (UTC)', ylab='Pressure (hPa)',
    main='A) MSLP at Broome Airport',
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM, xaxt='n')
axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)
abline(v=julian(explosion_start_time), col='red')
text(julian(explosion_start_time)-0.05, median(pg1$pressure, na.rm=TRUE) + 1, '   Explosion', srt=90, col=2, cex=1.3)

# Panel 2: Firstly plot main data invisibly, then add it in later
plot(pg1$juliant, pg1$resid2h, t='l', col='white',
    xlab='Time (UTC)', ylab='High-pass MSLP (hPa)',
    main='B) High-pass filtered MSLP (periods < 2 hours)', 
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM, ylim=c(-2.5,2.5), xaxt='n')
axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)

# Expected arrival times of wave travelling the short path from Tonga to the site
wave_1 = ( dist_tonga + c(0, 1, 2, 3, 4) * 2 * pi * EARTH_RADIUS)/LAMB_WAVE_SPEED
# Expected arrival times of wave travelling the long path from Tonga to the site
wave_2 = ( (2*pi*EARTH_RADIUS - dist_tonga) + c(0, 1, 2, 3, 4) * 2 * pi * EARTH_RADIUS)/LAMB_WAVE_SPEED

# Include rough theoretical lamb-wave arrivals
abline( v=julian(explosion_start_time + as.difftime(wave_1,units='secs')), col='blue', lty='dotted')
abline( v=julian(explosion_start_time + as.difftime(wave_2,units='secs')), col='blue', lty='dotted')

abline(v=julian(explosion_start_time), col='red')

# Add back in the main data 
points(pg1$juliant, pg1$resid2h, t='l')

# Panel 3: Like Panel 2 with a different xlim
plot(pg1$juliant, pg1$resid2h, t='l', col='white',
    xlab='Time (UTC)', ylab='High-pass MSLP (hPa)',
    main='C) Zoom of high-pass filtered MSLP on January 15 2022', 
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM_ZOOM, ylim=c(-2.5,2.5), xaxt='n')
axis(1, at=XLAB_JULIAN_ZOOM, labels=format(XLAB_JULIAN_ZOOM_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)

# Expected arrival times of wave travelling the short path from Tonga to the site
wave_1 = ( dist_tonga + c(0, 1, 2, 3, 4) * 2 * pi * EARTH_RADIUS)/LAMB_WAVE_SPEED
# Expected arrival times of wave travelling the long path from Tonga to the site
wave_2 = ( (2*pi*EARTH_RADIUS - dist_tonga) + c(0, 1, 2, 3, 4) * 2 * pi * EARTH_RADIUS)/LAMB_WAVE_SPEED

# Include rough theoretical lamb-wave arrivals
abline( v=julian(explosion_start_time + as.difftime(wave_1,units='secs')), col='blue', lty='dotted')
abline( v=julian(explosion_start_time + as.difftime(wave_2,units='secs')), col='blue', lty='dotted')

abline(v=julian(explosion_start_time), col='red')
# Add back in the main data 
points(pg1$juliant, pg1$resid2h, t='l')

#text(julian(explosion_start_time)-0.01, 1, 'HT Explosion', srt=90, col=2)

dev.off()

#
# Tide gauge time-series plot
#

# Crowdy-Head record
tg1 = read.csv(paste0(OUTPUT_TIDE_DIR, '/Crowdy_Head_Fishermans_Wharf_MHL.csv'), comment.char='#')

# Get lamb-wave arrival time at crowdy head
tide_gauge_metadata = read.csv(TIDEGAUGE_METADATA_TABLE_FILE)
ci = grep('rowdy', tide_gauge_metadata$name); stopifnot(length(ci) == 1)
ch_loc = matrix(c(tide_gauge_metadata$lon[ci], tide_gauge_metadata$lat[ci]), ncol=2)
lamb_wave_arrival_ch = explosion_start_time + as.difftime(
    (distHaversine(ch_loc, matrix(hunga_tonga_volcano, ncol=2))/LAMB_WAVE_SPEED), 
    units='secs')

XLIM = c(19006, 19012)
XLIM_ZOOM = c(19007., 19008.)
MAINSIZE = 1.7
AXISSIZE = 1.3
LABSIZE = 1.5
# We plot with Julian days, but label axes with a date/time format
XLAB_JULIAN = seq(XLIM[1], XLIM[2])
XLAB_JULIAN_label = as.difftime(XLAB_JULIAN, units='days') + R_JULIAN_ORIGIN
XLAB_JULIAN_ZOOM = seq(XLIM_ZOOM[1], XLIM_ZOOM[2], by=1/6)
XLAB_JULIAN_ZOOM_label = as.difftime(XLAB_JULIAN_ZOOM, units='days') +  R_JULIAN_ORIGIN

png(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_processing_example.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(4.5,6,3,2))

plot(tg1$juliant, tg1$stage, t='l',
    xlab='Time (UTC)', ylab='Sea level (m)',
    main='A) Crowdy Head tide gauge',
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM, xaxt='n')
axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)

abline(v=julian(explosion_start_time), col='red')
text(julian(explosion_start_time)-0.05, 1.3, '   Explosion', srt=90, col=2, cex=1.3)
plot(tg1$juliant, tg1$resid3h, t='l',
    xlab='Time (UTC)', ylab='High-pass sea level (m)',
    main='B) High-pass filtered sea level (periods < 3 hours)',
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM, xaxt='n')
axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)
abline(v=julian(explosion_start_time), col='red')

plot(tg1$juliant, tg1$resid3h, t='l',
    xlab='Time (UTC)', ylab='High-pass sea level (m)',
    main='C) Zoom of high-pass filtered sea level on January 15 2022',
    cex.main = MAINSIZE, cex.axis=AXISSIZE, cex.lab=LABSIZE,
    xlim=XLIM_ZOOM, xaxt='n')
abline(v=julian(explosion_start_time), col='red')
axis(1, at=XLAB_JULIAN_ZOOM, labels=format(XLAB_JULIAN_ZOOM_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)

# Add 'tsunami travel time' [gravity wave]
abline(v=julian(CROWDY_HEAD_TTT_ARRIVAL_TIME), lty='dashed', col='darkgreen')
text(julian(CROWDY_HEAD_TTT_ARRIVAL_TIME)+0.012, -1.0, '  TTT', srt=0, col='darkgreen', cex=1.3)

# Add Lamb wave arrival time
abline(v=julian(lamb_wave_arrival_ch), lty='dotted', col='blue')
text(julian(lamb_wave_arrival_ch)+0.012, -1.0, '  LW', srt=0, col='blue', cex=1.3)
dev.off()


#
# Pressure data with size of initial Lamb wave, and arrival time.
#
all_mslp = lapply(pressure_metadata$postprocessed_file, function(x) read.csv(x, comment.char='#'))

# Distance of all stations to the volcano
distance_mslp_tonga = distHaversine(cbind(pressure_metadata$lon, pressure_metadata$lat), 
    matrix(hunga_tonga_volcano, ncol=2, nrow=nrow(pressure_metadata), byrow=TRUE))

# Find index corresponding to the initial peak of the lamb-wave, by searching for the first
# peak following the expected arrival of the lamb wave. The peak will come after the arrival, but not too far. 
# The earlier graphical checks indicate that this should work.
expected_arrival_time = explosion_start_time + as.difftime( (distance_mslp_tonga/LAMB_WAVE_SPEED), units='secs')

# Search for the index of the highpass pressure peak following the expected
# arrival, searching 2h before and after (way too much, but this will help catch
# problems).
lamb_peak_index = rep(NA, length(all_mslp))
lamb_peak_jtime = rep(NA, length(all_mslp))
lamb_peak_size = rep(NA, length(all_mslp))
lamb_trough_index = rep(NA, length(all_mslp))
lamb_trough_jtime = rep(NA, length(all_mslp))
lamb_trough_size = rep(NA, length(all_mslp))

j_arrival = julian(expected_arrival_time)
for(i in 1:length(all_mslp)){
    k = which( (all_mslp[[i]]$juliant > (j_arrival[i] - 2/24)) &
               (all_mslp[[i]]$juliant < (j_arrival[i] + 2/24)))

    lpi = k[which.max(all_mslp[[i]]$resid2h[k])]
    lamb_peak_index[i] = lpi
    lamb_peak_jtime[i] = all_mslp[[i]]$juliant[lpi]
    lamb_peak_size[i] = all_mslp[[i]]$resid2h[lpi]

    lti = k[which.min(all_mslp[[i]]$resid2h[k])]   
    lamb_trough_index[i] = lti
    lamb_trough_jtime[i] = all_mslp[[i]]$juliant[lti]
    lamb_trough_size[i] = all_mslp[[i]]$resid2h[lti]
}

# Compare the 'time of peak' to the 'expected arrival time'.
time_difference_min = (as.numeric(24*60*(lamb_peak_jtime - j_arrival)))
#sort(time_difference_min, decreasing=TRUE)[1:5]
#[1] 45.13686 15.76457 14.23174 13.10642 12.31802

XLAB_JULIAN = seq(19007, 19008, by=1/12)
XLAB_JULIAN_label = as.difftime(XLAB_JULIAN, units='days') +  R_JULIAN_ORIGIN

png(paste0(OUTPUT_GRAPHICS_DIR, '/lamb_wave_arrival_and_pressure_maxima_arrival.png'), width=9, height=7, units='in', res=200) 
par(mar=c(5.1, 5.1, 2.1, 2.1))
plot(j_arrival, lamb_peak_jtime, 
    xlab = paste0('Theoretical arrival time of Lamb Wave'),
    ylab = 'Observed time of high-pass filtered MSLP maxima',
    cex.lab=1.5, asp=1, xaxt='n', yaxt='n')
axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)
axis(2, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)
#grid(col='orange')
abline(h=XLAB_JULIAN, col='orange', lty='dotted')
abline(v=XLAB_JULIAN, col='orange', lty='dotted')
abline(0, 1, col='red')
abline(0, 1/(LAMB_WAVE_SPEED), col='red')
dev.off()

png(paste0(OUTPUT_GRAPHICS_DIR, '/distance_tonga_and_pressure_maxima_arrival.png'), width=9, height=7, units='in', res=200) 
par(mar=c(5.1, 5.1, 2.1, 2.1))
plot(distance_mslp_tonga/1000, lamb_peak_jtime, 
    xlab = paste0('Distance from Hunga Volcano (km)'),
    ylab = 'Time of high-pass filtered MSLP maxima',
    cex.lab=1.5, yaxt='n')
axis(2, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=AXISSIZE)
grid(col='orange')
dev.off()

png(paste0(OUTPUT_GRAPHICS_DIR, '/distance_tonga_and_pressure_maxima_size.png'), width=7, height=7, units='in', res=200) 
plot(distance_mslp_tonga, lamb_peak_size, 
    xlab= 'Distance from Hunga Volcano (m)', 
    ylab= 'Size of high-pass filtered MSLP maxima (hPa)',
    cex.lab=1.5)
grid(col='orange')
dev.off()

png(paste0(OUTPUT_GRAPHICS_DIR, '/lamb_wave_peak_size.png'), width=8, height=5.5, units='in', res=300)

dc = 1000 # Discretize colors
Z_max = 4.5 # Upper limit of pressure perturbation for colour scheme
stopifnot(Z_max > max(lamb_peak_size))
library(RColorBrewer)
target_col = brewer.pal(n=9, name='PRGn')
target_col[5] = 'lightyellow' # Avoid white, as some reviewers thought it could imply zero.
mycols = colorRampPalette(rep(rev(target_col), each=20))(dc)
final_col = mycols[pmax(round((dc-1)*lamb_peak_size/Z_max), 1)]

par(mar=c(2.1, 4.1, 2.1, 5.1))
plot(pressure_metadata$lon, pressure_metadata$lat, col=final_col, pch=19, cex=1.4,
    xlim=c(90, 190), ylim=c(-55, -5),
    xlab="", ylab="", main='Maximum high-pass filtered MSLP within 2h of Lamb wave arrival (hPa)', cex.main=1.3,
    asp=1/cos(30/180*pi))
# Add outer circles around points for visual clarity
points(pressure_metadata$lon, pressure_metadata$lat, col='darkgrey', asp=1, pch=1, cex=1.4, lwd=0.2)
# Use image plot purely for the colorbar
library(fields)
image.plot(
    matrix(seq(0, Z_max, len=dc), ncol=dc, nrow=dc),
    zlim = c(0, Z_max),  
    col=mycols, legend.only=TRUE, horizontal=FALSE, legend.mar=3.1)
source('get_simple_world_map_data.R')
plot(wrld_simpl, add=TRUE)
points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)
dev.off()

png(paste0(OUTPUT_GRAPHICS_DIR, '/pressure_arrival_time_peak_and_trough.png'), width=7, height=7, units='in', res=200) 
theoretical_arrival_secs = (j_arrival - julian(explosion_start_time))*3600*24
observed_peak_arrival_secs = as.numeric(lamb_peak_jtime -julian(explosion_start_time))* 3600*24
observed_trough_arrival_secs = as.numeric(lamb_trough_jtime -julian(explosion_start_time))* 3600*24
plot(theoretical_arrival_secs, observed_peak_arrival_secs,
    xlab='Theoretical time to arrive (seconds)', 
    ylab='Observed time to arrive, peak or trough (seconds)',
    main='Time (seconds) for arrival of \n peak and trough of initial Lamb wave',
    ylim=c(0, max(observed_trough_arrival_secs)), asp=1, cex.main=1.5,
    cex.lab=1.5, cex.axis=1.4)
points(theoretical_arrival_secs, observed_trough_arrival_secs, col='blue')
grid(col='orange')
abline(0, 1, col='red')
legend('bottomright', c('Lamb peak', 'Lamb trough'), col=c('black', 'blue'), pch=c(1,1))
dev.off()

## Linear regression to estimate speeds
#> lm(observed_peak_arrival_secs ~ distance_mslp_tonga)
#
#Call:
#lm(formula = observed_peak_arrival_secs ~ distance_mslp_tonga)
#
#Coefficients:
#        (Intercept)  distance_mslp_tonga  
#          8.455e+02            3.111e-03  
#
#> 1/3.111e-03
#[1] 321.4401
#> 
#> lm(observed_trough_arrival_secs ~ distance_mslp_tonga)
#
#Call:
#lm(formula = observed_trough_arrival_secs ~ distance_mslp_tonga)
#
#Coefficients:
#        (Intercept)  distance_mslp_tonga  
#          1.518e+03            3.486e-03  



png(paste0(OUTPUT_GRAPHICS_DIR, '/lamb_wave_time_from_peak_to_trough.png'), width=8, height=5.5, units='in', res=300)
peak_to_trough_time = as.numeric(observed_trough_arrival_secs - observed_peak_arrival_secs)
dc = 1000 # Discretize colors
Z_max = max(peak_to_trough_time + 100) # Upper limit of pressure perturbation for colour scheme
library(RColorBrewer)
mycols = colorRampPalette(rep(rev(brewer.pal(n=9, name='PRGn')), each=20))(dc)
final_col = mycols[pmax(round((dc-1)*peak_to_trough_time/Z_max), 1)]
par(mar=c(2.1, 4.1, 2.1, 5.1))
plot(pressure_metadata$lon, pressure_metadata$lat, 
    cex=1, asp=1/cos(30/180*pi), pch=19, col=final_col,
    main='Time (seconds) between peak and trough of initial Lamb wave', cex.main=1.4)
# Add outer circles around points for visual clarity
points(pressure_metadata$lon, pressure_metadata$lat, col='darkgrey', asp=1, pch=1, cex=1.4, lwd=0.2)
# Use image plot purely for the colorbar
library(fields)
image.plot(
    matrix(seq(0, Z_max, len=dc), ncol=dc, nrow=dc),
    zlim = c(0, Z_max),  
    col=mycols, legend.only=TRUE, horizontal=FALSE, legend.mar=4.1)
source('get_simple_world_map_data.R')
plot(wrld_simpl, add=TRUE)
points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)
dev.off()

#
# Plots of tidal residual maxima on Jan 15-16
#
tg_metadata = read.csv(TIDEGAUGE_METADATA_TABLE_FILE)
all_tg = lapply(tg_metadata$postprocessed_file, function(x) read.csv(x, comment.char='#'))

tg_maxima = rep(NA, length(all_tg))
tg_minima = rep(NA, length(all_tg))
for(i in 1:length(all_tg)){
    # Jan 15 and 16
    k = which(all_tg[[i]]$juliant > 19007 & all_tg[[i]]$juliant < 19009)
    tg_maxima[i] = max(all_tg[[i]]$resid3h[k], na.rm=TRUE)
    tg_minima[i] = min(all_tg[[i]]$resid3h[k], na.rm=TRUE)
}

ZOOMBOX = c(148, 155, -37, -26)

png(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_peak_size.png'), width=7, height=4.5, units='in', res=300)

dc = 1000 # Discretize colors
Z_max = 1.4 # Upper limit of pressure perturbation for colour scheme
stopifnot(Z_max > max(tg_maxima))
library(RColorBrewer)
mycols = colorRampPalette(rep(rev(brewer.pal(n=11, name='Spectral')), each=1), bias=1.5)(dc)
final_col = mycols[pmax(round((dc-1)*tg_maxima/Z_max), 1)]
ptsize = sqrt(tg_maxima)*3
par(mar=c(3.1, 3.1, 2.1, 5.1))
plot(tg_metadata$lon, tg_metadata$lat, col=final_col, pch=19, cex=ptsize,
    xlim=c(95, 189), ylim=c(-55, -10),
    xlab="", ylab="", xaxs='i', yaxs='i',
    #main='A) Maximum residual sea level (m), January 15-16', cex.main=1.6,
    asp=1/cos(30/180*pi))
# Add outer circles around points for visual clarity
points(tg_metadata$lon, tg_metadata$lat, col='darkgrey', asp=1, pch=1, cex=ptsize, lwd=0.2)
# Use image plot purely for the colorbar
library(fields)
image.plot(
    matrix(seq(0, Z_max, len=dc), ncol=dc, nrow=dc),
    zlim = c(0, Z_max),  
    col=mycols, legend.only=TRUE, horizontal=FALSE, legend.mar=4.5)
source('get_simple_world_map_data.R')
plot(wrld_simpl, add=TRUE, border='lightblue')
points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)
rect(ZOOMBOX[1], ZOOMBOX[3], ZOOMBOX[2], ZOOMBOX[4], border='blue', lty='dotted')
text(100, -17, 'A)', col='blue', cex=2.7)
text(ZOOMBOX[1]+1.8, ZOOMBOX[4]-1.8, 'B)', col='blue', cex=1.3)
dev.off()

# As above zoomed around east coast, as suggested by Kaya

png(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_peak_size_inset.png'), width=2.7, height=4.5, units='in', res=300)
dc = 1000 # Discretize colors
Z_max = 1.4 # Upper limit of pressure perturbation for colour scheme
stopifnot(Z_max > max(tg_maxima))
library(RColorBrewer)
mycols = colorRampPalette(rep(rev(brewer.pal(n=11, name='Spectral')), each=1), bias=1.5)(dc)
final_col = mycols[pmax(round((dc-1)*tg_maxima/Z_max), 1)]
ptsize = sqrt(tg_maxima)*3
par(mar=c(2.8, 1.1, 2.1, 2.1))
plot(tg_metadata$lon, tg_metadata$lat, col=final_col, pch=19, cex=ptsize*1.5,
    xlim=c(148, 155), ylim=c(-37, -26),
    xlab="", ylab="", main='', cex.main=1.4, cex.axis=1.2, yaxt='n',
    asp=1/cos(30/180*pi))
axis(side=4, cex.axis=1.2)
# Add outer circles around points for visual clarity
points(tg_metadata$lon, tg_metadata$lat, col='darkgrey', asp=1, pch=1, cex=ptsize, lwd=0.2)
source('get_simple_world_map_data.R')
plot(wrld_simpl, add=TRUE, border='lightblue')
points(hunga_tonga_volcano[1], hunga_tonga_volcano[2], pch=17, col='red', cex=3)
rect(ZOOMBOX[1], ZOOMBOX[3], ZOOMBOX[2], ZOOMBOX[4], border='blue', lty='dashed')
text(ZOOMBOX[1]+1, ZOOMBOX[4]-1, 'B)', col='blue', cex=2.5)
dev.off()


#
# Comparisons of nearby sites
#

# 
tg_pairs = list(
    list(
        g1 = "Rosslyn_Bay_IOC",
        g2 = "Rosslyn_Bay_DES",
        title = 'Rossyln Bay, two nearby gauges'
    ),
    list(
        g1 = "Gold_Coast_Sand_Bypass_Jetty_IOC",
        g2 = "Gold_Coast_DES",
        title = 'Gold Coast Sand Bypass Jetty, two nearby gauges'
    ),
    list(
        g1 = "Eden_MHL",
        g2 = "Cruise_Wharf_Tide_PANSW",
        title = 'Eden Wharf, two nearby gauges'
    ),
    list(
        g1 = "Fort_Denison_Primary_Tide_PANSW",
        g2 = "Sydney_Fort_Denison_BOMPorts",
        title = 'Fort Denison, one gauge with 1min vs 6min sampling'
    )
)

png(paste0(OUTPUT_GRAPHICS_DIR, '/tide_highpass_comparisons.png'), width=8.5, height=10, units='in', res=300)
par(mfrow=c(4,1))
par(mar = c(2,4.3,2.5,1))
par(oma = c(2,0,0,0))
XLAB_JULIAN = seq(19007, 19008, by=1/6)
XLAB_JULIAN_label = as.difftime(XLAB_JULIAN, units='days') +  R_JULIAN_ORIGIN
for(i in 1:length(tg_pairs)){
    g1 = tg_pairs[[i]]$g1
    g2 = tg_pairs[[i]]$g2
    i1 = which(tg_metadata$name == g1)
    i2 = which(tg_metadata$name == g2)

    # Y-limits with a bit of fat for labels
    YLIM = range(c(all_tg[[i1]]$resid3h, all_tg[[i2]]$resid3d), na.rm=TRUE)*c(1., 1.35)

    plot(all_tg[[i1]]$juliant, all_tg[[i1]]$resid3h, t='l', xlim=c(19007.1, 19008), ylim=YLIM,
        xlab='Time (UTC)', ylab='High-pass sea level (m)', cex.lab=1.6, cex.axis=1.4, xaxt='n')
    axis(1, at=XLAB_JULIAN, labels=format(XLAB_JULIAN_label, format='%b%d %H:%M'), cex.axis=1.4)
    points(all_tg[[i2]]$juliant, all_tg[[i2]]$resid3h, t='l', col='red')
    #grid(col='orange')
    abline(h=seq(-2, 2, by=0.2), lty='dotted', col='orange')
    abline(v=XLAB_JULIAN, lty='dotted', col='orange')
    abline(v=julian(explosion_start_time), col='red', lwd=1, lty='dashed')
    #title(paste0(LETTERS[i], ') ', g1, ' compared to ', g2), cex.main=1.65)
    title(paste0(LETTERS[i], ') ', tg_pairs[[i]]$title), cex.main=1.65)

    if(i %in% c(1, 4)){
        # Standard case.
        legend('top', gsub("_", " ", c(g1, g2)), lty=c(1,1), col=c('black', 'red'), bty='n', cex=1.7, horiz=TRUE)
    }else if(i == 2){
        # Gold Coast legend looks weird, needs adjustment
        legend(x=19007.2, y=YLIM[2]*1.1, gsub("_", " ", c(g1, g2)), lty=c(1,1), col=c('black', 'red'), bty='n', cex=1.7, horiz=TRUE)
    }else if(i == 3){
        # Eden legend also needs adjustment
        legend(x=19007.2, y=YLIM[2]*1.1, gsub("_", " ", c(g1, g2)), lty=c(1,1), col=c('black', 'red'), bty='n', cex=1.7, horiz=TRUE)
    }
    #legend(x=19007.55, y=YLIM[2], c(g1, g2), lty=c(1,1), col=c('black', 'red'), bty='n', cex=1.7, horiz=TRUE)
    #legend('topleft', g1, lty=1, col='black', bty='n', cex=1.7)
    #legend('topright', g2, lty=1, col='red', bty='n', cex=1.7)

    if(i == 1){
        text(julian(explosion_start_time)+0.015, -0.14, 'Explosion', col='red', srt=90, cex=1.5)
    }

}
mtext('Time (UTC)', side=1, outer=TRUE, cex=1.4, line=0.5)
dev.off()
