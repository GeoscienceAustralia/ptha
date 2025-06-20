library(terra)
source('osm_backdrop.R')

source_zones = list(
    'New Hebrides' = 'newhebrides2',
    'Outer Rise New Hebrides' = 'outerrisenewhebrides',
    'Outer Rise Solomon' = 'outerrisesolomon',
    'Solomon' = 'solomon2',
    'Kermadec Tonga' = 'kermadectonga2',
    'Kurils Japan' = 'kurilsjapan',
    'Alaska Aleutians' = 'alaskaaleutians',
    'South America' = 'southamerica'
)
n_source_zones = length(source_zones)
n_rows = 2
n_cols = n_source_zones / n_rows
n_hours = 24
range=c(0, n_hours*3600)
lab_seq=seq(0, n_hours, by=4)
contour_seq=seq(0, n_hours, by=4)

plot_panel<-function(raster, title_str, xlim, ylim, add_cities=FALSE, title_line=0.5){
    if (length(xlim) == 1 && is.na(xlim)){
        xlim = ext(raster)[1:2]
    }
    if (length(ylim) == 1 && is.na(ylim)){
        ylim = ext(raster)[3:4]
    }
    plot(raster,
        col=hcl.colors(n_hours , 'Spectral'),
        range=range,
        xlim=xlim,
        ylim=ylim, #asp=cos(-32/180*pi), 
        plg=list(side=4, at=lab_seq[-1]*3600, labels=as.character(paste0(lab_seq[-1], 'h')), cex.lab=1.5),
        pax=list(cex.axis=2)
    )
    title(main=title_str, cex.main=2, line=title_line)
    contour(raster, levels=contour_seq*3600, add=TRUE, maxcells=1e+07, col='black', method='flattest',
        labels=as.character(paste0(contour_seq, 'h')), labcex=1.5)
    if(add_cities){
        plot_cities()
    }
}

# Mean arrival time plots
png('mean_arrival_times.png', width=14, height=6, units='in', res=300)
par(mfrow=c(n_rows, n_cols))
par(oma=c(0,0,3,0))
par(mar=c(0, 0, 0, 0))
for (sz in names(source_zones)){
    mean_arrival = rast(Sys.glob(paste0('../gladstone_arrival_time_min_and_scenario_average/', source_zones[[sz]], '/mean_arrival_time.vrt')))
    plot_panel(mean_arrival, sz, xlim=NA, ylim=NA)
}
mtext('Mean arrival times', side=3, line=1, cex=2, outer=TRUE)
dev.off()

# Minimum arrival time plots
png('minimum_arrival_times.png', width=14, height=6, units='in', res=300)
par(mfrow=c(n_rows, n_cols))
par(oma=c(0,0,3,0))
par(mar=c(0, 0, 0, 0))
for (sz in names(source_zones)){
    minimum_arrival = rast(
        Sys.glob(paste0('../gladstone_arrival_time_min_and_scenario_average/', source_zones[[sz]], '/minimum_arrival_time.vrt')))
    plot_panel(minimum_arrival, sz, xlim=NA, ylim=NA)
}
mtext('Minimum arrival times', side=3, line=1, cex=2, outer=TRUE)
dev.off()

# Zoom in on the local area
source('cities.R')
plot_cities <-function(){
    points(cities$lon, cities$lat, pch=19)
    text(cities$lon, cities$lat, cities$name, pos=1, cex=1.2)
}
xlim = c(150.5, 153)
ylim = c(-25.01, -23)
contour_seq=seq(0, n_hours, by=1)

# Mean arrival time plots
png('mean_arrival_times_zoom.png', width=14, height=6, units='in', res=300)
par(mfrow=c(n_rows, n_cols))
par(oma=c(0,0,3,0))
par(mar=c(0, 0, 0, 0))
for (sz in names(source_zones)){
    mean_arrival = rast(
        Sys.glob(paste0('../gladstone_arrival_time_min_and_scenario_average/', source_zones[[sz]], '/mean_arrival_time.vrt')))
    plot_panel(mean_arrival, sz, xlim, ylim, add_cities = TRUE)
}
mtext('Mean arrival times', side=3, line=1, cex=2, outer=TRUE)
dev.off()

# Minimum arrival time plots
png('minimum_arrival_times_zoom.png', width=14, height=6, units='in', res=300)
par(mfrow=c(n_rows, n_cols))
par(oma=c(0,0,3,0))
par(mar=c(0, 0, 0, 0))
for (sz in names(source_zones)){
    minimum_arrival = rast(
        Sys.glob(paste0('../gladstone_arrival_time_min_and_scenario_average/', source_zones[[sz]], '/minimum_arrival_time.vrt')))
    plot_panel(minimum_arrival, sz, xlim, ylim, add_cities = TRUE)
}
mtext('Minimum arrival times', side=3, line=1, cex=2, outer=TRUE)
dev.off()
