library(terra)
source('osm_backdrop.R')

minimum_arrival_sunda2 = rast(
    Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/arrival_time/sunda2/arrival_time_minimum/*.vrt'))
mean_arrival_sunda2 = rast(
    Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/arrival_time/sunda2/arrival_time_scenario_average/*.vrt'))
minimum_arrival_outerrisesunda = rast(
    Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/arrival_time/outerrisesunda/arrival_time_minimum/*.vrt'))
mean_arrival_outerrisesunda = rast(
    Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/arrival_time/outerrisesunda/arrival_time_scenario_average/*.vrt'))


#
# Sunda 2 plot
#
png('Sunda2_arrival_times.png', width=14, height=6, units='in', res=300)
par(mfrow=c(1,2))

plot(minimum_arrival_sunda2, col=hcl.colors(24, 'Spectral'), zlim=c(0, 24),
    plg=list(side=4, at=seq(1,24)*3600, labels=as.character(paste0(seq(1,24), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(minimum_arrival_sunda2, levels=seq(1, 24)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, 24), 'h')), labcex=1.5)
title(main='Minimum arrival time (sunda2)', cex.main=2)

plot(mean_arrival_sunda2, col=hcl.colors(24, 'Spectral'), zlim=c(0, 24),
    plg=list(side=4, at=seq(1,24)*3600, labels=as.character(paste0(seq(1,24), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(mean_arrival_sunda2, levels=seq(1, 24)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, 24), 'h')), labcex=1.5)
title(main='Average arrival time (sunda2)', cex.main=2)
dev.off()

#
# Outerrise Sunda plot
#
png('Outerrisesunda_arrival_times.png', width=14, height=6, units='in', res=300)
par(mfrow=c(1,2))

plot(minimum_arrival_outerrisesunda, col=hcl.colors(24, 'Spectral'), zlim=c(0, 24),
    plg=list(side=4, at=seq(1,24)*3600, labels=as.character(paste0(seq(1,24), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(minimum_arrival_outerrisesunda, levels=seq(1, 24)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, 24), 'h')), labcex=1.5)
title(main='Minimum arrival time (outerrisesunda)', cex.main=2)

plot(mean_arrival_outerrisesunda, col=hcl.colors(24, 'Spectral'), zlim=c(0, 24),
    plg=list(side=4, at=seq(1,24)*3600, labels=as.character(paste0(seq(1,24), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(mean_arrival_outerrisesunda, levels=seq(1, 24)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, 24), 'h')), labcex=1.5)
title(main='Average arrival time (outerrisesunda)', cex.main=2)
dev.off()


#
# Zoom-in plots 
#
source('cities.R')
add_cities<-function(){
    points(cities$lon, cities$lat, pch=19)
    text(cities$lon, cities$lat, cities$name, pos=1, cex=1.3)
}

plot_panel<-function(variable, variable_title, XLIM, YLIM, labseq){
    plot(variable, col=hcl.colors(length(labseq)-1, 'Spectral'), range=range(labseq)*3600,
        xlim=XLIM, ylim=YLIM, #asp=cos(-32/180*pi), 
        plg=list(side=4, at=labseq[-1]*3600, labels=as.character(paste0(labseq[-1], 'h')), cex.lab=2),
        pax=list(cex.axis=2))
    title(main=variable_title, cex.main=1.7, line=-1)
    contour(variable, levels=seq(1, 24, by=0.5)*3600, add=TRUE, maxcells=1e+07, col='grey', method='flattest',
        labels=as.character(paste0(seq(1, 24, by=0.5), 'h')), labcex=1.5)
    add_cities()
}

# Make a 4panel plot, given the XLIM/YLIM
plot_panels<-function(XLIM, YLIM, labseq){
    par(mfrow=c(1,4))
    plot_panel(minimum_arrival_sunda2, 'Minimum arrival time \n sunda2', XLIM, YLIM, labseq)
    plot_panel(mean_arrival_sunda2, 'Mean arrival time \n sunda2', XLIM, YLIM, labseq)
    plot_panel(minimum_arrival_outerrisesunda, 'Minimum arrival time \n outerrisesunda', XLIM, YLIM, labseq)
    plot_panel(mean_arrival_outerrisesunda, 'Mean arrival time \n outerrisesunda', XLIM, YLIM, labseq)
}

png('Local_arrival_times_zoom.png', width=12, height=6, units='in', res=300)
XLIM = c(113, 116.5)
YLIM = c(-35, -27)
labseq = seq(2, 6, by=0.5)
plot_panels(XLIM, YLIM, labseq)
dev.off()

png('Local_arrival_times_zoom_north.png', width=12, height=6, units='in', res=300)
XLIM = c(114, 115.5)
YLIM = c(-32, -28.5)
par(oma=c(0, 0, 0, 1))
labseq = seq(2, 6, by=0.25)
plot_panels(XLIM, YLIM, labseq)
dev.off()

png('Local_arrival_times_zoom_south.png', width=12, height=6, units='in', res=300)
XLIM = c(114.5, 116)
YLIM = c(-34.25, -30.75)
par(oma=c(0, 0, 0, 1))
labseq = seq(2, 6, by=0.25)
plot_panels(XLIM, YLIM, labseq)
dev.off()


#
# Zoom in close.
#
