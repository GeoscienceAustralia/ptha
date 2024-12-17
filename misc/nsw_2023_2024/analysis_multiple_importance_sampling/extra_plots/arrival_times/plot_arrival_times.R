library(terra)
source('osm_backdrop.R')

minimum_arrival_puysegur2 = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/puysegur2/minimum_arrival_time/*.vrt'))
mean_arrival_puysegur2 = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/puysegur2/mean_arrival_time/*.vrt'))

minimum_arrival_southamerica = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/southamerica/minimum_arrival_time/*.vrt'))
mean_arrival_southamerica = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/southamerica/mean_arrival_time/*.vrt'))

minimum_arrival_kermadectonga2 = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/kermadectonga2/minimum_arrival_time/*.vrt'))
mean_arrival_kermadectonga2 = rast(
    Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/kermadectonga2/mean_arrival_time/*.vrt'))

#
# Puysegur 2 plot
#
png('Puysegur2_arrival_times.png', width=14, height=6, units='in', res=300)
par(mfrow=c(1,2))
NH = 36
plot(minimum_arrival_puysegur2, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600, 
    xlim=c(75, 200), ylim=c(-80, 20),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(minimum_arrival_puysegur2, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Minimum arrival time (puysegur2)', cex.main=2)

plot(mean_arrival_puysegur2, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600,
    xlim=c(75, 200), ylim=c(-80, 20),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(mean_arrival_puysegur2, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Mean arrival time (puysegur2)', cex.main=2)
dev.off()

#
# Southamerica 2 plot
#
png('southamerica_arrival_times.png', width=14, height=4, units='in', res=300)
par(mfrow=c(1,2))
NH = 36
plot(minimum_arrival_southamerica, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600, 
    xlim=c(75, 320), ylim=c(-80, 80),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(minimum_arrival_southamerica, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Minimum arrival time (southamerica)', cex.main=2)

plot(mean_arrival_southamerica, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600,
    xlim=c(75, 320), ylim=c(-80, 80),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(mean_arrival_southamerica, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Mean arrival time (southamerica)', cex.main=2)
dev.off()

#
# Kermadectonga plot
#
png('kermadectonga2_arrival_times.png', width=14, height=4, units='in', res=300)
par(mfrow=c(1,2))
NH = 36
plot(minimum_arrival_kermadectonga2, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600, 
    xlim=c(75, 320), ylim=c(-80, 80),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(minimum_arrival_kermadectonga2, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Minimum arrival time (kermadectonga2)', cex.main=2)

plot(mean_arrival_kermadectonga2, col=hcl.colors(NH, 'Spectral'), range=c(0, NH)*3600,
    xlim=c(75, 320), ylim=c(-80, 80),
    plg=list(side=4, at=seq(1,NH, by=2)*3600, labels=as.character(paste0(seq(1,NH, by=2), 'h')), cex.lab=1.3),
    pax=list(cex.axis=1.4))
contour(mean_arrival_kermadectonga2, levels=seq(1, NH)*3600, add=TRUE, 
    labels=as.character(paste0(seq(1, NH), 'h')), labcex=1.5)
title(main='Mean arrival time (kermadectonga2)', cex.main=2)
dev.off()


#
# Zoom-in plots 
#
source('cities.R')
add_cities<-function(){
    points(cities$lon, cities$lat, pch=19)
    text(cities$lon, cities$lat, cities$name, pos=2, cex=0.9)
}

plot_panel<-function(variable, variable_title, XLIM, YLIM, labseq, maxcells, include_contour){
    plot(variable, col=hcl.colors(length(labseq)-1, 'Spectral'), range=range(labseq)*3600,
        xlim=XLIM, ylim=YLIM, #asp=cos(-32/180*pi), 
        plg=list(side=4, at=labseq[-seq(1, length(labseq), by=2)]*3600, labels=as.character(paste0(labseq[-seq(1, length(labseq), by=2)], 'h')), cex.lab=1.5),
        pax=list(cex.axis=1.5))

    # Subset to stop contouring using too much memory
    #cropped_variable = crop(
    #    variable, 
    #    rast(xmin=XLIM[1], xmax=XLIM[2], ymin=YLIM[1], ymax=YLIM[2], crs=crs(variable), nrows=1000, ncols=round(1000*diff(YLIM)/diff(XLIM))))

    if(include_contour) contour(variable, levels=seq(1, 24, by=0.5)*3600, add=TRUE, maxcells=maxcells, col='black', method='flattest',
        labels=as.character(paste0(seq(1, 24, by=0.5), 'h')), labcex=1, xlim=XLIM, ylim=YLIM)
    title(main=variable_title, cex.main=1.7, line=-1)
    add_cities()
}

# Make a 4panel plot, given the XLIM/YLIM
plot_panels<-function(XLIM, YLIM, labseq, maxcells=1e+07, include_contour=TRUE){
    par(mfrow=c(1,4))
    plot_panel(minimum_arrival_southamerica, 'Minimum arrival time \n southamerica', XLIM, YLIM, labseq, maxcells=maxcells, include_contour=include_contour)
    plot_panel(mean_arrival_southamerica, 'Mean arrival time \n southamerica', XLIM, YLIM, labseq, maxcells=maxcells, include_contour=include_contour)
    plot_panel(minimum_arrival_kermadectonga2, 'Minimum arrival time \n kermadectonga2', XLIM, YLIM, labseq, maxcells=maxcells, include_contour=include_contour)
    plot_panel(mean_arrival_kermadectonga2, 'Mean arrival time \n kermadectonga2', XLIM, YLIM, labseq, maxcells=maxcells, include_contour=include_contour)
}

png('Local_arrival_times_zoom.png', width=12, height=4.5, units='in', res=300)
XLIM = c(140, 170)
YLIM = c(-50, -10)
labseq = seq(0, 24, by=0.5)
plot_panels(XLIM, YLIM, labseq)
dev.off()

png('Local_arrival_times_zoom2.png', width=9, height=4.8, units='in', res=300)
XLIM = c(149.83, 149.95)
YLIM = c(-37.15, -36.9)
labseq = seq(0, 24, by=0.5)
plot_panels(XLIM, YLIM, labseq, include_contour=FALSE)
dev.off()

#png('Local_arrival_times_zoom_north.png', width=12, height=6, units='in', res=300)
#XLIM = c(114, 115.5)
#YLIM = c(-32, -28.5)
#par(oma=c(0, 0, 0, 1))
#labseq = seq(2, 6, by=0.25)
#plot_panels(XLIM, YLIM, labseq)
#dev.off()
#
#png('Local_arrival_times_zoom_south.png', width=12, height=6, units='in', res=300)
#XLIM = c(114.5, 116)
#YLIM = c(-34.25, -30.75)
#par(oma=c(0, 0, 0, 1))
#labseq = seq(2, 6, by=0.25)
#plot_panels(XLIM, YLIM, labseq)
#dev.off()


#
# Zoom in close.
#
