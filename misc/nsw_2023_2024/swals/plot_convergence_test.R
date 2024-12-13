library(terra)
library(cptcity)
library(basemaps)
library(sf)

contour_zoom   = vect('elevation_contour_level_0')
contour_lowres = vect('/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/zero_contour/zero_contour.shp')

r_standard = rast('OUTPUTS/run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1/RUN_20241112_155633335/all_max_stage.vrt')
r_highres  = rast('OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/all_max_stage_CONVERGENCE.vrt')

r_standard_speed = rast('OUTPUTS/run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1/RUN_20241112_155633335/all_max_speed.vrt')
r_highres_speed  = rast('OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/all_max_speed_CONVERGENCE.vrt')

output_dir = 'convergence_test_images'
dir.create(output_dir, showWarnings=FALSE)

plotcol = colorRampPalette(cpt("grass_precipitation_monthly", n=2500), bias=2.)(2500)
plotcol_max_stage = c(plotcol, rep(plotcol[length(plotcol)], length(plotcol)))
plotcol_max_speed = colorRampPalette(cpt("grass_precipitation_monthly", n=2500), bias=1.)(2500)
transparent = rgb(1,1,1, maxColorValue=1, alpha=0)
plotcol_max_speed[1] = transparent
plotcol_max_speed = c(plotcol_max_speed, rep(plotcol_max_speed[length(plotcol_max_speed)], length(plotcol_max_speed)))

add_textbox_upper_left<-function(string, location='topleft', ...){
    legend(location, legend=string, col=NA,
        bg=rgb(1,1,1,alpha=0.7,maxColorValue=1), box.col=NA, ...)
}

twoplots<-function(xlim=NULL, ylim=NULL, label_cex=1.8, mar=c(2,1,1,3)){

    # Choose a contour depending on the plot scale
    if(is.null(xlim)[1]){
        contour_l = contour_lowres
    }else if(diff(xlim) > 7){
        contour_l = contour_lowres
    }else{
        contour_l = contour_zoom
    }

    # First panel plot
    plot(r_standard, col=plotcol, xlim=xlim, ylim=ylim, range=c(1.1, 20), maxcell=1e+06, mar=mar, plg=list(title='max-stage (m)'))
    plot(contour_l, add=TRUE)
    add_textbox_upper_left('Standard', cex=label_cex)
    # Second panel plot
    plot(r_highres, col=plotcol, xlim=xlim, ylim=ylim, range=c(1.1, 20), maxcell=1e+06, mar=mar, plg=list(title='max-stage (m)'))
    plot(contour_l, add=TRUE)
    add_textbox_upper_left('High res', cex=label_cex)
}

twoplots_with_basemap<-function(xlim, ylim, label_cex=1.8, label_line=0, mar=c(2,1,1,3), map_service = 'esri', map_type='world_imagery', bg_buffer_frac=0.1, maxcell_scale=1, plotvar='max_stage'){

    # Choose a contour depending on the plot scale
    if(is.null(xlim)[1]){
        contour_l = contour_lowres
    }else if(diff(xlim) > 7){
        contour_l = contour_lowres
    }else{
        contour_l = contour_zoom
    }

    # Get a background map
    bmap_xlim = xlim + c(-1,1)*diff(xlim)/10*bg_buffer_frac
    bmap_ylim = ylim + c(-1,1)*diff(ylim)/10*bg_buffer_frac
    myext = st_bbox(ext(c(bmap_xlim, bmap_ylim)), crs='EPSG:4326')
    r0 = basemap_terra(myext, map_service=map_service, map_type=map_type, map_res=1)
    r0_4326 = project(r0, 'EPSG:4326')

    # Crop the main rasters of interest -- this ensures they have OK resolution when plotted with terra::plot(..., add=TRUE, ...)
    if(plotvar == 'max_stage'){
        r_standard_small = crop(r_standard, ext(r0_4326))
        r_highres_small = crop(r_highres, ext(r0_4326))
        plotcol = plotcol_max_stage
        plot_range = c(1.1, 20)
        legend_title = 'max stage (m)'
    }else if(plotvar == 'max_speed'){
        r_standard_small = crop(r_standard_speed, ext(r0_4326))
        r_highres_small = crop(r_highres_speed, ext(r0_4326))
        plotcol = plotcol_max_speed
        plot_range = c(0, 10)
        legend_title = 'max speed (m/s)'
    }else{
        stop(paste0('unknown plot_var: ', plot_var))
    }

    # First panel
    plot(r_standard_small, col=plotcol, xlim=xlim, ylim=ylim, range=plot_range, maxcell=1e+06*maxcell_scale, mar=mar, plg=list(title=legend_title))
    plot(r0_4326, add=TRUE)
    plot(r_standard_small, col=plotcol, range=plot_range, add=TRUE, maxcell=1e+06*maxcell_scale, legend=FALSE)
    plot(contour_l, add=TRUE)
    title('Standard', cex.main=label_cex, line=label_line)

    # Second panel
    plot(r_highres_small, col=plotcol, xlim=xlim, ylim=ylim, range=plot_range, maxcell=1e+06*maxcell_scale, mar=mar, plg=list(title=legend_title))
    plot(r0_4326, add=TRUE)
    plot(r_highres_small, col=plotcol, range=plot_range, add=TRUE, maxcell=1e+06*maxcell_scale, legend=FALSE)
    plot(contour_l, add=TRUE)
    title('High res', cex.main=label_cex, line=label_line)

}

#
# Global plot
#
png(paste0(output_dir, '/Global_scale_plot.png'), width=8.5, height=7, units='in', res=300)
par(mfrow=c(2,1))
twoplots()
dev.off()

#
# Regional plot
#
png(paste0(output_dir, '/Regional_scale_plot.png'), width=9, height=4, units='in', res=300)
par(mfrow=c(1,2))
twoplots(xlim=c(150,161.2), ylim=c(-37, -27))
dev.off()

#
# Port Stephens to Crescent Head
#
png(paste0(output_dir, '/Port_Stev_to_Crescent_Head.png'), width=10, height=5.2, units='in', res=400)
par(mfrow=c(1,2))
twoplots(xlim=c(151.8, 153.1), ylim=c(-32.8, -31.4))
dev.off()

#
# Plots with background imagery below
#
png(paste0(output_dir, '/Manning_river.png'), width=11, height=4.5, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(152.53, 152.76), ylim=c(-31.97, -31.82), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()

png(paste0(output_dir, '/Forster.png'), width=11, height=4.5, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(152.48, 152.55), ylim=c(-32.2, -32.15), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()

png(paste0(output_dir, '/Newcastle.png'), width=11, height=4.5, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(151.745, 151.814), ylim=c(-32.931, -32.890), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()

png(paste0(output_dir, '/Manly.png'), width=11, height=4.5, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(151.22, 151.31), ylim=c(-33.83, -33.77), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()

png(paste0(output_dir, '/Lake_Illawarra.png'), width=8, height=7.3, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(150.85, 150.9), ylim=c(-34.595, -34.5), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2, maxcell_scale=8)
dev.off()

png(paste0(output_dir, '/Batemans_Bay.png'), width=11, height=4.0, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(150.15, 150.3), ylim=c(-35.76, -35.68), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()


png(paste0(output_dir, '/Lord_Howe.png'), width=11, height=4.0, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(159.02, 159.12), ylim=c(-31.56, -31.5), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2)
dev.off()


#
# Max-speed
#

png(paste0(output_dir, '/Hawkesbury_max_speed.png'), width=11, height=4.5, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(151.14,151.5), ylim=c(-33.67,-33.42), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2, plotvar='max_speed')
dev.off()

png(paste0(output_dir, '/Sydney_max_speed.png'), width=11, height=7, units='in', res=400)
par(mfrow=c(1,2))
twoplots_with_basemap(xlim=c(151.08,151.33), ylim=c(-34.05,-33.77), bg_buffer_frac=0.3, mar=c(0.1,2,0.1,4.5), label_line=-2, plotvar='max_speed')
dev.off()

