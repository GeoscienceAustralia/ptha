#
# Plot all model types and data at all sites
#
output_dir = 'all_models_sites_vs_data'
dir.create(output_dir, showWarnings=FALSE)
model_base = c(
    'Tohoku2011_YamakaziEtAl2018-risetime_0-full-',
    'Tohoku2011_SatakeEtAl2013-risetime_0-full-',
    'Tohoku2011_RomanoEtAl2015-risetime_0-full-',
    'Sumatra2004_PiatanesiLorito2007-risetime_0-full-',
    'Sumatra2004_LoritoEtAl2010-risetime_0-full-',
    'Sumatra2004_FujiSatake2007-risetime_0-full-',
    'Chile2015_WilliamsonEtAl2017-risetime_0-full-',
    'Chile2015_RomanoEtAl2016-risetime_0-full-',
    'Chile2010_LoritoEtAl2011-risetime_0-full-',
    'Chile2010_FujiSatake2013-risetime_0-full-',
    'Chile1960_HoEtAl2019-risetime_0-full-',
    'Chile1960_FujiSatake2013-risetime_0-full-')
model_site = c(
    'NSW',
    'NSW',
    'NSW',
    'australia',
    'australia',
    'australia',
    'NSW',
    'NSW',
    'NSW',
    'NSW',
    'NSW',
    'NSW')
vertical_scale = c(0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.2, 0.2, 0.5, 0.5, 0.6, 0.6)

model_types = c('linear_with_no_friction', 'linear_with_manning', 'linear_with_linear_friction', 'linear_with_reduced_linear_friction', 'linear_with_delayed_linear_friction')

make_plot<-function(Fujii_Hillarys, plot_filename, plot_ylim, legend_xy=c(35, -0.25), extra_title = '', plot_xlim=c(0, 60)){

    obs_t = (Fujii_Hillarys[[1]]$event_data$obs$juliant - Fujii_Hillarys[[1]]$model_start_julian)*24
    obs_h = Fujii_Hillarys[[1]]$event_data$obs$resid

    png(plot_filename, , width=10, height=10, units='in', res=300)
    par(oma=c(4, 4, 0, 0))
    par(mfrow=c(5,1))
    par(mar=c(1,1,3,1))

    plot(Fujii_Hillarys[[1]]$model_t/3600, Fujii_Hillarys[[1]]$model_stage, t='l', col='red', lwd=0.7, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.axis=1.5, las=1)
    title(paste0('LSWE + frictionless on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=legend_xy[1], y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[2]]$model_t/3600, Fujii_Hillarys[[2]]$model_stage, t='l', col='red', lwd=1, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.7, cex.axis=1.5, las=1)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    title(paste0('LSWE + Manning-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[3]]$model_t/3600, Fujii_Hillarys[[3]]$model_stage, t='l', col='red', lwd=1.2, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.4, cex.axis=1.5, las=1)
    title(paste0('LSWE + Linear-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[4]]$model_t/3600, Fujii_Hillarys[[4]]$model_stage, t='l', col='red', lwd=1.2, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.4, cex.axis=1.5, las=1)
    title(paste0('LSWE + Reduced-Linear-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    plot(Fujii_Hillarys[[5]]$model_t/3600, Fujii_Hillarys[[5]]$model_stage, t='l', col='red', lwd=1.2, 
         ylim=plot_ylim, xlim=plot_xlim, 
         xlab='', ylab='', main='', cex.main=1.4, cex.axis=1.5, las=1)
    title(paste0('LSWE + Delayed-Linear-friction on global grid', extra_title), cex.main=1.7, line=0.5)
    points(obs_t, obs_h, t='o', col='black', cex=0.2)
    grid(col='orange')
    legend(x=35, y=legend_xy[2], c('Observed', 'Modelled'), lty=c(1, 1), pch=c(1, NA), 
           pt.cex=c(0.4, 0), col=c('black', 'red'), bty='n', horiz=TRUE, cex=2)

    mtext('Time (hours after earthquake)', side=1, outer=TRUE, cex=1.4, line=2.5)
    mtext('Detided water-surface (m)', side=2, outer=TRUE, cex=1.4, line=2)
    dev.off()
}


# For each inversion, make a plot for every gauge
for(source_ind in 1:length(model_base)){

    # Get data for all gauges and models
    all_RDS = sapply(model_types, f<-function(model_type){
        out = Sys.glob(paste0('../gauge_RDS_files/OUTPUTS/', 
                              model_base[source_ind], '*', 
                              model_type, '*', 
                              model_site[source_ind], 
                              '/RUN*/gauge*.RDS'))
        if(length(out) != 1) stop('Not finding the right number of files')
        return(out)
    })


    event_RDS = lapply(all_RDS, f<-function(x) readRDS(x))
    names(event_RDS) = all_RDS

    # Check there are the same number of gauges for each model
    n0 = length(event_RDS[[1]])
    nall = unlist(lapply(event_RDS, length))
    if(any(nall != n0)) stop('Inconsistent number of gauges in a file')
    site_names = names(event_RDS[[1]])

    # Loop over gauges, one plot for each
    for(h in 1:n0){
        event_site = lapply(event_RDS, f<-function(x) x[[h]])
        plot_ylim = c(-1, 1)*vertical_scale[source_ind]
        plot_filename = paste0(output_dir, '/', model_base[source_ind], site_names[h], '.png')
        make_plot(event_site, plot_filename, plot_ylim, extra_title = paste0(' @ ', site_names[h]) )
    }

    # Loop over gauges, first 10 hours after arrival, case-specific ylimit
    n_hours = 10
    for(h in 1:n0){
        event_site = lapply(event_RDS, f<-function(x) x[[h]])
        model_arrival_ind = min(which(abs(event_site[[1]]$model_stage) > 5e-04*max(abs(event_site[[1]]$model_stage))))
        model_arrival_time = event_site[[1]]$model_t[model_arrival_ind]
        include_inds = which((event_site[[1]]$model_t > model_arrival_time) & 
                             (event_site[[1]]$model_t < (model_arrival_time + n_hours*3600)))
        model_maxima = max(abs(event_site[[1]]$model_stage[include_inds]))

        plot_filename = paste0(output_dir, '/ZOOM_', model_base[source_ind], site_names[h], '.png')
        plot_ylim = c(-1, 1)*1.5 * model_maxima
        plot_xlim = model_arrival_time/3600 + c(0, n_hours) # 10 hours
        make_plot(event_site, plot_filename, plot_ylim, extra_title = paste0(' @ ', site_names[h]) , plot_xlim=plot_xlim)
    }
}

