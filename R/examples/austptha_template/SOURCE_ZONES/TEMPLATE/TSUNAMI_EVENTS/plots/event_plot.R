#
# Plot time-series of the 'optimal' events vs DART buoys
#
PLOT_DURATION_HOURS = as.numeric(commandArgs(trailingOnly=TRUE)[1])
if(is.na(PLOT_DURATION_HOURS)) stop('Must pass numeric number of hours for x axis')

library(rptha)

#'
#' R image files from gauge_summary_statistics.R
#'
all_Rdata = Sys.glob('*.Rdata')

#'
#' lon/lat gauge coordinates, to help make plotting interpretable
#'
fid = nc_open(Sys.glob('../all_uniform_slip_earthquake_events_tsunami*.nc')[1],
    readunlim=FALSE)
model_lonlat = cbind(ncvar_get(fid, 'lon'), ncvar_get(fid, 'lat'))
nc_close(fid)



#' Plot all gauge data and selected model results as time-series
#'
#' Note this function assumes the existence of a number of variables from the 
#' *.Rdata session referenced above. So for it to work, one of those sessions
#' needs to be loaded already.
#'
#' @param ui index of the uniform slip event in 
#'     uniform_slip_stats[[i]][[ui]]. Can be an integer (to plot one event),
#'     or a vector of length = n_gauges (to plot a different event at each gauge).
#' @param si index of the stochastic slip event in 
#'     stochastic_slip_stats[[i]][[ui]]. Can be integer or vector, as for ui.
#' @param vui index of the variable_uniform slip event in 
#'     variable_uniform_slip_stats[[i]][[ui]]. Can be integer or vector, as for ui.
#' @param png_name_stub character for the last part of the png filename
#' @return nothing but make a nice plot
#'
multi_gauge_time_series_plot<-function(ui, si, vui, png_name_stub, output_dir = '.',
    allow_time_offset=TRUE){

    n_gauges = length(uniform_slip_stats)

    if(length(ui) > 1){
        more_than_one_event = TRUE
        stopifnot(length(ui) == n_gauges)
        stopifnot(length(si) == n_gauges)
        stopifnot(length(vui) == n_gauges)
    }else{
        more_than_one_event = FALSE 
        stopifnot(length(ui) == 1) 
        stopifnot(length(si) == 1) 
        stopifnot(length(vui) == 1)

        ui = rep(ui, len = n_gauges)
        si = rep(si, len = n_gauges)
        vui = rep(vui, len = n_gauges)
    }

    
    # Split the plot into nrow/ncol rows and columns
    ncol = (n_gauges > 7) + 1
    nrow = ceiling((n_gauges+1)/ncol)

    # Extract the gauge names from the files
    gauge_names = unlist(lapply(
        strsplit(basename(names(uniform_slip_stats)), split='_'), 
        f<-function(x) substr(x[2], 1, 5)))

    output_png = paste0(output_dir, '/', output_name_base, '_', png_name_stub, '.png')
    png(output_png, width=7.5*ncol, height = max(nrow,5), units='in', 
        res=300)
    par(mfrow=c(nrow, ncol))
    par(mar=c(0.3,1.5,0.5,1.5))
    par(omd=c(0.05, 0.95, 0,1))

    # Make a panel plot for each gauge observation
    for(i in 1:n_gauges){

        uss = uniform_slip_stats[[i]][[ui[i]]]
        sss = stochastic_slip_stats[[i]][[si[i]]]
        svu = variable_uniform_slip_stats[[i]][[vui[i]]]

        # Choose the x-limits of the plot based on the stochastic slip x-limits
        xmin = min(c(sss$data_t, sss$model_t))
        xlim = c(xmin, xmin+PLOT_DURATION_HOURS*3600)

        ymin = min(c(uss$data_s, uss$model_s, sss$data_s, sss$model_s, 
            svu$model_s, svu$data_s))
        ymax = max(c(uss$data_s, uss$model_s, sss$data_s, sss$model_s, 
            svu$model_s, svu$data_s))
        ylim = c(ymin, ymax)

        plot(uss$data_t, uss$data_s, t='o', col='black', pch=19, 
            cex=0.4, xlim=xlim, ylim=ylim, axes=FALSE)
        points(sss$data_t, sss$data_s, t='o', col='black', pch=19, 
            cex=0.4)
        points(svu$data_t, svu$data_s, t='o', col='black', pch=19, 
            cex=0.4)
        points(uss$model_t - uss$model_time_offset*allow_time_offset, uss$model_s, 
            t='l', col='blue', lwd=2, lty='longdash')
        points(sss$model_t - sss$model_time_offset*allow_time_offset, sss$model_s, 
            t='l', col='red', lwd=1)
        points(svu$model_t - svu$model_time_offset*allow_time_offset, svu$model_s, t='l', 
            col='green', lwd=1)
        abline(v = seq(xmin + 1, xmin + 7*3600, len=7), col='tan2', xpd=TRUE)
        abline(h=0, col='tan2')

        tmp = pretty(ylim) 
        axis(side=4, at = c(min(tmp[tmp>ylim[1]]), max(tmp[tmp<=ylim[2]])), 
            las=1, line=-3, cex.axis=2.0)
        text(xlim[2]-1800, ylim[2]*0.8, labels=gauge_names[i], cex=3)
    }

    # Add a legend
    if(more_than_one_event){
        # Simple legend [do not show event numbers, since it takes too much space]
        plot(c(0,1), c(0,1), axes=FALSE, col=0)
        legend('center', 
            c('Data', 'Stoch. slip', 'Unif. slip ', 'VarU. slip'), 
            bty='n',
            lty=c('solid', 'solid', 'longdash', 'solid'),
            col=c('black', 'red', 'blue', 'green'), 
            pch=c(19, NA, NA, NA), cex= 2.0 + (nrow > 2),
            pt.cex=1,
            ncol=2) 
    }else{
        # Legend including event numbers
        plot(c(0,1), c(0,1), axes=FALSE, col=0)
        legend('center', 
            c('Data', paste0('Stoch. slip (', si[1], ')'), 
                paste0('Unif. slip (', ui[1], ')'), 
                paste0('VarU. slip (', vui[1], ')')), 
            bty='n',
            lty=c('solid', 'solid', 'longdash', 'solid'),
            col=c('black', 'red', 'blue', 'green'), 
            pch=c(19, NA, NA, NA), cex= 2.0 + (nrow > 2),
            pt.cex=1,
            ncol=2) 
    }

    dev.off()
}

#' Various plots comparing uniform and variable slip events with NGDC data
#'
#' The function assumes the existence of ....
#'
#' @param 
ngdc_comparison_plot<-function(ui, si, vui, png_name_stub, output_dir = '.'){

    output_png = paste0(output_dir, '/', output_name_base, '_', 
        png_name_stub, '.png')
    
    png(output_png, width=15, height = 10, units='in', 
        res=200)
    #par(mfrow=c(2,3))
    layout(matrix(c(1, 1, 4, 7,
                    2, 2, 5, 8,
                    3, 3, 6, 9),
        ncol=4, nrow=3, byrow=TRUE))

    par(mar=c(2,2,2,2))

    # A useful color scheme
    cols10 = c(c('hotpink', 'orange'), rev(rainbow(6))[3:6], terrain.colors(4))

    # Convenience function to make a map with the observed data
    # Use sqrt vertical scale
    map_panel_plot<-function(model_data, ind, titlewords){
        # First panel -- spatial plot of locations
        ngdc_lonlat = cbind(
            model_data$tsunami_obs$LONGITUDE,
            model_data$tsunami_obs$LATITUDE)

        # Ensure the longitudes are in the same range (e.g. -180-180 or 0-360)
        ngdc_lonlat = adjust_longitude_by_360_deg(ngdc_lonlat, 
            model_data$gauge_lonlat[,1:2])

        # Only compare points with 'distance to nearest model result < 20000'
        kp = which(model_data$distance_to_nearest < 20000)

        if(length(kp) == 0){
            plot(c(0,1), c(0,1), main='No points within distance')
            return(invisible())
        }
        
        # Apply sqrt transformation, to make it easier to see a range of scales
        ht = sqrt(model_data$tsunami_obs$WATER_HT[kp])
        ms  = sqrt(model_data$max_stage[ind,kp] * model_data$gcf_mat[ind,kp])

        # Maybe scale arrows to improve visibility
        tmp = max(diff(range(ngdc_lonlat[kp,1])), diff(range(ngdc_lonlat[kp,2])))
        # By default, largest vertical bar takes up 1/5 of vertical range
        vb = 0.2
        scaler = max(1, tmp * vb/max(c(max(ht), max(ms))) )
        
        plot_y_range = range(ngdc_lonlat[kp,2])
        plot_y_range[2] = plot_y_range[2] + vb * diff(plot_y_range)

        plot(ngdc_lonlat[kp,], asp=1, xlab="Lon", ylab="Lat", pch='.',
            ylim=plot_y_range)
        # Add all model points -- helps to see coast, etc
        points(model_lonlat[,1:2], pch='.')

        arrows(ngdc_lonlat[kp,1], ngdc_lonlat[kp,2], ngdc_lonlat[kp,1], 
            ngdc_lonlat[kp,2] + ht*scaler, length=0, lwd=2, col='black')
        # Make the model results thinner, and color by the measurement type
        arrows(ngdc_lonlat[kp,1],  ngdc_lonlat[kp,2], ngdc_lonlat[kp,1], 
            ngdc_lonlat[kp,2] + ms*scaler, length=0, lwd=0.5, 
            col=cols10[model_data$tsunami_obs$TYPE_MEASUREMENT_ID[kp]])
        title(titlewords, cex.main=1.5)

    }

    # Maps for all events
    map_panel_plot(NGDC_comparison$uniform, ui, 
        paste0('Uniform event ', ui))
    map_panel_plot(NGDC_comparison$stochastic, si, 
        paste0('Stochastic event ', si))
    map_panel_plot(NGDC_comparison$variable_uniform, vui, 
        paste0('Variable_uniform event ', vui))

    # Convenience function to plot predicted-vs-measured
    scatter_panel_plot<-function(model_data, ind, titlewords){
       
        # Only keep obs within 20km of model point 
        kp = which(model_data$distance_to_nearest < 20000)
    
        plot(c(0.01, 100), c(0.01, 100), log='xy', asp=1, xlab='Predicted',
            ylab='Measured', col=0)
        pred = model_data$max_stage[ind,kp]*model_data$gcf_mat[ind,kp]
        meas = model_data$tsunami_obs$WATER_HT[kp]
        points(pred, meas, pch=19, 
            col=cols10[model_data$tsunami_obs$TYPE_MEASUREMENT_ID[kp]])

        abline(h=10**(seq(-2,2)), col='orange', lty='dotted')
        abline(v=10**(seq(-2,2)), col='orange', lty='dotted')
        abline(0,1,col='red')

        # Only report statistics for tidal gauges and deep ocean gauges
        kk = which(model_data$tsunami_obs$TYPE_MEASUREMENT_ID[kp] %in% c(2,3))

        title(paste0('Gauge med(p/m): ', round(median(pred[kk]/meas[kk]), 2),
            ', med(|p-m|/m): ', round(median(abs(pred[kk]-meas[kk])/meas[kk]), 2)), 
            cex.main=1.5)

        legend('topleft', legend = as.character(1:10), pch=rep(19, 10), 
            col=cols10, bty='n')
    }

    scatter_panel_plot(NGDC_comparison$uniform, ui, 
        paste0('Uniform event ', ui))
    scatter_panel_plot(NGDC_comparison$stochastic, si, 
        paste0('Stochastic event ', si))
    scatter_panel_plot(NGDC_comparison$variable_uniform, vui, 
        paste0('Variable_uniform event ', vui))

    # Final plot -- slip rasters
    slip_rast1 = make_slip_raster(1, uniform_slip_stats[[1]][[ui]]$events_with_Mw, 
        unit_source_statistics)
    image(slip_rast1$slip_rast, xlim=slip_rast1$xlim, col=rev(terrain.colors(255)))
    slip_rast1 = make_slip_raster(1, stochastic_slip_stats[[1]][[si]]$events_with_Mw, 
        unit_source_statistics)
    image(slip_rast1$slip_rast, xlim=slip_rast1$xlim, col=rev(terrain.colors(255)))
    slip_rast1 = make_slip_raster(1, variable_uniform_slip_stats[[1]][[vui]]$events_with_Mw, 
        unit_source_statistics)
    image(slip_rast1$slip_rast, xlim=slip_rast1$xlim, col=rev(terrain.colors(255)))

    dev.off()
}

#
# Make lots of gauge comparison plots, and some summary statistics stuff.
#

for(RdataFile in all_Rdata){

    # Load the Rimage associated with gauge_summary_statistics.R
    attach(RdataFile)

    #'
    #' Files from NGDC comparison
    #'
    all_NGDC_comparison_RDS = Sys.glob(paste0('../', output_name_base, '*/event_NGDC_comparison.RDS'))

    # Order the NGDC files into 'uniform', 'stochastic', 'variable_uniform'
    order_NGDC_RDS = c(
        grep('uniform_uniform', all_NGDC_comparison_RDS),
        grep('stochastic_stochastic', all_NGDC_comparison_RDS),
        grep('variable_uniform_variable_uniform', all_NGDC_comparison_RDS)
        )
    stopifnot(
        (length(order_NGDC_RDS) == length(all_NGDC_comparison_RDS)) &
        (length(unique(order_NGDC_RDS)) == length(all_NGDC_comparison_RDS))
        )
    all_NGDC_comparison_RDS = all_NGDC_comparison_RDS[order_NGDC_RDS]
    # Read the NGDC RDS files into a list, with suitable names
    NGDC_comparison = lapply(as.list(all_NGDC_comparison_RDS), readRDS)
    names(NGDC_comparison) = c('uniform', 'stochastic', 'variable_uniform')

    # Make plots of 'optimal' model based on 3 different goodness of fit
    # measures
    for(stat_measure in 1:3){

        # Get the 'best' stochastic and uniform events
        if(stat_measure == 1){
            # Time-domain similarity statistic
            stat_name = 'time'
            stoc_score = do.call(cbind, similar_s_time)
            unif_score = do.call(cbind, similar_u_time)
            vu_score = do.call(cbind, similar_vu_time)
        }else if(stat_measure == 2){
            # Spectral domain similarity statistic
            stat_name = 'spec'
            stoc_score = do.call(cbind, similar_s_spec)
            unif_score = do.call(cbind, similar_u_spec)
            vu_score = do.call(cbind, similar_vu_spec)
        }else if(stat_measure == 3){
            # Combined Time-domain and spectral domain statistic
            stat_name = 'spectime'
            stoc_score = 0.8*do.call(cbind, similar_s_time) + 
                0.2*do.call(cbind, similar_s_spec)
            unif_score = 0.8*do.call(cbind, similar_u_time) + 
                0.2*do.call(cbind, similar_u_spec)
            vu_score = 0.8*do.call(cbind, similar_vu_time) + 
                0.2*do.call(cbind, similar_vu_spec)
        }

        # Find the median GoF statistic over all gauges
        stoc_score_median = apply(stoc_score, 1, median)
        unif_score_median = apply(unif_score, 1, median)
        vu_score_median = apply(vu_score, 1, median)

        # Find the 'best' modelled events, using the minimum(median of the
        # goodness-of-fit metric at all gauges) as the definition of 'best'
        si = which.min(stoc_score_median)
        ui = which.min(unif_score_median)
        vui = which.min(vu_score_median)
       
        #
        # Plot the single 'overall best' model event at all gauges.
        # 
        output_png = paste0(output_name_base, '_allmodel_goodness_fit_', 
            stat_name, '_plot.png')
        png(output_png, width=13, height = 5, units='in', res=300)
        par(mfrow=c(1,3))
        hist(stoc_score_median, freq=FALSE, 
            main=paste0('Stochastic slip model median GoF: ', stat_name), 
            xlim=c(-0.1,1.1), col='red', density=20)
        hist(unif_score_median, freq=FALSE, 
            main=paste0('Uniform slip model median GoF: ', stat_name), 
            xlim=c(-0.1,1.1), col='blue', density=20)
        hist(vu_score_median, freq=FALSE, 
            main=paste0('Variable uniform slip model median GoF: ', 
                stat_name), 
            xlim=c(-0.1,1.1), col='green', density=20)
        dev.off()

        # Plot up as gauges
        multi_gauge_time_series_plot(ui, si, vui, 
            png_name_stub=paste0('best_fit_', stat_name, '_gauges_plot'),
            allow_time_offset=(stat_name != 'spec')
            )
        ngdc_comparison_plot(ui, si, vui, 
            png_name_stub=paste0('best_fit_', stat_name, '_NGDC_plot')
            )

        #
        # Plot the 'gauge specific best' event at each gauge.  
        #
        si = apply(stoc_score, 2, which.min)
        ui = apply(unif_score, 2, which.min)
        vui = apply(vu_score, 2, which.min)
        multi_gauge_time_series_plot(ui, si, vui, 
            png_name_stub=paste0('LOCAL_best_fit_', stat_name, '_gauges_plot'),
            allow_time_offset=(stat_name != 'spec')
            )

    }

    #
    # Make gauges plots for EVERY event 
    #

    # Devise a mapping between the uniform slip events (no variation), and the
    # corresponding stochastic and variable_uniform events
    uniform_event_rows = event_metadata$events_with_Mw$uniform_event_row

    stopifnot(
        length(uniform_event_rows) == length(stochastic_slip_stats[[1]]))
    stopifnot(
        length(uniform_event_rows) == length(variable_uniform_slip_stats[[1]]))
    stopifnot(
        length(unique(uniform_event_rows)) == length(uniform_slip_stats[[1]]))

    variable_to_uniform = match(uniform_event_rows, unique(uniform_event_rows))

    output_dir = 'event_fig'
    dir.create(output_dir, showWarnings=FALSE)

    for(i in 1:length(variable_to_uniform)){
        # Indices of corresponding uniform, stochastic, and variable uniform
        ui = variable_to_uniform[i]
        si = i
        vui = i
        multi_gauge_time_series_plot(ui, si, vui, 
            png_name_stub=paste0('event_', ui, '_', i, '_gauges_plot'),
            output_dir=output_dir
            )
        ngdc_comparison_plot(ui, si, vui, 
            png_name_stub=paste0('event_', ui, '_', i, '_NGDC_plot'),
            output_dir=output_dir
            )

    }
    

    # Remove nearly all objects, except the ones controlling the outer loop
    rm(list=setdiff( ls(all=TRUE),
            c('RdataFile', 'all_Rdata', 'multi_gauge_time_series_plot', 'ngdc_comparison_plot',
                'model_lonlat')
        ))

}
