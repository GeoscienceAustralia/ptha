#
# Plot
#
library(rptha)

source('parse_gauge_outputs.R')

#
# Preliminaries. 
#

arrival_time_analysis<-function(){
    #
    # Look at the time of arrival/maxima for these events.
    #

    # Plot arrival time & time of maxima
    png('Times_of_arrivals_and_maxima.png', width=10, height=9, units='in', res=300)
    k = which(event_stats$good_nearshore & event_stats$is_inversion)
    par(mar=c(21, 7, 3, 2))
    boxplot(event_stats$model_time_to_arrive[k]/3600 ~ event_stats$site_and_event[k], 
        log='', ylim=c(0, 60), drop=TRUE, las=2, xlab="", ylab="Hours after earthquake",
        main='Time of maximum observed tsunami for each event/site combination & \n modelled tsunami arrival times',
        cex.main=1.5)
    boxplot(event_stats$obs_max_time_to_arrive[k]/3600 ~ event_stats$site_and_event[k],
        add=TRUE, col='red', border='red', las=2, xlab="", ylab="", names.arg="", ann=FALSE)
    abline(h=seq(10, 60, by=10), lty='dotted')
    legend('topleft', c('Time of observed tsunami maxima', 'Modelled tsunami arrival times'), 
           col=c('red', 'black'), lwd=c(2,2), bg='white')
    dev.off()

    # Check some key statistics at the 'same-site, same-event'
    # We will have different 'arrival times' because of the nominal way that was defined.
    arrival_time_range = aggregate(event_stats$model_time_to_arrive[k], 
        by=list(site_and_event=event_stats$site_and_event[k]), 
        f<-function(x) diff(range(x)))
    # Check the time of the observed maxima is identical for all inversions when
    # the same data is used.
    # (Clearly it should be -- but because the model was used to define the tsunami
    # "arrival time" it is possible that they are not, e.g. due to very noisy early
    # observation series, or error in model selection code, or strong differences in inversions)). 
    data_max_range = aggregate(event_stats$data_max[k],
        by=list(site_and_event=event_stats$site_and_event[k]), 
        f<-function(x) diff(range(x)))
    if(any(data_max_range$x > 0, na.rm=TRUE) ){
        stop('Observed maxima varies with the same data -- issue with model arrival times')
    }

    # Print some summaries of the arrival times that we use in the paper
    print('Median arrival time nearshore, inversions')
    print(median(event_stats$model_time_to_arrive[k]/3600))
    print('Range arrival time nearshore, inversions')
    print(range(event_stats$model_time_to_arrive[k]/3600))
}
arrival_time_analysis()
#[1] "Median arrival time nearshore, inversions"
#[1] 12.16667
#[1] "Range arrival time nearshore, inversions"
#[1]  5.77500 16.01667
    

initial_plots_not_used<-function(){
    # Generic plot of model-range vs data-range. Better plots are below.
    #
    # Note there is a tag for 'data we like to highlight in the plot'. TRUE values get a
    # bigger symbol.
    panel_plot<-function(model_range2, data_range2, tag, col, 
                         XYLIM=c(0.01, 1.3), xlab='model_range_stat', 
                         ylab='data_range_stat', title='Default title',
                         ...){
    
        plot(model_range2, data_range2, 
             ylim=XYLIM, xlim=XYLIM, pch=c('.', 'x')[tag+1], cex=c(0.1, 1)[tag+1],
             col=col, xlab='Model max', ylab='Data max',
             main=title, ... )
        grid()
        rng = seq(XYLIM[1], XYLIM[2], len=1000)
        points(rng, rng, t='l', col='red')
        points(rng, 2*rng, t='l', col='red', lty='dashed')
        points(rng, 0.5*rng, t='l', col='red', lty='dashed')
        points(rng, 1.5*rng, t='l', col='red', lty='dashed')
        points(rng, (1/1.5)*rng, t='l', col='red', lty='dashed')
    }
    
    png('Maxima_and_minima_over_time.png', width=12, height=8, units='in', res=300)
    # 
    # This plot shows the model-vs-data variability (MAX and MIN) for the full 60 hrs,
    # and the first 8 hours after modelled tsunami arrival. Note the latter duration
    # might be quite similar to the Adams and LeVeque study that compared with MOST.
    #
    # It is notable that the linear model is closer to the other models at shorter times.
    #
    tag = event_stats$good_nearshore & event_stats$is_inversion & (!is.na(event_stats$model_type))
    col = as.numeric(as.factor(model_type))
    par(mfrow=c(2,2))
    par(mar=c(4,4,2,1))
    panel_plot(event_stats$model_max, event_stats$data_max, tag=tag, col=col,
               xlab='Model MAX after arrival', ylab='Data MAX after model arrival',
               title = 'Tsunami MAX after model arrival', cex.main=1.5, log='xy')
    panel_plot(-event_stats$model_min, -event_stats$data_min, tag=tag, col=col,
               xlab='Model MIN after arrival', ylab='Data MIN after model arrival',
               title = 'Tsunami MIN after model arrival *(-1)', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_8hrs, event_stats$data_max_8hrs, tag=tag, col=col,
               xlab='Model MAX within 8hrs of model arrival ', ylab='Data MAX within 8hrs of model arrival',
               title = 'Tsunami MAX within 8hrs of model arrival ', cex.main=1.5, log='xy')
    panel_plot(-event_stats$model_min_8hrs, -event_stats$data_min_8hrs, tag=tag, col=col,
               xlab='Model MIN within 8hrs of model arrival ', ylab='Data MIN within 8hrs of model arrival',
               title = 'Tsunami MIN within 8hrs of model arrival *(-1)', cex.main=1.5, log='xy')
    dev.off()
    
    png('Maxima_over_time.png', width=12, height=8, units='in', res=300)
    # This plot shows the model-vs-data variability (MAX and MIN) for the full 60 hrs,
    # and the first 8 hours after modelled tsunami arrival. Note the latter duration
    # might be quite similar to the Adams and LeVeque study that compared with MOST.
    #
    # It is notable that the linear model is closer to the other models at shorter times.
    tag = event_stats$good_nearshore & event_stats$is_inversion & (!is.na(event_stats$model_type))
    col = as.numeric(as.factor(model_type))
    par(mfrow=c(2,2))
    par(mar=c(4,4,2,1))
    panel_plot(event_stats$model_max, event_stats$data_max, tag=tag, col=col,
               xlab='Model MAX after arrival', ylab='Data MAX after arrival',
               title = 'Tsunami MAX after model arrival', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_24, event_stats$data_max_24, 
               xlab='Model MAX 24 hr', ylab='Data MAX 24 hrs', tag=tag, col=col,
               title = 'Tsunami MAX after model arrival, 24 hr simulation', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_12hrs, event_stats$data_max_12hrs, 
               xlab='Model MAX 12hrs ', ylab='Data MAX 12hrs ',tag=tag, col=col,
               title = 'Tsunami MAX, 12hrs after model arrival ', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_8hrs, event_stats$data_max_8hrs, tag=tag, col=col,
               xlab='Model MAX 8hrs ', ylab='Data MAX 8hrs ',
               title = 'Tsunami MAX, 8hrs after model arrival ', cex.main=1.5, log='xy')
    dev.off()
    
    png('Maxima_over_time2.png', width=12, height=8, units='in', res=300)
    par(mfrow=c(2,2))
    par(mar=c(4,4,2,1))
    tag = event_stats$good_nearshore & event_stats$is_inversion & (!is.na(event_stats$model_type))
    col = as.numeric(as.factor(model_type))
    panel_plot(event_stats$model_max, event_stats$data_max, tag=tag, col=col,
               xlab='Model MAX 60 hrs', ylab='Data MAX 60 hrs',
               title = 'Tsunami MAX after model arrival, 60 hr simulation', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_36, event_stats$data_max_36, tag=tag, col=col,
               xlab='Model MAX after model arrival, 36 hr simulation', 
               ylab='Data MAX after model arrival, 36 hr simulation',
               title = 'Tsunami MAX after model arrival, 36 hr simulation', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_24, event_stats$data_max_24, tag=tag, col=col,
               xlab='Model MAX 24 hrs', ylab='Data MAX FIRST 24 hrs',
               title = 'Tsunami MAX after model arrival, 24 hr simulation', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_lastday, event_stats$data_max_lastday, tag=tag, col=col,
               xlab='Model MAX FINAL 24 hrs', ylab='Data MAX FINAL 24 hrs',
               title = 'Tsunami MAX 36-60 hours after earthquake', cex.main=1.5, log='xy')
    dev.off()
    
    
    #
    # Let's look at this kind of thing, for DARTS ONLY
    #
    png('Maxima_over_time2_DARTS.png', width=12, height=8, units='in', res=300)
    dart_tag = grepl('DART', event_stats$sites)
    col = as.numeric(as.factor(model_type))
    par(mfrow=c(2,2))
    par(mar=c(4,4,2,1))
    panel_plot(event_stats$model_max, event_stats$data_max, 
               tag=dart_tag, col=col, XYLIM=c(0.01, 2),
               xlab='Model MAX after model arrival', ylab='Data MAX after model arrival',
               title = 'Tsunami MAX after model arrival', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_36, event_stats$data_max_36, 
               tag=dart_tag, col=col, XYLIM=c(0.01, 2),
               xlab='Model MAX FIRST 12 hrs from arrival', 
               ylab='Data MAX FIRST 12 hrs from arrival',
               title = 'Tsunami MAX after model arrival, 36 hr sim', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_24, event_stats$data_max_24, 
               tag=dart_tag, col=col, XYLIM=c(0.01, 2),
               xlab='Model MAX FIRST 24 hrs', ylab='Data MAX FIRST 24 hrs',
               title = 'Tsunami MAX after model arrival, 24 hr sim', cex.main=1.5, log='xy')
    panel_plot(event_stats$model_max_lastday, event_stats$data_max_lastday, 
               tag=dart_tag, col=col, XYLIM=c(0.01, 1),
               xlab='Model MAX FINAL 24 hrs', ylab='Data MAX FINAL 24 hrs',
               title = 'Tsunami MAX 36-60 hours after earthquake', cex.main=1.5, log='xy')
    dev.off()
}

#
# Plot the maxima/minima results on a per-model basis, for good nearshore data,
# and separately for high-res DART data.
#

model_types = c('Frictionless', 'Manning0.035', 'LinearFriction', 'LinearReducedFriction', 
                'LinearDelayedFriction')
#plot_types = c('nearshore', 'DART')
plot_types = c('nearshore')

plot_model_vs_data<-function(model_types, plot_types){
    # Make many plots, and extract 'nearshore data' and 'DART data' for each model type
    # that will be useful for later analysis.

    for(pt in plot_types){

        # Depending on the plot-type, get a data-subset with either the nearshore 1min
        # data, or the DART highres data, for each model type
        model_stats = lapply(model_types,  
            f<-function(x, pt){
                if(pt == 'nearshore'){
                    # Get the good quality 1min data nearshore
                    keep = which(event_stats$model_type == x & 
                                 event_stats$good_nearshore & 
                                 event_stats$is_inversion)
                }else if(pt == 'DART'){
                    # Get the high-res dart data 1min data nearshore
                    # Should probably avoid sites where the waves are very very small
                    # because results are dominated by other stuff.
                    keep = which(event_stats$model_type == x & 
                                 event_stats$good_dart &
                                 event_stats$is_inversion)
                }else{
                    stop('pt unrecognized')
                }
                                    
                event_stats[keep, ]
                }, 
            pt=pt)
        names(model_stats) = model_types

        # Store the data
        if(pt == 'nearshore'){
            model_stats_nearshore = model_stats
        }else if(pt == 'DART'){
            model_stats_DART = model_stats
        }else{
            stop('unknown plot type')
        }

        # Make a bunch of model-vs-data plots. 
        stat_types = c('_max', '_max_36', '_max_8hrs', '_max_12hrs', '_max_24', '_max_lastday', 
                       '_min', '_min_36', '_min_8hrs', '_min_12hrs', '_min_24', '_min_lastday',
                       '_stgrng', '_stgrng_36', '_stgrng_8hrs', '_stgrng_12hrs', '_stgrng_24', '_stgrng_lastday'
                       )
        for(st in stat_types){

            if(pt == 'nearshore'){
                filename_prefix = 'Model_predictive_performance_nearshore'
            }else if(pt == 'DART'){
                filename_prefix = 'Model_predictive_performance_DARThighres'
            }else{
                stop('pt unrecognized')
            }

            png(paste0(filename_prefix, st, '.png'), width=13, height=8, units='in', res=300)
            par(mfrow=c(2,3))
            par(mar=c(4,4,4,1))
            for(mt in model_types){
                model = model_stats[[mt]][[ paste0('model', st) ]]
                data = model_stats[[mt]][[ paste0('data', st) ]]
                cols = model_stats[[mt]]$event_name_int
                event_and_source_int = model_stats[[mt]]$event_and_source_int
                pchs = event_and_source_int
                plot(model, data, main="", xlab='Model', ylab='Data', col=cols, pch=pchs)
                title(main=mt, line=2)
                ## Find the median of (inversion median gauge error)
                mre  = median(aggregate(    (model-data)/abs(data), by=list(event_and_source_int), 
                                        f<-function(x) median(x, na.rm=TRUE))$x)
                mare = median(aggregate( abs(model-data)/abs(data), by=list(event_and_source_int), 
                                        f<-function(x) median(x, na.rm=TRUE))$x)
                title(main= paste0('Median ( Median    (relative error per inversion) ): ', round(mre, 3)), 
                      line=1, col.main='blue')
                title(main= paste0('Median ( Median abs(relative error per inversion) ): ', round(mare, 3)), 
                      line=0, col.main='blue')
                grid()
                abline(0,1)
                abline(0, 1.5, lty='dashed')
                abline(0, 1/1.5, lty='dashed')
                legend('bottomright', unique(model_stats[[mt]]$event_name), 
                       col=unique(model_stats[[mt]]$event_name_int), 
                       pch=19, bty='n')
            }
            dev.off()

            png(paste0(filename_prefix, st, '_relativeErr.png'), width=13, height=8, units='in', res=300)
            par(mfrow=c(2,3))
            par(mar=c(4,4,4,1))
            for(mt in model_types){
                model = model_stats[[mt]][[ paste0('model', st) ]]
                data = model_stats[[mt]][[ paste0('data', st) ]]
                cols = model_stats[[mt]]$event_name_int
                event_and_source_int = model_stats[[mt]]$event_and_source_int
                pchs = event_and_source_int
                YLIM = c(-1,1)*min(3, max(2, max(abs(model-data)/data, na.rm=TRUE)))
                #print(paste(c('YLIM[1]: ', YLIM[1], ' ', mt)))
                plot(abs(data), (model-data)/data, main="", xlab='Data', ylab='(Model-Data)/Data', 
                     col=cols, pch=pchs, ylim=YLIM, log='x')
                title(main=mt, line=2)
                ## Find the median of (inversion median gauge error)
                mre  = median(aggregate(    (model-data)/abs(data), by=list(event_and_source_int), 
                                        f<-function(x) median(x, na.rm=TRUE))$x)
                mare = median(aggregate( abs(model-data)/abs(data), by=list(event_and_source_int), 
                                        f<-function(x) median(x, na.rm=TRUE))$x)
                title(main= paste0('Median ( Median    (relative error per inversion) ): ', round(mre, 3)), 
                      line=1, col.main='blue')
                title(main= paste0('Median ( Median abs(relative error per inversion) ): ', round(mare, 3)), 
                      line=0, col.main='blue')
                grid()
                abline(h=0)
                abline(h=0.5, lty='dashed')
                abline(h=-0.5, lty='dashed')
                legend('bottomright', unique(model_stats[[mt]]$event_name), 
                       col=unique(model_stats[[mt]]$event_name_int), 
                       pch=19, bty='n')
            }
            dev.off()
        }

    }
    #return(list(model_stats_nearshore=model_stats_nearshore, model_stats_DART=model_stats_DART))
    return(list(model_stats_nearshore=model_stats_nearshore, model_stats_DART=NULL))
}
output = plot_model_vs_data(model_types, plot_types)
model_stats_nearshore = output$model_stats_nearshore
#model_stats_DART = output$model_stats_DART
rm(output)

boxplots_and_alternatives<-function(model_stats_nearshore){
    # More graphical display of the model errors.
    #
    # For each inversion, the median error (over all gauges) gives an idea of the
    # bias, but the accuracy will depend on the number of gauge observations (varying per inversion).
    # Here we use statistics that attempt to weight each source-inversion equally.
    nearshore_median_max_err_by_event = lapply(model_stats_nearshore, 
           f<-function(x) aggregate(((x$model_max-x$data_max)/x$data_max), 
                                    by=list(x$event_and_source_int), f<-function(y) median(y, na.rm=TRUE))$x)
    nearshore_median_max_36_err_by_event = lapply(model_stats_nearshore, 
           f<-function(x) aggregate(((x$model_max_36-x$data_max_36)/x$data_max_36), 
                                    by=list(x$event_and_source_int), f<-function(y) median(y, na.rm=TRUE))$x)
    nearshore_median_max_24_err_by_event = lapply(model_stats_nearshore, 
           f<-function(x) aggregate(((x$model_max_24-x$data_max_24)/x$data_max_24), 
                                    by=list(x$event_and_source_int), f<-function(y) median(y, na.rm=TRUE))$x)
    png('Median_error_proportion_for_each_event_BOXPLOT.png', width=12, height=9, units='in', res=300)
    par(mfrow=c(3,1))
    boxplot(nearshore_median_max_err_by_event, 
            main='Tsunami maxima (60h): Median (error/obs), distribution of all inversions',
            notch=TRUE, cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    boxplot(nearshore_median_max_36_err_by_event, 
            main='Tsunami maxima (36h): Median (error/obs), distribution of all inversions',
            notch=TRUE, cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    boxplot(nearshore_median_max_24_err_by_event, 
            main='Tsunami maxima (24h): Median (error/obs), distribution of all inversions',
            notch=TRUE, cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    dev.off()
    png('Median_error_proportion_for_each_event_VIOPLOT.png', width=12, height=9, units='in', res=300)
    library(vioplot)
    par(mfrow=c(3,1))
    vioplot(nearshore_median_max_err_by_event, 
            main='Tsunami maxima (60h): Distribution of [Median (error/obs)] for all inversions',
            cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    vioplot(nearshore_median_max_36_err_by_event, 
            main='Tsunami maxima (36h): Distribution of [Median (error/obs)] for all inversions',
            cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    vioplot(nearshore_median_max_24_err_by_event, 
            main='Tsunami maxima (24h): Distribution of [Median (error/obs)] for all inversions',
            cex.main=1.9, cex.axis=1.5, ylim=c(-0.8, 0.8))
        grid(); abline(h=0)
    dev.off()
}
boxplots_and_alternatives(model_stats_nearshore)

#
# Statistics if we have "random inversion, and random gauge"
#

check_model_type_data_alignment<-function(model_stats){
    # Check that the model-specific data is ordered in the same way for all model types
    # For each model type, the 'event_name_and_source' and 'site_and_event' columns should be the same.
    # That way we can compare each model type without reordering the data 

    local_match = unlist(
        lapply(model_stats, f<-function(x){
            all(x$event_name_and_source == model_stats[[1]]$event_name_and_source) &
            all(x$event_and_source_int == model_stats[[1]]$event_and_source_int) &
            all(x$site_and_event == model_stats[[1]]$site_and_event) }
         ))
    print(local_match)
    if(! all(local_match)){

        stop('Model-specific data is not aligned in the same way, cannot proceed with comparison')
    }else{
        print('PASS')
    }
}
check_model_type_data_alignment(model_stats_nearshore)
#check_model_type_data_alignment(model_stats_DART)


#
#
#

source('quantile_CI.R')
alternative_summary_statistics<-function(model_stats_nearshore, stat_type = '_max'){
    # Summarise the "relative error" [ (model - data)/data ] at gauges in various ways:
    # - Treating all gauges as the same
    # - Weighting all scenarios the same, and accounting for the unequal number of gauges per scenario by:
    #     - Taking the median_per_scenario( median_per_gauge_for_each_scenario)
    #     - Creating a 'random sample of gauges' with equal representation for each scenario, and computing statistics on that
    # This function relies on all the entries of model_stats_nearshore having gauges ordered in the same way, which was checked above
    # by check_model_type_data_alignment

    # Statistics 
    st = stat_type # Choose the statistic 
    modstat = paste0('model', st)
    datstat = paste0('data',  st)

    # Make a large sample of 'random scenario with single random gauge', where
    # each scenario is equally likely to occur.

    mt = 'Frictionless' 
    random_event_and_source = sample(
        unique(model_stats_nearshore[[mt]]$event_and_source_int), 
        replace=TRUE, size=100000) # Large random sample of scenarios
    # The checks above ensure we can define random_event_and_source no matter the model type
    # (i.e mt='Manning0.035' or another valid value is equivalent).
    # Similarly we can pick a random gauge no matter the model type.
    random_gauge = sapply(random_event_and_source, 
       f<-function(x){
           # Find indices with valid data values, where event_and_source_int matches what we want.
           inds = which( (model_stats_nearshore[[mt]]$event_and_source_int == x) & 
                         (is.finite(model_stats_nearshore[[mt]][[datstat]])))
           # Pick a random gauge
           if(length(inds) == 1){
               # Work around 'sample' behaviour if the input has length 1
               output = inds
           }else{
               # Regular case
               output = sample( (1:nrow(model_stats_nearshore[[mt]]))[inds], size=1)
           }
           return(output)
       } )

    for(mt in model_types){
        # For the given statistic, we will print various summaries of
        #    median( (model-data)/data) 
        # and
        #    median( abs( (model-data)/data) )
        # for the random inverions / gauges sampled above

        # Relative model error
        stat = (model_stats_nearshore[[mt]][[modstat]] - model_stats_nearshore[[mt]][[datstat]])/
            model_stats_nearshore[[mt]][[datstat]]

        # Statistics that treat all gauges as the same (i.e. no adjustment for the differing number of gauges per event)
        median_stat = median(stat[random_gauge])
        abs_median_stat = median(abs(stat[random_gauge]))

        # Statistics that work with median_over_scenarios(median_over_gauges(...))
        median_stat_per_scenario = aggregate(stat,
            by=list(model_stats_nearshore[[mt]]$event_and_source_int), 
            f<-function(x) median(x, na.rm=TRUE))$x
        median_abs_stat_per_scenario = aggregate(abs(stat),
            by=list(model_stats_nearshore[[mt]]$event_and_source_int), 
            f<-function(x) median(x, na.rm=TRUE))$x
        median_median_per_scenario = median(median_stat_per_scenario)
        median_median_per_scenario_95CI = quantile_CI(median_stat_per_scenario, Quant=0.5, level=0.95)
        median_median_abs_per_scenario = median(median_abs_stat_per_scenario)
        median_median_abs_per_scenario_95CI = quantile_CI(median_abs_stat_per_scenario, Quant=0.5, level=0.95)

        # # Jacknife error estimate -- how does the result vary if we leave-one-scenario-out?
        # # BEWARE: This is a poor way to estimate uncertainty in the median and many other statistics,
        # # although in a small sample like this it seems like a pretty intuitive sensitivity test.
        # median_stat_jknf = rep(NA, 12)
        # abs_median_stat_jknf = rep(NA, 12)
        # for(i in 1:12){
        #     gauges_to_keep = setdiff(1:12, i) #sample(1:12, size=12, replace=TRUE)
        #     keep = which(model_stats_nearshore[[mt]]$event_and_source_int %in% gauges_to_keep)
        #     rg2 = random_gauge[which(!is.na(match(random_gauge, keep)))] # event_and_source_int is in the random sample of gauges
        #     median_stat_jknf[i] = median(stat[rg2])
        #     abs_median_stat_jknf[i] = median(abs(stat[rg2]))
        # }

        # Non parametric percentile bootstrap error estimate - how does the result vary if we
        # sample 12 scenarios randomly with replacement?
        N = 1000
        median_stat_boot = rep(NA, length=N)
        abs_median_stat_boot = rep(NA, length=N)
        for(i in 1:N){
            gauges_to_keep = sample(1:12, size=12, replace=TRUE)
            # Random gauges that match each 'gauge to keep'. Notice this is 'sampling with replacement'
            # in the event that gauges_to_keep has repeated values.
            rg2 = c(unlist(lapply(gauges_to_keep, 
                f<-function(x) random_gauge[which(model_stats_nearshore[[mt]]$event_and_source_int[random_gauge] == x)])))
            median_stat_boot[i] = median(stat[rg2])
            abs_median_stat_boot[i] = median(abs(stat[rg2]))
        }

        # Print the info
        print('')
        print('')
        print('@@@@@@@@@@@@@@@@@@@@@@@')
        print(paste0('@ ', mt))
        print(paste0('@ ', stat_type))
        print('@@@@@@@@@@@@@@@@@@@@@@@')
        print('Median relative error over gauges, without accounting for unequal numbers of gauges per inversion')
        print(median(stat, na.rm=TRUE))
        print('Median abs(relative error) over gauges, without accounting for unequal numbers of gauges per inversion')
        print(median(abs(stat), na.rm=TRUE))

        print('Median_over_scenarios(median_relative_error_over_gauges_for_each_scenario)')
        print(median_median_per_scenario)
        print('     95% CI, assuming a random sample of scenarios')
        print(median_median_per_scenario_95CI)

        print('Median_over_scenarios(median_abs_relative_error_over_gauges_for_each_scenario)')
        print(median_median_abs_per_scenario)
        print('     95% CI, assuming a random sample of scenarios')
        print(median_median_abs_per_scenario_95CI)

        print('median relative error for "random scenario with random gauge"')
        print(median_stat)
        #print(range(median_stat_jknf))
        print('     bootstrap 95% CI, assuming we have a random sample of inversions')
        print(quantile(median_stat_boot, p=c(0.025, 0.975), type=6))
        #print(quantile(median_stat_boot, p=c(0.05, 0.95), type=6))
        #print(quantile(median_stat_boot, p=c(0.1, 0.9), type=6))
        print('median abs(relative error) for "random scenario with random gauge"')
        print(abs_median_stat)
        #print(range(abs_median_stat_jknf))
        print('     bootstrap 95% CI, assuming we have a random sample of inversions')
        print(quantile(abs_median_stat_boot, p=c(0.025, 0.975), type=6))
        #print(quantile(abs_median_stat_boot, p=c(0.05, 0.95), type=6))
        #print(quantile(abs_median_stat_boot, p=c(0.1, 0.9), type=6))

        print('AIDA 78 Statistics (raw)')
        aida_K = exp(mean(log(stat+1), na.rm=TRUE))
        aida_k = exp(sqrt(mean( (log(stat+1)^2 - log(aida_K)^2), na.rm=TRUE)))
        print(c('    K = ', aida_K))
        print(c('    k = ', aida_k))

        print('AIDA 78 Statistics (random event, random gauge)')
        aida_K_rg = exp(mean(log(stat[random_gauge]+1), na.rm=TRUE))
        aida_k_rg = exp(sqrt(mean( (log(stat[random_gauge]+1)^2 - log(aida_K_rg)^2), na.rm=TRUE)))
        print(c('    K = ', aida_K_rg))
        print(c('    k = ', aida_k_rg))
    }
}
capture.output({alternative_summary_statistics(model_stats_nearshore, '_max')},
    file = 'model_stats_nearshore_summary_max.txt')
capture.output({alternative_summary_statistics(model_stats_nearshore, '_max_36')},
    file = 'model_stats_nearshore_summary_max_36.txt')
capture.output({alternative_summary_statistics(model_stats_nearshore, '_max_8hrs')},
    file = 'model_stats_nearshore_summary_max_8hrs.txt')
capture.output({alternative_summary_statistics(model_stats_nearshore, '_max_lastday')},
    file = 'model_stats_nearshore_summary_max_lastday.txt')


#
# Try a statistical model of the "relative model errors", grouped by source inversion
#
modelling_the_scenario_errors<-function(model_stats_nearshore, stat_type = '_max', 
    mt = 'Frictionless', YLIM=c(-1.5, 1.5)){
          
    modstat = paste0('model', stat_type) 
    datstat = paste0('data', stat_type)

    # Relative model error, log space
    stat = log(model_stats_nearshore[[mt]][[modstat]]/model_stats_nearshore[[mt]][[datstat]])
    # Factor for the scenario
    scenario = as.factor(model_stats_nearshore[[mt]]$event_and_source_int)

    #k = which(!is.na(stat) & model_stats_nearshore[[mt]]$event_name != 'Chile-1960')
    # Estimate the log(model/data) for each scenario. Force a zero intercept
    #model_fit = lm(stat[k] ~ scenario[k] + 0)
    
    model_fit = lm(stat ~ scenario + 0)
    print(summary(model_fit))

    par(mfrow=c(3,2))
    # To make the stats / boxplots align, we need to control the order of the input data. 
    k = order(model_stats_nearshore[[mt]]$event_and_source_int)
    plot(stat[k], col=model_stats_nearshore[[mt]]$event_and_source_int[k], 
         pch=model_stats_nearshore[[mt]]$event_and_source_int[k], 
         main=paste0(mt, '; ', stat_type, '; Median coef=', round(median(coef(model_fit)), 3)), 
         ylim=YLIM)
    grid(col='orange')
    abline(h=0, col='red')

    # For boxplot
    stat_agg = aggregate(stat[k], by=list(model_stats_nearshore[[mt]]$event_name_and_source[k]), f<-function(x) x)
    stat_var = stat_agg$x
    names(stat_var) = stat_agg[,1]
    oldmar = par('mar')
    par(mar=c(9,4, 4,1))
    boxplot(stat_var, names=names(stat_var), varwidth=TRUE, las=2,
            main=paste0(mt,'; ', stat_type, '; Median abs(coef)=', round(median(abs(coef(model_fit))), 3)), 
            ylim=YLIM)
    abline(h=0, col='red')
    grid(col='orange')
    par(mar=oldmar)

    plot(model_fit)
}
# Make plots for various statistics that model the scenario errors. This is
# complementary to the other statistical approaches
for(stat_type in c('_max', '_stgrng', '_max_36', '_stgrng_36')){
    pdf(paste0('Scenario_errors', stat_type, '.pdf'), width=8, height=10)
    modelling_the_scenario_errors(model_stats_nearshore, stat_type=stat_type, mt='Frictionless')
    modelling_the_scenario_errors(model_stats_nearshore, stat_type=stat_type, mt='Manning0.035')
    modelling_the_scenario_errors(model_stats_nearshore, stat_type=stat_type, mt='LinearDelayedFriction')
    modelling_the_scenario_errors(model_stats_nearshore, stat_type=stat_type, mt='LinearFriction')
    modelling_the_scenario_errors(model_stats_nearshore, stat_type=stat_type, mt='LinearReducedFriction')
    dev.off()
}

modelling_scenario_errors_combined<-function(model_stats_nearshore, stat_type = '_max'){
    # Model log(model/obs) as a function of the factors "scenario" and "model type" simultaneously
    # Note this will be different to separate models, because the former will use the same 'sigma' term
    # for each group.

    stat = unlist(lapply(model_stats_nearshore, f<-function(x){
        modstat = paste0('model', stat_type)
        datstat = paste0('data', stat_type)
        stat = log(x[[modstat]]/x[[datstat]])
        return(stat)
        }))
    model_type = as.factor(unlist(lapply(model_stats_nearshore, f<-function(x) x$model_type)))
    scenario = as.factor(unlist(lapply(model_stats_nearshore, f<-function(x) x$event_name_and_source)))

    stat_model_all_in_one = lm(stat ~ (scenario + model_type + 0))
    # The above model has 'Frictionless' integrated into the 'scenario' coefficients -- that's how R's lm works.
    # So positive bias in 'Frictionless' is shown by positive bias in the scenario coefficients.
    # To compute the scenario bias for the other models, add the model constant term to the scenario coefficients.
    stat_models_per_model_type = vector(mode='list', length=length(model_stats_nearshore))
    for(i in 1:length(model_stats_nearshore)){
        k = which(model_type == model_stats_nearshore[[i]]$model_type[1])
        stat_models_per_model_type[[i]] = lm(stat[k] ~ scenario[k] + 0)
    }

    return(environment())
}
sm = modelling_scenario_errors_combined(model_stats_nearshore)

combined_plot_model_vs_data<-function(model_stats_nearshore, model_names, stat_types, 
    stat_names, score_names, XYLIM = c(0.01, 1), modcol=c('blue', 'darkred', 'darkgreen')){
    #
    # Plot model-vs-data for each model type, with summary statistics.
    #

    nmodels = length(model_stats_nearshore)
    nstats = length(stat_types)

    par(mfrow=c(nstats, nmodels))
    par(mar=c(3,4,4.3,1.2))
    par(oma = c(4,4,0,0.3))

    #modcol = c('blue', 'darkred', 'darkgreen')
    for(j in 1:nstats){

        modstat = paste0('model', stat_types[j]) 
        datstat = paste0('data' , stat_types[j])

        for(i in 1:nmodels){

            mod = model_stats_nearshore[[i]][[modstat]]
            dat = model_stats_nearshore[[i]][[datstat]]

            plot(mod, dat, xlab='Model', ylab = 'Data', main="", 
                 asp=1, log='xy', xlim=XYLIM, ylim=XYLIM, ann=FALSE, axes=FALSE, 
                 frame.plot=TRUE, cex.lab=1.8, col=modcol[i], 
                 #pch=letters[model_stats_nearshore[[i]]$event_and_source_int],
                 pch=19,
                 cex=1.4)
            if(model_names[i] == 'Frictionless'){
                title(paste0(model_names[i], ' on global grid \n', stat_names[j]), cex.main=1.8)
            }else{
                title(paste0(model_names[i], '-friction on global grid \n', stat_names[j]), cex.main=1.8)
            }
            add_log_axis_ticks(side=1)
            add_log_axis_ticks(side=2)
            axis(side=1, at = c(0.05, 0.1, 0.5, 1), cex.axis=1.8, las=1)
            axis(side=2, at = c(0.05, 0.1, 0.5, 1), cex.axis=1.8, las=1)
            abline(h=c(0.05, 0.1, 0.5, 1), lty='dotted', col='orange')
            abline(v=c(0.05, 0.1, 0.5, 1), lty='dotted', col='orange')
            abline(0, 1, col='red')
            abline(0, 1.5, col='red', lty='dashed', untf=TRUE)
            abline(0, 1/1.5, col='red', lty='dashed', untf=TRUE)

            mt = model_names[i]
            local_stat_GI = aggregate( (mod - dat)/dat, 
                by=list(model_stats_nearshore[[i]]$event_and_source_int),
                f<-function(x) median(x, na.rm=TRUE))$x
            star = '**'
            heuristic_significance = (sum(local_stat_GI < 0 ) < 3) | (sum(local_stat_GI < 0) > 9) # Heuristic statistical significance -- true median not inside 95% confidence interval for median, which spans 3rd-9th values inclusive
            local_stat_G = median(local_stat_GI)
            if(heuristic_significance){
                text(0.08, 0.82, bquote(G^.(mt) == .(round(local_stat_G, 2))^.(star)), cex=2.5)
            }else{
                text(0.08, 0.82, bquote(G^.(mt) == .(round(local_stat_G, 2))), cex=2.5)
            }

            # As above for abs(relative error) -- note the statistical test
            local_stat_abs_GI = aggregate( abs( (mod - dat)/dat ), 
                by=list(model_stats_nearshore[[i]]$event_and_source_int),
                f<-function(x) median(x, na.rm=TRUE))$x
            local_stat_abs_G = median(local_stat_abs_GI)
            text(0.25, 0.025, bquote(abs(G)^.(mt) == .(round(local_stat_abs_G, 2))), cex=2.5)
        }
    }
    mtext('Modelled tsunami maxima (m)', side=1, cex=2, line=1, outer=TRUE)
    mtext('Observed tsunami maxima (m)', side=2, cex=2, line=1, outer=TRUE)

}
png('combined_model_vs_data.png', width=12, height=13, units='in', res=300)
# Scatterplots of model-vs-predicted over different simulation time intervals
combined_plot_model_vs_data(model_stats_nearshore[1:3], 
                            model_names = c('Frictionless', 'Manning', 'Linear'),
                            stat_types = c('_max_8hrs', '_max_36', '_max', '_max_lastday'), 
                            stat_names = c('(0-8 hours after arrival)', '(0-36 hours after earthquake)', 
                                           '(0-60 hours after earthquake)', '(36-60 hours after earthquake)'),
                            XYLIM=c(0.02, 1),
                            modcol = c('blue', 'darkred', 'darkgreen'))
dev.off()
png('combined_model_vs_data_linear.png', width=12, height=13, units='in', res=300)
combined_plot_model_vs_data(model_stats_nearshore[3:5], 
                            model_names = c('Linear', 'Reduce-Lin', 'Delay-Lin'),
                            stat_types = c('_max_8hrs', '_max_36', '_max', '_max_lastday'), 
                            stat_names = c('(0-8 hours after arrival)', '(0-36 hours after earthquake)', 
                                           '(0-60 hours after earthquake)', '(36-60 hours after earthquake)'),
                            XYLIM=c(0.02, 1),
                            modcol = c('darkgreen', 'wheat4', 'purple')
                            )
dev.off()

boxplot_relative_errors<-function(model_stats_nearshore, stat_type = '_max', 
                                  title_extra = '0-60 hours post-earthquake', YLIM=c(0.4, 2.5),
                                  plot_height=5){
    #
    # Make a data.frame with the:
    #   - model/observed
    #   - source inversions
    #   - model type
    # and make a boxplot
    all_stat_types = paste0(stat_type, collapse="")

    png(paste0('boxplot_relative_errors_Manning_delayedLinear_', all_stat_types, '.png'), 
        width=10, height=plot_height, units='in', res=300)

    #par(mfrow=c(length(stat_type), 1))
    hts = rep(1, length(stat_type))
    hts[length(stat_type)] = 1.5
    layout(matrix(1:length(stat_type), ncol=1), widths=1, heights=hts)

    for(k in 1:length(stat_type)){

        if(k == length(stat_type)){
            final_plot = TRUE
            par(mar=c(10, 4, 3, 2))
        }else{
            final_plot = FALSE
            par(mar=c(0.5, 4, 3, 2))
        }

        data_we_want<-function(ms){
            modstat = paste0('model', stat_type[k])
            datstat = paste0('data', stat_type[k])

            #stat = (ms$model_max - ms$data_max)/ms$data_max
            stat = (ms[[modstat]]/ms[[datstat]])
            sourceinv = ms$event_name_and_source
            model_type = ms$model_type

            sourceinv_split = strsplit(sourceinv, ' ')
            sourceinv_start = unlist(lapply(sourceinv_split, f<-function(x) x[1]))
            sourceinv_end = unlist(lapply(sourceinv_split, f<-function(x){
                firstLetter = substring(x[2], 1, 1)
                lastYear = substring(x[2], nchar(x[2])-1, nchar(x[2]))
                return(paste0(firstLetter, lastYear))
                                }))
            # Now fix some of the labels where I used the wrong year for the paper:
            #   F06 --> F07
            #   P08 --> P07
            sourceinv_end = gsub("F06", "F07", sourceinv_end)
            sourceinv_end = gsub("P08", "P07", sourceinv_end)
            sourceinv_for_label = paste0(sourceinv_start, ' ', sourceinv_end)


            unique_sourceinv_for_label = unique(sourceinv_for_label)
            sourceinv_for_label = factor(sourceinv_for_label, levels = unique_sourceinv_for_label) # Control the boxplot order

            return(data.frame(rel_err = stat, sourceinv = sourceinv_for_label, model_type = model_type))
        }

        all_data = lapply(model_stats_nearshore, data_we_want)
        names(all_data) = names(model_stats_nearshore)

        # Plot with just Manning and Delayed Linear -- these are the best models,
        # and the others have quite a bit of clutter.
        n = 12
        if(final_plot){
            bp_names = levels(all_data$Manning0.035$sourceinv)
        }else{
            bp_names = FALSE
        }
        boxplot(rel_err ~ sourceinv + model_type, data=all_data$Manning0.035, log='y', col='red', border='darkred',
                las=2, boxwex=0.3, at=(1:n), ylim=YLIM, ylab = 'Model / Observed', varwidth=FALSE, 
                names=bp_names, xlab="", cex.axis=1.2, cex.lab=1.2, 
                main=paste0('Modelled-maxima / Observed-maxima,', title_extra[k]), cex.main=1.3)
        boxplot(rel_err ~ sourceinv + model_type, data=all_data$LinearDelayedFriction, col='violet', border='purple',
                boxwex=0.3, add=TRUE, at=(1:n)+0.3, names=FALSE, axes=FALSE, varwidth=FALSE)
        abline(h=c(1/1.5, 1.5), lty='dotted', col='brown')
        abline(h=1, lty='solid', col='brown')
        if(final_plot){
            legend('bottomleft', c('Manning-friction on global grid', 'Delayed-linear-friction on global grid'), 
                   fill=c('red', 'violet'), border = c('darkred', 'purple'), 
                   horiz=TRUE, cex=1.2, bty='n')
        }
        text( (1:n) + 0.15, YLIM[2] - 0.1, 
             as.character(table(all_data$Manning0.035$sourceinv[!is.na(all_data$Manning0.035$rel_err)])),
             cex=1.1)
    }

    dev.off()
    
}
boxplot_relative_errors(model_stats_nearshore, '_max', ' 0-60 hours after earthquake', YLIM = c(1/4, 4))
boxplot_relative_errors(model_stats_nearshore, '_max_8hrs', ' 0-8 hours after tsunami arrival', YLIM=c(1/4, 4))
boxplot_relative_errors(model_stats_nearshore, '_max_12hrs', ' 0-12 hours after tsunami arrival', YLIM=c(1/4, 4))
boxplot_relative_errors(model_stats_nearshore, '_max_24', ' 0-24 hours post-earthquake', YLIM=c(1/4, 4))
boxplot_relative_errors(model_stats_nearshore, '_max_36', ' 0-36 hours post-earthquake', YLIM = c(1/4, 4))
boxplot_relative_errors(model_stats_nearshore, '_max_lastday', ' 36-60 hours post-earthquake', YLIM = c(1/4, 4))

boxplot_relative_errors(model_stats_nearshore, c('_max_24', '_max'), 
                        c(' 0-24 hours post earthquake', ' 0-60 hours post earthquake'), 
                        YLIM = c(1/3, 3), plot_height=8)


compare_with_other_studies<-function(model_stats_nearshore){
    #
    # Compare my 'model-max' -vs- 'data-max' with other studies.
    # Add results from Allen and Greenslade (2016) (Port Kembla) and 
    # Adams and Leveque (2017) (various mostly USA sites).
    #
    source('../other_studies/parse_other_study_data.R', chdir=TRUE)

    filename = paste0('Model-max-vs-data-max.png')
    png(filename, width=10, height=6.66, units='in', res=300)
    par(mfrow=c(2,3))
    par(mar=c(4,4,3,1))
    for(i in 1:length(model_stats_nearshore)){
        plot(model_stats_nearshore[[i]]$data_max_36, 
             model_stats_nearshore[[i]]$model_max_36,
             asp=1, ylim=c(0, 1), xlim=c(0, 1), pch=4, cex=1.5,
             xlab='Observed maxima (m)', ylab= 'Modelled maxima (m)',
             cex.axis=1.5, cex.lab=1.5)
        title(names(model_stats_nearshore)[i], cex.main=1.5)
        grid(col='orange')
        abline(0, 1)
        abline(0, 1.5, col='red', lty='dashed')
        abline(0, 1/1.5, col='red', lty='dashed')
        points(adams_leveque$data_max, adams_leveque$model_max_GEOCLAW, 
               col='purple', pch=19, cex=1.3)
        points(adams_leveque$data_max, adams_leveque$model_max_MOST, 
               col='green', pch=19, cex=1.3)
        points(allen_greenslade$data_max, allen_greenslade$model_max_BOM, 
               col='red', pch=19, cex=1.3)
    }
    par(mar=c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann=FALSE, axes=FALSE, col='white')
    legend('center', 
           c('This study (36 hour duration)', 'Allen & Greenslade (2016), MOST', 
             'Adams & LeVeque (2017), MOST', 'Adams & LeVeque (2017), GEOCLAW',
             'y=x', 'y=x/1.5', 'y=1.5x'),
           pch=c(4, 19, 19, 19, NA, NA, NA), 
           col=c('black', 'red', 'green', 'purple', 'black', 'red', 'red'), 
           pt.cex = c(1.5, 1.3, 1.3, 1.3, NA, NA, NA), 
           lty=c(NA, NA, NA, NA, 'solid', 'dashed', 'dashed'),
           bg=rgb(1,1,1, alpha=0.8), 
           bty='n',
           cex=1.5)
    dev.off()


    print('## Adams and LeVeque, MOST')
    print(summary(-(adams_leveque$data_max - adams_leveque$model_max_MOST)/adams_leveque$data_max))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # -0.67391 -0.23684 -0.01220  0.02679  0.25819  1.40000 
    print(summary(abs(adams_leveque$data_max - adams_leveque$model_max_MOST)/adams_leveque$data_max))
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #0.01220 0.08623 0.33333 0.34275 0.55162 1.40000 

    # Statistics aggregated by inversion, like I do in the paper to avoid events with many gauges having more weight
    MOST_G_m = aggregate(
        (adams_leveque$model_max_MOST - adams_leveque$data_max)/adams_leveque$data_max, 
        by=list(adams_leveque$Inversion), median)
    print(median(MOST_G_m$x))
    # [1] -0.01219512
    MOST_absG_m = aggregate(
        abs(adams_leveque$model_max_MOST - adams_leveque$data_max)/adams_leveque$data_max, 
        by=list(adams_leveque$Inversion), median)
    print(median(MOST_absG_m$x))
    # [1] 0.2138158

    print('## Adams and LeVeque, GEOCLAW')
    print(summary(-(adams_leveque$data_max - adams_leveque$model_max_GEOCLAW)/adams_leveque$data_max))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # -0.71739 -0.35227 -0.02033 -0.03333  0.24662  1.05263 
    print(summary(abs(adams_leveque$data_max - adams_leveque$model_max_GEOCLAW)/adams_leveque$data_max))
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #0.02033 0.17449 0.25000 0.36508 0.54096 1.05263 

    # Statistics aggregated by inversion, like I do in the paper to avoid events with many gauges having more weight
    GEOCLAW_G_m = aggregate(
        (adams_leveque$model_max_GEOCLAW - adams_leveque$data_max)/adams_leveque$data_max, 
        by=list(adams_leveque$Inversion), median)
    print(median(GEOCLAW_G_m$x))
    # [1] -0.1158565
    GEOCLAW_absG_m = aggregate(
        abs(adams_leveque$model_max_GEOCLAW - adams_leveque$data_max)/adams_leveque$data_max, 
        by=list(adams_leveque$Inversion), median)
    print(median(GEOCLAW_absG_m$x))
    # [1] 0.375

    print('## Allen and Greenslade')
    print(summary(-(allen_greenslade$data_max - allen_greenslade$model_max_BOM)/allen_greenslade$data_max))
    #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    #-0.65000 -0.43750 -0.10000 -0.07639  0.25000  0.56667 
    print(summary(abs(allen_greenslade$data_max - allen_greenslade$model_max_BOM)/allen_greenslade$data_max))
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #0.05556 0.25000 0.36111 0.35046 0.46667 0.65000 
    ## This is a single gauge, so the above calculation is equivalent to the |G|^m statistic.


    for(i in 1:length(model_stats_nearshore)){
        print(names(model_stats_nearshore)[i])
        print(summary(    -(model_stats_nearshore[[i]]$data_max_36 - model_stats_nearshore[[i]]$model_max_36)/model_stats_nearshore[[i]]$data_max_36))
        print(summary( abs(model_stats_nearshore[[i]]$data_max_36 - model_stats_nearshore[[i]]$model_max_36)/model_stats_nearshore[[i]]$data_max_36))
    }


    #[1] "Frictionless"
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    #-0.44057 -0.09881  0.05410  0.21183  0.44040  1.41147        2 
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    #0.00095 0.09421 0.21891 0.35086 0.44883 1.41147       2 
    #[1] "Manning0.035"
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    #-0.4467 -0.2366 -0.0951  0.0429  0.2179  1.1126       2 
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    #0.003243 0.120978 0.233915 0.288264 0.345229 1.112636        2 
    #[1] "LinearFriction"
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    #-0.67171 -0.41171 -0.23895 -0.15915  0.04607  0.96944        2 
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    #0.01456 0.18359 0.31211 0.32509 0.44546 0.96944       2 
    #[1] "LinearReducedFriction"
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    #-0.62852 -0.34649 -0.20360 -0.08851  0.10104  1.06169        2 
    #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    #0.009557 0.174881 0.313314 0.313264 0.419874 1.061689        2 
    #[1] "LinearDelayedFriction"
    #     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
    #-0.587568 -0.274591 -0.108557  0.008524  0.215206  1.252448         2 
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    #0.01779 0.14215 0.27130 0.31818 0.43141 1.25245       2 

}
compare_with_other_studies(model_stats_nearshore)


## Save some info
write.csv(model_stats_nearshore, 'model_stats_nearshore_snapshot.csv', row.names=FALSE)
saveRDS(model_stats_nearshore, 'model_stats_nearshore_snapshot.RDS')
#write.csv(model_stats_DART, 'model_stats_DART_snapshot.csv', row.names=FALSE)
#saveRDS(model_stats_DART, 'model_stats_DART_snapshot.RDS')


model_vs_model_maxima_change<-function(){
    #
    # Have a quick look at model-vs-model ratios with a gauge-by-gauge approach.
    #
    png('Tsunami_maxima_model_change_2.png', width=10, height=7, units='in', res=300)
    par(mfrow=c(2,3))
    #par(oma = c(0, 4, 3, 0))
    par(oma = c(0, 4, 0, 0))
    YLIM = c(0.55, 1/0.55)
    colz = 'blue' #rainbow(15)[model_stats_nearshore$Frictionless$event_and_source_int]
    pchz = 19 #model_stats_nearshore$Frictionless$event_and_source_int
    titles = rev(c('60 hour simulation', '36 hour simulation', '24 hour simulation'))
    stats = rev( c('model_max'         , 'model_max_36'      , 'model_max_24'))
    #titles = rev(c('60 hour simulation', '36 hour simulation', '0-8 hours post arrival'))
    #stats = rev( c('model_max'         , 'model_max_36'      , 'model_max_8hrs'))
    model_names = c('Frictionless', 'Delayed-linear')
    model_types = c('Frictionless', 'LinearDelayedFriction')
    for(i in 1:length(model_types)){
        mt = model_types[i]
        for(j in 1:3){
            stat = stats[j]
            #plot(model_stats_nearshore$Manning0.035[[stat]], 
            plot(model_stats_nearshore$Manning0.035$model_time_to_arrive/3600, 
                model_stats_nearshore[[mt]][[stat]]/model_stats_nearshore$Manning0.035[[stat]], 
                #model_stats_nearshore[[mt]][[stat]]/model_stats_nearshore$Manning0.035[[gsub('model', 'data', stat)]], 
                #log='xy', ylim=YLIM, col=colz, pch=pchz, xlab="Modelled maxima (Manning)", 
                log='y', ylim=YLIM, col=colz, pch=pchz, xlab="Modelled arrival time (hours)", 
                ylab=paste0(model_names[i], ' / Manning '),
                main=titles[j], cex.main=2, cex.lab=1.6, cex.axis=1.5, las=1)
            #grid(col='brown')
            abline(h=seq(0.6, 1.8, by=0.2), lty='dotted')
            if(i == 2){
                xs = seq( min(model_stats_nearshore$Manning0.035$model_time_to_arrive/3600) - 1,
                          max(model_stats_nearshore$Manning0.035$model_time_to_arrive/3600) + 1,
                          len=100)
                #ys = exp(-1e-05/2 * xs * 3600) / exp(-1e-05/2 * 12 * 3600) #
                ys = exp(-1e-05/2 * (xs - 12) * 3600)
                points(xs, ys, t='l', col='red', lty='twodash', lwd=1.5)
                #legend('topright', legend=bquote(y == exp(-ft/2) / exp(-f*t[12]/2)), 
                #       lty='twodash', lwd=1.5, col='red', bty='n')
            }
        }
    }
    mtext(side=2, 'Ratio of modelled tsunami maxima ', cex=1.7, outer=TRUE, line=1.5)
    dev.off()
}
model_vs_model_maxima_change()

#
# A few checks to support the discussion.
#

# For the linear-with-delayed-linear-friction model, only 2 of the 68 model-vs-observed series
# do not have BOTH modelled-max < 14hours AND modelled-max > 3hours-post-arrival
late_model_maxima = (model_stats_nearshore$LinearDelayedFriction$model_time_of_max > 14*3600) & 
                  # Modelled maxima occurs well after 12 hours
              (model_stats_nearshore$LinearDelayedFriction$model_time_of_max - 
               model_stats_nearshore$LinearDelayedFriction$model_time_to_arrive > 3*3600)
                  # Modelled maxima is well after arrival time.

summary(late_model_maxima)
#> summary(late_model_maxima)
#   Mode   FALSE    TRUE 
#logical       2      66 
late_observed_maxima = model_stats_nearshore$LinearDelayedFriction$obs_max_time_to_arrive/3600 < 18
model_stats_nearshore$LinearDelayedFriction$site_and_event[which(late_observed_maxima)]
#> model_stats_nearshore$LinearDelayedFriction$site_and_event[which(late_observed_maxima)]
#[1] "Sumatra-2004_Hillarys_BOM_1min_2004" "Sumatra-2004_Hillarys_BOM_1min_2004"
#[3] "Sumatra-2004_Hillarys_BOM_1min_2004"

