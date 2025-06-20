source('analysis_data_and_functions.R')

# Function to randomly sample a single tide gauge at random for each event, and
# then compute the fraction of model scenarios below the data (and summary
# statistics thereof). This is a simple way to remove the possible of
# dependence between tide-gauges for a single event.
source('stats_with_random_site_per_event.R')
#source('stats_median_fraction_below.R')
source('stats_under_null_hypothesis.R')

library(rptha)
# Start a cluster to help with some Monte Carlo calculations.
# This is only used to compute the mean and sd of "fraction of model scenarios below the observations,
# sampling 1 gauge per event at random".
library(parallel)
MC_CORES = 16
MC_SAMPLES_PER_CORE = 100 # 4
cl = makeCluster(MC_CORES)
# Make random numbers reproducible yet distinct on each thread. Uses RNGkind("L'Ecuyer-CMRG")
clusterSetRNGStream(cl, iseed=1234)

for(use_downsampled_model in c(FALSE, TRUE)){

    # Get indices of random scenarios for every event, expanding out scenarios
    # sampled more than once (for statistical testing).
    sites_compared = lapply(unique(event_stats$event_name), function(event_name){
        get_table_indices_to_process(event_stats, event_name, 
            #good_nearshore_and_close_and_highres=FALSE, stagerange_threshold_15min=0.4, 
            good_nearshore_and_close_and_highres = TRUE,
            expand_multi_counts=TRUE, run_type='random_like_historic',
            rigidity_type='constant', exclude_batch2_doubleups=FALSE)
        })
    all_inds = unlist(sites_compared) 

    # Control order of events in plot
    event_order =  c('chile1960', 'sumatra2004', 'sumatra2005', 'java2006',  
                     'solomon2007', 'sumatra2007', 'puysegur2009', 
                     'chile2010', 'tohoku2011', 'southamerica2014', 'southamerica2015', 
                     'newhebrides2021', 'kermadec2021', 'sandwich2021')

    # Plot a few different variables, and present the data in a few ways for each
    for(plot_variable in c(
        'max_during_obs', 'stgrng_during_obs', 'rms_during_obs', 
        'max_36_during_obs', 'max_24_during_obs', 'max_12hrs_during_obs', 
        'max_8hrs_during_obs', 'max_24hrs_during_obs', 'max_lastday_during_obs', 
        'max_36_during_obs', 
        'stgrng_36_during_obs', 'stgrng_lastday_during_obs', 'stgrng_12hrs_during_obs', 
        'stgrng_8hrs_during_obs', 'stgrng_24_during_obs', 'stgrng_24hrs_during_obs',
        'stgrng_36_during_obs')){

        if(use_downsampled_model){
            # Data for plot, downsampled stats
            downsampled_stats_flag = "_downsampled"
            es = event_stats_downsampled_model
        }else{
            downsampled_stats_flag = ''
            es = event_stats
        }

        # Get key variables for plots
        slip_type = as.factor(es$slip_type[all_inds])
        site = as.factor(es$sites[all_inds])
        event = as.factor(es$event_name[all_inds])
        event_int = es$event_name_int[all_inds]
        obs   = es[[paste0('data_' , plot_variable)]][ all_inds]
        model = es[[paste0('model_', plot_variable)]][ all_inds]
        var_label = plot_variable
    
        obs_range = range(obs, na.rm=TRUE)
        model_range = range(model, na.rm=TRUE)
        obs_model_range = range(c(obs_range, model_range))
    
        # For the loess regression, use weights to ensure every event counts equally
        # (in practice, the weights are inversely proportional to the number of tide gauges)
        event_weight_inversely_by_num_sites = 1/(sapply(event, function(x) length(unique(site[which(x == event)]))))
        #event_weight_inversely_by_num_sites = 1 + 0/(sapply(event, function(x) length(unique(site[which(x == event)]))))

        # Use better labelling for cases where higher quality figures are needed
        label_extras = list(
            'max_during_obs' = 'maxima ratio (< 60 hours post-earthquake)',
            'max_8hrs_during_obs' = 'maxima ratio (8 hours post-arrival)',
            'max_12hrs_during_obs' = 'maxima ratio (12 hours post-arrival)',
            'max_lastday_during_obs' = 'maxima ratio (36-60 hours post-earthquake)',
            'stgrng_during_obs' = 'stage range ratio (< 60 hours post-earthquake)',
            'stgrng_8hrs_during_obs' = 'stage range ratio (8 hours post-arrival)',
            'stgrng_12hrs_during_obs' = 'stage range ratio (12 hours post-arrival)',
            'stgrng_lastday_during_obs' = 'stage range ratio (36-60 hours post-earthquake)')
        if(!is.null(label_extras[[plot_variable]])){
            var_label = label_extras[[plot_variable]]
        }


        # Make a scatterplot, a boxplot, a combined version, or a
        # "fraction_below" plot that includes some statistics
        for(plot_type in c('scatter', 'boxplot', 'combined', 'fraction_below')){

            if(plot_type == 'scatter'){

                png(paste0(FIG_OUTPUT_BASEDIR, '/Model_obs_ratio_', plot_variable, '_', plot_type, downsampled_stats_flag, '.png'), 
                    width=9, height=2.7, units='in', res=300)
                par(mfrow=c(1,3))

            }else if(plot_type == 'fraction_below'){

                png(paste0(FIG_OUTPUT_BASEDIR, '/Model_obs_ratio_', plot_variable, '_', plot_type, downsampled_stats_flag, '.png'), 
                    width=9, height=2.7, units='in', res=300)
                par(mfrow=c(1,3))

            }else if(plot_type == 'boxplot'){

                png(paste0(FIG_OUTPUT_BASEDIR, '/Model_obs_ratio_', plot_variable, '_', plot_type, downsampled_stats_flag, '.png'), 
                    width=9, height=5.0, units='in', res=300)
                par(oma=c(1.4, 0, 0, 0))
                par(mfrow=c(3,1))

            }else if(plot_type == 'combined'){

                png(paste0(FIG_OUTPUT_BASEDIR, '/Model_obs_ratio_', plot_variable, '_', plot_type, downsampled_stats_flag, '.png'), 
                    width=12, height=6, units='in', res=300)
                layout(matrix(c(1, 2, 3, 4, 5, 6), ncol=2, nrow=3, byrow=TRUE),
                       widths = c(0.25, 0.75), heights=c(1, 1, 1))

            }else{
                stop('unknown plot_type')
            }

            # Useful parameters for all plots
            mar_bp = c(1.2, 5.5, 2, 1)
            mar_sp = c(4, 5.5, 2, 1)
            unique_slip_types = c('FAUS', 'HS', 'VAUS')
            panel_labels = c('A) ', 'B) ', 'C) ')
            boxcols = c('lightblue', 'red', 'green')
            bordercols = rgb(r=c(0, 1, 0), g=c(0, 0, 1), b=c(1, 0, 0), alpha=0.1, maxColorValue=1)

            for(i in 1:length(unique_slip_types)){
                #
                # Get model/obs ratio results for each slip type, split by site and event
                #

                boxcol = boxcols[i]
                bordercol = bordercols[i]
                k = which(slip_type == unique_slip_types[i])
                # To enable good control of boxplot plotting order and labels, it
                # helps to split data into a list
                bp_model_on_data = aggregate(
                    model[k]/obs[k], 
                    # Include both event_int and event in the grouping.
                    # This is redundant, but it's convenient to have both variables in the output
                    by=list(site=site[k], event=event[k], event_int=event_int[k]),
                    function(x) list(x))
                bp_obs = aggregate(
                    obs[k], 
                    by=list(site=site[k], event=event[k], event_int=event_int[k]), 
                    function(x) list(x))
                bp_model = aggregate(
                    model[k], 
                    by=list(site=site[k], event=event[k], event_int=event_int[k]), 
                    function(x) list(x))
                # Ensure the ordering is identical
                stopifnot(all(bp_model_on_data$event == bp_obs$event & bp_model_on_data$site == bp_obs$site & 
                              bp_model_on_data$event == bp_model$event & bp_model_on_data$site == bp_model$site ))

                # Re-order the rows so the event-order is as desired
                desired_order = order(match(bp_model_on_data$event, event_order))
                bp_model_on_data = bp_model_on_data[desired_order,]
                bp_obs = bp_obs[desired_order,]
                bp_model = bp_model[desired_order,]


                if(plot_type %in% c('scatter', 'combined')){
                    #
                    # Scatterplot
                    #
                    par(mar=mar_sp)

                    plot(jitter(obs[k]), model[k], 
                        log='xy', col=boxcol, pch=19, cex=0.1, 
                        xlab="", ylab="", las=1, cex.axis=1.4,
                        xlim=obs_range, ylim=obs_model_range)

                    median_model = unlist(lapply(bp_model$x, median))
                    median_obs = unlist(lapply(bp_obs$x, median))
                    points(median_obs, median_model, col='black', cex=0.5, pch=19)

                    # Consider making your own moving average of these medians, where the average uses a window,
                    # and weights are inversely related to the number of tide gauges for each event in the window

                    # Fit a smoother to the medians
                    #median_wts = sapply(bp_model$event, function(x) 1/sum(x==bp_model$event))
                    #smoother = loess(log(median_model) ~ log(median_obs), weights=median_wts, span=0.3, degree=1)
                    ### FIXME: -- make unique
                    #sm_inds = order(smoother$x) # For plotting as a line
                    #points(exp(smoother$x[sm_inds]), exp(smoother$fitted[sm_inds]), t='l', 
                    #       col='black', lwd=2, lty='solid')

                    # Add a smoother
                    #smooth_in_log_space = TRUE
                    #if(smooth_in_log_space){
                    #    smoother = loess(log(model[k]) ~ log(obs[k]), weights=event_weight_inversely_by_num_sites[k], span=0.3, degree=1)
                    #    # FIXME: -- make unique
                    #    sm_inds = order(smoother$x) # For plotting as a line
                    #    points(exp(smoother$x[sm_inds]), exp(smoother$fitted[sm_inds]), t='l', 
                    #           col='black', lwd=2, lty='solid')
                    #}else{
                    #    smoother = loess((model[k]) ~ (obs[k]), weights=event_weight_inversely_by_num_sites[k], span=0.3, degree=1)
                    #    sm_inds = order(smoother$x) # For plotting as a line
                    #    points((smoother$x[sm_inds]), (smoother$fitted[sm_inds]), t='l', 
                    #           col='black', lwd=2, lty='solid')
                    #}

                    #browser()
                    ##
                    ## Try to make a confidence interval for the smoother. But this fails because the resampling method
                    ## tends to make the range of teh "observations" too small, and then extrapolation is capped by approx.
                    ## Might be able to make a better method?
                    #smoother_CI<-function(ignored, event, model, obs){
                    #    tmp = swap_data_and_random_scenario_from_model_or_data(
                    #        event, model, obs)
                    #    fake_model_vals = unlist(tmp$model_values_list)
                    #    fake_data_vals = unlist(tmp$data_values_list)
                    #    weights_list = vector(mode='list', length=length(tmp$model_values_list))
                    #    for(i in 1:length(weights_list)){
                    #        weights_list[[i]] = 1/rep(sum(event[i] == event), length(tmp$model_values_list[[i]]))
                    #    }
                    #    fake_weights = unlist(weights_list)
                    #    
                    #    
                    #    lfit = loess(log(fake_model_vals) ~ log(fake_data_vals), weights=fake_weights, span=0.3, degree=1)
                    #    sm_inds = order(lfit$x)
                    #    # FIXME: -- make x values unique, and extrapolate better
                    #    xout = seq(0.01, 1, by=0.01)
                    #    yout = approx((lfit$x[sm_inds]), (lfit$fitted[sm_inds]), xout=log(xout), rule=2)$y
                    #    return(list(x=xout, y=exp(yout)))
                    #}
                    #clusterExport(cl, varlist=c('smoother_CI', 'swap_data_and_random_scenario_from_model_or_data'))
                    #all_smooths = parLapply(cl, 1:1000, smoother_CI, event=bp_model$event, model=bp_model$x, obs=bp_obs$x)
                    #all_ys = do.call(rbind, lapply(all_smooths, function(x) x$y))
                    #tlower = apply(all_ys, MARGIN=2, FUN=function(x) quantile(x, 0.05))
                    #tupper = apply(all_ys, MARGIN=2, FUN=function(x) quantile(x, 0.95))
                    #tmed = apply(all_ys, MARGIN=2, FUN=function(x) quantile(x, 0.95))
                    #tx = all_smooths[[1]]$x 
                    #points(tx, tlower, t='l', col='purple')
                    #points(tx, tupper, t='l', col='purple')
                    #points(tx, tmed, t='l', col='purple')


                    mtext('Observed (m)', side=1, line=2.4, cex=1.0)
                    # For 'lastday' variables, need to move the label to avoid axis overlap
                    mtext('Modelled (m)', side=2, line=3.4 + 0.8*grepl('36-60', var_label), cex=1.0)

                    add_log_axis_ticks(side=1, lwd=0.5)
                    add_log_axis_ticks(side=2, lwd=0.5)
                    # Shorten the title
                    title_here = gsub('earthquake', 'EQ', 
                        #gsub('maxima', 'max', 
                            gsub(' hours', 'h', 
                                gsub(" ratio", "", var_label)))#)
                    title_cex = ifelse(grepl('stage', var_label), 1.3, 1.4)
                    title(main=paste0(unique_slip_types[i], ' ', title_here),
                        cex.main=title_cex)
                    abline(0, 1, lwd=1, col='black', lty='dashed')
                    abline(0, 2, untf=TRUE, lty='dotted', col='black')
                    abline(0, 1/2, untf=TRUE, lty='dotted', col='black')

                    # Add a legend for the scatterplot-only case
                    if(plot_type == 'scatter'){
                        legend('bottomright', 
                               c('y=x', 'y=0.5x; y=2x'),
                               lty=c('dashed', 'dotted'),
                               col=c('black', 'black'),
                               lwd=c(1    ,   1 ),
                               bty='n', cex=1.3)
                        #legend('topleft', c('smoother'), lty=c('solid'), 
                        #       col=c('black'), lwd=c(   2 ), bty='n', cex=1.3)
                        legend('topleft', paste0(unique_slip_types[i], ' median'), pch=19, col='black', bty='n', cex=1.3, pt.cex=0.5*1.3)

                    }
                }

                if(plot_type == 'fraction_below'){
                    # 
                    # Plot the fraction of models below the data, against the data statistic
                    #
                    par(mar=mar_sp)
                    frac_below = unlist(lapply(bp_model_on_data$x, function(x) mean(x < 1)))
                    obs_value = unlist(lapply(bp_obs$x, function(x) median(x))) # Usually identical for all values
                    stopifnot(length(obs_value) == nrow(bp_obs))
                    stopifnot(length(frac_below) == nrow(bp_obs))
                    plot(obs_value, frac_below, log='x', col=boxcol, pch=bp_obs$event_int, xlab="", ylab="", ylim=c(0,1))
                    grid(col='orange')
                    mtext('Observed (m)', side=1, line=2.4, cex=1.2)
                    mtext('Fraction below observed', side=2, line=3.4, cex=1.2)
                    title(main=paste0(unique_slip_types[i], ' tsunami ', var_label),
                        cex.main=1.1)
   
                    #
                    # Monte Carlo calculation [number of samples = MC_SAMPLES_PER_CORE * MC_CORES]
                    # - Many times, take one random tide gauge for each event
                    #   and compute the fraction of random scenarios below the
                    #   observation at this gauge.
                    # - Return the mean and sd of these fractions, as well as the number of tide gauges with non-missing data.
                    #
                    clusterExport(cl, 
                        varlist=c('bp_model', 'bp_obs',
                            'MC_CORES', 'MC_SAMPLES_PER_CORE',
                            'mean_and_sd_of_fraction_below_at_one_random_site_per_event', 
                            'simulate_random_data_by_sampling_model_at_tide_gauges'))

                    stats_parallel = clusterApply(cl, 1:MC_CORES, function(ignored){
                        replicate(MC_SAMPLES_PER_CORE, 
                            mean_and_sd_of_fraction_below_at_one_random_site_per_event(
                                bp_model$event,
                                bp_model$x,
                                bp_obs$x)
                        )})
                    stats_parallel = do.call(cbind, stats_parallel) # Combine parallel runs
                    # At this point the number of samples is (MC_SAMPLES_PER_CORE * MC_CORES)
                    stat_values = rowMeans(stats_parallel) # This is the statistic

                    mtext(paste0(
                            'Mean=', signif(stat_values[1], 3), 
                            ', sd=', signif(stat_values[2], 3), 
                            ', n=' , signif(stat_values[3], 3)), 
                        side=3, line=-2, cex=0.9)

                    if(FALSE){
                        # This is expensive. Instead we can use the statistical test from Davies (2019).
                        # While a different test, it provides a similar indication of the model performance, 
                        # and is also much faster to calculate.

                        # To compute a p-value for stat_values above, we'll repeatedly compute a quantity
                        # like "stat_values" with random data meeting the null hypothesis
                        stat_values_fake_data<-function(ignored){
                            # Make random data with the same structure as the input data
                            fake_data_NH = simulate_random_data_by_sampling_model_at_tide_gauges(
                                bp_model$event, bp_model$x, bp_obs$x, jitter_scale=1e-06)
                            # Compute statistic like "stat_values". Since this is happening inside a single core,
                            # we need all samples on this core (MC_SAMPLES_PER_CORE*MC_CORES)
                            stats_parallel = replicate(MC_SAMPLES_PER_CORE*MC_CORES, 
                                mean_and_sd_of_fraction_below_at_one_random_site_per_event(
                                    fake_data_NH$event_id,
                                    fake_data_NH$model_values_list,
                                    fake_data_NH$data_values_list)
                                )
                            stat_values = rowMeans(stats_parallel)
                            }
                        clusterExport(cl, varlist='stat_values_fake_data')
                        # Repeat the computation a large number of times (doesn't have to be the same as the
                        # number of samples used to assess the statistics with each fake dataset, although that
                        # is fine)
                        stat_values_NH = clusterApply(cl, 1:(MC_SAMPLES_PER_CORE*MC_CORES), stat_values_fake_data)
                        stat_values_NH = do.call(cbind, stat_values_NH)

                        mtext(paste0(
                                '(Mean_p=', signif(mean(stat_values_NH[1,]<stat_values[1]), 3), 
                                ', sd_p=', signif(mean(stat_values_NH[2,]<stat_values[2]), 3), ')'),  
                            side=3, line=-4, cex=0.9)
                    }
                    #
                    # End Monte Carlo calculations
                    # 

                    #
                    # Statistic from Davies (2019)
                    #

                    R_me = compute_Rme(bp_model$event, bp_model$x, bp_obs$x, return_stats=TRUE) #, multi_gauge_summary_function=function(x){mean(x, na.rm=TRUE)})


                    mtext(paste0('AD.test p=', signif(R_me$adtest$p.value,3), 
                        ' , sd= ', signif(R_me$sd_R_me_data, 3), 
                        ' ( ', signif(R_me$sd_test_fraction_of_sds_below, 3), ')'),
                        side=3, line=-5, cex=0.9)
                }

                if(plot_type %in% c('boxplot', 'combined')){

                    #
                    # Boxplot
                    # Note a few extra parameters described here: https://r-charts.com/distribution/boxplot-function/
                    par(mar=mar_bp)
                    box_scale = (unlist(lapply(bp_obs$x, median)))**0.5 # Larger waves --> wider boxes
                    boxplot(bp_model_on_data$x, 
                        col=boxcol, ylim=exp(c(-2,2)), log='y', las=1, cex.axis=1.3, pars=list(xaxs='i'), 
                        #boxcol=bordercol, boxlwd=0.5, medcol='black', medlwd=1, 
                        names=rep("", length(bp_model_on_data$x)), 
                        # Ensure box is visible even if the data is NA for some reason
                        width = ifelse(is.na(box_scale), 0.01, box_scale),
                        boxwex = 1.2, staplewex=0.5,
                        range=0 # This makes the whiskers extend to the most extreme data point
                        )
                    axis(side=1, at=seq(1, nrow(bp_obs), 2), 
                        #labels=as.character(seq(1, nrow(bp_obs), 2)), 
                        labels=NA,
                        cex.axis=1.3)
                    add_log_axis_ticks(side=2, lwd.ticks=0.5)
                    # Add a label to the bottom outer margin to cover all plots
                    if(i == length(unique_slip_types)){
                        mtext(side=1, 'Tide gauges', outer=TRUE, line=0.3)
                        mtext(side=1, 'F - models fail to envelope data', adj=0.9, outer=TRUE, line=0.3, cex=0.8, col='darkred')
                        #mtext(side=1, '                              - Failed to envelope data', outer=TRUE, line=0.3, cex=0.8)
                    }

                    # Find entries where the box does not cover 1, and add a mark to highlight those cases
                    not_enveloped = which(unlist(lapply(bp_model_on_data$x, function(y) (all(y > 1)|all(y<1)))))
                    if(length(not_enveloped) > 0){
                        # Since plot uses log scale, par('usr') needs to be transformed by 10^
                        points(not_enveloped, 10^rep(par('usr')[3]-0.17, length(not_enveloped)), pch='F', xpd=TRUE, col='darkred')
                        #points(not_enveloped, rep(4, length(not_enveloped)), pch=4, xpd=TRUE)
                    }

                    ## Plot Annotations
                    plot_title = paste0(panel_labels[i], unique_slip_types[i], 
                        ' model/observed tsunami ', var_label)
                    title(main=plot_title, cex.main=1.7)

                    # Find columns where the event changes
                    NN = length(bp_model_on_data$event)
                    event_sep = which(bp_model_on_data$event[2:NN] != bp_model_on_data$event[1:(NN-1)]) + 1
                    # Find the x mid-point of each event data (the "0" helps with the
                    # first column label, because the x-axis extends a bit to the left.
                    # Same for the "+2" on the last column)
                    mid_event_sep = filter(
                       c(0, event_sep, length(bp_model_on_data$event)+2), 
                       c(1,1)/2)[1:(length(event_sep)+1)]

                    # Add event labels
                    for(ename in 1:length(mid_event_sep)){
                        event_name = as.character(unique(bp_model_on_data$event))[ename]
                        short_event_name = paste0(
                            substring(event_name, 1,1), 
                            substring(event_name, nchar(event_name)-1, nchar(event_name)))
                        # Replace 'southamerica2014' and 'southamerica2015' with 'chile2014' and 'chile2015'
                        if(event_name %in% c('southamerica2014', 'southamerica2015')){
                            short_event_name = gsub('s', 'c', short_event_name)
                        }
                        short_event_name = toupper(short_event_name)
                        # For Sumatra, add a 'u' to disambiguate with Solomon
                        if(grepl('umatra', event_name)) short_event_name = gsub('S', 'Su', short_event_name, fixed=TRUE)
                        if(grepl('andwich', event_name)) short_event_name = gsub('S', 'Sa', short_event_name, fixed=TRUE)

                        #if(event_name == 'solomon2007'){
                        #    # Nudge location slightly
                        #    text(mid_event_sep[ename], exp(1.9), short_event_name,
                        #         adj=c(0.5,0.7), cex=1.)
                        #}else if(event_name == 'sandwich2021'){
                        if(event_name == 'sandwich2021'){
                            # Nudge location slightly
                            text(mid_event_sep[ename]+1, exp(1.9), short_event_name,
                                 adj=c(0.7,0.7), cex=1.)
                        }else{
                            text(mid_event_sep[ename], exp(1.9), short_event_name,
                                 adj=c(0.7,0.7), cex=1.)
                        }
                    }
                    abline(v=event_sep-1/2)
                    abline(h=exp(0), col='orange', lwd=2)
                    #abline(h=c(1/5, 1/2, 1, 2, 5), lty='dotted', col='orange')
                    abline(h=c(1/2, 1, 2), lty='dotted', col='orange')

                    mtext('Model / Obs ', side=2, cex=1.2, line=2.8)

                }

            }

            dev.off()
        }
    }

}
# Clean up parallel
stopCluster(cl)

#library(rptha)
#png(paste0(FIG_OUTPUT_BASEDIR, '/Energy_in_scenarios.png'), width=20, height=7, units='in', res=300)
#print('FIXME BOX COLOURS')
#par(mar=c(14, 4, 3, 2))
#boxplot(event_stats$energy_start ~ event_stats$slip_type + event_stats$event_name, log='y', las=3, col=c('blue', 'red', 'black', 'green'))
#add_log_axis_ticks(side=2)
#dev.off()
