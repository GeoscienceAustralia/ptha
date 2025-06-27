#
# Find the "best" modelled scenarios according to a variety of GOF statistics.
# Plot them and store statistics for later analysis.
#


# Get base datasets and functions for analysis -- these were split out to keep the
# file size managable
source('analysis_data_and_functions.R')

for(scenario_GOF_metric in 1:10){
    #
    # Plot the best_scenarios for each slip type and event
    #

    # Try a range of goodness of fit metrics -- these will select different scenarios, to some extent.
    if(scenario_GOF_metric == 1){
        # median(abs(log(stage_range_ratio)))
        plot_dir = './best_scenarios_good_nearshore_median_abs_log_stage_range_ratio/'
        GOFFUN<-function(x) median(abs(log(x$data_stgrng_during_obs/x$model_stgrng_during_obs)))
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 2){
        # median(abs(log(stage_range_ratio_first_12_hrs)))
        plot_dir = './best_scenarios_good_nearshore_median_abs_log_stage_range_12hrs_ratio/'
        GOFFUN<-function(x) median(abs(log(x$data_stgrng_12hrs_during_obs/x$model_stgrng_12hrs_during_obs)))
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 3){
        # median(abs(stage_range_ratio - 1)))
        plot_dir = './best_scenarios_good_nearshore_median_abs_stage_range_ratio_less_1/'
        GOFFUN<-function(x) median(abs(x$data_stgrng_during_obs/x$model_stgrng_during_obs-1), na.rm=TRUE)
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 4){
        # weighted.mean(log(stage-range-ratio))
        plot_dir = './best_scenarios_good_nearshore_weighted_mean_abs_log_stage_range_ratio/'
        GOFFUN<-function(x){
            weighted.mean(abs(log(x$data_stgrng_during_obs/x$model_stgrng_during_obs)), 
                          w=x$data_stgrng_during_obs)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 5){
        # Statistic similar to Davies (2019) first 10 hours from arrival
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm_median/'
        GOFFUN<-function(x){
            median(x$time_domain_hybrid_norm_stat)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 6){
        # Statistic similar to Davies (2019) first 10 hours from arrival,
        # but using the average instead of the median
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm_average/'
        GOFFUN<-function(x){
            mean(x$time_domain_hybrid_norm_stat)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 7){
        # Statistic similar to Davies (2019) first 10 hours from arrival,
        # but ONLY at the tide gauge with largest observed maxima
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm_biggestwave/'
        GOFFUN<-function(x){
            k = which.max(x$data_max_during_obs)
            x$time_domain_hybrid_norm_stat[k]
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset'
    }else if(scenario_GOF_metric == 8){
        # Alternative to time-varying-hybrid-norm
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm2_median/'
        GOFFUN<-function(x){
            median(x$time_domain_hybrid_norm_stat2)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset2'
    }else if(scenario_GOF_metric == 9){
        # Alternative to time-varying-hybrid-norm
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm3_median/'
        GOFFUN<-function(x){
            median(x$time_domain_hybrid_norm_stat3)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset3'
    }else if(scenario_GOF_metric == 10){
        # Alternative to time-varying-hybrid-norm
        plot_dir = './best_scenarios_good_nearshore_time_varying_hybrid_norm4_median/'
        GOFFUN<-function(x){
            median(x$time_domain_hybrid_norm_stat4)
            }
        time_offset_varname = 'time_domain_hybrid_norm_time_offset4'
    }else{
        stop(paste0('unrecognized scenario_GOF_metric', scenario_GOF_metric))
    }

    #
    # Apply the GOFFUN to PTHA18 random scenarios.
    #

    # First get the relevant indices in event_stats. In this case it isn't
    # useful to include double-counts, since we are just searching for good scenarios,
    # not looking at the statistics
    unique_event_names = unique(event_stats$event_name)
    # Indices for random PTHA scenarios
    scenario_sites_compared = lapply(unique_event_names, function(event_name){
        get_table_indices_to_process(event_stats, event_name, 
            #good_nearshore_and_close_and_highres=FALSE, stagerange_threshold_15min=0.4, 
            good_nearshore_and_close_and_highres=TRUE,
            expand_multi_counts=FALSE, run_type='random_like_historic',
            rigidity_type='constant', exclude_batch2_doubleups=TRUE)
        })
    names(scenario_sites_compared) = unique_event_names

    # Now get the equivalent of "scenario_sites_compared" for the source inversions
    # There are not double-ups
    ssc_source_inv = lapply(unique_event_names, function(event_name){
        get_table_indices_to_process(event_stats, event_name, 
            #good_nearshore_and_close_and_highres=FALSE, stagerange_threshold_15min=0.4, 
            good_nearshore_and_close_and_highres=TRUE,
            expand_multi_counts=FALSE, run_type='source_inversion',
            rigidity_type='source_inversion', exclude_batch2_doubleups=FALSE)
        })
    names(ssc_source_inv) = unique_event_names

    # Append the source inversions to the scenario_sites_compared data
    for(uen in unique_event_names){
        scenario_sites_compared[[uen]] = c(scenario_sites_compared[[uen]], ssc_source_inv[[uen]])
    }


    # For each unique_event, apply the statistic to each scenario, and store
    # the top few scenarios for each slip type
    best_scenarios = vector(mode='list', length=length(scenario_sites_compared))
    names(best_scenarios) = names(scenario_sites_compared)
    scenario_scores_store = best_scenarios # Also store the scores of every scenario
    for(i in 1:length(scenario_sites_compared)){

        # Subset the event statistics for the sites 
        #local_stats = event_stats_downsampled_model[scenario_sites_compared[[i]],]
        local_stats = event_stats[scenario_sites_compared[[i]],]

        # Compute GOFFUN for each modelled scenario (so we are aggregating over
        # tide gauges for a single scenario). I also group by slip_type as it is
        # convenient to have in the output table, although each run_name has a unique
        # slip type so this does not otherwise affect the grouping.
        library(plyr)
        scenario_score = ddply(local_stats, .(run_name, slip_type),
            .fun<-function(x) GOFFUN(x)
            )
        # Replace name 'V1' with name 'x' for consistency with 'aggregate'
        names(scenario_score)[3] = 'x'
        scenario_scores_store[[i]] = scenario_score # Store for later analysis

        # For each slip type, get the best few scenaros
        slip_types = unique(local_stats$slip_type) #c('FAUS', 'HS', 'VAUS', 'source_inversion')
        best_scenarios[[i]] = lapply(slip_types, function(x){
            score_for_slip_type = scenario_score[which(scenario_score$slip_type == x), ]
            score_sort = sort(score_for_slip_type$x, decreasing=FALSE, index.return=TRUE)
            #threshold = score_sort$x[best_N]
            best_N = min(3 , length(score_sort$ix)) # Source inversions may have less than 3 scenarios
            best_scenarios = score_for_slip_type[score_sort$ix[1:best_N],]
            return(best_scenarios)
        })
        names(best_scenarios[[i]]) = slip_types
    }

    dir.create(plot_dir, showWarnings=FALSE)

    # Save the scenario_scores_store
    output_RDS_file = paste0(plot_dir, '/all_scenarios_GOF.RDS')
    saveRDS(scenario_scores_store, output_RDS_file)
    # Save the 'best scenarios'
    output_RDS_file = paste0(plot_dir, '/best_scenarios_GOF.RDS')
    saveRDS(best_scenarios, output_RDS_file)

    # Make plots of the top few scenarios for each event and slip type
    for(ue in unique_event_names){
        for(slip_type in names(best_scenarios[[ue]])){
            top_scenarios = best_scenarios[[ue]][[slip_type]]
            for(i in 1:nrow(top_scenarios)){
                run_i = top_scenarios$run_name[i] # Name of scenario

                # Get event_stats rows for this scenario, for "good_nearshore" sites in highres domains
                esr = which(event_stats$run_name == run_i &
                            event_stats$good_nearshore & 
                            event_stats$is_gauge_in_highres_domain &
                            event_stats$distance_to_gauge < GOOD_NEARSHORE_DISTANCE_THRESHOLD_M)
                event_i = event_stats$event_name[esr[1]] # Name of event
                stopifnot(all(event_stats$event_name[esr] == event_i))

                # Find gauges to be plotted
                gauge_IDS = event_stats$sites[esr]
                stopifnot(all(length(gauge_IDS) == length(unique(gauge_IDS))))

                plot_ylim_scale = max(event_stats$data_max[esr], -event_stats$data_min[esr],
                    event_stats$model_max[esr], -event_stats$model_min[esr], na.rm=TRUE)

                # Make plot
                score = as.character(signif(top_scenarios$x[i], 3))
                plot_type = 'simple' #'original'

                if(plot_type == 'original'){
                    #
                    # One type of plot
                    #
                    output_png = paste0(plot_dir,
                        gsub('.RDS', paste0('_', i, '_score_', score, '_.png'), basename(run_i), fixed=TRUE)
                        )
                    # When there are > 7 gauges, use multiple columns
                    plot_cols = ceiling(length(gauge_IDS)/7)
                    plot_rows = ceiling(length(gauge_IDS)/plot_cols)
                    png(output_png, width=10*plot_cols, height=2 + plot_rows, res=200, units='in')
                    par(mfrow=c(plot_rows, plot_cols))
                    par(mar = c(2,2,2,1))
                    for(j in 1:length(gauge_IDS)){
                        plot_model_and_data_at_event_row(event_stats, row_index = esr[j],
                            downsample_model_like_data = FALSE, add_energy=FALSE,
                            xlim_as_hours_post_arrival = c(0, 18),
                            YLIM=c(-1,1)*plot_ylim_scale, use_optimal_time_offset=TRUE,
                            time_offset_varname = time_offset_varname)
                    }
                    dev.off()
                }else if(plot_type == 'simple'){
                    #
                    # Another plot that better uses the space
                    #

                    # When there are > 12 gauges, use multiple columns
                    plot_cols = ceiling(length(gauge_IDS)/12)
                    plot_rows = ceiling(length(gauge_IDS)/plot_cols)

                    # Repeat the plot with/without the time offset
                    for(use_time_offset in c(FALSE, TRUE)){

                        if(use_time_offset){
                            local_tag = '_simple_plot_time_offset'
                        }else{
                            local_tag = '_simple_plot'
                        }

                        output_png = paste0(plot_dir,
                            gsub('.RDS', paste0('_', i, '_score_', score, local_tag, '.png'), basename(run_i), fixed=TRUE)
                            )
                        png(output_png, width=10*plot_cols, height=0.666*(2 + plot_rows), res=200, units='in')
                        par(oma=c(0, 0, 0.1, 0))
                        par(mfrow=c(plot_rows, plot_cols))
                        par(mar=c(0,0,0.1,0))
                        if(plot_rows > 2){
                            title_cex = 2.4 #plot_cols
                        }else{
                            title_cex = 1.7 #plot_cols
                        }

                        # Order gauges as
                        #   - West coast before east coast
                        #   - North before south
                        gauge_order = order(event_stats$obs_lon[esr] < 130, event_stats$obs_lat[esr], decreasing=TRUE)

                        for(j0 in 1:length(gauge_IDS)){
                            j = gauge_order[j0] # To enforce a gauge order
                            vscale = max(1, floor(plot_ylim_scale*10))/10
                            plot_model_and_data_at_event_row(event_stats, row_index = esr[j],
                                downsample_model_like_data = FALSE, add_energy=FALSE,
                                xlim_as_hours_post_arrival = c(0, 18),
                                YLIM=c(-1,1)*plot_ylim_scale, 
                                use_optimal_time_offset=use_time_offset, time_offset_varname = time_offset_varname,
                                axes=FALSE, frame.plot=FALSE, 
                                horizontal_line_spacing=999999, #vscale/2,
                                title_location='topleft', title_cex=title_cex, include_run_name=FALSE,
                                title_type='coordinate', add_10h_line=FALSE)
                            axis(side=2, at=c(0, vscale), las=1, cex.axis=1.8, line=-3.5)
                            axis(side=1, at=c(0, 1)/24 + event_stats$model_vs_gauge_start_time[esr[j]], 
                                labels=c(0,'1h'), las=1, cex.axis=1.8, line=-2.0) 
                            abline(v=0, col='orange')
                        }
                        dev.off()
                    }
                }
            }
        }
    }

    ## # Plot the goodness-of-fit summary statistics
    ## # BEWARE: THIS IS IGNORING REPEATED COUNTS.
    ## output_png = paste0(plot_dir, '/GOF_summary_', gsub('best_scenarios_good_nearshore_', '', basename(plot_dir)), '.png')
    ## png(output_png, width=12, height=9, units='in', res=300)
    ## x = scenario_scores_store
    ## par(mfrow=c(5,3))
    ## for(ii in 1:length(x)){ 
    ##     if(!any(x[[i]]$slip_type == 'source_inversion')){
    ##         boxcol = c('blue', 'red', 'green')
    ##     }else{
    ##         boxcol = c('blue', 'red', 'black', 'green')
    ##     }
    ##     boxplot(x[[ii]]$x ~ x[[ii]]$slip_type, col=boxcol, 
    ##         main=names(x)[ii], cex.main=2, ylab='GOF score', xlab='Slip type', 
    ##         cex.lab=1.4, cex.axis=1.3, range=0)
    ##     abline(h=seq(0, 10, by=0.1), lty='dotted')
    ## }
    ## dev.off()
}

