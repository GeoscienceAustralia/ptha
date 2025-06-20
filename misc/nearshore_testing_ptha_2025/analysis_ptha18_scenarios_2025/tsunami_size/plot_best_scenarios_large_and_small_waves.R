source('analysis_data_and_functions.R')


for(GOF_statistic in c('time_domain_hybrid_norm_stat', 'time_domain_hybrid_norm_stat2', 'time_domain_hybrid_norm_stat3', 'time_domain_hybrid_norm_stat4')){
    # For each event, find the data with the largest waves, and plot the 'best'
    # model there
    gauge_size_statistic = 'data_max'
    #GOF_statistic = 'time_domain_hybrid_norm_stat'

    # Get unique events and order them by time
    unique_events = unique(event_stats$event_name)
    tmp = order(event_stats$model_start[match(unique_events, event_stats$event_name)])
    unique_events = unique_events[tmp]

    unique_slip_types = unique(event_stats$slip_type)

    for(plot_type in c('largest_gauge', 'second_largest_gauge', 'smallest_gauge')){

        output_png = paste0(FIG_OUTPUT_BASEDIR, '/best_scenario_at_', plot_type, '_', GOF_statistic, '.png')
        png(output_png, width = 16, height=0.6666*(2+28), res=300, units='in')

        #par(mfrow=c(28, 2))
        # For each event, we have a title panel [taking up 2 rows] as well as 4 tide gauge plots
        s2 = rep(NA, length=length(unique_events)*6)
        for(i in 1:14){
            i0 = (i-1)*5 + 1
            s2[((i-1)*6 + 1):(i*6)] = (i-1)*5 + c(1, 2, 3, 1, 4, 5)
        }
        s2 = matrix(s2, ncol=3, byrow=TRUE)
        layout(s2, widths=c(0.13, 0.435, 0.435), heights=rep(1, nrow(s2)))

        par(mar=c(0,0,0,0))
        par(oma=c(2,2,0.5,0.5))

        for(event in unique_events){

            # Get stats for event
            k = which((event_stats$event_name == event) & 
                event_stats$good_nearshore & 
                event_stats$is_gauge_in_highres_domain &
                (event_stats$distance_to_gauge < GOOD_NEARSHORE_DISTANCE_THRESHOLD_M))
            event_stats_k = event_stats[k,]

            # Find the gauge with the largest waves for this event. This assumes that the
            # statistics do not change when the same data is compared to different models, 
            # and we check that holds
            unique_gauges = unique(event_stats_k$site_and_event)
            unique_gauges_i = match(unique_gauges, event_stats_k$site_and_event)
            gauge_size_order = order(event_stats_k[[gauge_size_statistic]][unique_gauges_i], decreasing=TRUE)

            if(plot_type == 'largest_gauge'){
                # Choose the largest gauge
                chosen_gauge = unique_gauges[gauge_size_order[1]]
            }else if(plot_type == 'second_largest_gauge'){
                chosen_gauge = unique_gauges[gauge_size_order[min(2, length(gauge_size_order))]]
            }else if(plot_type == 'smallest_gauge'){
                chosen_gauge = unique_gauges[gauge_size_order[length(gauge_size_order)]]
            }else{
                stop('unknown plot_type')
            }
           
            # Check the statistic is always the same at the largest gauge for this event 
            # It should be, also for other gauges..... unless I forgot to apply a
            # consistent comparison start time.
            tmp = which(event_stats_k$site_and_event == chosen_gauge)
            tmp_stat_range = range(event_stats_k[[gauge_size_statistic]][tmp])
            if(abs(diff(tmp_stat_range)) > 1e-05*abs(max(tmp_stat_range))){
                print(tmp_stat_range)
                stop('Unequal values of gauge_size_statistic for a single gauge')
            }

            # For each slip type, find the row index with the 'best' scenario according
            # to the GOF stat   
            best_scenario_index = sapply(unique_slip_types, function(ust){
                i0 = which(event_stats_k$slip_type == ust & event_stats_k$site_and_event == chosen_gauge)
                return(i0[which.min(event_stats_k[[GOF_statistic]][i0])])
                })

            # 
            # Plotting
            #

            plot_ylim_scale = max(c(
                event_stats_k$data_max[best_scenario_index],
                -event_stats_k$data_min[best_scenario_index],
                event_stats_k$model_max[best_scenario_index],
                -event_stats_k$model_min[best_scenario_index]))
        
            vscale = floor(plot_ylim_scale*20)/20 # Round-ish number for axis

            obs_lonlat = c(event_stats_k$obs_lon[best_scenario_index[1]], event_stats_k$obs_lat[best_scenario_index[1]])

            # Plot the event name and gauge location
            plot(c(0, 1), c(0, 1), axes=FALSE, ann=FALSE, col='white')
            source_year_name = gsub('southamerica', 'chile', event)
            text(0.4, 0.6, paste0(
                substring(source_year_name, 1, nchar(source_year_name)-4), '\n',
                substring(source_year_name, nchar(source_year_name)-3, nchar(source_year_name)), '\n',
                round(obs_lonlat[1], 3), 'E \n', 
                round(-obs_lonlat[2], 3), 'S'), 
                cex=2.5)
            abline(h=-0.02)


            # Plot model vs data for each
            for(i in 1:length(best_scenario_index)){

                ki = k[best_scenario_index[i]]

                # Label placement depending on event [to reduce overlap with important waves]
                if(event %in% c('chile1960', 'sumatra2005', 'chile2010', 'tohoku2011', 'southamerica2014', 'southamerica2015', 'sandwich2021')){
                    title_loc = 'topleft'
                }else{
                    title_loc = 'topright'
                }

                plot_model_and_data_at_event_row(event_stats, row_index = ki ,
                    downsample_model_like_data = FALSE, add_energy=FALSE,
                    xlim_as_hours_post_arrival = c(0, 18),
                    YLIM=c(-1,1)*plot_ylim_scale, 
                    use_optimal_time_offset=TRUE, time_offset_varname = gsub('_stat', '_time_offset', GOF_statistic),
                    axes=FALSE, frame.plot=FALSE, 
                    horizontal_line_spacing=999999, #vscale/2,
                    title_location=title_loc, title_cex=2, include_run_name=FALSE,
                    title_type = 'slip_type', add_10h_line=FALSE)
                abline(v=0+event_stats$model_vs_gauge_start_time[ki], col='orange')
                if(i == 3) axis(side=2, at=c(0, vscale) + 0*vscale/2, labels = c(0, vscale), las=1, cex.axis=1.8, line=0.2)
                if(i == 1 | i == 2){
                    axis(side=1, at=c(0, 1)/24 + event_stats$model_vs_gauge_start_time[ki], labels=c('0','1h'), las=1, cex.axis=1.8, line=-0.3) 

                    if(GOF_statistic %in% c('time_domain_hybrid_norm_stat3', 'time_domain_hybrid_norm_stat4')){
                        # Show when the comparison with data finishes
                        comparison_end = min(event_stats$model_vs_gauge_start_time[ki] + 10/24, event_stats$data_time_of_max[ki] + 3/24)
                        axis(side=1, at=c(0,as.numeric(comparison_end - event_stats$model_vs_gauge_start_time[ki])) + 
                            event_stats$model_vs_gauge_start_time[ki], labels=NA, line=-0.3)
                    }
                }

                # Add lines to separate events
                if(i == 3 | i == 4) abline(h=-plot_ylim_scale)
            } 

        }
        dev.off()
    }


    #
    # Plot the best GOF by tide gauge
    #
    k = which(event_stats$good_nearshore & 
        event_stats$is_gauge_in_highres_domain &
        (event_stats$distance_to_gauge < GOOD_NEARSHORE_DISTANCE_THRESHOLD_M))
    event_stats_k = event_stats[k,]
    # Get the best GOF at each tide-gauge/event for each slip type. Append a
    # bunch of additional variables to the output for convenience when plotting.
    best_gof = aggregate(event_stats_k[[GOF_statistic]], 
        by=list(slip_type = event_stats_k$slip_type, 
                site_and_event=event_stats_k$site_and_event, 
                event=event_stats_k$event_name, 
                obs_lon=event_stats_k$obs_lon, 
                obs_lat=event_stats_k$obs_lat,
                data_max = event_stats_k$data_max), 
        min)

    #
    # For each gauge and event, rank the 4 model types from 1 (best) to 4 (worst) and summarise the result
    #
    rank_per_gauge = rep(NA, nrow(best_gof))
    gauges_per_event = rep(NA, nrow(best_gof))
    for(i in 1:length(rank_per_gauge)){
        ei = which(best_gof$site_and_event == best_gof$site_and_event[i])
        rank_per_gauge[i] = sum(best_gof$x[ei] <= best_gof$x[i]) # Rank from 1 (best) to 4 (worst)
        gauges_per_event[i] = sum(best_gof$slip_type == best_gof$slip_type[i] & 
                                  best_gof$event == best_gof$event[i]) # Count of gauges for each event
    }
    best_gof = cbind(best_gof, data.frame(rank_per_gauge = rank_per_gauge, gauges_per_event=gauges_per_event))

    # Compute a weighted average of the ranks
    mean_rank_per_event = aggregate(best_gof$rank_per_gauge, by=list(slip_type=best_gof$slip_type, event=best_gof$event), mean)
    mean_rank_per_slip_model = aggregate(mean_rank_per_event$x, by=list(slip_type=mean_rank_per_event$slip_type), mean)
    print('##################################')
    print('#')
    print('##################################')
    print(paste0('GOF_statistic: ', GOF_statistic))
    print('#  For each slip model and event, mean_over_tide_gauges(rank) ')
    print(mean_rank_per_event)
    print('#  For each slip model, mean_over_events(mean_over_tide_gauges(rank)): ')
    print(mean_rank_per_slip_model)
    # Same as a weighted average per source type, like
    #ip = which(best_gof$slip_type == 'FAUS')
    #weighted.mean(best_gof$rank_per_gauge[ip], w=1/best_gof$gauges_per_event[ip])

    ##
    ## Scenario based GOF statistic, median over gauges
    ##
    #gof_median_over_gauges = aggregate(event_stats_k[[GOF_statistic]],
    #    by=list(slip_type = event_stats_k$slip_type,
    #            event = event_stats_k$event_name,
    #            run_name = event_stats_k$run_name),
    #    median)


    # Get the order right for plotting
    i0  = order(match(best_gof$event, unique_events), 
        #best_gof$site_and_event, 
        1000*(best_gof$obs_lon < 130) - best_gof$obs_lat, # WA before eastern Australia, North before south
        gsub('source', 'A', best_gof$slip_type))
    best_gof = best_gof[i0,]

    #
    # Plot the minimum GOF
    #

    # Plot colours
    colbase = c('purple', 'red', 'green', 'blue')
    plotcol = colbase[match(best_gof$slip_type, c('source_inversion', 'HS', 'VAUS', 'FAUS'))]
    #plotcol[plotcol == 1] = 6
    pchbase = c(1, 15, 17, 19)
    plotPCH = pchbase[match(best_gof$slip_type, c('source_inversion', 'HS', 'VAUS', 'FAUS'))] #15 + plotcol
    plot_cex = sqrt(best_gof$data_max) + 0.3

    output_png = paste0(FIG_OUTPUT_BASEDIR, '/best_scenario_GOF_at_gauges_', GOF_statistic, '.png')
    png(output_png, width=8, height=12, units='in', res=250)
    par(mar=c(5.1, 3.1, 1.5, 1))
    spacing_4 = c(0.5, 0.15, -0.15, -0.5)/0.8 # To cluster the bars into groups
    maxX = ceiling(max(best_gof$x))
    plot(best_gof$x, 
        1 + length(best_gof$x) - ( seq(1, length(best_gof$x)) + spacing_4 ), 
        xlim=c(0, maxX), xlab="", ylab="",
        ylim=c(1, length(best_gof$x)),
        col=plotcol, cex = plot_cex,
        type='p', pch=plotPCH, 
        frame.plot=FALSE, axes=FALSE)
    for(i in 1:length(best_gof$x)){
        points(c(0, best_gof$x[i]),
            1 + length(best_gof$x) - rep(i + spacing_4[(i-1)%%4 + 1],2),
            type='l', lwd=0.5, col=plotcol[i])

        if((i-4)%%4 == 1){
            points(c(0, 0), 1 + length(best_gof$x) - ( c(i+spacing_4[1], i+3+spacing_4[4])), t='l')
        }
    }
    #points(seq(1, length(best_gof$x)) + c(0.5, 0.15, -0.15, -0.5), 
    #    best_gof$x, 
    #    col=match(best_gof$slip_type, c('source_inversion', 'HS', 'VAUS', 'FAUS')), 
    #    type='p', pch=19, cex=0.5)

    axis(side=1, cex.axis=1.5, line=-1)

    mtext(side=2, "Tide gauges with 4 models per site", cex=2)
    mtext(side=1, 'Goodness of fit (smaller is better)', cex=2, line=2.)
    mtext(side=1, 'Larger points where observed tsunami is larger', cex=1.3, line=3.)
    #text(0, seq(1, length(best_gof$x), by=4) + 1, length(best_gof$x)/4 + 1 - seq(1, length(best_gof$x)/4))

    if(GOF_statistic %in% c('time_domain_hybrid_norm_stat', 'time_domain_hybrid_norm_stat3')){
        legend_x = 0.35*maxX
    }else{
        legend_x = 0.4*maxX
    }
    legend(legend_x, length(best_gof$x)*(61/80), 
    #legend('top', 
        c('Source inversion', 'HS', 'VAUS', 'FAUS'), 
        col=colbase, pch=pchbase, 
        cex=1.5)# , horiz=TRUE)


    #axis(side=3)
    event_change_rle = rle(best_gof$event)
    event_change_index = 1 + length(best_gof$x) - cumsum(event_change_rle$lengths)
    event_labels = gsub('southamerica', 'chile', event_change_rle$values)
    event_midpoints = event_change_index + event_change_rle$lengths/2

    #abline(h=event_change_index + 1/2, lty='dashed', col='darkgrey', lwd=1.5)
    for(i in 1:(length(event_change_index)-1)){
        points(c(0, maxX), event_change_index[i] - c(1,1)/2, type='l', lty='dashed', col='darkgrey', lwd=1.5)
    }
    text(0.89*maxX, 
        event_midpoints-1/2 - 3*(event_labels=='sandwich2021') + 2*(event_labels == 'chile1960'), 
        event_labels, cex=1.5)
    dev.off()
    #library(magick)
    #im0 = image_read(output_png)
    #im1 = image_rotate(im0, 90)
    #image_write(im1, path=gsub('gauges', 'gauges_rotated', output_png))
}
