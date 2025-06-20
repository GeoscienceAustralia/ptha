# Get an object with the time-series (data and model)
target_RDS_data = readRDS('target_RDS_data.RDS')
#print('FIXME: Not reading target_RDS_data to save time, uncomment this to make everything work!')

# Get event statistics, including "raw model results" and data.
event_stats = readRDS('event_stats.RDS')

# Get alternative event statistics, computed from model results that have been
# smoothed and/or downsampled to mimic the observation process at the tide gauges.
# Typically this involves a 1-min smooth + decimating to the data frequency.
event_stats_downsampled_model = readRDS('event_stats_DOWNSAMPLED.RDS')

# Get the function "downsample_and_smooth_model_at_nearshore_tide_gauge" which
# adjusts the model such that the sampling is "like the data". 
source('downsample_and_smooth_model_at_nearshore_tide_gauge.R')

# When sampling "GOOD_NEARSHORE" gauges, skip those where the model point is > 200m
# from the gauge. This happened in a few instances where I originally had the wrong
# tide gauge coordinate, with distances up to 600 m (although in practice the
# model still works well at those sites).
GOOD_NEARSHORE_DISTANCE_THRESHOLD_M = 200

# Location for figures
FIG_OUTPUT_BASEDIR = 'FIG/'
dir.create(FIG_OUTPUT_BASEDIR, showWarnings=FALSE)

#' Get event observation indices for further processing
#' 
#' For analysis, we won't use all the data in event_stats -- some data is not
#' 'good enough', and sometimes we have multiple measurements at the same site.
#' This function finds indices of the data we will use. 
#'
#' @param event_stats -- data.frame as per above
#' @param event_name -- we get table indices for this event 
#' @param good_nearshore_and_close_and_highres -- if TRUE, we only use 1min tide-gauge data, or
#' 5min for WA, or whatever we have for Chile1960, and only in high-res model
#' regions when gauges are within GOOD_NEARSHORE_DISTANCE_THRESHOLD of the
#' nearest stored model gauge. If FALSE then we allow 15min data if its
#' stage-range exceeds stagerange_threshold_15min -- some 15min data clearly
#' shows a significant tsunami although may not well capture it.
#' @param stagerange_threshold_15min threshold stage-range to keep 15min data
#' @param expand_multi_counts. If TRUE, then repeat indices as many times as
#' the scenario was sampled [noting we sample scenarios at random with
#' replacement, so some scenarios were repeated, but we only simulate the tsunami
#' once]. This expands repeated samples into repeated values, which is often
#' most appropriate for statistical analysis.
#' @param run_type -- the category of model run, either 'random_like_historic'
#' or 'nonrandom_like_historic'
#' @param rigidity_type -- the model rigidity type
#' @param exclude_batch2_doubleups -- Drop scenarios from "batch2" if they appear
#' in the first batch. Recall the first set of "random simulations like historic events"
#' were run in the first two "batches" (30 scenarios each), whereas later I ran 60 scenarios
#' in each batch. If you really only want unique scenarios [e.g. searching for 
#' unique 'good' scenarios'], then this can be set to TRUE, but for most analysis 
#' it should be FALSE
#' @return indices in event_stats matching the criteria 
#'
get_table_indices_to_process<-function(event_stats, event_name, 
    good_nearshore_and_close_and_highres=TRUE, stagerange_threshold_15min=0.4, 
    expand_multi_counts=FALSE, run_type='random_like_historic',
    rigidity_type='constant', exclude_batch2_doubleups=FALSE){

    if(good_nearshore_and_close_and_highres){
        # Use sites where both the data and model are relatively high quality,
        # and the model gauge is close to the data gauge.
        other_use_criteria = event_stats$good_nearshore & 
            event_stats$is_gauge_in_highres_domain &
            (event_stats$distance_to_gauge < GOOD_NEARSHORE_DISTANCE_THRESHOLD_M)

    }else{
        # Accept 'not so good' nearshore data, but avoid the DARTS and the 6
        # minute data (which is redundant - we have nearby 1min data) and
        # anything including Cocos Islands (as our model doesn't have good
        # elevation data), and any 'not big' 15min data.
        other_use_criteria = !(
            grepl('DART', event_stats$sites) | 
            grepl('_6min_', event_stats$sites) |
            grepl('Cocos', event_stats$sites) |
            (grepl('_15min_', event_stats$sites) & 
             (event_stats$data_stgrng < stagerange_threshold_15min)) |
            (!event_stats$is_gauge_in_highres_domain)
            )
    }

    k = which(event_stats$event_name == event_name & 
              event_stats$run_type == run_type & 
              other_use_criteria &
              event_stats$rigidity == rigidity_type)

    if(exclude_batch2_doubleups){
        # Batch1 and Batch 2 cover the same historical events with 30 random scenarios each
        # The other batches cover distinct historical events with 60 randoms scenarios each.

        if(expand_multi_counts){
            stop('Does not make sense to exclude batch-2 double-ups, but expand multi-counts, because batch2 doubleups ARE multi-counts')
        }

        in_batch_1 = which(event_stats$batch_number[k] == 1)
        in_batch_2 = which(event_stats$batch_number[k] == 2)
        if(length(in_batch_1) > 0 & length(in_batch_2) > 0){
            # Define a unique scenario based on the slip type and ID. This works because
            # we already have a unique event and rigidity and run type
            useful_flag = paste0(event_stats$slip_type[k], '_', event_stats$scenario_ID[k])
            b2_matches_b1 = match(useful_flag[in_batch_2], useful_flag[in_batch_1])
            # The unique b2 scenarios will not have a match
            to_remove = in_batch_2[!is.na(b2_matches_b1)]
            k = setdiff(k, k[to_remove])
        }
    }

    # For statistical analysis it is most often easiest to repeat scenario results 
    # based on their counts. This is only done if (!exclude_batch2_doubleups), otherwise
    # we'd hit an error above.
    if(expand_multi_counts & any(event_stats$scenario_count[k] > 1)){ 
        k = rep(k, times=event_stats$scenario_count[k])
    }


    # Quick check
    stopifnot(all(event_stats$event_name[k] == event_name &
                  event_stats$rigidity[k] == rigidity_type &
                  other_use_criteria[k] &
                  event_stats$run_type[k] == run_type))

    return(k)

}

##' Make a boxplot of "model_stat_type"/"data_stat_type" for each random model
##' type
##'
##' @param event_stats -- data.frame with model and data statistics -- either
##' 'event_stats' or 'event_stats_downsampled_model'
##' @param event_name -- the value of event_stats$event_name for the plotted
##' scenarios
##' @param run_type -- the value of event_stats$run_type for the plotted
##' scenarios
##' @param stat_type -- the statistic plotted is
##' event_stats[[paste0('model_', stat_type)]]/event_stats[[paste0('data_', stat_type)]]
##' @param YLIM plot y limits on the log scale
##' @param good_nearshore_and_close_and_highres -- if TRUE then only include nearshore gauges
##' with (event_stats$good_nearshore==TRUE) that are in high-res regions and
##' where the gauge is close to the model., which allows focus on the
##' 'high-quality' data.  If FALSE, then include other nearshore data too [e.g.
##' 15min] which requires care.
##' @return Invisibly return the row-indices in event_stats that are involved in
##' the plot.  This is useful because the function eliminates some
##' 'double-up-observations' (at Eden, when we've got both BOM and DPIE
##' tide-gauge data at the same site,  I preference the latter as by inspection
##' it is generally better).
##'
##'
#boxplot_sites<-function(event_stats, event_name='chile1960', 
#    run_type='random_like_historic', stat_type = '_max', YLIM=c(0.25, 4),
#    good_nearshore_and_close_and_highres=TRUE){
#
#    par(mfrow=c(3,1))
#
#    slip_types = c('HS', 'VAUS', 'FAUS')
#    cols = c('red', 'green', 'blue')
#    borders = c('darkred', 'darkgreen', 'darkblue')
#
#    indices_to_use = get_table_indices_to_process(event_stats, event_name,
#        good_nearshore_and_close_and_highres=good_nearshore_and_close_and_highres, 
#        expand_multi_counts=TRUE, run_type=run_type)
#    
#    kstore = c()
#    for(i in 1:length(slip_types)){
#
#        # Get indices with the desired slip type that are also in 'indices_to_use'
#        slptp_inds = which(event_stats$slip_type[indices_to_use] == slip_types[i])
#        k = indices_to_use[slptp_inds]
#
#        # Store the indices for later.....
#        kstore = c(kstore, k)
#
#        if(length(k) == 0){
#            print(paste0('Skipping : ', slip_types[i]))
#            next
#        }
#
#        # Ratio of the model stat to the corresponding data stat
#        model_stat = paste0('model', stat_type)
#        data_stat = paste0('data', stat_type)
#        model_on_data = event_stats[[model_stat]][k] / 
#                        event_stats[[data_stat ]][k]
#
#        # Make the plot                    
#        sites = event_stats$sites[k]
#        boxplot(model_on_data ~ sites, log='y', border=borders[i], 
#                col=cols[i], ylim=YLIM, pars=list(las=3, cex.axis=0.8))
#        grid()
#        abline(h=1, col='orange')
#        abline(h=c(1/1.5, 1.5), lty='dotted', col='orange')
#        title(paste0(event_name, ' ', stat_type, ' ..... ', 
#                     round(median(model_on_data, na.rm=TRUE), 3)))
#    }
#
#    # Optionally get the corresponding row indices in event_stats.
#    # This can help with subsequent study of data
#    return(invisible(kstore))
#}

#
# Make a time-series plot of the model/data for a given row in the event
# statistics table
#
plot_model_and_data_at_event_row<-function(event_stats, row_index, 
    downsample_model_like_data = FALSE, add_energy=TRUE, 
    xlim_as_hours_post_arrival = NA, YLIM = NA, 
    use_optimal_time_offset = FALSE, time_offset_varname = 'time_domain_hybrid_norm_time_offset',
    title_location = 'regular', horizontal_line_spacing = NA, title_cex = NULL,
    include_run_name = TRUE, add_10h_line=TRUE, title_type='site_and_event',
    ...){

    k = match(event_stats$run_name[row_index], names(target_RDS_data))
    site = event_stats$sites[row_index]

    # Get julian time, stage for the model and data
    data = target_RDS_data[[k]][[site]]$event_data$obs[,c(2,4)]
    if(downsample_model_like_data){
        # At tide gauges, smooth and downsample the model to be like the data
        # By plotting this we can check that everything looks OK.
        model = data.frame(
            juliant = target_RDS_data[[k]][[site]]$model_start_julian + 
                      target_RDS_data[[k]][[site]]$model_time/(3600*24),
            resid = target_RDS_data[[k]][[site]]$model_stage)
        gauge_name = event_stats$sites[row_index]

        if(!grepl('DART', gauge_name)){
            model = downsample_and_smooth_model_at_nearshore_tide_gauge(
                model, gauge_name = gauge_name)
        }

        # Convert back to a matrix
        model = cbind(model$juliant, model$resid)

    }else{
        model = cbind(target_RDS_data[[k]][[site]]$model_start_julian + 
                        target_RDS_data[[k]][[site]]$model_time/(3600*24),
                      target_RDS_data[[k]][[site]]$model_stage)
    }

    if(use_optimal_time_offset){
        model[,1] = model[,1] - event_stats[[time_offset_varname]][row_index]/(3600*24)
    }

    # Optionally set YLIM
    if(all(is.na(YLIM))){
        YLIM = max(abs(c(range(data[,2], na.rm=TRUE), range(model[,2], na.rm=TRUE))))
        YLIM = c(-1,1)*YLIM
    }

    # Optionally set xlim
    if(!all(is.na(xlim_as_hours_post_arrival))){
        t0 = event_stats$model_vs_gauge_start_time[row_index] + xlim_as_hours_post_arrival[1]/24
        t1 = event_stats$model_vs_gauge_start_time[row_index] + xlim_as_hours_post_arrival[2]/24
        XLIM = c(t0, t1)
    }else{
        XLIM = NULL
    }

    plot(data[,1], data[,2], t='o', ylim=YLIM, xlim=XLIM, cex=0.3,
         xlab='Julian day', ylab='Stage residual (m)', ...)

    # Use a different colour for each slip type
    model_col_chooser = list('FAUS' = 'blue', 
        'VAUS' = 'darkgreen', 
        'HS' = 'red', 
        'source_inversion' = 'purple')
    model_col = model_col_chooser[[event_stats$slip_type[row_index]]]
    points(model[,1], model[,2], t='l', col=model_col)
    #grid(col='orange')
    #abline(h=seq(-3,3,by=0.2), col='orange', lty='dotted')
    dy = horizontal_line_spacing
    if(is.na(dy)) dy = max(round(min(abs(YLIM)/2), 1), 0.1) # Grid plotting interval to nearest 10 cm
    abline(h=seq(-10,10)*dy, col='orange', lty='dotted')
    # Hourly vertical lines
    abline(v=event_stats$model_vs_gauge_start_time[row_index] + seq(0, 2, by=1/24), col='orange', lty='dotted')
    if(add_10h_line) abline(v=event_stats$model_vs_gauge_start_time[row_index] + 10/24, col='orange', lty='solid') # Add in 10 hr line for relevance to GOF stat

    if(title_type == 'coordinate'){
        title_start = paste0(round(event_stats$obs_lon[row_index], 3), 'E, ', -round(event_stats$obs_lat[row_index], 3), 'S')
    }else if(title_type == 'site_and_event'){
        title_start = event_stats$site_and_event[row_index]
    }else if(title_type == 'slip_type'){
        title_start = event_stats$slip_type[row_index]
    }else{
        title_start = ""
    }

    if(use_optimal_time_offset){
        time_offset = event_stats[[time_offset_varname]][row_index]
        GOF_stat = gsub('_time_offset', '_stat', time_offset_varname)
        # Title the time offset to show whether the model has been left or right shifted
        sign_char = c('R+', 'L+')[1 + (time_offset > 0)]
        title_word = paste0(
            title_start, ' (', sign_char, abs(round(time_offset)), 's)')
            #title_start, ' (', sign_char, abs(round(time_offset)), 's)', round(event_stats[[GOF_stat]][row_index], 3)) # Include GOF score
    }else{
        title_word = title_start
    }

    if(title_location == 'regular'){
        title(title_word, cex.main=title_cex)
    }else if(title_location == 'topleft'){
        text(XLIM[1]+2/24, YLIM[2]*1.0, title_word, cex=title_cex, adj=c(0.0,0.9)) 
    }else if(title_location == 'topright'){
        text(XLIM[2]-2/24, YLIM[2]*1.0, title_word, cex=title_cex, adj=c(1.0,0.9)) 
    }

    if(add_energy) title(format(event_stats$energy_start[row_index]), line=0)
    if(include_run_name){
        title(sub=substring(basename(dirname(dirname(event_stats$run_name[row_index]))), 1, 120), 
            cex.sub=0.7)
    }
}


#plot_model_and_data_multiple_rows<-function(event_stats, row_indices,
#    downsample_model_like_data = FALSE, add_energy=TRUE, 
#    xlim_as_hours_post_arrival = NA, YLIM = NA, use_optimal_time_offset = FALSE,
#    ...){
#
#    # Index into target_RDS_data
#    kmatch = match(event_stats$run_name[row_indices], names(target_RDS_data))
#
#    # Use a different colour for each slip type
#    model_col_chooser = list(
#        'FAUS' = 'blue', 
#        'VAUS' = 'darkgreen', 
#        'HS' = 'red', 
#        'source_inversion' = 'purple')
#
#    # Pack the data to be plotted in a convenient form
#    plot_data = vector(mode='list', length=length(row_indices))
#    model_range = c(0, 0)
#    data_range = c(0,0)
#    panel_range = c(0,0)
#
#    for(i in 1:length(row_indices)){
#        ri = row_indices[i]
#        tk = kmatch[i] # Index into target_RDS_data
#        site = event_stats$sites[ri]
#        run_type = event_stats$run_type[ri]
#        model_col = model_col_chooser[[run_type]]
#
#        # Get julian time, stage for the model and data
#        data = target_RDS_data[[tk]][[site]]$event_data$obs[,c(2,4)]
#        model = data.frame(
#            juliant=target_RDS_data[[tk]][[site]]$model_start_julian + 
#                    target_RDS_data[[tk]][[site]]$model_time/(3600*24),
#            resid=target_RDS_data[[tk]][[site]]$model_stage)
#
#        # Update the model and data ranges
#        #model_range = range(c(model_range, range(model[,2], na.rm=TRUE)))
#        #data_range  = range(c(data_range,  range(data[,2] , na.rm=TRUE)))
#        data_model_range = range(c(range(model[,2], na.rm=TRUE), range(data[,2], na.rm=TRUE)))
#        panel_offset = panel_range[2] + diff(model_data_range)*0.5
#
#        time_to_subtract_from_model = 0
#        if(use_optimal_time_offset){
#            time_to_subtract_from_model = event_stats$time_domain_hybrid_norm_time_offset[ri]/(3600*24)
#        }
#        
#        # Optionally set YLIM
#        if(all(is.na(YLIM))){
#            YLIM = max(abs(c(range(data[,2], na.rm=TRUE), range(model[,2], na.rm=TRUE))))
#            YLIM = c(-1,1)*YLIM
#        }
#
#        # Optionally set XLIM
#        if(!all(is.na(xlim_as_hours_post_arrival))){
#            t0 = event_stats$model_vs_gauge_start_time[ri] + xlim_as_hours_post_arrival[1]/24
#            t1 = event_stats$model_vs_gauge_start_time[ri] + xlim_as_hours_post_arrival[2]/24
#            XLIM = c(t0, t1)
#        }else{
#            XLIM = NULL
#        }
#
#        # Pack data for later plotting
#        plot_data[[i]] = list(ri=ri, tk=tk, site=site, run_type=run_type, model_col=model_col, data=data, model=model,
#            panel_offset, XLIM=XLIM, YLIM=YLIM, time_to_subtract_from_model) 
#
#    }
#
#}
