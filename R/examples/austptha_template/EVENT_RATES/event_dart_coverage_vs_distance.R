library(rptha)

# INPUT PARAMETER
variable_mu = FALSE # Manually change from TRUE/FALSE to treat each case
# END INPUT

# Get the *.Rdata files that compare models and data
source('config_DART_test_files.R', local=TRUE)

# Get the peak_slip_limit_factor
source('config_peak_slip_limit_factor.R', local=TRUE)

#' Extract statistics from the "corresponding family of model scenarios"
#'
#' This computes some summary statistics from objects made in
#' gauge_summary_statistics.R, that are lists (per dart buoy) containing lists
#' (per model scenario) with summary statistics
#'
#' Returns a data.frame with the goodness of fit statistic and other columns summarising
#' the rupture properties
#' 
#' @param gauge_stats object like 'stochastic_slip_stats' or 'uniform_slip_stats', etc, as created
#' by the script gauge_summary_statistics.R {in e.g. SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/plots/ }
#' for sources where we did DART buoy comparisons
#' @param unit_source_statistics the unit source statistics
#'
family_stats<-function(gauge_stats, unit_source_statistics, peak_slip_limit_factor=Inf){

    # Get time goodness-of-fit statistic for each model scenario
    gf_mat = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) x$model_data_similarity_time))
    # Convert from list of lists to matrix
    for(i in 1:length(gf_mat)) gf_mat[[i]] = unlist(gf_mat[[i]])
    gf_mat = matrix(unlist(gf_mat), ncol=length(gf_mat))
    # Use the 'median GF over all dart buoys' as our GOF value
    gf_median = apply(gf_mat, 1, median)

    # Get the stage range (median over all darts). This is a crude indicator 
    # of the tsunami size
    stage_range_mat = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) diff(x$model_range)))
    # Convert from list of lists to matrix
    for(i in 1:length(stage_range_mat)) stage_range_mat[[i]] = unlist(stage_range_mat[[i]])
    stage_range_mat = matrix(unlist(stage_range_mat), ncol=length(stage_range_mat))
    stage_range_median = apply(stage_range_mat, 1, median)

    #
    dart_names = basename(names(gauge_stats))
    stage_range_data_mat = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) diff(x$data_range)))
    site_time = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) x$data_t[1]))


    # Useful to have the reference Mw in any case
    reference_Mw = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$Mw))

    # Get the peak slip for each model scenario. The earthquake metadata is
    # stored repeatedly for each DART, so just pull it out of the first one
    if('event_slip_string' %in% names(gauge_stats[[1]][[1]]$events_with_Mw)){
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) max(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
        # Mean slip 
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) mean(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))

        # Sum slip
        sum_slip_sum = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) sum(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
    }else{
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
    }

    if(peak_slip_limit_factor < Inf){
        # FIXME: Currently this is not treating scaling relation variability or shear modulus variability
        # This would have to happen if we were to treat normal faults, or source-zones which use other
        # scaling relations
        keep = which(peak_slip_sum < (peak_slip_limit_factor * slip_from_Mw(reference_Mw)))
    }

    return(environment())
}


# Store the statistics in a list (one entry per historical event)
uniform_stat = vector(mode='list', length=length(all_Rdata))
stochastic_stat = vector(mode='list', length=length(all_Rdata))
variable_uniform_stat = vector(mode='list', length=length(all_Rdata))


for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    # Load the R session associated with the gauge_summary_statistics.R script
    # for the ith event
    load(all_Rdata[i], envir=event_env)

    # Main computation here
    stochastic_stat[[i]] = family_stats(event_env$stochastic_slip_stats, 
        event_env$unit_source_statistics, peak_slip_limit_factor)
    uniform_stat[[i]] = family_stats(event_env$uniform_slip_stats, 
        event_env$unit_source_statistics, peak_slip_limit_factor)
    variable_uniform_stat[[i]] = family_stats(event_env$variable_uniform_slip_stats, 
        event_env$unit_source_statistics, peak_slip_limit_factor)
}
event_env = new.env() # Clear the memory

# Store the results
if(variable_mu){
    save.image('event_dart_coverage_vs_distance_session_varyMu.Rdata')
}else{
    save.image('event_dart_coverage_vs_distance_session.Rdata')
}

#
# Nice plotting code here
#

get_results<-function(i, model_type_stat = stochastic_stat){
    site_time = lapply(model_type_stat[[i]]$site_time, f<-function(x) unlist(x))
    site_time = matrix(unlist(site_time), ncol=length(site_time))

    stage_range_data_mat = lapply(model_type_stat[[i]]$stage_range_data_mat, 
        f<-function(x) unlist(x))
    stage_range_data_mat = matrix(unlist(stage_range_data_mat), 
        ncol=length(stage_range_data_mat))

    stage_range_model_mat = model_type_stat[[i]]$stage_range_mat

    # The above results do not discard events that exceed the peak-slip-limit-factor.
    # Let's do that
    keep = model_type_stat[[i]]$keep
    stage_range_data_mat = stage_range_data_mat[keep,,drop=FALSE]
    stage_range_model_mat = stage_range_model_mat[keep,,drop=FALSE]

    # Useful for plotting
    percentile_values = apply(stage_range_model_mat < stage_range_data_mat, 2, mean)
    travel_time = apply(site_time, 2, median)

    return(environment())
}


plotter<-function(model_type_stat, mytitle){
    #
    # Heterogeneous-slip plot
    #

    event_results = lapply(1:length(model_type_stat), get_results, model_type_stat=model_type_stat)

    # Enough colours to make the plot
    #mycol = rainbow(20)
    mycol = colorRampPalette(c('blue', 'orange', 'green', 'purple', 'red'))(18)

    par(oma = c(0, 3,0,0))
    layout(matrix(c(1, 2), nrow=1))
    #par(mar=c(5,10,4,2))
    options(scipen=5) # Suppress scientific notation
    par(mar=c(5,4,3,0))
    plot(c(180, 3600*30)/3600, c(0, 1), log='x', col='white', xlab="", ylab="", cex.axis=1.5,
        main=mytitle, cex.main=2)
    mtext('Time for tsunami to reach DART (hours)', side=1, cex=1.5, line=2.7) 
    mtext('Observed stage-range as percentile of models',
        side=2, cex=1.5, line=3.5)
    mtext('Small obs. <--------------------> Large obs.',
        side=2, cex=1.2, line=2.2)
    for(i in 1:18){
        points(event_results[[i]]$travel_time/3600, event_results[[i]]$percentile_values, 
            pch=i, col=mycol[i], cex=2)
    }
    add_log_axis_ticks(side=1, lwd.ticks=0.5)

    event_names = basename(all_Rdata)
    event_names = gsub('gauge_summary_stats_session_', '', event_names)
    event_names = gsub('.Rdata', '', event_names)
    grid(col='brown')
    plot(c(0,1), c(0,1), col='white', ann=FALSE, frame.plot=FALSE, axes=FALSE)
    legend('left', event_names, pch=1:18, col=mycol[1:18],
        bg='white', cex=1.3, box.col='white', pt.cex=2)

}

#
# Make the plot
# 
if(variable_mu){
    pdf('observation_percentiles_model_types_variable_mu.pdf', width=17, height=7)
}else{
    pdf('observation_percentiles_model_types.pdf', width=17, height=7)
}

plotter(stochastic_stat, mytitle='Heterogeneous-slip')
plotter(variable_uniform_stat, mytitle='Variable-area-uniform-slip')
plotter(uniform_stat, mytitle='Fixed-area-uniform-slip')

dev.off()

