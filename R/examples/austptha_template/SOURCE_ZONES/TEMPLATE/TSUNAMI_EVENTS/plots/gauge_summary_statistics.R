# Produce summary statistics to compare gauges + models
library(rptha)
source('time_domain_hybrid_norm.R')

# Find part of the time-series on which we apply the model data comparison, and
# make a plot. Although originally mainly for plotting, this function actually
# extracts data in a useful way, and may be used just for that purpose with
# make_plot=FALSE
#
# Idea: Plot comparison when the time is:
# A) After the earthquake occurs
# B) The time is not more than half-an-hour before the modelled |stage| exceeds
#    0.5% of its maximum absolute value
# C) The start-time/end-time do not exceed the start/end of the DART sampling
#    with at least 60 seconds frequency
# D) The entire window is <= 12 hours long (so we can focus on the most significant waves)
#   - But summary statistics are computed over a potentially shorter period -- ending
#     at most 3 hours after the data maxima occurs
#
plot_model_gauge_vs_data_gauge<-function(
    model_index, 
    event_data, 
    event_metadata, 
    unit_source_statistics, 
    time_window_hrs=12, 
    title_extra='',
    max_model_time_shift_min = 15, 
    make_plot=TRUE){

    # Start obs time = 30 min before model absolute value exceeds 0.5% of its
    # maxima, but start is never < 0 minutes post-event 
    peak_abs_modelled_stage = max(abs(event_data$model_events[[model_index]][1,,1]))
    tmp = min(which(abs(event_data$model_events[[model_index]][1,,1]) > 
        5.0e-03*peak_abs_modelled_stage))
    obs_time_zero = max(event_data$model_times[tmp] - 0.5 * 3600, 0)

    # The data has a column 'allowed', which is FALSE before a manually
    # identified start time [used to skip over Rayleigh waves], and TRUE after.
    # It goes FALSE again beyond the graphical window we used to extract DART data,
    # which is often designed to finish and avoid some artefact
    min_allowed_data_index = min(which(event_data$gauge_obs$allowed & (event_data$gauge_obs$zoom == 1)))
    max_allowed_data_index = max(which(event_data$gauge_obs$allowed & (event_data$gauge_obs$zoom == 1)))

    # Do not allow the start of the model data comparison to occur AFTER the peak
    # stage of the observation.
    i1i2 = min_allowed_data_index:max_allowed_data_index
    peak_stage_index = min_allowed_data_index -1 + which.max(event_data$gauge_obs$resid[i1i2])
    peak_stage_time = event_data$gauge_obs_times[peak_stage_index]
    # Do not start comparison later than "10min before peak data stage". May be overruled by the
    # 'min_allowed_data_index' later, or if the DART sampling frequency is too low.
    latest_comparison_start_time = max(peak_stage_time - 10*60, 0) 
    obs_time_zero = min(obs_time_zero, latest_comparison_start_time)

    # Adjust the beginning of the model-data comparison time to be within the
    # allowed range of the data
    obs_time_zero = max(obs_time_zero, 
        event_data$gauge_obs_times[min_allowed_data_index])
    obs_time_max = obs_time_zero + 3600 * (time_window_hrs) 
    obs_time_max = min(obs_time_max, 
        event_data$gauge_obs_times[max_allowed_data_index])

    # Identify the time region bounding where the DART sampling rate is <= 1
    # minutes
    dt_threshold = 61
    obs_dt = c(dt_threshold, diff(event_data$gauge_obs_times, lag=2)/2, 
        dt_threshold)
    obs_keep = which(obs_dt < dt_threshold & 
        event_data$gauge_obs_times > obs_time_zero & 
        event_data$gauge_obs_times < obs_time_max)
    if(length(obs_keep) == 0){
        print(c('No obs to keep for model_index ', model_index))
        print(range(obs_dt))
        print(obs_time_zero)
        print(obs_time_max)
    }
    time_range = as.numeric(c(event_data$gauge_obs_times[min(obs_keep)], 
        event_data$gauge_obs_times[max(obs_keep)]))

    # Extract the desired subset of the data
    data_inds = which(event_data$gauge_obs_times > time_range[1] & 
        event_data$gauge_obs_times < time_range[2])
    data_t = as.numeric(event_data$gauge_obs_times[data_inds])
    data_s = event_data$gauge_obs$resid[data_inds]
    # Get the data range [previously used filtering, but since using better
    # DART extraction this is no-longer required]
    data_filtered_stats = gauge_range_filtered(data_t, data_s, 
        filter_freq=1/(0.001*60), interp_dt=15, detailed=TRUE)
    data_range = range(data_s) #data_filtered_stats$range
    data_peak_freq = data_filtered_stats$mintomax_peak_frequency

    # Extract the desired subset of the model
    model_inds = which( (event_data$model_times > 
        (time_range[1] - 60*max_model_time_shift_min)) & 
        (event_data$model_times < (time_range[2] + 60*max_model_time_shift_min)))
    model_s = event_data$model_events[[model_index]][1,model_inds,1]
    model_t = event_data$model_times[model_inds]

    # Get the model range [we previously used this for filtering, but since
    # using better DART extraction this is no-longer required]
    model_filtered_stats = gauge_range_filtered(model_t, model_s, 
        filter_freq=1/(0.001*60), interp_dt=15, detailed=TRUE)
    model_range = range(model_s) #model_filtered_stats$range
    model_peak_freq = model_filtered_stats$mintomax_peak_frequency

    #
    # Model-vs-data similarity statistic  
    # 

    # Allow the model to be time-shifted (delayed) by up to a few minutes to
    # get the 'best' goodness of fit statistic. The maximum allowed delay is
    # at least 2 min, at most 15 min (for far away stations), but [1/50 *
    # (model_t[1])] if the latter is between 2min and 15 min. See Watada et al
    # (2014) for example of time delays for Tohoku and Chile 2010 -- attributed
    # to loading, etc.
    allowed_model_time_delay_min = pmax(-15, pmin(-2, -1/50 * (model_t[1]/60)))
    stopifnot( (allowed_model_time_delay_min >= -15) & 
               (allowed_model_time_delay_min <= -2) )

    # Try a model-data similarity statistic
    post_peak_time_range_comparison_hours = 3
    data_max_time = data_t[which.max(abs(data_s))]

    model_data_similarity_time_detailed = gauge_similarity_time_domain(
        data_t, data_s, 
        model_t, model_s,
        interp_dt = 15, 
        allowed_lag_minutes=c(allowed_model_time_delay_min, 0), 
        time_range = c(time_range[1], 
            min(time_range[2], data_max_time + post_peak_time_range_comparison_hours*3600)),
        detailed=TRUE)

    model_data_similarity_time = model_data_similarity_time_detailed$objective
    model_time_offset = model_data_similarity_time_detailed$minimum

    # Energy in bands
    # Bin boundaries of 2min, 6min, 20min, 60min, 180 min, along the lines of
    # Rabinovich's analyses of tsunami spectra
    bin_divisors = c(1/(2*60), 1/(6*60), 1/(20*60), 1/(60*60), 1/(3*60*60))
    energy_data = gauge_energy_banding(data_t, data_s, interp_dt=15, 
        bin_divisors = bin_divisors)
    energy_model = gauge_energy_banding(model_t, model_s, interp_dt=15, 
        bin_divisors = bin_divisors)

    # Compute a spectral measure of model similarity, using the 6-20, 20-60,
    # and 60-180 min bands:
    # 1 - 2*(energy_data .dot. energy_model)/( |energy_data|^2 + |energy_model|^2)
    model_data_similarity_spec = 1 - 2*sum(energy_data[3:5]*energy_model[3:5])/
        (sum(energy_data[3:5]**2) + sum(energy_model[3:5]**2))


    if(make_plot){
        # Plot
        par(mfrow=c(4,1))
        par(mar=c(3,3,2,1))
        ylim1 = range(c(data_range, model_range))

        # Data
        plot(data_t, data_s, t='l', ylim=ylim1, 
            main=paste0('Data: ', title_extra), col='red')
        points(model_t - model_time_offset, model_s, t='l', col='grey', 
            lty='dashed')
        grid()
        abline(v=seq(min(data_t), max(data_t), by=3600), col='green', 
            lty='dashed')
        abline(v=range(data_t), col='orange', lty='dashed')
        abline(h=data_range, col='brown', lty='dashed')

        # Model
        plot(model_t - model_time_offset, model_s, t='l', ylim=ylim1,
            main=paste0('Model: ', model_index, ' , Stat: ',
                round(model_data_similarity_time, 2),
                ' ', round(model_data_similarity_spec, 2)),
            col='blue')
        points(data_t, data_s, t='l', col='grey', lty='dashed')
        grid()
        abline(v=seq(min(data_t), max(data_t), by=3600), col='green', 
            lty='dashed')
        abline(h=model_range, col='brown', lty='dashed')

        # Slip raster
        slip_rast_list = make_slip_raster(model_index,
            event_metadata$events_with_Mw, 
            unit_source_statistics)
        plot(slip_rast_list$slip_rast, asp=1, xlim=slip_rast_list$xlim)

        # Banded energy-in-model vs energy-in-data
        plot(1:6, energy_model/energy_data, log='y', type='h', lend=2, lwd=4, 
            main = 'Energy model / Energy data', axes=FALSE, 
            ylim=c(1e-02, 1e+02))
        grid(col='orange')
        abline(h=1, col='red')
        axis(side=2)
        axis(side=1, at=1:6, 
            labels=c('<2 min', '2-6 min', '6-20', '20 - 60', '60 - 180', '>180'))
    }

    return(list(
        bin_divisors = bin_divisors, 
        energy_data=energy_data, 
        energy_model=energy_model, 
        data_t=data_t, data_s=data_s, 
        model_t=model_t, model_s = model_s,
        model_range=model_range, data_range=data_range,
        model_peak_freq = model_peak_freq, data_peak_freq = data_peak_freq,
        model_data_similarity_time=model_data_similarity_time,
        model_data_similarity_spec = model_data_similarity_spec,
        model_time_offset = model_time_offset,
        events_with_Mw = event_metadata$events_with_Mw[model_index,])
    )
}


#' Make a raster with the stochastic or uniform or variable_uniform slip values
#'
#' @param event_index index of stochastic slip values for table
#' @param event_metadata data.frame with the event metadata, for a stochastic
#' or uniform slip model
#' @param unit_source_statistics data.frame with the summary statistics for the
#' unit sources
#' @return a raster with the slip values
#'
make_slip_raster<-function(event_index, event_metadata, unit_source_statistics){

    # Construct a matrix to hold slip anywhere on the source-zone
    slip_mat = matrix(0, ncol=max(unit_source_statistics$alongstrike_number), 
        nrow=max(unit_source_statistics$downdip_number))

    # Read the slip indices & values
    event_meta = event_metadata[event_index,]
    slip_indices = scan(text=gsub("-", " ", event_meta$event_index_string), 
        quiet=TRUE)
    if('event_slip_string' %in% names(event_meta)){
        slip_vals = scan(text=gsub("_", " ", event_meta$event_slip_string), 
            quiet=TRUE)
    }else{
        slip_vals = rep(event_meta$slip, length=length(slip_indices))
    }

    # Set the slip values
    slip_mat[slip_indices] = slip_vals

    # Make a raster with the slip values. Let dx/dy be based on the dx/dy at
    # the peak slip location. This is inexact but will give the idea.
    max_slip_loc = which.max(slip_vals)
    dx = mean(unit_source_statistics$length[slip_indices])
    dy = mean(unit_source_statistics$width[slip_indices])

    slip_rast = raster(slip_mat, xmn=0, xmx = dx*ncol(slip_mat), ymx=0, 
        ymn=-dy*nrow(slip_mat))

    # Limit output to regions where slip is non-zero
    alongstrike_keep = range(
        unit_source_statistics$alongstrike_number[slip_indices])
    downdip_keep = range(
        unit_source_statistics$downdip_number[slip_indices])

    # Return peak slip location for variable slip [or centriod for constant
    # slip]
    if(diff(range(slip_vals)) > 1.0e-02){
        
        peak_slip_loc = c(unit_source_statistics$lon_c[slip_indices[max_slip_loc]],
            unit_source_statistics$lat_c[slip_indices[max_slip_loc]])

        peak_slip_alongstrike = unit_source_statistics$alongstrike_number[slip_indices[max_slip_loc]]
        peak_slip_downdip = unit_source_statistics$downdip_number[slip_indices[max_slip_loc]]

    }else{

        peak_slip_loc = c(mean_angle(unit_source_statistics$lon_c[slip_indices]),
            mean_angle(unit_source_statistics$lat_c[slip_indices]))
        # Make sure longitude has same domain as original points
        peak_slip_loc = adjust_longitude_by_360_deg(peak_slip_loc, 
            cbind(unit_source_statistics$lon_c[1], unit_source_statistics$lat_c[1]))

        peak_slip_alongstrike = mean(unit_source_statistics$alongstrike_number[slip_indices])
        peak_slip_downdip = mean(unit_source_statistics$downdip_number[slip_indices])
    }

    return(list(
        slip_rast=slip_rast, 
        xlim=dx*c(alongstrike_keep[1] - 1, alongstrike_keep[2]), 
        alongstrike_keep = alongstrike_keep,
        downdip_keep = downdip_keep,
        dx = dx, 
        dy = dy,
        approx_L = dx*(diff(alongstrike_keep)+1),
        approx_W = dy*(diff(downdip_keep)+1),
        peak_slip_alongstrike=peak_slip_alongstrike,
        peak_slip_downdip=peak_slip_downdip)
        )
}


###############################################################
#
# Main code here
#
###############################################################

#
# Identify data directories [made with 'check_dart.R' in '..']
#
source_name = basename(dirname(dirname(getwd())))
event_basedirs_uniform = dirname(
    Sys.glob(paste0('../*uniform_uniform/event_metadata.RDS')))
event_basedirs_stochastic = gsub('uniform', 'stochastic', 
    event_basedirs_uniform)
event_basedirs_variable_uniform = gsub('uniform', 'variable_uniform', 
    event_basedirs_uniform)

unit_source_statistics = read_table_from_netcdf(
    Sys.glob('../unit_source_statistics*.nc'))

# Loop over all events
for(dir_ind in 1:length(event_basedirs_uniform)){

    event_basedir_uniform = event_basedirs_uniform[dir_ind]
    event_basedir_stochastic = event_basedirs_stochastic[dir_ind]
    event_basedir_variable_uniform = event_basedirs_variable_uniform[dir_ind]

    # Apply the plot for each gauge, uniform slip
    event_metadata = readRDS(paste0(event_basedir_uniform, '/event_metadata.RDS'))
    event_data_files = Sys.glob(paste0(event_basedir_uniform, '/gauge*.RDS'))
    
    uniform_slip_stats = list()

    # Loop over all gauges
    for(event_data_file in event_data_files){
    
        uniform_slip_stats[[event_data_file]] = list()
        event_data = readRDS(event_data_file)
        event_gauge_name = basename(event_data_file)
        for(i in 1:length(event_data$model_events)){
            energies = plot_model_gauge_vs_data_gauge(i, event_data, 
                event_metadata, unit_source_statistics, time_window_hrs=12, 
                event_gauge_name, make_plot=FALSE)
            uniform_slip_stats[[event_data_file]][[i]] = energies
            # Also record the earthquake peak slip -- here just equal to the
            # uniform slip
            peak_slip = event_metadata$events_with_Mw$slip[i]
            uniform_slip_stats[[event_data_file]][[i]]$peak_slip = peak_slip
            
        }
    
    }
    
    # Apply the plot for each gauge, stochastic slip
    event_metadata = readRDS(
        paste0(event_basedir_stochastic, '/event_metadata.RDS'))
    event_data_files = Sys.glob(
        paste0(event_basedir_stochastic, '/gauge*.RDS'))
    
    stochastic_slip_stats = list()
    # Loop over all gauges
    for(event_data_file in event_data_files){
    
        stochastic_slip_stats[[event_data_file]] = list()
        event_data = readRDS(event_data_file)
        event_gauge_name = basename(event_data_file)

        # Loop over all events
        for(i in 1:length(event_data$model_events)){

            energies = plot_model_gauge_vs_data_gauge(i, event_data, 
                event_metadata, unit_source_statistics, time_window_hrs=12, 
                event_gauge_name, make_plot=FALSE)
            stochastic_slip_stats[[event_data_file]][[i]] = energies

            # Also record the earthquake peak slip 
            peak_slip = max(scan(text=gsub("_", " ", 
                event_metadata$events_with_Mw$event_slip_string[i]), 
                    quiet=TRUE))
            stochastic_slip_stats[[event_data_file]][[i]]$peak_slip = peak_slip
            stochastic_slip_stats[[event_data_file]][[i]]$peak_slip_downdip_ind = 
                event_metadata$events_with_Mw$peak_slip_downdip_ind[i]
        }
    }

    # Apply the plot for each gauge, variable_uniform slip
    event_metadata = readRDS(
        paste0(event_basedir_variable_uniform, '/event_metadata.RDS'))
    event_data_files = Sys.glob(
        paste0(event_basedir_variable_uniform, '/gauge*.RDS'))
    
    variable_uniform_slip_stats = list()
    # Loop over all gauges
    for(event_data_file in event_data_files){
    
        variable_uniform_slip_stats[[event_data_file]] = list()
        event_data = readRDS(event_data_file)
        event_gauge_name = basename(event_data_file)

        # Loop over all events
        for(i in 1:length(event_data$model_events)){

            energies = plot_model_gauge_vs_data_gauge(i, event_data, 
                event_metadata, unit_source_statistics, time_window_hrs=12, 
                event_gauge_name, make_plot=FALSE)
            variable_uniform_slip_stats[[event_data_file]][[i]] = energies

            # Also record the earthquake peak slip 
            peak_slip = max(scan(text=gsub("_", " ", 
                event_metadata$events_with_Mw$event_slip_string[i]), 
                    quiet=TRUE))
            variable_uniform_slip_stats[[event_data_file]][[i]]$peak_slip = 
                peak_slip
            variable_uniform_slip_stats[[event_data_file]][[i]]$peak_slip_downdip_ind = 
                event_metadata$events_with_Mw$peak_slip_downdip_ind[i]
        }
    }

    # Save the uniform, stochastic, and variable_uniform outputs to RDS
    output_name_base = gsub('_uniform', '', basename(event_basedir_uniform))
    saveRDS(list(
            uniform_slip_stats = uniform_slip_stats, 
            stochastic_slip_stats = stochastic_slip_stats,
            variable_uniform_slip_stats=variable_uniform_slip_stats),
        file=paste0(output_name_base, '_energies.RDS'))

    # Compute some summary statistics -- ratios of banded energy in model +
    # data
    pdf(paste0(output_name_base, '_energy_band_ratios.pdf'), width=20, 
        height=15)
    par(mfrow=c(6,5))
    par(mar=c(2,2,3,2))
    
    nstats = length(uniform_slip_stats[[1]][[1]]$energy_data)
    
    for(gauge_ind in 1:length(uniform_slip_stats)){

        # Investigate -- uniform slip
        energiesU_data = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_data)), 
            byrow=TRUE, ncol=nstats)
        # The energy of the data varies (slightly) because of changes in the
        # model start time.
        energiesU_model = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_model)), 
            byrow=TRUE, ncol=nstats)

        # Key information on the model performance, stochastic
        energiesS_data = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_data)), 
            byrow=TRUE, ncol=nstats)
        # The energy of the data varies (slightly) because of changes in the
        # model start time.
        energiesS_model = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_model)), 
            byrow=TRUE, ncol=nstats)

        # Key information on the model performance, variable_uniform
        energiesVU_data = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_data)), 
            byrow=TRUE, ncol=nstats)
        # The energy of the data varies (slightly) because of changes in the
        # model start time.
        energiesVU_model = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$energy_model)), 
            byrow=TRUE, ncol=nstats)


        # Key information on the model performance
        #summary(energiesS_model/energiesS_data)
        bp_widths = colMeans(energiesS_model) + 1.0e-010

        boxplot(energiesS_model/energiesS_data, width=bp_widths, 
            col=rgb(1.0, 0.0, 0.0, alpha=0.3), log='y', border='red', 
            names=c('<2 min', '2-6', '6-20', '20-60', '60-180', '>180'), las=2)

        boxplot(energiesVU_model/energiesVU_data, width=bp_widths, 
            col=rgb(0.0, 1.0, 0.0, alpha=0.3), log='y', border='green', 
            add=TRUE, axes=FALSE)

        boxplot(energiesU_model/energiesU_data, width=bp_widths,
            col=rgb(0, 0, 1, alpha=0.3), add=TRUE, border='blue',
            axes=FALSE)

        title(main = basename(names(uniform_slip_stats)[gauge_ind]))
        abline(h=1, col='red')
    }

    #
    # Boxplot of stage range, vs data, for all model types
    #

    par(mfrow=c(6,5))
    for(gauge_ind in 1:length(uniform_slip_stats)){
    
        # Uniform
        stageU_data = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)),
            byrow=TRUE, ncol=2)
        stageU_model = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
       
        # Stochastic 
        stageS_data = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        stageS_model = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)

        # Variable uniform
        stageVU_data = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        stageVU_model = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)

        dU = stageU_model[,2] - stageU_model[,1]
        dS = stageS_model[,2] - stageS_model[,1]
        dVU = stageVU_model[,2] - stageVU_model[,1]

        dS_obs = stageS_data[,2] - stageS_data[,1]
        dU_obs = stageU_data[,2] - stageU_data[,1]
        dVU_obs = stageVU_data[,2] - stageVU_data[,1]

        ldU = length(dU)
        ldS = length(dS)
        ldVU = length(dVU)

        boxplot(c(dU, dS, dVU) ~ c(rep('Uniform', ldU), rep('Stochastic', ldS), 
                rep('VarUnif', ldVU)), 
            horizontal = TRUE, density=20, 
            col=c(rgb(1,0,0, alpha=0.3), rgb(0,0,1, alpha=0.3), 
                rgb(0,1,0,alpha=0.3)),
            border=c('red', 'blue', 'green'),
            log='x')
        title(main=paste0('Stage range: ', 
            basename(names(uniform_slip_stats)[gauge_ind]),
            ' \n Random (red), Uniform (blue), VarUnif (green), Obs (black)'))
        abline(v=dS_obs, col='black', lty='solid', lwd=2)
        abline(v=dU_obs, col='black', lty='solid', lwd=2)
        abline(v=dVU_obs, col='black', lty='solid', lwd=2)
    }
   
    # 
    # Boxplot of model-data similarity statistic, for all model types
    #

    par(mfrow=c(6,5))
    similar_s_time = list()
    similar_u_time = list()
    similar_vu_time = list()
    similar_s_spec = list()
    similar_u_spec = list()
    similar_vu_spec = list()
    for(gauge_ind in 1:length(uniform_slip_stats)){

        similar_S = unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_time))
        similar_VU = unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_time))
        similar_U = unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_time))

        similar_s_time[[gauge_ind]] = similar_S
        similar_u_time[[gauge_ind]] = similar_U
        similar_vu_time[[gauge_ind]] = similar_VU

        similar_s_spec[[gauge_ind]] = unlist(
            lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_spec))
        similar_u_spec[[gauge_ind]] = unlist(
            lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_spec))
        similar_vu_spec[[gauge_ind]] = unlist(
            lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_data_similarity_spec))

        ldU = length(similar_U)
        ldS = length(similar_S)
        ldVU = length(similar_VU)

        boxplot(c(similar_U, similar_S, similar_VU) ~ c(rep('Uniform', ldU), 
                rep('Stochastic', ldS), rep('VarUnif', ldVU)), 
            horizontal = TRUE, density=20, 
            col=c(rgb(1,0,0, alpha=0.3), rgb(0,0,1, alpha=0.3), 
                rgb(0,1,0, alpha=0.3)),
            border=c('red', 'blue', 'green'),
            log='x')
        title(main=paste0('Similarity statistic: ', 
            basename(names(uniform_slip_stats)[gauge_ind]),
            ' \n Random (red), Uniform (blue), VariUnif (green), lower is better'))
    }

    #
    # Plot stage range vs peak slip for all models
    #

    par(mfrow=c(6,5))
    for(gauge_ind in 1:length(uniform_slip_stats)){

        # Uniform
        stageU_data = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        stageU_model = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        peak_slip_U = unlist(lapply(uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$peak_slip))
       
        # Stochastic 
        stageS_data = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)),
            byrow=TRUE, ncol=2)
        stageS_model = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        peak_slip_S = unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
            f<-function(x) x$peak_slip))
        peak_slip_dd = unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
            f<-function(x) x$peak_slip_downdip_ind))

        # Variable uniform
        stageVU_data = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)),
            byrow=TRUE, ncol=2)
        stageVU_model = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        peak_slip_VU = unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$peak_slip))


        dU = stageU_model[,2] - stageU_model[,1]
        dS = stageS_model[,2] - stageS_model[,1]
        dVU = stageVU_model[,2] - stageVU_model[,1]
        dS_obs = stageS_data[,2] - stageS_data[,1]
        dU_obs = stageU_data[,2] - stageU_data[,1]
        dVU_obs = stageVU_data[,2] - stageVU_data[,1]

        xlim1=range(c(peak_slip_S, peak_slip_U, peak_slip_VU))
        ylim1=range(c(dS, dU, dVU))
        plot(peak_slip_S, dS, log='xy', xlab='Peak slip', ylab='Stage range', 
            col='red', pch=as.character(peak_slip_dd), xlim=xlim1, ylim=ylim1)
        points(jitter(peak_slip_U), dU, col='blue', pch=19)
        points(jitter(peak_slip_VU), dVU, col='green', pch=19)
        title(main=paste0('Peak slip vs stage range: ', 
            basename(names(uniform_slip_stats)[gauge_ind]),
            ' \n Random (red), Uniform (blue), VariUniform (green), Obs (black)'))

    }

    #
    # Plot stage range vs peak spectral period, for models and data
    #

    par(mfrow=c(6,5))
    for(gauge_ind in 1:length(uniform_slip_stats)){

        # Uniform
        stageU_data = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        freqU_data = unlist(lapply(uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$data_peak_freq))
        stageU_model = matrix(
            unlist(lapply(uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        freqU_model = unlist(lapply(uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$model_peak_freq))
       
        # Stochastic 
        stageS_data = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        freqS_data = unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
            f<-function(x) x$data_peak_freq))
        stageS_model = matrix(
            unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        freqS_model = unlist(lapply(stochastic_slip_stats[[gauge_ind]], 
            f<-function(x) x$model_peak_freq))

        # Variable uniform
        stageVU_data = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$data_range)), 
            byrow=TRUE, ncol=2)
        freqVU_data = unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$data_peak_freq))
        stageVU_model = matrix(
            unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) x$model_range)),
            byrow=TRUE, ncol=2)
        freqVU_model = unlist(lapply(variable_uniform_slip_stats[[gauge_ind]], 
            f<-function(x) x$model_peak_freq))


        dU = stageU_model[,2] - stageU_model[,1]
        dU_freq = freqU_model
        dS = stageS_model[,2] - stageS_model[,1]
        dS_freq = freqS_model
        dVU = stageVU_model[,2] - stageVU_model[,1]
        dVU_freq = freqVU_model

        dS_obs = stageS_data[,2] - stageS_data[,1]
        dS_obs_freq = freqS_data
        dVU_obs = stageVU_data[,2] - stageVU_data[,1]
        dVU_obs_freq = freqVU_data
        dU_obs = stageU_data[,2] - stageU_data[,1]
        dU_obs_freq = freqU_data

        xlim1=range(c(dU_freq, dS_freq, dVU_freq, dS_obs_freq, dU_obs_freq, 
            dVU_obs_freq))
        ylim1=range(c(dS, dU, dU_obs, dS_obs, dVU_obs))
        plot(dS_freq, dS, log='xy', xlab='Peak frequency', ylab='Stage range', 
            col='red', pch=as.character(peak_slip_dd), xlim=xlim1, ylim=ylim1)
        points(jitter(dU_freq), dU, col='blue', pch=19)
        points(jitter(dVU_freq), dVU, col='green', pch=19)
        points(dU_obs_freq, dU_obs, col='black', pch=15)
        points(dS_obs_freq, dS_obs, col='black', pch=15)
        points(dVU_obs_freq, dVU_obs, col='black', pch=15)
        title(main=paste0('Spectral peak frequency vs stage range: ', 
            basename(names(uniform_slip_stats)[gauge_ind]),
            ' \n Random (red), Uniform (blue), VariUnif (green), Obs (black)'))

    }


    #
    # Image plot of model-data similarity statistic
    #
    
    gauge_time_stat_S = matrix(NA, ncol=length(stochastic_slip_stats), 
        nrow=length(stochastic_slip_stats[[1]]))
    for(i in 1:length(stochastic_slip_stats) ){
        for(j in 1:length(stochastic_slip_stats[[1]]) ){
            gauge_time_stat_S[j,i] = 
                stochastic_slip_stats[[i]][[j]]$model_data_similarity_time
        }
    }

    gauge_time_stat_VU = matrix(NA, ncol=length(variable_uniform_slip_stats), 
        nrow=length(variable_uniform_slip_stats[[1]]))
    for(i in 1:length(variable_uniform_slip_stats) ){
        for(j in 1:length(variable_uniform_slip_stats[[1]]) ){
            gauge_time_stat_VU[j,i] = 
                variable_uniform_slip_stats[[i]][[j]]$model_data_similarity_time
        }
    }

    gauge_time_stat_U = matrix(NA, ncol=length(uniform_slip_stats), 
        nrow=length(uniform_slip_stats[[1]]))
    for(i in 1:length(uniform_slip_stats) ){
        for(j in 1:length(uniform_slip_stats[[1]]) ){
            gauge_time_stat_U[j,i] = 
                uniform_slip_stats[[i]][[j]]$model_data_similarity_time
        }
    }

    # Get the median over gauges, for each event
    med_gts_S = apply(gauge_time_stat_S, 1, median)
    med_gts_VU = apply(gauge_time_stat_VU, 1, median)
    med_gts_U = apply(gauge_time_stat_U, 1, median)

    zlims = range(c(range(gauge_time_stat_S), range(gauge_time_stat_U), 
        range(gauge_time_stat_VU)))

    par(mfrow=c(3,1))
    par(mar=c(4,4,2,2))
    image(gauge_time_stat_S, col=terrain.colors(255), xlab='Events', ylab='Gauges', 
        main='Stochastic slip time-similarity statistic (points give event median)',
        cex.main=1.5, zlim=zlims)
    points(seq(0,1, len=length(med_gts_S)), med_gts_S)
    grid(col='purple'); abline(h=0.5, col='red')
    image(gauge_time_stat_VU, col=terrain.colors(255), xlab='Events', ylab='Gauges',
        main='Variable uniform slip time-similarity statistic (points give event median)',
        cex.main=1.5, zlim=zlims)
    points(seq(0,1, len=length(med_gts_VU)), med_gts_VU)
    grid(col='purple'); abline(h=0.5, col='red')
    image(gauge_time_stat_U, col=terrain.colors(255), xlab='Events', ylab='Gauges',
        main='Uniform slip time-similarity statistic (points give event median)',
        cex.main=1.5, zlim=zlims)
    points(seq(0,1, len=length(med_gts_U)), med_gts_U)
    grid(col='purple'); abline(h=0.5, col='red')

    #
    # Plots comparing the distributions of the 'median goodness of fit' stat for
    # each Gauge.
    #

    par(mfrow=c(2,2))
    qqplot(med_gts_S, med_gts_VU, 
        main='Stochastic vs variable-uniform median GOF',
        xlab='Stochastic', ylab='Variable uniform', asp=1)
    abline(0,1,col='red')
    grid()
    abline(h=0.5, v=0.5, col='orange')

    qqplot(med_gts_S, med_gts_U, 
        main='Stochastic vs uniform median GOF',
        xlab='Stochastic', ylab='Uniform', asp=1)
    abline(0,1,col='red')
    grid()
    abline(h=0.5, v=0.5, col='orange')

    qqplot(med_gts_S, med_gts_U, 
        main='Variable-uniform vs uniform median GOF',
        xlab='Variable-uniform', ylab='Uniform',
        asp=1)
    abline(0,1,col='red')
    grid()
    abline(h=0.5, v=0.5, col='orange')

    #boxplot(c(med_gts_S, med_gts_VU, med_gts_U) ~ 
    #    c(rep('Stoc.', length(med_gts_S)), rep('VarUnif.', length(med_gts_VU)),
    #        rep('Unif', length(med_gts_U))), col=c('red', 'green', 'blue'),
    #    horizontal=TRUE)

    plot(ecdf(med_gts_S), col='red', main='', xlim=zlims)
    plot(ecdf(med_gts_VU), col='green', add=TRUE, main='', xlim=zlims)
    plot(ecdf(med_gts_U), col='blue', add=TRUE, main='', xlim=zlims)
    abline(v=0.5, col='orange')
    title('Median GOF statistic, ECDF ')
    legend('bottomright', c('Stoch.', 'VarUnif.', 'Unif'), 
        col=c('red', 'green', 'blue'), pch=19)

    #
    # Source slip properties vs goodness of fit
    #
    stochastic_slip_source = lapply(stochastic_slip_stats[[1]], 
        f<-function(x){
            make_slip_raster(1, x$events_with_Mw,
                unit_source_statistics = unit_source_statistics)}
        )
    variable_uniform_slip_source = lapply(variable_uniform_slip_stats[[1]], 
        f<-function(x){
            make_slip_raster(1, x$events_with_Mw,
                unit_source_statistics = unit_source_statistics)}
        )
    uniform_slip_source = lapply(uniform_slip_stats[[1]], 
        f<-function(x){
            make_slip_raster(1, x$events_with_Mw,
                unit_source_statistics = unit_source_statistics)}
        )
 
    #
    # Plot goodness of fit vs location
    #
    panel_plot<-function(slip_source, GOF, titleword='', ...){

        alongstrike_locs = unlist(lapply(slip_source, 
            f<-function(x) x$peak_slip_alongstrike))
        downdip_locs = unlist(lapply(slip_source, 
            f<-function(x) x$peak_slip_downdip))

        approx_area = unlist(lapply(slip_source, 
            f<-function(x) x$approx_L*x$approx_W))

        plot(alongstrike_locs, GOF, 
            cex= 1.0*sqrt(approx_area/min(approx_area)),
            col = downdip_locs,
            xlab='Alongstrike location', ylab='GOF', ...)
        title(titleword, cex.main=1.5)
    }

    par(mfrow=c(2,2))
    zlims = range(c(med_gts_S, med_gts_VU, med_gts_U))
    panel_plot(stochastic_slip_source, med_gts_S, 
        'Stochastic slip GOF spatial (point-size ~ area)', ylim=zlims)
    panel_plot(variable_uniform_slip_source, med_gts_VU, 
        'Variable_uniform slip GOF spatial (point-size ~ area)', ylim=zlims)
    panel_plot(uniform_slip_source, med_gts_U, 
        'Uniform slip GOF spatial (point-size ~ area)', ylim=zlims)

    dev.off()


    #
    #
    #
    save.image(paste0('gauge_summary_stats_session_', output_name_base, '.Rdata'))
}
