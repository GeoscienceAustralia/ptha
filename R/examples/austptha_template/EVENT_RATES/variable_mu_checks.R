library(rptha)

# Read key user variables
config = new.env()
source('config.R', local=config)

#
ptha = new.env()
#source('../../CODE/ptha/ptha_access/get_PTHA_results.R', local=ptha, chdir=TRUE)
source('../../../CODE/ptha/ptha/ptha_access/get_PTHA_results.R', local=ptha, chdir=TRUE)

# Read key sourcezone parameters
sourcezone_parameter_file = config$sourcezone_parameter_file 
sourcezone_parameters = read.csv(sourcezone_parameter_file, stringsAsFactors=FALSE)


# Plot tsunami stage-vs-exceedance_rate curve at these points
nearby_points = rbind(c(151.42, -34.05), # Sydney @ 100m depth
                      c(186.97, -33.01), # Dart straight offshore of Kermadec
                      c(153.5592, -24.6837) ) # Roberts point

#
# Load the event_properties sessions, which provide some bias correction functions
#
ep_fixed_mu = new.env()
load('event_properties_and_GOF_session_end.Rdata', envir=ep_fixed_mu)
ep_variable_mu = new.env()
load('event_properties_and_GOF_session_varyMu_end.Rdata', envir=ep_variable_mu)

#
# Utility function to search a vector of file paths ('files') for
# a single file with a particular base-name ('file')
#
match_file<-function(files, file){
    ind = grep(file, basename(files), fixed=TRUE)
    if(length(ind) != 1) stop(paste0('Could not match unique file like ', file))
    return(files[ind])
}

#
# Plot function
#
plot_sourcezone_rate_curve_with_fixed_and_variable_mu<-function(sourcezone, slip_type='stochastic'){

    sourcezone_events = list()

    # 
    # Get the events (without tsunami) 
    #
    if(slip_type == 'stochastic'){
        # Stochastic slip
        stochastic_eq_file = match_file(
            gsub('_tsunami', '', config$all_source_stochastic_slip_tsunami, fixed=TRUE),
            paste0('all_stochastic_slip_earthquake_events_', sourcezone, '.nc'))
        sourcezone_events$events = read_table_from_netcdf(stochastic_eq_file)

        tsunami_file = match_file(
            config$all_source_stochastic_slip_tsunami,
            paste0('all_stochastic_slip_earthquake_events_tsunami_', sourcezone, '.nc'))

        # Functions for bias adjustment of rates
        bias_adjust_peak_slip_fixed_mu = ep_fixed_mu$stochastic_quantile_adjust$peak_slip$bias_adjustment_factor
        bias_adjust_peak_slip_variable_mu = ep_variable_mu$stochastic_quantile_adjust$peak_slip$bias_adjustment_factor

        bias_adjust_area_including_zeros_fixed_mu = ep_fixed_mu$stochastic_quantile_adjust$area_including_zeros$bias_adjustment_factor
        bias_adjust_area_including_zeros_variable_mu = ep_variable_mu$stochastic_quantile_adjust$area_including_zeros$bias_adjustment_factor

    }else if(slip_type == 'variable_uniform'){
        # Variable uniform slip
        variable_uniform_eq_file = match_file(
            gsub('_tsunami', '', config$all_source_variable_uniform_slip_tsunami, fixed=TRUE),
            paste0('all_variable_uniform_slip_earthquake_events_', sourcezone, '.nc'))
        sourcezone_events$events = read_table_from_netcdf(variable_uniform_eq_file)

        tsunami_file = match_file(
            config$all_source_variable_uniform_slip_tsunami,
            paste0('all_variable_uniform_slip_earthquake_events_tsunami_', sourcezone, '.nc'))

        # Functions for bias adjustment of rates
        bias_adjust_peak_slip_fixed_mu = ep_fixed_mu$variable_uniform_quantile_adjust$peak_slip$bias_adjustment_factor
        bias_adjust_peak_slip_variable_mu = ep_variable_mu$variable_uniform_quantile_adjust$peak_slip$bias_adjustment_factor

        bias_adjust_area_including_zeros_fixed_mu = ep_fixed_mu$variable_uniform_quantile_adjust$area_including_zeros$bias_adjustment_factor
        bias_adjust_area_including_zeros_variable_mu = ep_variable_mu$variable_uniform_quantile_adjust$area_including_zeros$bias_adjustment_factor

    }else{

        stop('Unknown slip type')

    }


    # Get the unit source statistics
    uss_file = match_file(config$unit_source_statistics_netcdf_files,
        paste0('unit_source_statistics_', sourcezone, '.nc'))
    sourcezone_events$unit_source_statistics = read_table_from_netcdf(uss_file)

    #
    # Extract event properties from the file
    #
    # Event inds/slip
    event_inds = lapply(sourcezone_events$events$event_index_string, 
        f<-function(x) as.numeric(strsplit(x, '-')[[1]]) )
    event_slip = lapply(sourcezone_events$events$event_slip_string, 
        f<-function(x) as.numeric(strsplit(x, '_')[[1]]) )
    # Rates
    event_rate = sourcezone_events$events$rate_annual
    event_rate_lower = sourcezone_events$events$rate_annual_lower_ci
    event_rate_upper = sourcezone_events$events$rate_annual_upper_ci
    depth = sourcezone_events$unit_source_statistics$depth
    # Shear modulus
    mu_constant = depth * 0 + 3e+10
    mu_variable = shear_modulus_depth(depth)
    # Unit source area
    area = sourcezone_events$unit_source_statistics$length * 
        sourcezone_events$unit_source_statistics$width
    #
    # Pre-allocate memory and populate in loop
    #
    moment_fixed_mu = rep(0, length(event_slip)) 
    moment_variable_mu = rep(0, length(event_slip))
    nonzero_rupture_area = rep(0, length(event_slip))
    full_rupture_area = rep(0, length(event_slip))
    peak_slip = rep(0, length(event_slip))
    # Loop over all events
    for(i in 1:length(moment_fixed_mu)){

        inds = event_inds[[i]]
        slip = event_slip[[i]]
        moment_fixed_mu[i] = sum(slip * area[inds] * mu_constant[inds] * 1e+06)
        moment_variable_mu[i] = sum(slip * area[inds] * mu_variable[inds] * 1e+06)
        nonzero_rupture_area[i] = sum(area[inds])

        ddi = sourcezone_events$unit_source_statistics$downdip_number[inds]
        asi = sourcezone_events$unit_source_statistics$alongstrike_number[inds]
        full_inds = which(
            (sourcezone_events$unit_source_statistics$downdip_number >= min(ddi) ) &
            (sourcezone_events$unit_source_statistics$downdip_number <= max(ddi) ) &
            (sourcezone_events$unit_source_statistics$alongstrike_number >= min(asi) ) &
            (sourcezone_events$unit_source_statistics$alongstrike_number <= max(asi) ) )

        full_rupture_area[i] = sum(area[full_inds])
        peak_slip[i] = max(slip)
    }

    # Magnitude
    # If mu is variable, we can consider it a relabeling of magnitude.
    Mw_mu_fixed = round(M0_2_Mw(moment_fixed_mu), 1) # Rounding is 'almost' not needed -- just accounting for numerical imprecision in file storage 
    Mw_mu_vary = M0_2_Mw(moment_variable_mu)


    #
    # Compute various bias-corrections, by applying the bias correction functions to
    # sets of events defined by their corresponding_uniform_row
    #
    rate_adjustment_peak_slip_fixed_mu = rep(0, length(event_slip))
    rate_adjustment_peak_slip_variable_mu = rep(0, length(event_slip))
    rate_adjustment_area_including_zeros_fixed_mu = rep(0, length(event_slip))
    rate_adjustment_area_including_zeros_variable_mu = rep(0, length(event_slip))
    corresponding_uniform_row = sourcezone_events$events$uniform_event_row
    unique_uniform_row = unique(corresponding_uniform_row)
    for(ur in unique_uniform_row){

        # Find events which 'correspond' to this uniform event row.
        k = which(corresponding_uniform_row == ur)
        ps = peak_slip[k]
        nza = full_rupture_area[k]
        # The events are random, so we can break ties by 'first' and it is still random
        slip_empirical_quantile = rank(ps, ties='first')/(length(ps)+1)
        nonzero_area_empirical_quantile = rank(nza, ties='first')/(length(nza) + 1)

        #
        # Bias correction factors
        # 

        # Variable mu, peak slip
        tmp = bias_adjust_peak_slip_variable_mu(slip_empirical_quantile)
        rate_adjustment_peak_slip_variable_mu[k] = tmp/sum(tmp)*length(tmp)
        # Fixed mu, peak slip
        tmp = bias_adjust_peak_slip_fixed_mu(slip_empirical_quantile)
        rate_adjustment_peak_slip_fixed_mu[k] = tmp/sum(tmp)*length(tmp)

        # Variable mu, full area
        tmp = bias_adjust_area_including_zeros_variable_mu(nonzero_area_empirical_quantile)
        rate_adjustment_area_including_zeros_variable_mu[k] = tmp/sum(tmp)*length(tmp)
        # Fixed mu, full area
        tmp = bias_adjust_area_including_zeros_fixed_mu(nonzero_area_empirical_quantile)
        rate_adjustment_area_including_zeros_fixed_mu[k] = tmp/sum(tmp)*length(tmp)

    }


    # Here we compute the exceedance rate
    # Note our fixed-mu Mws range 7.2, 7.3, ...
    # These rates can be interpreted as representing a 'bin' of events
    Mw_bin_size = 0.1
    mws = seq(7.15, 9.95, by=Mw_bin_size) # Lower bin boundary
    rate_mu_fixed = sapply(mws, f<-function(x) sum(event_rate * (Mw_mu_fixed >= x)))
    rate_mu_fixed_lower = sapply(mws, f<-function(x) sum(event_rate_lower * (Mw_mu_fixed >= x)))
    rate_mu_fixed_upper = sapply(mws, f<-function(x) sum(event_rate_upper * (Mw_mu_fixed >= x)))
    rate_mu_vary = sapply(mws, f<-function(x) sum(event_rate * (Mw_mu_vary >= x)))

    # Mw-vs-rate as a single page plot
    par(mfrow=c(1,1))
    plot(mws, pmax(rate_mu_fixed, 1.0e-100), t='o', log='y', ylim=c(1.0e-06, 1), 
        xlab="Mw", ylab='Exceedance rate')
    points(mws, pmax(rate_mu_vary, 1.0e-100), t='l', col='red')
    points(mws, rate_mu_fixed_lower, t='l', col='orange', lty='dashed')
    points(mws, rate_mu_fixed_upper, t='l', col='orange', lty='dashed')
    grid()
    title(paste0('Rates with fixed and variable mu, ', slip_type, ', ', sourcezone))

    abline(h=c(1/500, 1/1000, 1/2500, 1/10000), col='brown', lty='dotted')

    #
    # Plot dart stages, one site per panel
    #
    par(mfrow=c(3,1))
    oldmar = par('mar')
    par(mar=c(2,3,3,2))
    stage_rate_list = vector(mode='list', length=nrow(nearby_points))
    for(i in 1:nrow(nearby_points)){

        # Get the peak stage
        rates_dart = ptha$get_stage_exceedance_rate_curve_at_hazard_point(
            target_point = nearby_points[i,])
        max_stage_dart = ptha$get_peak_stage_at_point_for_each_event(
            target_index = rates_dart$target_index, all_source_names=sourcezone)

        # Compute stage-vs-rate with various rate adjustment factors
        stages = rates_dart$stage
        event_stages = max_stage_dart[[1]]$max_stage
        rates_raw = sapply(stages, f<-function(x) sum(event_rate * (event_stages > x)))
        rates_1 = sapply(stages, f<-function(x) sum(event_rate * rate_adjustment_peak_slip_fixed_mu * (event_stages > x)))
        rates_2 = sapply(stages, f<-function(x) sum(event_rate * rate_adjustment_peak_slip_variable_mu * (event_stages > x)))
        rates_3 = sapply(stages, f<-function(x) sum(event_rate * rate_adjustment_area_including_zeros_fixed_mu * (event_stages > x)))
        rates_4 = sapply(stages, f<-function(x) sum(event_rate * rate_adjustment_area_including_zeros_variable_mu * (event_stages > x)))

        plot(stages, rates_raw, log='xy', ylim=c(1.0e-05, 1), xlim=c(0.01, 10), t='l')
        points(stages, rates_1, t='l', col='red')
        points(stages, rates_2, t='l', col='green')
        points(stages, rates_3, t='l', col='blue')
        points(stages, rates_4, t='l', col='orange')
        grid()
        abline(h=c(1/500, 1/1000, 1/2500, 1/10000), col='pink', lty='dashed')
        legend('topright', 
            c('Raw', 'slip_fixed_mu', 'slip_v_mu', 'area_fixed_mu', 'area_v_mu'),
            col=c('black', 'red', 'green', 'blue', 'orange'),
            pch=rep(NA, 5), lty=rep(1,5))

        title(paste0(slip_type, ' ', round(nearby_points[i,], 1)))

        stage_rate_list[[i]] = data.frame(stages=stages, 
            rates_raw = rates_raw, rates_slip_const_mu = rates_1, 
            rates_slip_var_mu = rates_2, rates_area_const_mu = rates_3, 
            rates_area_var_mu = rates_4)
    }
    par(mar=oldmar)
    # Return the key data as output
    output1 = data.frame(mws=mws, rate_mu_fixed=rate_mu_fixed, rate_mu_vary=rate_mu_vary)

    output = list('mw_rates' = output1, 'stage_rates' = stage_rate_list)
    return(output)
}


#
# Main run code
#
all_source_zones = sort(unique(sourcezone_parameters$sourcename))
rate_curves_store = vector(mode='list', length=length(all_source_zones))
names(rate_curves_store) = all_source_zones
slip_type = 'stochastic'

pdf(paste0('Mw_rate_curves_with_and_without_mu_variation_', slip_type, '.pdf'), width=10, height=12)
for(i in 1:length(all_source_zones)){

    # Only work on subduction and thrust type zones, because we use a different
    # shear modulus treatment elsewhere.
    k = min(which(sourcezone_parameters$sourcename == all_source_zones[i]))
    if(sourcezone_parameters$scaling_relation[k] != 'Strasser') next

    print(all_source_zones[i])

    # Make the panel plot, and return the mw-rate results
    tmp = plot_sourcezone_rate_curve_with_fixed_and_variable_mu(all_source_zones[i], slip_type=slip_type)

    # Store a globally integrated rate
    if(i == 1){
        global_rates = tmp$mw_rates
    }else{
        global_rates[,2:3] = global_rates[,2:3] + tmp$mw_rates[,2:3]
    }

    # Store the results
    rate_curves_store[[i]] = tmp
}
#
# Plot the globally integrated rate
#
plot(global_rates[,1], pmax(global_rates[,2], 1e-100), t='o', log='y', 
    ylim=c(1e-04,10), xlab="Mw", ylab='Exceedance rate')
points(global_rates[,1], pmax(global_rates[,3], 1e-100), t='l', col='red')
grid()
title('Globally integrated rates')

dev.off()

save.image(paste0('variable_mu_checks_session_', slip_type, '.Rdata'))

# Look at how much mw changes at different return periods as a result of varying mu
# We generally see a small shift that doesn't change much with return period, but
# varies by source-zone.
mw_change = lapply(rate_curves_store, 
    f<-function(y){ 
        x = y$mw_rates
        try(
        # Note that x[,2] and x[,3] are rates, which are decreasing but might
        # finish with a collection of zeros. In that case, we want to use the 
        # minimum magnitude corresponding to rate=0 for interpolation. The ties='min'
        # argument takes care of that.
        approx(x[,2], x[,1], xout=c(1/100, 1/500, 1/1000, 1/2500), ties='min')$y - 
        approx(x[,3], x[,1], xout=c(1/100, 1/500, 1/1000, 1/2500), ties='min')$y
    )}
)


