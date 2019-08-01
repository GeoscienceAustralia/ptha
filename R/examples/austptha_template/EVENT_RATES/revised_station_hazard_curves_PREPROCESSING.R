#
# Extract data from the 'compute_rates_all_sources_session.RData' session
# to enable an improved stage-vs-rate-percentile calculation.
#
#

out_dir = 'preprocessed_source_rate_revised_stage_percentiles/'
dir.create(out_dir, showWarnings=FALSE)

library(rptha)
rc = new.env()
load('compute_rates_all_sources_session.RData', envir=rc)

# All source names (full sources only)
source_names = unique(rc$sourcezone_parameters$sourcename)

# Stages at which we compute exceedance-rates
stage_seq = rc$config$stage_seq

for(source_name in source_names){

    all_segments_on_source = which(unlist(lapply(rc$source_envs, f<-function(x) x$source_name)) == source_name)

    # Skip source-zones with zero weight (typically ones that were replaced with SLAB2)
    if(rc$source_envs[[all_segments_on_source[1]]]$sourcepar$sourcezone_parameters_row$row_weight == 0){
        next
    }

    # Loop over segments on source_name, and store Mw-frequency information for all logic-tree branches
    # for each segment (and the unsegmented branch).
    seg = vector(mode='list', length=length(all_segments_on_source))
    for(i in 1:length(all_segments_on_source)){

        seg[[i]] = list() # Store key info for this segment
        seg_ind = all_segments_on_source[i]



        #
        # Get all the logic-tree data
        #
        seg_mw_rate_curve_data = rc$source_envs[[seg_ind]]$mw_rate_function(NA, 
            return_all_logic_tree_branches=TRUE)
        # Mw values at which logic-tree curves are tabulated
        seg[[i]]$Mw_seq = seg_mw_rate_curve_data$Mw_seq
        # Tabulated values for all logic-tree curves
        seg[[i]]$all_rate_matrix = seg_mw_rate_curve_data$all_rate_matrix
        # Curve weights (fixed mu)
        seg[[i]]$logic_tree_weights = seg_mw_rate_curve_data$all_par_prob
        # Curve weights (variable mu) -- recall we treat this as though it is an Mw error
        seg[[i]]$logic_tree_weights_varyMu = seg_mw_rate_curve_data$all_par_prob_with_Mw_error
        seg[[i]]$source_name = rc$source_envs[[seg_ind]]$source_name
        seg[[i]]$segment_name = rc$source_envs[[seg_ind]]$segment_name
        seg[[i]]$is_a_segment = rc$source_envs[[seg_ind]]$is_a_segment
        seg[[i]]$row_weight = rc$source_envs[[seg_ind]]$sourcezone_parameters_row$row_weight


        # Open files with event information
        nc_events = list()
        nc_events$uniform = nc_open(paste0('../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_',
                source_name, '.nc'), readunlim=FALSE)
        nc_events$stochastic = nc_open(paste0('../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_',
                source_name, '.nc'), readunlim=FALSE)
        nc_events$variable_uniform = nc_open(paste0('../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_variable_uniform_slip_earthquake_events_',
                source_name, '.nc'), readunlim=FALSE)

        #
        # Get the event 'fixed-mu-magnitudes'. This depends only on the slip-type
        #
        # When reading magnitudes, use rounding to ensure that finite-precision storage did not
        # introduce minor differences in 'the same' magnitudes. Since our fixed-mu magnitudes have
        # rates of 7.2, 7.3, .... 9.8, rounding to 2-decimal-places is "more than
        # sufficient" to achieve this.
        event_Mw = list()
        for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){
            event_Mw[[slip_type]] = round( ncvar_get(nc_events[[slip_type]], 'Mw'), 2)
        }

        # 
        # Get the event conditional probabilities.
        event_conditional_probability = list()
        # This depends on the slip type and the rigidity type (except for uniform slip, where rigidity doesn't matter).
        # Recall different bias-adjustment functions were derived for stochastic/variable uniform slip, depending on rigidity. 
        #
        event_conditional_probability$uniform = list()
        event_conditional_probability$uniform$fixed_mu    = rc$source_envs[[seg_ind]]$event_conditional_probabilities
        event_conditional_probability$uniform$variable_mu = rc$source_envs[[seg_ind]]$event_conditional_probabilities
        #
        # For stochastic and variable_uniform, the uniform rates are spread over
        # the corresponding events (often 15 of them), using a conditional probability model based on maximum slip.
        #
        # To compute their conditional probabilities, the following variables are useful 
        event_rate_mean = list()
        for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){
            event_rate_mean[[slip_type]] = list()
        }
        event_rate_mean$uniform$fixed_mu = ncvar_get(nc_events$uniform, 'rate_annual')
        event_rate_mean$uniform$variable_mu = ncvar_get(nc_events$uniform, 'variable_mu_rate_annual')
        event_rate_mean$stochastic$fixed_mu = ncvar_get(nc_events$stochastic, 'rate_annual')
        event_rate_mean$stochastic$variable_mu = ncvar_get(nc_events$stochastic, 'variable_mu_rate_annual')
        event_rate_mean$variable_uniform$fixed_mu = ncvar_get(nc_events$variable_uniform, 'rate_annual')
        event_rate_mean$variable_uniform$variable_mu = ncvar_get(nc_events$variable_uniform, 'variable_mu_rate_annual')

        # This function is useful for double-checking our calculation of the stochastic/variable_uniform
        # bias adjustment (used to 'spread' the conditional probability of the uniform event over the multiple corresponding
        # stochastic or variable_uniform events)
        check_bias_adjuster<-function(bias_adjuster, uniform_event_row, uniform_event_rate){
            # When we sum the 'bias_adjuster' in groups with the same uniform_event_row, it should
            # give 1, unless all values are zero (in which case the rate is zero, so no problem)
            check = aggregate(bias_adjuster, by=list(uniform_event_row), sum)
            if(!all((abs(check$x  - 1) < 1.0e-06) | (check$x == 0 & uniform_event_rate==0))) stop('Problem with bias adjuster')
        }

        # A very small number for zero-divide protection
        TINY = (.Machine$double.xmin)*10

        # We need to match between stochastic and uniform events (and similarly, variable_uniform and uniform events)
        uniform_event_row = list()
        uniform_event_row$stochastic = ncvar_get(nc_events$stochastic, 'uniform_event_row')
        uniform_event_row$variable_uniform = ncvar_get(nc_events$variable_uniform, 'uniform_event_row')

        #
        # Stochastic and variable_uniform conditional probabilities
        #
        for(slip_type in c('stochastic', 'variable_uniform')){
            for(mu_type in c('fixed_mu', 'variable_mu')){
                # The accounts for the 'max-slip' related bias adjustment
                bias_adjuster = event_rate_mean[[slip_type]][[mu_type]] / 
                    pmax(event_rate_mean$uniform[[mu_type]][uniform_event_row[[slip_type]]], TINY)
                # Double-check the result is sensible!
                check_bias_adjuster(bias_adjuster, uniform_event_row[[slip_type]], event_rate_mean$uniform[[mu_type]])
                # Conditional probability --> spread the uniform-slip conditional probability over child events
                event_conditional_probability[[slip_type]][[mu_type]] = bias_adjuster *
                    event_conditional_probability$uniform[[mu_type]][uniform_event_row[[slip_type]]]
            }
        }

        seg[[i]]$event_conditional_probability = event_conditional_probability
        seg[[i]]$event_Mw = event_Mw

        # Close event information files
        for(ii in 1:length(nc_events)) nc_close(nc_events[[ii]])

    }

    output_list = list(seg = seg, source_name=source_name, stage_seq=stage_seq)
    output_file = paste0(out_dir, 'preprocessed_rate_info_', source_name, '.RDS')
    saveRDS(output_list, output_file)
    rm(output_list)

}
