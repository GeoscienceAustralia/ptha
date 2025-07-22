#
# Implement the F^m_e and R^m_e statistics from Davies (2019), adapted to the situation
# where we have a random sample of PTHA18 scenarios matching the bulk earthquake properties
# (rather than every scenario in PTHA18 as in Davies (2019)).
#
source('analysis_data_and_functions.R')

#' Default function used at a single tide gauge for a single event to summarise the models and data.
#' Computes the fraction of models below the data (using a quantile equation that avoids 0 or 1)
default_F_m_e_d_function = function(data_at_gauge, all_models_at_gauge){
    N = length(all_models_at_gauge)
    (sum(data_at_gauge > all_models_at_gauge) + 0.5)/(N+1)
}

#' Default function that is used when there are multiple tide gauges for a single event, to provide
#' a single number from the multiple values.
default_multi_gauge_summary_function = function(x){median(x, na.rm=TRUE)}


#' Statistic similar to F^m_e in Davies (2019)
#'
#' Compute the fraction of models below the data at every tide gauge for every
#' event (or some other F_m_e_d_function), then for each unique event take the
#' median over gauges (or some other multi_gauge_summary_function). Return the
#' latter statistic for every unique event.
#' 
#' @param event_id vector of character event IDs, same length as model_values_list
#' @param model_values_list list (one value for every unique "event & tide
#' gauge" combination).  Each entry contains a vector of model statistics (one
#' value per model scenario for that event).
#' @param data_values_list list with the same length as model_values_list (i.e.
#' one value for every unique "event & tide gauge" combination).  Each entry is
#' EITHER a single value with the data statistic for the given "event & tide
#' gauge", or a vector with the same length as the corresponding index of
#' model_data_list. (The latter potentially allows different data for each
#' model, which normally we wouldn't want - but could happen if e.g.  the
#' statistics were computed in a way that depended on the modelled tsunami
#' arrival time, different for each model scenario).
#' @param F_m_e_d_function A function that is applied for each individual tide-gauge and event combination,
#' to summarise the data and model values. In Davies (2019) this gives the quantile of the observation 
#' relative to the model, which is also implemented here (using a quantile equation that avoids 0 or 1).
#' @param multi_gauge_summary_function Applied over tide gauges, for each event separately.
#' Collapses the value at each tide gauge to a single value (important because
#' results at different tide gauges are not usually independent).
#' @return vector with one entry for each unique(event_id), giving the
#' median_over_tide_gauges(fraction_of_models_below_data) for each event (or
#' the result of the user-specified multi_gauge_summary_function, if not the median).
#'
median_over_tide_gauges_of_fraction_models_below_data_per_event<-function(
    event_id, 
    model_values_list, 
    data_values_list,
    F_m_e_d_function=default_F_m_e_d_function,
    multi_gauge_summary_function=default_multi_gauge_summary_function){

    unique_event_ids = unique(event_id)
    median_frac_below = rep(NA, length(unique_event_ids))
    data_is_NA = unlist(lapply(data_values_list, function(x) any(is.na(x))))

    # Sanity check the length and size of inputs
    ndata = unlist(lapply(data_values_list, length))
    nmodels = unlist(lapply(model_values_list, length))
    if( (length(ndata) != length(nmodels)) | 
        (any((ndata > 1) & (ndata != nmodels))) ){
        stop('For each tide gauge, the data_values_list must either have a single entry, or one entry per model scenario')
    }

    for(i in 1:length(unique_event_ids)){
         # Find all tide gauges with non-missing data for the i'th event ID 
        tg_inds = which(event_id == unique_event_ids[i] & !data_is_NA)
        
        if(length(tg_inds) == 0){
            # Issues for some statistics when the data finishes before the model
            median_frac_below[i] = NA
        }else{
            # At each tide gauge with non-missing data, compute the fraction of
            # models below the data (or some other user-specified statistic)
            frac_below = rep(NA, length(tg_inds))
            for(j in 1:length(tg_inds)){
                tgj = tg_inds[j]
    
                # Fraction of models below the data at this tide gauge (or some other statistic)
                frac_below[j] = F_m_e_d_function(data_values_list[[tgj]], model_values_list[[tgj]])
            }
            # Record the median over all non-missing tide gauges (or some other summary statistic)
            median_frac_below[i] = multi_gauge_summary_function(frac_below)
        }
    }
    names(median_frac_below) = unique_event_ids
    return(median_frac_below)
}

.basic_example_median_over_tide_gauges_of_fraction_models_below_data_per_event<-function(){
    # Suppose we have 3 events with
    #   - event 1: 2 tide gauges, 4 scenarios
    #   - event 2: 1 tide gauge, 3 scenarios
    #   - event 3: 2 tide gauges, 4 scenarios
    # The input data could be
    event_id = c('e1', 'e1', 'e2', 'e3', 'e3')
    # The model data could be (
    model_values_list = list(
        c(1.3, 1.2, 1.5, 1.4), # event 1, site 1, 4 scenarios
        c(0.5, 0.7, 0.1, 0.2), # event 1, site 2, 4 scenarios
        c(2.2, 2.7, 2.1     ), # event 2, site 1, 3 scenarios
        c(3.3, 3.2, 3.5, 3.4), # event 3, site 1, 4 scenarios
        c(4.5, 4.7, 4.1, 4.2)) # event 3, site 2, 4 scenarios
    # Usually the data will be constant for each model scenario, but in principle it could vary
    # (e.g. for statistics that rely on a property of the model scenario, such
    # as its arrival time, unless we are careful to define that in the same way
    # for every scenario)
    data_values_list = list(
        c(1.25, 1.25, 1.25, 1.25), # event 1, site 1, 4 scenarios
        c( 0.6,  0.6,  0.6,  0.6), # event 1, site 2, 4 scenarios
        c( 2.8,  2.8,  2.8      ), # event 2, site 1, 3 scenarios
        c( 3.1,  3.1,  3.1,  3.1), # event 3, site 1, 4 scenarios
        c( 4.4,  4.4,  4.4,  4.4)) # event 3, site 2, 4 scenarios

    # The results will be:
    # - Event 1, 
    #   - site1: [1/4 model is lower] - (1+0.5)/(4+1)
    #   - site2[3/4 models are lower] - (3+0.5)/(4+1)
    #    median = 1/2
    # - Event 2
    #   - site 1: [3/3 models are lower] - (3+0.5)/(3+1)
    #   median = 0.875
    # - Event 3
    #   - site 1: 0/4 models are lower - (0+0.5)/(4+1)
    #   - site 2: 2/4 models are lower - (2 + 0.5)/(4+1)
    #   median = 0.3
    #
    tmp = median_over_tide_gauges_of_fraction_models_below_data_per_event(event_id, model_values_list, data_values_list)
    print(tmp)
    stopifnot( all(abs(tmp - c(0.5, 0.875, 0.3)) < 1e-08)  )

    #
    # Make another example with ties in the models 
    #

    model_values_list_with_ties = list(
        c(1.3, 1.2, 1.5, 1.4), # event 1, site 1, 4 scenarios
        c(0.5, 0.7, 0.1, 0.2), # event 1, site 2, 4 scenarios
        c(2.2, 2.7, 2.1     ), # event 2, site 1, 3 scenarios
        c(3.3, 3.2, 3.2, 3.4), # event 3, site 1, 4 scenarios
        c(4.5, 4.1, 4.1, 4.2)) # event 3, site 2, 4 scenarios

    data_values_list_with_ties = list(
        c(1.25, 1.25, 1.25, 1.25), # event 1, site 1, 4 scenarios
        c( 0.6,  0.6,  0.6,  0.6), # event 1, site 2, 4 scenarios
        c( 2.8,  2.8,  2.8      ), # event 2, site 1, 3 scenarios
        c( 3.2,  3.2,  3.2,  3.2), # event 3, site 1, 4 scenarios
        c( 4.1,  4.1,  4.1,  4.1)) # event 3, site 2, 4 scenarios

    a = median_over_tide_gauges_of_fraction_models_below_data_per_event(event_id, model_values_list_with_ties, data_values_list_with_ties)
    print(a)
}


#' Swap the data with a random sample from the set of "all model scenarois and the data" for each event.
#' This means the "real data" is likely to be included with the random models, and in those case the 
#' "fake data" will be one of the random models. This is useful for simulating statistics under the 
#' null hypothesis that the data is just another model simulation (since then we would be just as likely
#' to have observed any of these swapped model/data results)
#'
#' If the data is NA then don't swap (as these cases should anyway be dropped in summary statistics).
#'
#' @param event_id vector of character event IDs, same length as model_values_list
#' @param model_values_list list (one value for every unique "event & tide
#' gauge" combination).  Each entry contains a vector of model statistics (one
#' value per model scenario for that event).
#' @param data_values_list list with the same length as model_values_list (i.e.
#' one value for every unique "event & tide gauge" combination).  Each entry is
#' EITHER a single value with the data statistic for the given "event & tide
#' gauge", or a vector with the same length as the corresponding index of
#' model_data_list. (The latter potentially allows different data for each
#' model, which normally we wouldn't want - but could happen if e.g.  the
#' statistics were computed in a way that depended on the modelled tsunami
#' arrival time, different for each model scenario).
#'
swap_data_and_random_scenario_from_model_or_data<-function(
    event_id, 
    model_values_list, 
    data_values_list){

    N_model_all = unlist(lapply(model_values_list, length))
    N_data_all  = unlist(lapply(data_values_list, length))
    unique_event_id = unique(event_id)

    # Make data that will be updated below
    new_data_values_list = data_values_list
    new_model_values_list = model_values_list

    swap_store = rep(NA, length(unique_event_id))
    for(i in 1:length(unique_event_id)){
        # Find rows matching this event
        k = which(event_id == unique_event_id[i])

        # Should be the same number of model scenarios for each event
        stopifnot(all(N_model_all[k] == N_model_all[k[1]]))

        # Each tide gauge should either have a single data value, or one per scenario.
        # (if one per scenario then anyway they are typically identical)
        stopifnot(all(N_data_all[k] == 1 | N_data_all[k] == N_model_all[k]))

        # Should always have > 1 model scenario
        stopifnot(all(N_model_all[k] > 1))

        # Choose 1 model scenario at random (or 0 == use data as is)
        ind = sample(0:N_model_all[k[1]], size=1)
        swap_store[i] = ind

        for(j in k){
            # Don't swap if the data is NA. In these cases the gauge will be
            # dropped from later statistical calculations.
            data_is_NA = any(is.na(data_values_list[[j]]))
            if(!data_is_NA){
                Nd = length(data_values_list[[j]])
                if(Nd > 1){
                    # There is a separate data value for each model. Check they are all identical, as the
                    # logic of resampling won't hold otherwise
                    if(!all(data_values_list[[j]] == data_values_list[[j]][1])){
                        stop('Not all data values are identical for a single site and event, the code logic does not treat this case')
                    }
                }

                if(ind > 0){
                    # Swap the data and model

                    model_value_to_make_data = model_values_list[[j]][ind]
                    # Replace the data with the model
                    new_data_values_list[[j]] = rep(model_value_to_make_data, Nd)
                    # Replace the model with the data (first entry OK since they are all identical as checked above)
                    new_model_values_list[[j]][ind] = data_values_list[[j]][1] 
                }else{
                    # We sampled the data, no swapping needed
                }
            }
        }
    }
    # Even though we only modified the data, here we return all arguments for convenience
    return(list(event_id = event_id, model_values_list=new_model_values_list, data_values_list=new_data_values_list, swapped_index = swap_store))
}

#' For each event, compute a "null hypothesis" distribution of the event statistic by
#' - For each scenario
#'    Swap the data with the model scenario
#'    Compute the event statistic
#' This gives a set of event statistics that can be compared with the "real" statistic. The
#' null hypothesis is that the observation is just another sample from the model, so under that
#' hypothesis we are just as likely to have observed any one of these synthetic statistics. 
#'
#' No swapping is applied tide gauges where the data is NA, since the statistic should anyway
#' take care of these cases.
#'
#' @param event_id vector of character event IDs, same length as model_values_list
#' @param model_values_list list (one value for every unique "event & tide
#' gauge" combination).  Each entry contains a vector of model statistics (one
#' value per model scenario for that event).
#' @param data_values_list list with the same length as model_values_list (i.e.
#' one value for every unique "event & tide gauge" combination).  Each entry is
#' EITHER a single value with the data statistic for the given "event & tide
#' gauge", or a vector with the same length as the corresponding index of
#' model_data_list. (The latter potentially allows different data for each
#' model, which normally we wouldn't want - but could happen if e.g.  the
#' statistics were computed in a way that depended on the modelled tsunami
#' arrival time, different for each model scenario).
#' @param include_no_swap_case If TRUE, the output empirical distribution will also include
#' the 'real' event statistic (where we don't swap the model and data). 
#' @param F_m_e_d_function A function that is applied for each individual tide-gauge and event combination,
#' to summarise the data and model values. In Davies (2019) this gives the quantile of the observation 
#' relative to the model, which is also implemented here (using a quantile equation that avoids 0 or 1).
#' @param multi_gauge_summary_function Applied over tide gauges, for each event separately.
#' Collapses the value at each tide gauge to a single value (important because
#' results at different tide gauges are not usually independent).
#'
compute_distribution_of_event_statistic_when_data_is_swaped_with_each_scenario<-function(
    event_id,
    model_values_list,
    data_values_list,
    include_no_swap_case = TRUE,
    F_m_e_d_function=default_F_m_e_d_function,
    multi_gauge_summary_function=default_multi_gauge_summary_function){

    unique_event_id = unique(event_id)
    # Number of model scenarios for each event.
    N_model_all = unlist(lapply(model_values_list, length))

    # Store outputs
    median_over_tg_of_frac_below_for_each_scenario_and_event = vector(mode='list', length=length(unique_event_id))
    names(median_over_tg_of_frac_below_for_each_scenario_and_event) = unique_event_id

    for(i in 1:length(unique_event_id)){
        k = which(event_id == unique_event_id[i])

        # For the i'th event, at every site there should be the same number of
        # model scenarios
        N_model = N_model_all[k[1]]
        stopifnot(all(N_model_all[k] == N_model))

        # Number of tide gauges
        N_gauges = length(k)


        if(include_no_swap_case){
            j0 = 0
        }else{
            j0 = 1
        }
        frac_below_for_each_model = rep(NA, N_model + include_no_swap_case)
        for(j in j0:N_model){
            # Model values for the i'th event 
            model_values_this_event = model_values_list[k]
            data_values_this_event = data_values_list[k]
            frac_below = rep(NA, N_gauges)

            for(tgi in 1:N_gauges){
                data_1 = data_values_this_event[[tgi]][1]
                if(!is.na(data_1) & length(data_values_this_event[[tgi]]) > 1){
                    # If there are multiple data values at the gauge (one per model scenario) then they
                    # should all be identical for the resampling logic to work. This is the usual situation.
                    stopifnot(all(data_values_this_event[[tgi]] == data_1))
                }

                if(!is.na(data_1)){
                    new_model_values_this_event_and_gauge = model_values_this_event[[tgi]]
                    new_data_value_this_event_and_gauge = data_1

                    # Replace the j'th model value with the data, and the data value with the j'th model value
                    if(j > 0){
                        new_model_values_this_event_and_gauge[j] = data_1
                        new_data_value_this_event_and_gauge = model_values_this_event[[tgi]][j]
                    }

                    frac_below[tgi] = F_m_e_d_function(new_data_value_this_event_and_gauge, new_model_values_this_event_and_gauge)
                }
            }
            # Use the median to collapse the fraction over multiple gauges
            frac_below_for_each_model[j+include_no_swap_case] = multi_gauge_summary_function(frac_below)
        }


        median_over_tg_of_frac_below_for_each_scenario_and_event[[i]] = frac_below_for_each_model
    }
    return(median_over_tg_of_frac_below_for_each_scenario_and_event)
}


#' Compute the R^m_e statistic from Davies (2019), adapted to our situation
#' where we only have a sample of model scenarios (not every model scenario)
#'
#' This first computes F^m_e from the data (one value per event), then computes
#' the distribution of F^m_e if we swaped each model scenario with the data (leading to one empirical distribution
#' per event), then uses the latter to convert the former to a number R^m_e in
#' (0, 1) which is approximately uniformly distributed under the null hypothesis that the data behaves
#' like a random sample from the model (approximate due to discretization caused by a
#' finite number of model scenarios, ties, etc). Optionally return some
#' statistcs that check for uniformity.
#'
#' @param event_id vector of character event IDs, same length as model_values_list
#' @param model_values_list list (one value for every unique "event & tide gauge" combination). 
#' Each entry contains a vector of model statistics (one value per model scenario for that event).
#' @param data_values_list list with the same length as model_values_list (i.e.
#' one value for every unique "event & tide gauge" combination).  Each entry is
#' EITHER a single value with the data statistic for the given "event & tide
#' gauge", or a vector with the same length as the corresponding index of
#' model_data_list. (The latter potentially allows different data for each
#' model, which normally we wouldn't want - but could happen if e.g.  the
#' statistics were computed in a way that depended on the modelled tsunami
#' arrival time, different for each model scenario).
#' @param return_stats If FALSE, return R^m_e. If TRUE, return a list with the latter as well as some statistics
#' that can be used to check whether the distribution is close to uniform.
#' @param multi_gauge_summary_function The function used in computing F^m_e, to
#' collapse results at multiple tide gauges for a single event.
#' @return list with one entry for each unique(event_id), an array giving the
#' F^m_e value if the data is the model scenario. The number of rows is equal
#' to N_repeats. These values define the sampling distribution of F^m_e for the event.
#'
compute_Rme<-function(
    event_id, 
    model_values_list, 
    data_values_list,
    return_stats = FALSE,
    F_m_e_d_function=default_F_m_e_d_function,
    multi_gauge_summary_function= default_multi_gauge_summary_function){

    F_me_data = median_over_tide_gauges_of_fraction_models_below_data_per_event(
        event_id, 
        model_values_list, 
        data_values_list,
        F_m_e_d_function=F_m_e_d_function,
        multi_gauge_summary_function=multi_gauge_summary_function) 

    F_me_null = compute_distribution_of_event_statistic_when_data_is_swaped_with_each_scenario(
        event_id=event_id,
        model_values_list=model_values_list,
        data_values_list=data_values_list,
        include_no_swap_case = FALSE, # Don't include the F_me_data stat -- since our definition of R_me_data should not have this
        F_m_e_d_function=F_m_e_d_function,
        multi_gauge_summary_function=multi_gauge_summary_function)

    stopifnot(all(names(F_me_data) == names(F_me_null)))

    R_me_data = F_me_data * 0
    for(i in 1:length(F_me_data)){
        ## fraction = (count_below + 1/2)/(count + 1) -- avoid 0,1, which might mess with the ad.test
        R_me_data[i] = (sum(F_me_data[i] > F_me_null[[i]]) + 0.5)/(length(F_me_null[[i]]) + 1)
    }

    if(!return_stats){
        return(R_me_data)
    }else{
        library(ADGofTest)
        sd_R_me_data = sd(R_me_data, na.rm=TRUE)
        # Compute a 'fraction below' compared with samples from a uniform distribution
        sd_test_N = sum(!is.na(R_me_data))
        sd_test_fraction_of_sds_below = mean(sd_R_me_data > replicate(1e+04, sd(runif(sd_test_N))) )
        
        outputs = list(R_me_data = R_me_data,
            F_me_data = F_me_data,
            adtest = ad.test(R_me_data, punif),
            sd_R_me_data = sd_R_me_data,
            sd_test_fraction_of_sds_below=sd_test_fraction_of_sds_below)
        return(outputs)
    }
}

.basic_example_Rme<-function(){

    # The input data could be
    event_id = c('e1', 'e1', 'e2', 'e3', 'e3')
    # The model data could be (
    model_values_list = list(
        c(1.3, 1.2, 0.1, 1.4), # event 1, site 1, 4 scenarios
        c(0.5, 0.7, 0.1, 0.2), # event 1, site 2, 4 scenarios
        c(2.2, 2.7, 2.1     ), # event 2, site 1, 3 scenarios
        c(3.3, 3.2, 3.5, 3.4), # event 3, site 1, 4 scenarios
        c(4.5, 4.7, 4.1, 4.2)) # event 3, site 2, 4 scenarios
    # Usually the data will be constant for each model scenario, but in principle it could vary
    # (e.g. for statistics that rely on a property of the model scenario, such
    # as its arrival time, unless we are careful to define that in the same way
    # for every scenario)
    data_values_list = list(
        c(1.25, 1.25, 1.25, 1.25), # event 1, site 1, 4 scenarios
        c( 0.6,  0.6,  0.6,  0.6), # event 1, site 2, 4 scenarios
        c( 2.8,  2.8,  2.8      ), # event 2, site 1, 3 scenarios
        c( 3.1,  3.1,  3.1,  3.1), # event 3, site 1, 4 scenarios
        c( 4.4,  4.4,  4.4,  4.4)) # event 3, site 2, 4 scenarios

    F_me = median_over_tide_gauges_of_fraction_models_below_data_per_event(
        event_id, model_values_list, data_values_list)

    null_distr = compute_distribution_of_event_statistic_when_data_is_swaped_with_each_scenario(
        event_id, model_values_list, data_values_list)

    R_me = compute_Rme(event_id, model_values_list, data_values_list)
}

#' Helper function to get data that is useful for several tests
get_data_for_testing<-function(test_variable, test_slip_type, use_downsampled_model=FALSE){

    # 
    # Get some model simulations for testing.
    # It shouldn't matter exactly which variables are used. We will generate
    # fake data by randomly sampling from the model, the point being to check
    # the distribution of various R_me statistics under the null hypothesis
    # (which is known theoretically under some reasonable approximations, e.g.,
    # approximating uniform discrete distributions with 60 different values
    # with uniform distributions). 
    # The code here helps to confirm that the tests are correctly implemented
    # and do not suffer too much from approximations (such as the limited
    # number random scenarios used in this study, 60 per event and model).
    #

    if(use_downsampled_model){
        es = event_stats_downsampled_model
    }else{
        es = event_stats
    }

    sites_compared = lapply(unique(es$event_name), function(event_name){
        get_table_indices_to_process(es, event_name, 
            good_nearshore_and_close_and_highres = TRUE,
            expand_multi_counts=TRUE, run_type='random_like_historic',
            rigidity_type='constant', exclude_batch2_doubleups=FALSE)
        })
    all_inds = unlist(sites_compared) 

    #
    # Get key variables
    #
    slip_type = as.factor(es$slip_type[all_inds])
    site = as.factor(es$sites[all_inds])
    event = as.factor(es$event_name[all_inds])
    event_int = es$event_name_int[all_inds]
    obs   = es[[paste0('data_' , test_variable)]][ all_inds]
    model = es[[paste0('model_', test_variable)]][ all_inds]

    k = which(slip_type == test_slip_type)

    #
    # Split data into a list with expected structure
    #
    bp_obs = aggregate(
        obs[k], 
        by=list(site=site[k], event=event[k], event_int=event_int[k]), 
        function(x) list(x))
    bp_model = aggregate(
        model[k], 
        by=list(site=site[k], event=event[k], event_int=event_int[k]), 
        function(x) list(x))
    # Ensure the ordering is identical
    stopifnot(all((bp_obs$event == bp_model$event) & (bp_obs$site == bp_model$site)))

    return(list(bp_model=bp_model, bp_obs=bp_obs))
}

#' Check that the R^m_e statistics are distributed as expected under the null hypothesis.
#'
#' This was run under several cases, e.g. 
#'
#' # HS model
#' tmp = test_compute_Rme(test_slip_type='HS') 
#'
#' # FAUS model. This has more ties than HS, so good to see the distribution reasonably close to uniform despite that
#' tmp = test_compute_Rme(test_slip_type='FAUS')
#' 
#' # FAUS model with lastday values. This has NA values for tide gauges where the last day was truncated,
#' # so good to see it still works in that case (it should -- the underlying functions do treat this case)
#' tmp = test_compute_Rme(test_variable='max_lastday_during_obs', test_slip_type='FAUS')
#'
test_compute_Rme<-function(test_variable='max_during_obs', test_slip_type='HS', use_downsampled_model=FALSE, N_sims=1000, iseed=123){

    # Get some test data.
    tmp = get_data_for_testing(test_variable, test_slip_type, use_downsampled_model)
    bp_model = tmp$bp_model
    bp_obs = tmp$bp_obs
    rm(tmp)

    test_function<-function(ignored){
        source('stats_under_null_hypothesis.R')
        # Use the test data to make synthetic data
        fake_data_list = swap_data_and_random_scenario_from_model_or_data(
            event_id = bp_model$event, 
            model_values_list = bp_model$x, 
            data_values_list = bp_obs$x
            )

        Rme_stats = compute_Rme(
            event_id = fake_data_list$event_id,
            model_values_list = fake_data_list$model_values_list,
            data_values_list = fake_data_list$data_values_list,
            return_stats=TRUE)
        return(Rme_stats)
    }

    MC_CORES = 16
    library(parallel)
    cl = makeCluster(MC_CORES)
    # Make random numbers reproducible yet distinct on each thread. Uses RNGkind("L'Ecuyer-CMRG")
    clusterSetRNGStream(cl, iseed=iseed)
    # Get the variables on the cluster
    clusterExport(cl, varlist=c('bp_model', 'bp_obs', 'test_function'), envir=environment())

    # Compute the statistic with random data a given number of times
    repeated_sims = parLapply(cl, 1:N_sims, test_function)

    # The following variables should be nearly uniformly distributed if everything is working OK
    adtest_p_values = unlist(lapply(repeated_sims, function(x) x$adtest$p.value))
    sd_fractions_below = unlist(lapply(repeated_sims, function(x) x$sd_test_fraction_of_sds_below))

    par(mfrow=c(1,2))
    plot(ecdf(adtest_p_values)); abline(0,1,col='red'); title(sub='Should be nearly uniformly distributed')
    plot(ecdf(sd_fractions_below)); abline(0,1,col='red'); title(sub='Should be nearly uniformly distributed')

    # Are approximately 5% of the adtest values below 0.05? (as expected, since we simulated data that
    # satisfies the null hypothesis)
    print(binom.test(sum(adtest_p_values < 0.05), n=length(adtest_p_values), p = 0.05))

    # Are approximately 5% of the sd_fractions_below below 0.05? (as expected, since we simulated data that
    # satisfies the null hypothesis)
    print(binom.test(sum(sd_fractions_below < 0.05), n=length(sd_fractions_below), p = 0.05))
    # Are approximately 95% of the sd_fractions_below below 0.95? (as expected, since we simulated data that
    # satisfies the null hypothesis)
    print(binom.test(sum(sd_fractions_below < 0.95), n=length(sd_fractions_below), p = 0.95))

    print('Chance of ADGoF p < 0.05')
    print(mean(adtest_p_values < 0.05))
    print('quantiles of p value ')
    print(quantile(adtest_p_values, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)))
    print('Chance of sd_fractions_below < 0.05')
    print(mean(sd_fractions_below < 0.05))
    print('Chance of sd_fractions_below > 0.95')
    print(mean(sd_fractions_below > 0.95))
    print('Chance of ADGoF p < 0.05 OR sd_fractions_below < 0.05')
    print(mean(adtest_p_values < 0.05 | sd_fractions_below < 0.05))
    print('Chance of ADGoF p < 0.05 OR sd_fractions_below > 0.95')
    print(mean(adtest_p_values < 0.05 | sd_fractions_below > 0.95))

    stopCluster(cl)
    return(invisible(repeated_sims))


}

compute_expected_frequency_of_exceedances_under_null_hypothesis<-function(test_variable='max_during_obs', 
    test_slip_type='HS', definition_failure='under_estimate_or_equal', use_downsampled_model = FALSE, N_sims=1e+05){
    #' Estimate how often the model would fail to envelope the data under the null hypothesis.
    #' We define enveloping the data as "having model values both smaller and larger than the data".
    #' It's worth commenting about the case where the model is equal to the data
    #' - This will never happen for 'real data'. But it can happen with synthetic data as some scenarios are sampled more than once.
    #' - If we accept "equal to the data" as enveloping, then the expected rate
    #'   of failure from synthetic data will be smaller than if we require "both greater
    #'   and less", not just "greater than or equal" or "less than or equal".
    #'   - Thus, the "both greater and less" approach will increase the rate that we expect to see failures according to the null hypothesis.
    #'     - If we include "equal" as success, then for some events we could never get failure if the data is replaced with synthetic data from model scenarios. This
    #'       would happen if the maxima are achieved by scenarios that were sampled more than once. It seems undesirable.
    #'   - By requiring "both greater and less" it is more difficult for us to say that a real model is failing too often. 
    #'     - That seems a good thing -- we don't want our treatment of ties to
    #'       make it easier for us to get "statistically significant" results -- we'd
    #'       prefer to make it more difficult, so we are less likely to falsely reject the
    #'       null hypothesis. 

    tmp = get_data_for_testing(test_variable, test_slip_type, use_downsampled_model)
    bp_model = tmp$bp_model
    bp_obs = tmp$bp_obs
    rm(tmp)

    #
    # Estimate how likely we would be to see a given number of events with at
    # least 1 site where the model does not envelope the data.
    #

    # Alternative to F_m_e_d -- return 1 if the "data" is at an extreme
    is_data_extreme_or_equal=function(data_at_gauge, all_models_at_gauge){
        # NB: data_at_gauge is a constant value (even if it is a vector with length > 1)
        1*(any(data_at_gauge <= min(all_models_at_gauge)) | any(data_at_gauge >= max(all_models_at_gauge)))
    }
    is_data_extreme=function(data_at_gauge, all_models_at_gauge){
        # NB: data_at_gauge is a constant value (even if it is a vector with length > 1)
        1*(any(data_at_gauge < min(all_models_at_gauge)) | any(data_at_gauge > max(all_models_at_gauge)))
    }
    is_data_larger_or_equal=function(data_at_gauge, all_models_at_gauge){
        1*(any(data_at_gauge >= max(all_models_at_gauge)))
    }
    is_data_larger=function(data_at_gauge, all_models_at_gauge){
        1*(any(data_at_gauge > max(all_models_at_gauge)))
    }

    if(definition_failure == 'under_estimate_or_equal'){
        single_gauge_fun = is_data_larger_or_equal
    }else if(definition_failure == 'not_enveloped_or_equal'){
        single_gauge_fun = is_data_extreme_or_equal
    }else if(definition_failure == 'under_estimate'){
        single_gauge_fun = is_data_larger
    }else if(definition_failure == 'not_enveloped'){
        single_gauge_fun = is_data_extreme
    }else{
        stop('unknown definition_failure')
    }

    # Alternative multi gauge summary -- return 1 if any gauge 'fails' by the above standards
    multi_gauge_any_failures=function(x){1*any(x!=0, na.rm=TRUE)}

    # Using the functions above, this call is looking at non-enveloping (not the median)
    stat_value = median_over_tide_gauges_of_fraction_models_below_data_per_event(
        event_id = bp_model$event,
        model_values_list = bp_model$x,
        data_values_list = bp_obs$x, 
        F_m_e_d_function=single_gauge_fun,
        multi_gauge_summary_function=multi_gauge_any_failures)

    print('Statistic: ')
    ord = order(as.numeric(gsub('[a-z]', '', names(stat_value))))
    print(stat_value[ord])
    print(paste0('Statistic sum: ', sum(stat_value, na.rm=TRUE)))

    # Store various outputs in a list
    output_store = list()
    output_store$stat_by_event = stat_value[ord]
    output_store$stat_sum = sum(stat_value, na.rm=TRUE)

    # Pool the real scenario and the model scenarios, then for each of these, compute the
    # event_statistic if it is swaped with the data 
    sampling_distribution_failure_at_any_gauge = compute_distribution_of_event_statistic_when_data_is_swaped_with_each_scenario(
        event_id = bp_model$event,
        model_values_list = bp_model$x,
        data_values_list = bp_obs$x, 
        include_no_swap_case = TRUE, # This means the distibution includes the 'real' value
        F_m_e_d_function=single_gauge_fun,
        multi_gauge_summary_function=multi_gauge_any_failures)

    simulate_random_event_failure_any_gauge<-function(){
        # For each event, pick a random scenario. Value will be 1 if any of
        # the gauges was not enveloped.
        failed_at_gauge = unlist(lapply(sampling_distribution_failure_at_any_gauge, function(x) sample(x, size=1)))
        # Note we can have events where every statistics has an NA values
        sum(failed_at_gauge, na.rm=TRUE)
    }

    print('Chance of failure at any gauge by event')
    results_by_event = unlist(lapply(sampling_distribution_failure_at_any_gauge, function(x) mean(x)))
    ord = order(as.numeric(gsub('[a-z]', '', names(results_by_event))))
    print(results_by_event[ord])

    output_store$prob_failure_any_event = results_by_event[ord]

    distribution_of_non_enveloped_events=replicate(N_sims, simulate_random_event_failure_any_gauge())
    # Get the probabilities
    fails = seq(0, length(results_by_event))
    N = length(fails)
    cum_prob_fails = sapply(fails, function(x) mean(distribution_of_non_enveloped_events <= x, na.rm=TRUE))
    print('Probability of N events failing:')

    prob_failure_df = data.frame(fails=fails, prob_fails=diff(c(0, cum_prob_fails)), prob_LE_num_fails=cum_prob_fails, 
        prob_GE_num_fails=1-c(0, cum_prob_fails[-N]))
    output_store$prob_N_events_failing = prob_failure_df
    print(prob_failure_df)

    print('Additional F^m_e and R^m_e statistics below (which do not involve the "definition_failure" above)')

    Rme_stat = compute_Rme(
        bp_model$event,
        bp_model$x,
        bp_obs$x, 
        return_stats=TRUE)
    print('F^m_e')
    print(Rme_stat$F_me_data)
    print('median(F^m_e)')
    print(median(Rme_stat$F_me_data, na.rm=TRUE))
    print('ADTest')
    print(Rme_stat$adtest)
    print('SD')
    print(Rme_stat$sd_R_me_data)
    if(Rme_stat$sd_test_fraction_of_sds_below < 0.025 |
       Rme_stat$sd_test_fraction_of_sds_below > 0.975){
        print('OUTSIDE 95% CONFIDENCE INTERVAL')
    }else{
        print('(inside 95% confidence interval)')
    }
    print('Fraction synthetic scenarios with smaller sd')
    print(Rme_stat$sd_test_fraction_of_sds_below)

    output_store$Rme_stat = Rme_stat

    return(invisible(output_store))
}

stats_for_paper<-function(){

    sink('stats_under_null_for_paper.txt')

    for(use_downsampled_model in c(FALSE, TRUE)){

        # Prepare storage
        slip_types = c('FAUS', 'VAUS', 'HS')
        results = vector(mode='list', length=length(slip_types))
        names(results) = slip_types 

        for(slip_type in slip_types){

            # Prepare storage
            stats_to_store = c('max_during_obs', 'stgrng_during_obs', 'max_8hrs_during_obs', 'stgrng_8hrs_during_obs', 'max_lastday_during_obs', 'stgrng_lastday_during_obs')
            results[[slip_type]] = vector(mode='list', length=length(stats_to_store))
            names(results[[slip_type]]) = stats_to_store

            for(test_variable in stats_to_store){

                # Prepare storage
                def_failures = c('under_estimate_or_equal', 'not_enveloped_or_equal')
                results[[slip_type]][[test_variable]] = vector(mode='list', length=length(def_failures))
                names(results[[slip_type]][[test_variable]]) = def_failures

                for(definition_failure in def_failures){
                    print("")
                    print("")
                    print('##############################################')
                    print(paste0('# ', slip_type, ' ', ifelse(use_downsampled_model, 'downsampled', 'regular'), ' ', test_variable, ' ', definition_failure))
                    print('##############################################')
                    stat = compute_expected_frequency_of_exceedances_under_null_hypothesis(
                        test_variable=test_variable, 
                        test_slip_type=slip_type, 
                        definition_failure=definition_failure,
                        use_downsampled_model=use_downsampled_model)

                    results[[slip_type]][[test_variable]][[definition_failure]] = stat
                }
            }
        }

        if(use_downsampled_model){
            saveRDS(results, 'stats_under_null_for_paper_downsampled.RDS')
        }else{
            saveRDS(results, 'stats_under_null_for_paper_regular.RDS')
        }
    }
    sink()
}

convert_stats_to_latex_table<-function(tsunami_size_type = 'tsunami_maxima', downsample=FALSE){

    if(tsunami_size_type == 'tsunami_maxima'){
        stat_starts_with = 'max_' # tsunami maxima   
    }else if(tsunami_size_type == 'stage_range'){
        stat_starts_with = 'stgrng_' # tsunami maxima   
    }else{
        stop('unknown tsunami_size_type')
    }

    #downsample = FALSE
    slip_types = c('FAUS', 'HS', 'VAUS')


    if(downsample){
        stats = readRDS('stats_under_null_for_paper_downsampled.RDS')
    }else{
        stats = readRDS('stats_under_null_for_paper_regular.RDS')
    }

    # Get the names of the stats (for either the tsunami maxima, or the stage range)
    k = grep(stat_starts_with, names(stats[[1]]))
    stat_types = names(stats[[1]])[k]

    all_df_col = list()

        for(stat_type in stat_types){
    for(slip_type in slip_types){

            getstat<-function(stats_mt_st){
                failed_to_envelope_sum = stats_mt_st$stat_sum
                df = stats_mt_st$prob_N_events_failing
                i0 = which(df[,1] == failed_to_envelope_sum)
                failed_to_envelope_sum_prob_gt = df$prob_GE_num_fails[i0]
                output = paste0(failed_to_envelope_sum, ' (', signif(failed_to_envelope_sum_prob_gt, 3), ')')
                return(output)
            }

            # Failure to envelope
            failed_to_envelope = getstat(stats[[slip_type]][[stat_type]]$not_enveloped_or_equal)

            # Understimate and failure to envelope
            low_fail_to_envelope = getstat(stats[[slip_type]][[stat_type]]$under_estimate_or_equal)

            # Median value of F_me (same in under_estimate_or_equal and not_enveloped_all_equal)
            stats_mtn_Rme = stats[[slip_type]][[stat_type]]$not_enveloped_or_equal$Rme_stat
            median_Fme = signif(median(stats_mtn_Rme$F_me_data, na.rm=TRUE), 3)
            median_Rme = signif(median(stats_mtn_Rme$R_me_data, na.rm=TRUE), 3)
            adGOF_p = signif(stats_mtn_Rme$adtest$p.value, 3)
            sd_frac_below = signif(stats_mtn_Rme$sd_test_fraction_of_sds_below, 3)

            if(grepl('8', stat_type)){
                time_window = 'early'
            }else if(grepl('lastday', stat_type)){
                time_window = 'late'
            }else{
                time_window = 'full'
            }

            df_col = data.frame( 
                time_window = time_window, slip_type = slip_type, failed_to_envelope=failed_to_envelope, 
                low_fail_to_envelope=low_fail_to_envelope, median_Fme=as.character(median_Fme), 
                #median_Rme=as.character(median_Rme), 
                adGOF_p = as.character(adGOF_p), sd_frac_below=as.character(sd_frac_below))

            all_df_col = c(all_df_col, list(df_col))
            
        }
    }

    library(xtable)
    final_table = do.call(rbind, all_df_col)

    output_table_filename = paste0('stats_under_null_table_', tsunami_size_type, ifelse(downsample, '_downsample', '_regular'), '.tex')
    output_caption = paste0('Tsunami size statistics based on the ', ifelse(tsunami_size_type == 'tsunami_maxima', 'tsunami maxima', 'stage range'), '.', 
        ifelse(downsample, ' The model was downsampled to match the tide gauge data sampling rate prior to computing the statistics.', ''))
    latex_table = xtable(final_table, caption=output_caption)
    print.xtable(latex_table, file=output_table_filename, include.rownames=FALSE)
    return(latex_table)
}
