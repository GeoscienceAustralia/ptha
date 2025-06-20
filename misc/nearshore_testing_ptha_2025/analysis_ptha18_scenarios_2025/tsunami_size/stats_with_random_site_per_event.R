#' Sample a single tide gauge randomly from each event, compute the fraction of
#' model scenarios less than the data at that tide gauge, and return either the
#' random fractions for each event, or summary statistics (mean, standard
#' deviation, number of non-missing values)
#'
#' @param event_id event_name (or integer), with one entry for every site that
#' has tide gauge data. Typically there are multiple tide gauges per event, so
#' event_id will have repeated values (see example below). 
#' @param model_values_list a list with length equal to length(event_id), where
#' each entry is a vector of stochastic model values for that particular tide
#' gauge and event (see example below).
#' @param data_values_list a list with length equal to length(event_id), where
#' each entry is a vector of data values for that particular tide
#' gauge and event (see example below). 
#' @param return_summary_stats if TRUE return the mean,sd, and sum(!is.na()) of the 
#' random fractions of model values below the data. Otherwise return the random
#' fractions of model values below the data (i.e. one value per unique event_id)
#' 
#' @example
#' # See the test functions in this file
#' 
#' # Suppose we have 3 events with
#' #   - event 1: 2 tide gauges, 4 scenarios
#' #   - event 2: 1 tide gauge, 3 scenarios
#' #   - event 3: 2 tide gauges, 4 scenarios
#' # The input data could be
#' event_id = c('e1', 'e1', 'e2', 'e3', 'e3')
#' # The model data could be (
#' model_values_list = list(
#'     c(1.3, 1.2, 1.5, 1.4), # event 1, site 1, 4 scenarios
#'     c(0.5, 0.7, 0.1, 0.2), # event 1, site 2, 4 scenarios
#'     c(2.2, 2.7, 2.1     ), # event 2, site 1, 3 scenarios
#'     c(3.3, 3.2, 3.5, 3.4), # event 3, site 1, 4 scenarios
#'     c(4.5, 4.7, 4.1, 4.2)) # event 3, site 2, 4 scenarios
#' # Usually the data will be constant for each model scenario, but in principle it could vary
#' # (e.g. for statistics that rely on a property of the model scenario, such as its arrival time)
#' data_values_list = list(
#'     c(1.25, 1.25, 1.25, 1.25), # event 1, site 1, 4 scenarios
#'     c( 0.6,  0.6,  0.6,  0.6), # event 1, site 2, 4 scenarios
#'     c( 2.8,  2.8,  2.8      ), # event 2, site 1, 3 scenarios
#'     c( 3.1,  3.1,  3.1,  3.1), # event 3, site 1, 4 scenarios
#'     c( 4.4,  4.4,  4.4,  4.4)) # event 3, site 2, 4 scenarios
#' mean_and_sd_of_fraction_below_at_one_random_site_per_event(event_id, model_values_list, data_values_list)
#'
mean_and_sd_of_fraction_below_at_one_random_site_per_event<-function(
    event_id, 
    model_values_list, 
    data_values_list, 
    return_summary_stats=TRUE){

    # For each event, select one tide gauge at random.
    # Store the fraction of model values below the data value at that tide gauge
    unique_event_ids = unique(event_id)
    frac_values = rep(NA, length(unique_event_ids))
    data_is_NA = unlist(lapply(data_values_list, function(x) any(is.na(x))))

    # Check the entries make sense
    ndata = unlist(lapply(data_values_list, length))
    nmodels = unlist(lapply(model_values_list, length))
    if( (length(ndata) != length(nmodels)) | 
        (any((ndata > 1) & (ndata != nmodels))) ){
        stop('For each tide gauge, the data_values_list must either have a single entry, or one entry per model scenario')
    }

    for(i in 1:length(unique_event_ids)){
         # Find all tide gauges with non-missing data for the i'th event ID 
        possible_inds = which(event_id == unique_event_ids[i] & !data_is_NA)
        
        if(length(possible_inds) == 0){
            # Issues for some statistics when the data finishes before the model
            frac_values[i] = NA
        }else{
            if(length(possible_inds) == 1){
                # There is only 1 tide gauge that can be sampled.
                # Avoid calling sample with length(possible_inds)==1, as it can get
                # confused, see ?sample
                ki = possible_inds
            }else if(length(possible_inds) > 1){
                # Randomly choose the fraction at one of the tide gauges
                # This is the typical case
                ki = sample(possible_inds, size=1)
            } 
            frac_values[i] = mean(model_values_list[[ki]] < data_values_list[[ki]])
        }
    }

    if(!return_summary_stats){
        return(frac_values)
    }else{
        output = c(mean(frac_values, na.rm=TRUE), sd(frac_values, na.rm=TRUE), sum(!is.na(frac_values)))
        return(output)
    }
}

# Create synthetic data matching the structure of the input data
# with the same event_id, number of gauges, number of synthetic scenarios.
# The event_id remains the same, but model and data values are just draws from
# a uniform distribution, unless the data is NA (in which case the random data is NA too)
simulate_uniform_random_data_and_models_at_tide_gauges<-function(
    event_id, 
    model_values_list, 
    data_values_list){

    # Replace the model and data with synthetic results, 
    # assuming model and data come from the same distribution (runif)
    for(i in 1:length(model_values_list)){
        # Each 'i' corresponds to a single tide gauge

        N_model = length(model_values_list[[i]])
        N_data = length(data_values_list[[i]])
        stopifnot(N_data == 1 | N_data == N_model)

        model_values_list[[i]] = runif(N_model)
        if( any(is.na(data_values_list[[i]])) ){
            data_values_list[[i]] = rep(NA, N_data)
        }else{
            # Single data value at the tide gauge, repeated for each model
            data_values_list[[i]] = rep(runif(1), N_data)
        }
    }
    return(list(event_id = event_id, model_values_list=model_values_list, data_values_list=data_values_list))
}

simulate_random_data_by_sampling_model_at_tide_gauges<-function(
    event_id, 
    model_values_list, 
    data_values_list,
    jitter_scale=0){
    # Replace the data with a random model scenario (the same scenario for each event)

    N_model_all = unlist(lapply(model_values_list, length))
    N_data_all  = unlist(lapply(data_values_list, length))
    unique_event_id = unique(event_id)

    for(i in 1:length(unique_event_id)){
        # Find rows matching this event
        k = which(event_id == unique_event_id[i])

        # Should be the same number of model scenarios for each event
        stopifnot(all(N_model_all[k] == N_model_all[k[1]]))

        # Each tide gauge should either have a single data value, or one per scenario.
        stopifnot(all(N_data_all[k] == 1 | N_data_all[k] == N_model_all[k]))

        # Should always have > 1 model scenario
        stopifnot(all(N_model_all[k] > 1))

        # Choose 1 model scenario at random
        ind = sample(1:N_model_all[k[1]], size=1)

        # Populate the data with the random scenario, unless we have missing data
        for(j in k){
            if(!any(is.na(data_values_list[[j]]))){
                data_values_list[[j]] = rep(model_values_list[[j]][ind] + jitter_scale*rnorm(1), N_data_all[j])
            }
        }
    }
    return(list(event_id = event_id, model_values_list=model_values_list, data_values_list=data_values_list))
}

.basic_example_mean_and_sd_of_fraction_below_at_one_random_site_per_event<-function(){
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
    mean_and_sd_of_fraction_below_at_one_random_site_per_event(event_id, model_values_list, data_values_list)
    mean_and_sd_of_fraction_below_at_one_random_site_per_event(event_id, model_values_list, data_values_list, return_fraction_below=TRUE)
}

.test_mean_and_sd_of_fraction_below_at_one_random_site_per_event<-function(){
    # Suppose we have 3 events with
    #   - event 1: 2 tide gauges, 5 scenarios
    #   - event 2: 1 tide gauge , 2 scenarios
    #   - event 3: 2 tide gauges, 7 scenarios

    event_id   = c('e1', 'e1', 'e2', 'e3', 'e3')
    nscenarios = c(5   ,    5,    2,    7,    7) 

    # For each event and gauge, let the model be a random variable
    simulate_model_or_data<-function(nscenarios){
        random_mu = c(1, 2, 3, 4, 5) # Mean of a normal distribution
        output = vector(mode='list', length=length(random_mu))
        for(i in 1:length(output)){
            # Random sample, different distribution for each event/gauge pair
            output[[i]] = rnorm(nscenarios[i], mean=random_mu[i])
        }
        return(output)
    }

    # Make some random model values, and random data from the same distribution, and apply the function
    compute_once<-function(){
        model_values_list = simulate_model_or_data(nscenarios)
        data_values_list  = simulate_model_or_data((nscenarios*0 + 1))
        mean_and_sd_of_fraction_below_at_one_random_site_per_event(event_id, model_values_list, data_values_list)
    }

    # Used to compute the theoretical distribution under the null hypothesis
    compute_once_null_hypothesis<-function(){
        model_values_list = simulate_model_or_data(nscenarios)
        data_values_list  = simulate_model_or_data((nscenarios*0 + 1))
        fake_data = simulate_uniform_random_data_and_models_at_tide_gauges(
            event_id,
            model_values_list,
            data_values_list)
        mean_and_sd_of_fraction_below_at_one_random_site_per_event(fake_data$event_id, fake_data$model_values_list, fake_data$data_values_list)
    }

    bigN = 10000
    many_solutions = replicate(bigN, compute_once())
    # True mean should be 0.5
    true_mean = 0.5
    unif_dist_variance = 1/12
    NUNIQUE = length(unique(event_id))
    true_variance_of_mean = unif_dist_variance * 1/NUNIQUE * 1/bigN 
    tol_as_standard_errors = 5 # unlikely to fail
    if(abs(mean(many_solutions[1,]) - true_mean) < sqrt(true_variance_of_mean)*tol_as_standard_errors){
        print('PASS')
    }else{
        print('FAIL (or very unlikely stochastic event)')
    }

    # Theoretically the distribution should be pretty close to NUNIQUE samples from a uniform distribution
    comparison_fun<-function(){
        data = runif(NUNIQUE)
        return(c(mean(data), sd(data), length(data)))
    }
    comparison_data = replicate(bigN, comparison_fun())

    ## Theoretically this should be EXACTLY the distribution
    #better_comparison_fun<-function(){
    #    data = rep(NA, NUNIQUE)
    #    N = nscenarios[1]
    #    for(i in 1:NUNIQUE){
    #        data[i] = mean(runif(N) < runif(1))
    #    }
    #    return(c(mean(data), sd(data), length(data)))
    #}
    #better_comparison_data = replicate(bigN, better_comparison_fun())

    # This should also be the same as above
    better_comparison_data = replicate(bigN, compute_once_null_hypothesis())

    par(mfrow=c(2,2))
    qqplot(many_solutions[1,], better_comparison_data[1,])
    abline(0, 1, col='red')
    grid()
    title('many_solutions mean vs better_comparison_data mean \n (theoretically drawn from the same distribution)')

    qqplot(many_solutions[2,], better_comparison_data[2,])
    abline(0, 1, col='red')
    grid()
    title('many_solutions sd vs better_comparison_data sd \n (theoretically drawn from the same distribution)')

    qqplot(many_solutions[1,], comparison_data[1,])
    abline(0, 1, col='red')
    grid()
    title('many_solutions mean vs comparison_data mean \n (a slightly different distribution unless nscenarios --> Infty)')

    qqplot(many_solutions[2,], comparison_data[2,])
    abline(0, 1, col='red')
    grid()
    title('many_solutions sd vs comparison_data sd \n (a slightly different distribution unless nscenarios --> Infty)')

    ks1 = suppressWarnings( ks.test(many_solutions[2,], better_comparison_data[2,]) )# Same distribution
    ks2 = suppressWarnings( ks.test(many_solutions[2,], comparison_data[2,]) )# Different (but pretty similar)

    if(ks1$p.value > 1e-06){
        print('PASS')
    }else{
        print('FAIL (or very unlikely stochastic event)')
    }
     
}
