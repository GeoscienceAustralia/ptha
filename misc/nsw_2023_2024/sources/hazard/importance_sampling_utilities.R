#' Compute basic importance sampling weights, given the event rate and sampling probabilities
basic_importance_sampling_weights<-function(event_rate, sampling_prob){

    stopifnot(isTRUE(all.equal(sum(sampling_prob), 1)))

    TOTAL_RATE = sum(event_rate)
    # Basic importance sampling weights
    bis_wts = (event_rate/TOTAL_RATE) / sampling_prob
    bis_wts[sampling_prob == 0] = 0
    return(bis_wts)
}

#' Analytical exceedance rate and variance of Monte Carlo exceedance rates
#'
#' Given a site in the offshore PTHA (where all event_stage values are known), compute
#' the rate of exceedance of some threshold, and the variance of the Monte Carlo estimate.
#'
#' @param event_rate vector with occurrence rates for each scenario
#' @param event_stage vector with peak stage (or some other intensity measure) for each scenario
#' @param sampling_prob vector with probabilities used for non-uniform sampling of scenarios with replacement
#' @param bis_wts basic importance sampling weights
#' @param stage_threshold the threshold (we compute the rate with which event_stage>stage_threshold)
#' @param N_MC Number of monte carlo samples that would be taken (only used for MC_variance calculation)
#' @return Vector of length 2 with the (exact) exceedance rate, and variance of
#' Monte Carlo exceedance rates with N_MC samples.
analytical_exrate_and_MC_variance<-function(event_rate, event_stage, stage_threshold, sampling_prob, N_MC){

    # The following use ALL PTHA scenarios -- they are exact results, not Monte Carlo estimates
    TOTAL_RATE = sum(event_rate)
    prob_stage_gt_threshold = sum(event_rate * (event_stage > stage_threshold))/sum(event_rate)
    output_exrate = TOTAL_RATE * prob_stage_gt_threshold

    bis_wts = basic_importance_sampling_weights(event_rate, sampling_prob)

    # Analytical variance of Monte Carlo sampling with N_MC scenarios and the
    # given basic-importance-sampling weights and sample weights.
    output_var = (TOTAL_RATE**2 / N_MC)  * 
        sum((bis_wts*(event_stage > stage_threshold) - prob_stage_gt_threshold)**2 * sampling_prob)

    return(c(output_exrate, output_var))
}

sample_scenarios_monte_carlo<-function(event_rate, event_stage, stage_thresholds, sampling_prob, N_MC, myseed=NULL){

    # Option to set the random seed
    if(is.null(myseed)){
        # Set myseed to .Random.seed (initialise the latter if needed)
        myseed = rptha::get_random_seed()
    }else{
       .Random.seed = myseed
    }

    # Use Monte Carlo sampling with N_MC scenarios to estimate exceedance-rates at stage_thresholds
    si = sample(seq(1, length(sampling_prob)), size=N_MC, replace=TRUE, prob=sampling_prob)

    TOTAL_RATE = sum(event_rate)
    bis_wts = basic_importance_sampling_weights(event_rate, sampling_prob)

    sampled_stages = event_stage[si]
    sampled_bis_wts = bis_wts[si]

    # Empirical exceedance-rates with basic important sampling
    exrates = sapply(stage_thresholds, 
        # Monte Carlo exceedance-rate estimate
        function(x) TOTAL_RATE*sum(sampled_bis_wts * (sampled_stages > x)/N_MC), USE.NAMES=FALSE)

    output = list(si = si, sampled_stages = sampled_stages, sampled_bis_wts = sampled_bis_wts, 
        stage_thresholds = stage_thresholds, exrates = exrates, myseed=myseed, TOTAL_RATE=TOTAL_RATE)

    return(output)
}

#' Estimate exceedance-rates and their monte carlo variance from sampled scenarios
#'
#' Here the input vectors refer ONLY to the sampled scenarios (except for the
#' one relating to the thresholds of interest), which were sampled using basic
#' importance sampling.
#'
#' @param rates_sampled_scenarios. The rate r_i(e) of each of the sampled
#' scenarios (according to rate model r_i)
#' @param sampling_prob_sampled_scenarios. The chance of originally sampling
#' each of the sampled scenarios
#' @param stage_sampled_scenarios. The maximum stage (or other intensity
#' measure) for the sampled scenarios
#' @param desired_stage_thresholds. A vector containing one or more stage
#' thresholds, for which we will estimate exceedance-rates and their uncertainties.
#' @return A data.frame with the estimated exceedance-rate, and variance of the
#' exceedance rate, for all desired_stage_thresholds
estimate_exrate_and_variance_sampled_scenarios<-function(
    rates_sampled_scenarios,
    sampling_prob_sampled_scenarios,
    stage_sampled_scenarios,
    desired_stage_thresholds){

    N = length(rates_sampled_scenarios) # Number sampled scenarios
    stopifnot(N == length(stage_sampled_scenarios))
    stopifnot(N == length(sampling_prob_sampled_scenarios))

    NT = length(desired_stage_thresholds) # Number thresholds

    output = data.frame(
        exrate = rep(NA, length=NT),
        exrate_var = rep(NA, length=NT),
        threshold = desired_stage_thresholds)

    for(i in 1:NT){
        # Estimate the exceedance-rate and its variance for threshold i
        ti = desired_stage_thresholds[i]
        exceeded = (stage_sampled_scenarios > ti)
        r_on_w = rates_sampled_scenarios/sampling_prob_sampled_scenarios
        # Estimate of the exceedance-rate
        exrate_est = 1/N * sum(r_on_w*exceeded)
        # Estimate of the exceedance-rate variance
        exrate_var_est = 1/N * sum( (r_on_w*exceeded - exrate_est )**2 / N)

        output$exrate[i] = exrate_est
        output$exrate_var[i] = exrate_var_est
    }

    return(output)
}

