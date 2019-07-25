#
# Given a set of events & their max-stage values at a point, as well as all
# Mw-vs-exceedance-rate logic-tree curves on the source-zone, here we calculate
# a family of stage-vs-exceedance-rate curves (one per logic-tree branch).
#

#'
#' Slow R-implementation of the fortran version that is wrapped below. Useful for regression testing
#' 
R_convert_Mw_vs_exceedance_rate_2_stage_vs_exceedance_rate<-function(
    logic_tree_rate_curve_Mw,
    logic_tree_rate_curve_exceedance_rate,
    unique_Mw,
    event_Mw, 
    event_conditional_probability,
    event_Mw_to_unique_Mw_match,
    event_max_stage_sorted_decreasing,
    event_max_stage_sorted_decreasing_inds,
    output_stages){

    max_mw_max = max(logic_tree_rate_curve_Mw)
    dMw = unique_Mw[2] - unique_Mw[1]

    # Get the 'incremental rate' for events with magnitude = unique_Mw
    f1 = approxfun(logic_tree_rate_curve_Mw, logic_tree_rate_curve_exceedance_rate, rule=2, ties='ordered')
    unique_rates = 
        (f1(unique_Mw - dMw/2) * (unique_Mw - dMw/2 <= max_mw_max) - 
         f1(unique_Mw + dMw/2) * (unique_Mw + dMw/2 <= max_mw_max) )
   
    # Partition the 'unique_rates' among all scenarios 
    scenario_rates = event_conditional_probability * unique_rates[event_Mw_to_unique_Mw_match]

    # Interpolate to get the stage exceedance rate
    sorted_stages_with_cap = c(.Machine$double.xmax, 
                               event_max_stage_sorted_decreasing[1]+1.0e-03, 
                               event_max_stage_sorted_decreasing)
    sorted_exrates_with_cap = cumsum( c(0, 0, scenario_rates[event_max_stage_sorted_decreasing_inds]) )
    stage_exceedance_rates = approx(sorted_stages_with_cap, sorted_exrates_with_cap, xout=output_stages, rule=2:1, ties=max)$y

    return(stage_exceedance_rates)

}

#SEXP make_stage_vs_rate_curve_site_c(SEXP rate_curve_Mw, SEXP rate_curve_exrate, SEXP N1,
#    SEXP dMw, SEXP unique_Mw, SEXP N2, SEXP event_Mw_to_unique_Mw_match, SEXP event_conditional_prob,
#    SEXP N3, SEXP sorted_event_stages, SEXP sorted_stage_inds, SEXP output_stages, SEXP output_stage_exrates,
#    SEXP N4){
fortran_convert_Mw_vs_exceedance_rate_2_stage_vs_exceedance_rate<-function(
    logic_tree_rate_curve_Mw,
    logic_tree_rate_curve_exceedance_rate,
    unique_Mw,
    event_Mw, 
    event_conditional_probability,
    event_Mw_to_unique_Mw_match,
    event_max_stage_sorted_decreasing,
    event_max_stage_sorted_decreasing_inds,
    output_stages){
    
    dMw = as.double(unique_Mw[2] - unique_Mw[1])

    output_exrates = output_stages * 0.0

    N1 = as.integer(length(logic_tree_rate_curve_Mw))
    N2 = as.integer(length(unique_Mw))
    N3 = as.integer(length(event_Mw_to_unique_Mw_match))
    N4 = as.integer(length(output_stages))

    .Call('make_stage_vs_rate_curve_site_c',
          logic_tree_rate_curve_Mw, logic_tree_rate_curve_exceedance_rate, N1,
          dMw, unique_Mw, N2, event_Mw_to_unique_Mw_match, 
          event_conditional_probability, N3, event_max_stage_sorted_decreasing, 
          event_max_stage_sorted_decreasing_inds, output_stages, 
          output_exrates, # This is the variable that is updated
          N4)

    return(output_exrates)

}

#' Get a set of stage-vs-exceedance-rate curves (one for each magnitude-exceedance rate logic-tree branch)
#'
#' Suppose we have a source zone containing a set of Mw-vs-exceedance rate curves, and a corresponding set
#' of events, each of which lead to a given max-stage at a particular site. This routine computes a
#' stage-vs-exceedance-rate curve for each Mw-vs-exceedance rate curve.
#'
#' @param logic_tree_rate_curve_Mw vector with magnitudes at which the logic-tree exceedance-rate curves are tabulated
#' @param logic_tree_rate_curves_exceedance_rate matrix with one row per logic-tree branch, and one column 
#' per logic_tree_rate_curve_Mw, giving the exceedance rate at each magnitude.
#' @param event_Mw vector of magnitudes for each event. These should be drawn from an evenly spaced set (e.g. 7.2, 7.3, ... 9.8)
#' @param event_conditional_probability vector with teh same length as event_Mw giving the conditional probability that 
#' the event occurs, given that an event with the same Mw occurred.
#' @param event_max_stage The max stage for each event at the site of interest
#' @param output_stages the stages at which we tabulate outputs
#' @param use_Fortran if TRUE do the heavy computations using compiled Fortran code. Otherwise do them in R (slower)
#' @return A matrix with one row for each logic-tree branch, and one column for each output-stage, giving the exceedance-rates
#'
#' @export
#'
#' @examples
#' 
#' logic_tree_rate_curve_Mw = seq(7.2, 9.8, by=0.1)
#' # 20 logic-tree branches
#' NB = 20
#' logic_tree_rate_curves_exceedance_rate = matrix(0, nrow=NB, ncol=length(logic_tree_rate_curve_Mw))
#' bs = seq(0.7, 1.3, len=NB)
#' for(i in 1:NB){
#'     logic_tree_rate_curves_exceedance_rate[i,] = 10**(-bs[i]*(logic_tree_rate_curve_Mw - 7))
#' }
#' 
#' # Set of events -- 30 per magnitude
#' N_per_mag = 30
#' NE = N_per_mag*length(logic_tree_rate_curve_Mw)
#' event_Mw = rep(logic_tree_rate_curve_Mw, times=N_per_mag)
#' event_conditional_probability = rep(1/N_per_mag, length(event_Mw))
#' # Make up some max-stage values
#' event_max_stage = runif(length(event_Mw)) - 0.5 + 0.5*(event_Mw - 6)
#' N_output_stages = 100
#' output_stages = 10**seq(-2, 1.5, len=N_output_stages)
#' 
#' # Typical version
#' results = convert_Mw_vs_exceedance_rates_2_stage_vs_exceedance_rates(
#'     logic_tree_rate_curve_Mw,
#'     logic_tree_rate_curves_exceedance_rate,
#'     event_Mw,
#'     event_conditional_probability,
#'     event_max_stage,
#'     output_stages)
#' 
#' ## This is in pure R code
#' results_R = convert_Mw_vs_exceedance_rates_2_stage_vs_exceedance_rates(
#'     logic_tree_rate_curve_Mw,
#'     logic_tree_rate_curves_exceedance_rate,
#'     event_Mw,
#'     event_conditional_probability,
#'     event_max_stage,
#'     output_stages,
#'     use_Fortran=FALSE)
#' 
#' # They should be 'the same' up to floating point reordering
#' ## Small absolute error
#' stopifnot(all(abs(results - results_R) < 1.0e-15))
#' ## Small relative error (1e-20 zero divide protection)
#' stopifnot(all( ((results + 1e-20)/(results_R+1e-20) - 1) < 1.0e-12))
#'
convert_Mw_vs_exceedance_rates_2_stage_vs_exceedance_rates<-function(
    logic_tree_rate_curve_Mw,
    logic_tree_rate_curves_exceedance_rate,
    event_Mw,
    event_conditional_probability,
    event_max_stage,
    output_stages,
    use_Fortran=TRUE){

    #
    # Check inputs.
    #
    stopifnot(length(logic_tree_rate_curve_Mw) == ncol(logic_tree_rate_curves_exceedance_rate))
    stopifnot(length(event_Mw) == length(event_conditional_probability))
    stopifnot(length(event_Mw) == length(event_max_stage))

    #
    # Get the unique magnitudes, and ensure they satisfy our assumptions
    #
    unique_Mw = sort(unique(event_Mw))
    dMw = unique_Mw[2] - unique_Mw[1]
    # Check that unique_Mw is evenly spaced 
    diff_unique_Mw = diff(unique_Mw)
    tol = 1.0e-06
    if(any(abs(diff_unique_Mw - dMw) > tol*dMw)){
        msg = paste0('event_Mw must contain evenly spaced magnitudes, and evenly spaced up to one part in ', tol, '\n',
                     'This can be broken (e.g. due to finite precision file storage), in which case event_Mw should be rounded')
        stop(msg)
    }

    # Match event_Mw to unique_Mw
    event_Mw_to_unique_Mw_match = match(event_Mw, unique_Mw)

    # Sort the event_stages
    sorted_max_stage = sort(event_max_stage, decreasing=TRUE, index.return=TRUE)

    event_max_stage_sorted_decreasing = as.double(sorted_max_stage$x)
    event_max_stage_sorted_decreasing_inds = as.integer(sorted_max_stage$ix)

    # This will hold the exceedance-rates 
    output_stage_exrates = matrix(0.0, ncol = length(output_stages), nrow=dim(logic_tree_rate_curves_exceedance_rate)[1])

    if(use_Fortran){
        # Get the stage exceedance-rates for each Mw-frequency curve
        for(i in 1:(nrow(output_stage_exrates))){
            output_stage_exrates[i,] = fortran_convert_Mw_vs_exceedance_rate_2_stage_vs_exceedance_rate(
                logic_tree_rate_curve_Mw,
                logic_tree_rate_curves_exceedance_rate[i,],
                unique_Mw,
                event_Mw,
                event_conditional_probability,
                event_Mw_to_unique_Mw_match,
                event_max_stage_sorted_decreasing,
                event_max_stage_sorted_decreasing_inds,
                output_stages)
        }
    }else{
        # This is slower but offers some regression-test
        # Get the stage exceedance-rates for each Mw-frequency curve
        for(i in 1:(nrow(output_stage_exrates))){
            output_stage_exrates[i,] = R_convert_Mw_vs_exceedance_rate_2_stage_vs_exceedance_rate(
                logic_tree_rate_curve_Mw,
                logic_tree_rate_curves_exceedance_rate[i,],
                unique_Mw,
                event_Mw,
                event_conditional_probability,
                event_Mw_to_unique_Mw_match,
                event_max_stage_sorted_decreasing,
                event_max_stage_sorted_decreasing_inds,
                output_stages)
        }

    }

    return(output_stage_exrates)

}


#' Get a percentile from an empirical distribution
#'
#' Suppose a distribution is defined empirically by a set of values (vals), each
#' having a weight 'weight' which defines the probability mass function of each vals.
#' This function returns the smallest value of X in 'vals' such that PR(vals <= X) >= p.
#'
#' @param vals vector of values
#' @param weights vector of non-negative weights (one for each value). These will be
#' normalised to weights/sum(weights) inside the function
#' @param p vector, with length 1 or more, containing values for 'p' as defined above. We require 0<=p<=1
#' @return A vector 'X' with the same length as 'p', as defined above
#'
#' @examples
#'    vals = c(-10, 5, 3, -4, 5, 6)
#'    weights = c(0.5, 0.1, 0.1, 0.1, 0.1, 0.1)
#'
#'    test_p = c(0.0, 0.5, 0.501, 0.6, 0.601, 0.7, 0.701, 0.8, 0.9, 0.901, 1.0)
#'    empirical_q = weighted_percentile(vals, weights, test_p)
#'
#'    expected_vals = c(-10, -10, -4, -4, 3, 3, 5, 5, 5, 6, 6)
#'    stopifnot(all(empirical_q == expected_vals))
#'
#' @export
#'
weighted_percentile<-function(vals, weights, p){

    stopifnot( all(p >= 0 & p <= 1) )
    stopifnot( all(weights >= 0) )

    weights = weights/sum(weights)
    sorted_vals = sort(vals, index.return=TRUE)
    sorted_weights = weights[sorted_vals$ix]
    cum_sorted_weights = cumsum(sorted_weights)

    # Look up multiple indices at once
    ind = p * 0
    for(i in 1:length(p)){
        ind[i] = sum(cum_sorted_weights < p[i]) + 1
    }
    return(sorted_vals$x[ind])
}
