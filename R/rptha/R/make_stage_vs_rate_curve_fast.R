#
# Given a set of events & their max-stage values at a point, as well as all
# Mw-vs-exceedance-rate logic-tree curves on the source-zone, here we calculate
# a family of stage-vs-exceedance-rate curves (one per logic-tree branch).
#

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
#' @return A matrix with one row for each logic-tree branch, and one column for each output-stage, giving the exceedance-rates
#'
#' @export
convert_Mw_vs_exceedance_rates_2_stage_vs_exceedance_rates<-function(
    logic_tree_rate_curve_Mw,
    logic_tree_rate_curves_exceedance_rate,
    event_Mw,
    event_conditional_probability,
    event_max_stage,
    output_stages){

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

    return(output_stage_exrates)

}
