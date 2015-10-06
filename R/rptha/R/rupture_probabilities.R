#' Given Mw, compute parameters which characterise the 'probability distribution' of slip
#'
#' FIXME: This function is not currently exported -- do we want to keep it?
#' Seems a good idea.
#' We assume mu is constant, and log10(Area) has a normal distribution. Since 
#' the seismic moment = mu*slip*Area, this implies that log10(Slip) has a normal
#' distribution as well. So the mean and standard deviation of the log10 variable
#' completely describe it. As well as returning those values, we return the mean slip
#' (mean of untransformed slip) since its calculation is slightly inconvenient.
#'
#' @param Mw the moment magnitude
#' @param mu Shear Modulus (Pascals)
#' @param constant Value of constant passed to \code{M0_2_Mw}
#' @return A list containing 1) The log10(Slip) mean and stanard deviation and
#' 2) The associated mean slip (mean of the untransformed slip)
#'
compute_slip_density_parameters<-function(Mw, mu=3e+10, constant=9.05){

    M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

    size_scaling = Mw_2_rupture_size(Mw, detailed=TRUE, CI_sd=1)

    # km
    mu_log10_area = log10(size_scaling$values['area'])
    sd_log10_area = size_scaling$log10_sigmas['area']

    # To multiply area by 1e+06 (i.e. km^2 --> m^2) we just add log10(1e+06) to mu_log10_area
    mu_log10_area = mu_log10_area + log10(1.0e+06)

    mu_log10_slip = (-mu_log10_area + log10(M0) - log10(mu))
    sd_log10_slip = sd_log10_area

    log10_slip_par = c(mu_log10_slip, sd_log10_slip)
    names(log10_slip_par) = c('log10_mean', 'log10_sd')

    # If log10(slip) is normally distributed, then this is the mean_slip for
    # ruptures of size Mw
    mean_slip = exp( mu_log10_slip*log(10) + 0.5*( (sd_log10_slip*log(10))**2))
    # [NOTE: To derive this formula, consider a variable y ~ exp(x) where x is
    # normal with given mean and sigma. We know mean(y) = exp(mean + sigma^2/2).
    # Now consider x' = log10(y). We know this has mean' and sigma' which are the
    # same as x **but divided by ln(10)**
    # We know mean(y) = exp(mean + sigma^2/s), but obviously this is the same as
    #                   exp([mean' x ln(10)] + [sigma' x ln(10)]^2 /2)
    # Hence if y = 10**(x'), then mean(y) can be computed in this way

    return(list(log10_slip_par = log10_slip_par, mean_slip = mean_slip))
}

#' Assign conditional probabilities to earthquakes with the same Mw
#'
#' Consider all the earthquake events on a discrete source with magnitude Mw.
#' We wish to assign conditional probabilities to these events (conditional on
#' the fact that an event of size Mw occurred). This allows us to compute the rate
#' of any individual event as 'event_conditional_probability * rate_of_events_of_size_Mw' \cr
#' One way to do this is assume that for a given magnitude, all event
#' conditional probabilities are equal ('all_equal').
#' This might be reasonable if the events are identical in area and slip, and
#' just differ in location. \cr 
#' Another approach is to assume the conditional
#' probabilities vary inversely with the average slip of the event ('inverse_slip'). This might
#' be reasonable if the individual events are almost the same, but there are slight
#' differences in the event area and slip, and we want each event to make the same 
#' contribution to the long-term average slip rate. 
#'
#' @param event_table A data.frame containing information on all earthquake
#' events on the source-zone. It's the output of
#' \code{get_all_earthquake_events}
#' @param conditional_probability_model character. Supported values are
#' 'all_equal' and 'inverse_slip', corresponding to cases discussed above.
#' @return A numeric vector with length = nrow(event_table) giving the conditional probabilities. 
#' @export
get_event_probabilities_conditional_on_Mw<-function(
    event_table,
    conditional_probability_model = 'all_equal'){

    # Define functions corresponding to supported conditional_probability models
    # These take as input a subset of the event table containing events with identical Mw
    if(conditional_probability_model == 'all_equal'){
        pfun<-function(event_table_Mw){
            N = length(event_table_Mw[,1])
            return( rep(1, N)/N)
        }

    }else if(conditional_probability_model == 'inverse_slip'){
        pfun<-function(event_table_Mw){
            slip_inv = 1/event_table_Mw$slip
            return(slip_inv/sum(slip_inv))
        }

    }else{
        stop('conditional_probability_model value not recognized')
    }

    # Find unique Mw values in event table
    event_table_Mws = sort(unique(event_table$Mw))

    # Compute conditional_probabilities
    event_conditional_probabilities = rep(NA, length(event_table[,1]))
    for(Mw in event_table_Mws){
        ii = which(event_table$Mw == Mw)
        event_conditional_probabilities[ii] = pfun(event_table[ii,])

        # Ensure they sum to 1.0
        stopifnot(isTRUE(all.equal(sum(event_conditional_probabilities[ii]), 1.0)))
    }

    # Final check
    if(any(is.na(event_conditional_probabilities))){
        warnings('NA values among event_conditional_probabilities')
    }

    return(event_conditional_probabilities)
}


#' Get a function to compute the rate of earthquakes with magnitude > Mw
#'
#' The input parameters are based on a truncated Gutenberg Richter model: \cr
#' The parameter 'a' is adjusted to account for constant rescaling, and here it
#' is computed to match the long-term seismogenic slip rate (m/year). Details
#' are explained below. All other parameters are provided by the user. \cr
#' To account for uncertainties in those parameters, the user may specify each
#' as a vector of one-or-more parameters with associated probabilities (which of
#' course must sum to 1). The code then treats all combinations of those
#' parameters as having probability equal to the product of the individual
#' parameter probabilities. It then integrates the associated Gutenberg Richter
#' models to produce a single Rate function, which is returned. \cr
#' Next we clarify how the 'a' parameter is constrained to match the long-term seismogenic slip
#' rate. Given particular values of all the pararameters except a, and a long-term seisomogenic slip
#' rate (= long_term_slip_rate*coupling_coefficient), the long-term moment rate
#' should equal the rate integrated from the individual events: \cr
#' \deqn{ \mu  long_term_seismogenic_slip_rate  sourcezone_area = 
#' \sum_{i \in events} \big{(} \mu event_slip_{i} event_area_{i} Pr(event_{i} | event_of_size_Mw = Mw_{i}) rate_of_events_of_size_Mw_{i}\big{)} }\cr
#' where Pr(event_{i} | event_of_size_Mw_{i}) gives the conditional probability of each
#' event with the same Mw in the event table [i.e. if we take all events with a
#' given Mw from the table, their conditional probabilities should sum to 1]. 
#' The meaning of other variables should be clear. 
#' Note that a factor of 10^(a) will appear in the RHS term
#' rate_of_events_of_size_Mw_{i}, and nowhere else. Hence to compute a, we just
#' compute the LHS and RHS of the above equations, tentatively assuming a=0 in the latter
#' case, and then it follows that a = log10(LHS/RHS). 
#' 
#'
#' @param slip_rate numeric vector with 1 or more long-term seismogenic slip rates for
#' the source-zone (m/year). This gives the long-term average slip on the fault
#' due to seismic events (creep is ignored, so these rates must account for the
#' rate of coupling)
#' @param slip_rate_prob numeric vector with probabilities for the slip_rate parameters
#' @param b numeric vector with 1 or more 'b' parameters
#' @param b_prob numeric vector with probabilities for the b parameters
#' @param Mw_min numeric vector with 1 or more 'Mw_min' parameters
#' @param Mw_min_prob numeric vector with probabilities for the Mw_min parameters
#' @param Mw_max numeric vector with 1 or more 'Mw_max' parameters
#' @param Mw_max_prob numeric vector with probabilities for the Mw_max parameters
#' @param sourcezone_total_area The total area of the sourcezone (km^2)
#' @param event_table A data.frame containing information on all events in the
#' source-zone. Output of \code{get_all_earthquake_events}
#' @param event_conditional_probabilities A vector of conditional probabilities for
#' all events in the source-zone (conditional on the fact than an event with
#' their magnitude did occur). Output of
#' \code{get_event_probabilities_conditional_on_Mw}
#' @param computational_increment For each final branch of the logic tree, the
#' rate function with those parameters is computed at points in a sequence from
#' min(Mw_min) to max(Mw_max), with spacing computational_increment (adjusted if
#' required so that the extremes are included). The function value is then
#' evaluated via lookup. 
#' @return function of a numeric vector Mw, which returns the rate of events
#' with magnitude > Mw
#' @export
rate_of_earthquakes_greater_than_Mw_function<-function(
    slip_rate, slip_rate_prob, b, b_prob, 
    Mw_min, Mw_min_prob, Mw_max, Mw_max_prob, 
    sourcezone_total_area, event_table, 
    event_conditional_probabilities,
    computational_increment=0.01){

    # Check data
    stopifnot(length(slip_rate) == length(slip_rate_prob))
    stopifnot(length(b) == length(b_prob))
    stopifnot(length(Mw_min) == length(Mw_min_prob))
    stopifnot(length(Mw_max) == length(Mw_max_prob))

    # Check (conservative) on the input Mw_min/Mw_max
    stopifnot(min(Mw_max) > max(Mw_min))

    # Check probabilities sums to 1 (within floating point roundoff)
    stopifnot(isTRUE(all.equal(sum(slip_rate_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(b_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_min_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_max_prob), 1)))

    # Check that the range of Mw in the table covers the range of Mw_min/Mw_max)
    stopifnot(min(event_table$Mw) <= min(Mw_min))
    stopifnot(max(event_table$Mw) >= max(Mw_max))


    # Get all combinations of the parameters, and all combinations of the
    # probabilities
    all_par_combo = expand.grid(slip_rate = slip_rate, b = b, Mw_min = Mw_min, 
        Mw_max = Mw_max)
    all_par_combo_prob = expand.grid(slip_rate_prob = slip_rate_prob,
        b_prob = b_prob, Mw_min_prob = Mw_min_prob, Mw_max_prob = Mw_max_prob)

    # Collapse the probabilities to a single one for each set of table
    # parameters
    all_par_prob = rep(NA, nrow(all_par_combo))
    for(i in 1:nrow(all_par_combo)){
        all_par_prob[i] = prod(all_par_combo_prob[i,])    
    }
  
    # Check the summed probabilities still = 1 
    stopifnot(isTRUE(all.equal(sum(all_par_prob), 1.0)))

    # Make sequence at which function evaluations occur
    min_Mw_min = min(Mw_min)
    max_Mw_max = max(Mw_max)
    npts = round((max_Mw_max - min_Mw_min)/computational_increment) + 1
    Mw_seq = seq(min_Mw_min, max_Mw_max, len=npts)

    # Get all Mw values in the eq table, and find the spacing between unique values
    table_Mw_values = sort(unique(event_table$Mw))
    table_Mw_increment = diff(table_Mw_values)
    range_of_table_increments = diff(range(table_Mw_increment))
    
    if(!isTRUE(all.equal(range_of_table_increments, 0.))){
        print('Earthquake table has these Mw values:')
        print(table_Mw_values)
        stop('The sorted unique Mw values must be an evenly-spaced sequence')
    }
    table_Mw_increment = table_Mw_increment[1]


    # Shorthand notation
    eq_Mw = event_table$Mw
    eq_slip = event_table$slip
    eq_area = event_table$area

    # Check that all event_conditional_probabilities
    # sum to 1.0, for fixed Mw
    for(Mw in table_Mw_values){
        ii = which(eq_Mw == Mw)
        conditional_prob_sum = 
            sum(event_conditional_probabilities[ii])
        if(!isTRUE(all.equal(conditional_prob_sum, 1.0))){
            print(paste0('Error: conditional probabilities for Mw = ', Mw, 
                ' do not sum to 1.0'))
            print(paste0('They sum to: ', conditional_prob_sum))
            print('The values are:')
            print(event_conditional_probabilities[ii])
            stop('Invalid conditional probabilities')
        }
    }
    
    # Compute the rates for each final branch of the logic tree
    all_rate_matrix = matrix(NA, ncol = length(Mw_seq), nrow=nrow(all_par_combo))
    for(i in 1:nrow(all_rate_matrix)){
        # Get parameter vector
        par = all_par_combo[i,] 

        # Evaluate LHS of balance equation
        LHS = (sourcezone_total_area * 1e+06) * par$slip_rate

        lower_Mw = eq_Mw - table_Mw_increment/2
        upper_Mw = eq_Mw + table_Mw_increment/2

        rate_of_events_of_size_Mw_aeq0 = 
            Mw_exceedence_rate_truncated_gutenberg_richter(lower_Mw, a = 0, 
                b = par$b, Mw_min=par$Mw_min, Mw_max = par$Mw_max) -
            Mw_exceedence_rate_truncated_gutenberg_richter(upper_Mw, a = 0, 
                b = par$b, Mw_min=par$Mw_min, Mw_max = par$Mw_max)

        RHS = sum(eq_slip * (eq_area*1e+06) * rate_of_events_of_size_Mw_aeq0 * 
            event_conditional_probabilities)

        # Back-calculate 'a'
        # 10^a = LHS/RHS
        a_parameter = log10(LHS/RHS)

        # Compute rates for a truncated Gutenberg Richter model, assuming a = 0
        # We subsequently correct a to match the slip
        all_rate_vec = Mw_exceedence_rate_truncated_gutenberg_richter(Mw_seq, 
            a = a_parameter, b = par$b, Mw_min = par$Mw_min, Mw_max = par$Mw_max)

        all_rate_matrix[i,] = all_rate_vec
    }

    final_rates = rep(0, length(Mw_seq))
    for(i in 1:nrow(all_rate_matrix)){
        final_rates = final_rates + all_par_prob[i]*all_rate_matrix[i,]
    }

    # For the output function, if we exceed the bounds on the LHS we can return
    # the most extreme value, and on the RHS we should return zero.
    # We take care of the LHS in approxfun, and the RHS below
    output_function = approxfun(Mw_seq, final_rates, rule=2)
    output_function2<-function(Mw){
        output = (Mw <= (max_Mw_max))*output_function(Mw)
        return(output)
    }

    return(output_function2)

}

#' Evaluate the rate of earthquakes with magnitude > Mw using a doubly
#' truncated Gutenberg Richter distribution.
#'
#' The Gutenberg Richter distribution gives the number of earthquakes with magnitude > Mw as:\cr
#' N_{GR}(x > Mw) = 10^(a-bMw) \cr
#' where a and b are constants. Note 10^(a) is the rate of earthquakes with Mw > 0.
#' By differentiating N_{GR} we can estimate the number of events with 
#' (Mw - dx/2 <= magnitude <= Mw +dx/2) as:\cr
#' N(Mw -dx/2 <= x <= Mw + dx/2) ~= dx * n(Mw) = dx * ([10^(a-bMw)] * bln(10)) \cr
#' where n(Mw) is the negative of the derivative of N(x>Mw). \cr
#' n(Mw) is like a scaled probability-density-function of the number of
#' earthquakes. Unlike a pdf, n(Mw) integrates to the rate of earthquakes,
#' instead of 1.\cr
#' For the truncated Gutenberg Richter distribution, n(Mw) is truncated between
#' lower and upper Mw limits (i.e. set to zero outside these limits). \cr
#' We then have the equivalent of N_{GR} for the truncated distribution as \cr:
#' N_{TGR}(x > Mw) = 10^(-max(Mw, Mw_min)*b + a) - 10^(-Mw_max*b + a) \cr
#' Notice that now, 10^(a) is no-longer the rate of earthquakes with Mw > 0 -- instead
#' that rate is \cr:
#' N_{TGR}(x > Mw_min) = 10^(a - b*Mw_min) - 10^(a - b*Mw_max) \cr
#'
#' @param Mw Moment magnitude
#' @param a The a parameter
#' @param b The b parameter
#' @param Mw_min the lower truncated moment magnitude
#' @param Mw_max the upper truncated moment magnitude
#' @return The rate of events with magnitude > Mw
#' @export
Mw_exceedence_rate_truncated_gutenberg_richter<-function(Mw, a, b, Mw_min, Mw_max){
   
    N_mag_gt_Mw = 10^(a - b*pmax(Mw, Mw_min)) - 10^(a - b*Mw_max)

    output = N_mag_gt_Mw*(Mw <= Mw_max)
    
    return(output)
}
