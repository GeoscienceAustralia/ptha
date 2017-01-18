#'
#' Given Mw, compute parameters which characterise the 'probability
#' distribution' of slip
#'
#' FIXME: This function is not currently exported -- do we want to keep it?
#' Seems a good idea. \cr
#' We assume mu is constant, and log10(Area) has a normal distribution. Since
#' the seismic moment = mu*slip*Area, this implies that log10(Slip) has a normal
#' distribution as well. So the mean and standard deviation of the log10
#' variable completely describe it. As well as returning those values, we return
#' the mean slip (mean of untransformed slip) since its calculation is slightly
#' inconvenient.
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

    # To multiply area by 1e+06 (i.e. km^2 --> m^2) we just add log10(1e+06) to
    # mu_log10_area
    mu_log10_area = mu_log10_area + log10(1.0e+06)

    mu_log10_slip = (-mu_log10_area + log10(M0) - log10(mu))
    sd_log10_slip = sd_log10_area

    log10_slip_par = c(mu_log10_slip, sd_log10_slip)
    names(log10_slip_par) = c('log10_mean', 'log10_sd')

    # If log10(slip) is normally distributed, then this is the mean_slip for
    # ruptures of size Mw:
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

#'
#' Compute the fraction of seismic moment associated with earthquake > a
#' threshold
#'
#' For long-term slip rate computations we need to know how much slip is caused
#' by events with Mw > mwmin. This function helps compute that, though it
#' requires appropriately designed inputs
#'
#' @param Mws numeric vector with evenly spaced Mw values, increasing, densely
#' spaced enough for integration
#' @param rate_Mws numeric vector with an occurrence rate for each value in
#' Mws. THIS IS NOT AN EXCEEDANCE RATE, IT IS A RATE FOR EVENTS WITH Mw = Mws[i].
#' The value will thus depend on how closely spaced Mws is.
#' @param moment_Mws vector with the seismic moment for each Mws [follows from
#' standard formulae]
#' @param mwmin numeric We want to know the fraction of the moment associated
#' with events of this size or larger
#' @return number between 0-1, giving the fraction of moment associated with
#' events with Mw > mwmin
compute_moment_fraction_from_events_greater_or_equal_than_mwmin<-function(
    Mws, rate_Mws, moment_Mws, mwmin){

    dMw = Mws[2] - Mws[1]
    if( !isTRUE(all.equal(diff(Mws), rep(dMw,len=length(Mws)-1))) ){
        stop("Mws must be evenly spaced")
    }

    if(mwmin < Mws[1] | mwmin > max(Mws)){
        stop("mwmin not compatible with Mws")
    }

    events_high = which(Mws >= mwmin)

    fraction = sum( rate_Mws[events_high] * moment_Mws[events_high] ) / 
        sum(rate_Mws * moment_Mws)

    return(fraction)
}


#' Assign conditional probabilities to earthquakes with the same Mw
#'
#' Consider all the earthquake events on a discrete source with magnitude Mw.
#' We wish to assign conditional probabilities to these events (conditional on
#' the fact that an event of size Mw occurred). This allows us to compute the rate
#' of any individual event as: \cr
#' 'event_conditional_probability * rate_of_events_of_size_Mw' \cr
#' One way to do this is assume that for a given magnitude, all event
#' conditional probabilities are equal ('all_equal'). \cr
#' This might be reasonable if the events are identical in area and slip, and
#' just differ in location. \cr 
#' Another approach is to assume the conditional probabilities vary
#' inversely with the average slip of the event ('inverse_slip'). This might be
#' reasonable if the individual events are almost the same, but there are slight
#' differences in the event area and slip, and we want each event to make the
#' same contribution to the long-term average slip rate. \cr
#' Other methods could be proposed, and this routine also supports a
#' user-defined function
#'
#' @param event_table A data.frame containing information on all earthquake
#' events on the source-zone. It's the output of
#' \code{get_all_earthquake_events}
#' @param conditional_probability_model character or function. Supported
#' character values are 'all_equal' and 'inverse_slip', corresponding to cases
#' discussed above. \cr
#' If a function, it must take a single argument, which is a subset of the
#' event_table containing all events with fixed Mw, and return the conditional
#' probability of each event.
#' @return A numeric vector with length = nrow(event_table) giving the
#' conditional probabilities. 
#' @export
#'
get_event_probabilities_conditional_on_Mw<-function(
    event_table,
    conditional_probability_model = 'all_equal'){

    # Define functions corresponding to supported conditional_probability
    # models. 
    # These take as input a subset of the event table containing events with
    # identical Mw
    if(is.character(conditional_probability_model) && 
        (conditional_probability_model == 'all_equal')){
        # 'all_equal' probability function
        pfun<-function(event_table_Mw){
            N = length(event_table_Mw[,1])
            return( rep(1, N)/N)
        }

    }else if(is.character(conditional_probability_model) &&
        (conditional_probability_model == 'inverse_slip')){

        # 'inverse_slip' probability function
        pfun<-function(event_table_Mw){
            slip_inv = 1/event_table_Mw$slip
            return(slip_inv/sum(slip_inv))
        }

    }else if( class(conditional_probability_model) == 'function'){

        # Here the user must provide a function which can convert a subset of
        # the event table, containing all events with fixed Mw, to a
        # conditional probability

        pfun = conditional_probability_model

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
        if(!isTRUE(all.equal(sum(event_conditional_probabilities[ii]), 1.0))){
            msg = paste0('ERROR: event conditional probabilities with Mw = ',
                Mw, ' do not seem to sum to 1.0. This suggests an error in ',
                'the conditional_probability_model')
            stop(msg)
        }
    }

    # Final check
    if(any(is.na(event_conditional_probabilities))){
        warnings('NA values among event_conditional_probabilities')
    }

    return(event_conditional_probabilities)
}

#'
#' Get a function to compute the rate of earthquakes with magnitude > Mw
#'
#' Create a function f(Mw) which returns the rate of events with 
#' magnitude > Mw. Here Mw is a vector of earthquake magnitudes. \cr
#' The user provides parameters for a variant of the truncated Gutenberg Richter
#' model, although the parameter 'a' is adjusted to match the long-term
#' seismogenic slip rate (m/year). Details are explained below. \cr
#' To account for uncertainties in the GR model parameters, the user may specify
#' each as a vector of one-or-more parameters with associated probabilities
#' (which must sum to 1, of course). The code then treats all combinations of
#' those parameters as having probability equal to the product of the individual
#' parameter probabilities. It then integrates over the chosen Gutenberg Richter
#' models to produce a single exceedance rate function, which is returned. \cr
#' The output rate function also supports a number of optional arguments, which
#' provide various details on the rate function on each branch of the logic
#' tree. See details below \cr
#'
#' @note The 'a' parameter in the GR model is constrained to match the long-term
#' seismogenic slip rate.  Given particular values of all the parameters except
#' a, and a long-term seisomogenic slip rate (= long_term_slip_rate *
#' coupling_coefficient), the long-term moment rate should equal the rate
#' integrated from the individual events. Mathematically: \cr
#' \cr
#' (definitions) \cr
#' S = long term semsmogenic slip rate \cr
#' A = source-zone area \cr
#' \eqn{ s_{i} }= slip for event i
#' \eqn{ a_{i} } = area for event i
#' \eqn{rMw_{i} } = rate of events with the same magnitude Mw as event i
#' \eqn{ Pr(event_{i} | event with size Mw_{i}) } = the conditional probability
#' of each event with the same Mw in the event table [i.e. if we take all events
#' with a given Mw from the table, their conditional probabilities should sum to
#' 1]. \cr
#' (main equation) \cr
#' \deqn{ \mu  S  A = 
#' \sum_{i \in events} ( \mu s_{i} a_{i} Pr(event_{i} | event_of_size_Mw_{i}) rMw_{i}) }\cr
#' \cr
#' The form of the GR relation implies that a factor of 10^(a) will appear in
#' the RHS term rate_of_events_of_size_Mw_{i}, and nowhere else. Hence to
#' compute 'a', we just compute the LHS and RHS of the above equations,
#' tentatively assuming a=0 in the latter case, and then it follows that a =
#' log10(LHS/RHS). \cr
#' The code also optionally allows us to account for the fraction of seismic
#' moment released by earthquakes with magnitude smaller than Mw_min. This
#' reduces the LHS term in the equation above.
#'
#' @param slip_rate numeric vector with 1 or more long-term seismogenic slip
#' rates for the source-zone (m/year). This gives the long-term average slip on
#' the fault due to seismic events (creep is ignored, so these rates must
#' account for the rate of coupling)
#' @param slip_rate_prob numeric vector with weights for the slip_rate
#' parameters
#' @param b numeric vector with 1 or more 'b' parameters
#' @param b_prob numeric vector with weights for the b parameters
#' @param Mw_min numeric vector with 1 or more 'Mw_min' parameters
#' @param Mw_min_prob numeric vector with weights for the Mw_min parameters
#' @param Mw_max numeric vector with 1 or more 'Mw_max' parameters
#' @param Mw_max_prob numeric vector with weights for the Mw_max parameters
#' @param sourcezone_total_area The total area of the sourcezone (km^2)
#' @param event_table A data.frame containing information on all events in the
#' source-zone. Output of \code{get_all_earthquake_events}. This is needed to
#' evaluate the RHS of the seismic moment balance equation.
#' @param event_conditional_probabilities A vector of conditional probabilities
#' for all events in the source-zone (conditional on the fact than an event with
#' their magnitude did occur). Output of
#' \code{get_event_probabilities_conditional_on_Mw}. This is needed to evaluate
#' the RHS of the seismic moment balance equation.
#' @param computational_increment For each parameter combination in the logic
#' tree, the rate function is computed at points in a sequence from min(Mw_min)
#' to max(Mw_max), with spacing computational_increment (adjusted if required so
#' that the extremes are included). The function value is then evaluated via
#' lookup. 
#' @param Mw_frequency_distribution character vector giving the variant on the
#' Gutenberg Richter model used (if more than 1 they are passed to the logic
#' tree).  Either 'truncated_gutenberg_richter' or
#' 'characteristic_gutenberg_richter'. 
#' @param Mw_frequency_distribution_prob numeric vector. Logic tree weights for
#' each Mw_frequency_distribution type
#' @param update_logic_tree_weights_with_data logical. If TRUE, then a value for
#' Mw_count_duration must be provided. The weights for each parameter
#' combination in the tree are then updated with Bayes-theorem, based on the
#' probability that they would produce Mw_count_duration[2] events which exceed
#' Mw_count_duration[1] within time duration Mw_count_duration[3]. This assumes
#' earthquake timings follow a Poisson-process.
#' @param Mw_count_duration numeric vector of length 3 providing some data. The
#' first entry is an Mw value, the second is the number of events exceeding Mw
#' at the site, and the third is the duration of observation (in years).  If
#' \code{update_logic_tree_weights_with_data=TRUE}, this is used to re-weight
#' parameter combinations in the logic tree, with more weight given to those
#' which predict similar event frequencies as the data.
#' @param Mw_2_M0 function which takes an earthquake magnitude and returns the
#' corresponding seismic moment M0
#' @param account_for_moment_below_mwmin logical. If FALSE assume all
#' seismically coupled slip occurs due to earthquakes with Mw > Mwmin. If TRUE,
#' estimate the fraction of seismically coupled slip that occurs due to
#' earthquakes with Mw > Mwmin, and use that when computing the rates.
#' @return Function f(Mw) of a numeric vector Mw, which by default returns the
#' rate of events with magnitude > Mw. It has a number of optional arguments.
#' \cr
#' If the optional argument 'bounds=TRUE' is passed,
#' then the function returns upper/lower/median rates of the logic tree as well
#' as the probability weighted mean rate. \cr
#' If the optional argument \code{return_all_logic_tree_branches=TRUE} is
#' passed, then the function returns a list containing a data.frame with all
#' combinations of parameters, a vector with their respective
#' probabilities(weights), a vector with the Mw sequence at which the function
#' is tabulated, and a matrix with the tabulated values for every branch of the
#' logic tree. The latter option can be used to scrutinize the results from each
#' branch of the logic tree. \cr
#' If the optional argument \code{quantiles} = a vector of probabilities is
#' passed, then the function evaluates the rate function at each given quantile
#' of the logic tree (e.g. quantiles = c(0.1, 0.5, 0.9) would return the lower
#' 10\%, median, and upper 90\% rates from the logic tree at each Mw) \cr
#' If the optional argument \code{return_random_curve=TRUE} is passed, then the
#' function will randomly pick a parameter combination from the logic tree (with
#' probability equal to the probability of each branch), and evaluate the
#' exceedance rate of Mw for that parameter combination. This can be useful for
#' some monte-carlo computations. \cr
#' Only one of these optional arguments can be passed at a time.
#' @export
#'
rate_of_earthquakes_greater_than_Mw_function<-function(
    slip_rate, 
    slip_rate_prob, 
    b, 
    b_prob, 
    Mw_min, 
    Mw_min_prob, 
    Mw_max, 
    Mw_max_prob, 
    sourcezone_total_area, 
    event_table, 
    event_conditional_probabilities,
    computational_increment=0.01,
    Mw_frequency_distribution='truncated_gutenberg_richter',
    Mw_frequency_distribution_prob = 1,
    update_logic_tree_weights_with_data = FALSE,
    Mw_count_duration = c(NA, NA, NA),
    Mw_2_M0 = function(x) M0_2_Mw(x, inverse=TRUE),
    account_for_moment_below_mwmin=FALSE){

    # Check data
    stopifnot(length(slip_rate) == length(slip_rate_prob))
    stopifnot(length(b) == length(b_prob))
    stopifnot(length(Mw_min) == length(Mw_min_prob))
    stopifnot(length(Mw_max) == length(Mw_max_prob))
    stopifnot(length(Mw_frequency_distribution) == 
        length(Mw_frequency_distribution_prob))

    # Check (conservative) on the input Mw_min/Mw_max
    stopifnot(min(Mw_max) > max(Mw_min))

    # Check probabilities sums to 1 (within floating point roundoff)
    stopifnot(isTRUE(all.equal(sum(slip_rate_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(b_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_min_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_max_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_frequency_distribution_prob), 1)))

    # Check that the range of Mw in the table covers the range of
    # Mw_min/Mw_max)
    stopifnot(min(event_table$Mw) <= min(Mw_min))
    stopifnot(max(event_table$Mw) >= max(Mw_max))


    # Get all combinations of the parameters, and all combinations of the
    # probabilities
    all_par_combo = expand.grid(
        slip_rate = slip_rate, 
        b = b, 
        Mw_min = Mw_min, 
        Mw_max = Mw_max, 
        Mw_frequency_distribution=Mw_frequency_distribution)
    all_par_combo_prob = expand.grid(
        slip_rate_prob = slip_rate_prob,
        b_prob = b_prob, 
        Mw_min_prob = Mw_min_prob, 
        Mw_max_prob = Mw_max_prob,
        Mw_frequency_distribution_prob = Mw_frequency_distribution_prob)

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


    # Get all Mw values in the eq table, and find the spacing between unique
    # values
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
    all_rate_matrix = matrix(NA, ncol = length(Mw_seq), 
        nrow=nrow(all_par_combo))
    for(i in 1:nrow(all_rate_matrix)){
        # Get parameter vector
        par = all_par_combo[i,] 
        
        if(par$Mw_frequency_distribution == 'truncated_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_truncated_gutenberg_richter
        }else if(par$Mw_frequency_distribution == 'characteristic_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_characteristic_gutenberg_richter
        }else{
            stop(paste0(" Value of Mw_frequency_distribution: ", 
                par$Mw_frequency_distribution, " not recognized"))
        }

        if(account_for_moment_below_mwmin){
            # Compute the fraction of seismic moment associated with
            # earthquakes between mw_min and mw_max.

            # Finely spaced sequence used for numerical integration
            broad_Mw_seq = seq(min(par$Mw_min, 0), par$Mw_max, len = 1e+04)
            dmw0 = broad_Mw_seq[2] - broad_Mw_seq[1]
            broad_rates = Mfd(broad_Mw_seq - dmw0/2, a = 0, b = par$b, 
                              Mw_min=broad_Mw_seq[1], Mw_max = par$Mw_max) -
                          Mfd(broad_Mw_seq + dmw0/2, a = 0, b = par$b, 
                              Mw_min=broad_Mw_seq[1], Mw_max = par$Mw_max) 
            broad_seismic_moment = Mw_2_M0(broad_Mw_seq)

            moment_fraction = 
                compute_moment_fraction_from_events_greater_or_equal_than_mwmin(
                    broad_Mw_seq, broad_rates, broad_seismic_moment, par$Mw_min)
           
            stopifnot(moment_fraction >= 0 & moment_fraction <= 1)
        }else{
            moment_fraction = 1
        }

        # Evaluate LHS of seismic moment balance equation
        LHS = (sourcezone_total_area * 1e+06) * par$slip_rate * moment_fraction


        lower_Mw = eq_Mw - table_Mw_increment/2
        upper_Mw = eq_Mw + table_Mw_increment/2

        # We need to solve for 'a', given the long-term slip rate.
        # This can be done by initially setting 'a' to zero, then
        # back-calculating the value required to match the long-term slip
        rate_of_events_of_size_Mw_aeq0 = 
            Mfd(lower_Mw, a = 0, b = par$b, Mw_min=par$Mw_min, 
                Mw_max = par$Mw_max) -
            Mfd(upper_Mw, a = 0, b = par$b, Mw_min=par$Mw_min, 
                Mw_max = par$Mw_max)

        RHS = sum(eq_slip * (eq_area*1e+06) * rate_of_events_of_size_Mw_aeq0 * 
            event_conditional_probabilities)

        # Back-calculate 'a'
        # 10^a = LHS/RHS
        # NOTE: This assumes that 'a' appears in Mfd as 10^(a) {only}
        a_parameter = log10(LHS/RHS)

        # Compute rates of exceedance for a truncated Gutenberg Richter model.
        all_rate_vec = Mfd(Mw_seq, a = a_parameter, b = par$b, 
            Mw_min = par$Mw_min, Mw_max = par$Mw_max)

        all_rate_matrix[i,] = all_rate_vec
    }


    # Save original probabilities in case we update them
    all_par_prob_prior = all_par_prob   
 
    if(update_logic_tree_weights_with_data){
        # Adjust the weights of logic-tree branches based on the data

        if(sum(is.na(Mw_count_duration)) > 0){
            msg = paste0('Must provide Mw_count_duration data if ', 
                'update_logic_tree_weights_with_data=TRUE')
            stop(msg)
        }

        data_thresh = Mw_count_duration[1]
        data_count = Mw_count_duration[2]
        data_observation_duration = Mw_count_duration[3]
        
        if( (data_thresh < min(Mw_seq)) | (data_thresh > max(Mw_seq)) ){
            stop('The Mw threshold of the data must be within [Mw_min, Mw_max]')
        }
        
        model_rates = all_par_prob*0
        for(i in 1:length(model_rates)){
            model_rates[i] = approx(Mw_seq, all_rate_matrix[i,], 
                xout=data_thresh)$y
        }
        pr_data_given_model = dpois(data_count, 
            lambda = (model_rates*data_observation_duration))

        sum_pr_data_given_model = sum(pr_data_given_model)
        if(sum_pr_data_given_model == 0){
            stop('Mw_count_duration data is impossible under every model')
        }
        
        all_par_prob =  all_par_prob * pr_data_given_model / 
            sum(all_par_prob * pr_data_given_model)
    }else{

        all_par_prob = all_par_prob_prior

    }

    # Compute the average rate
    final_rates = rep(0, length(Mw_seq))
    for(i in 1:nrow(all_rate_matrix)){
        final_rates = final_rates + all_par_prob[i] * all_rate_matrix[i,]
    }

    # Store the max/min rates from all branches of the logic tree,
    # to help quantify uncertainty later.
    # Consider producing a fuller picture from all_rate_matrix
    upper_rate = apply(all_rate_matrix, 2, max)
    lower_rate = apply(all_rate_matrix, 2, min)

    quantile_rate_fun<-function(rates, p){
        sorted_rates = sort(rates, index.return=TRUE)
        sorted_prob = all_par_prob[sorted_rates$ix]
        cumulative_sorted_prob = cumsum(sorted_prob)
        # The median occurs when the cumulative probability >= 0.5
        ind = min(which(cumulative_sorted_prob >= p))
        return(sorted_rates$x[ind])
    }

    # Compute the 50th percentile
    median_rate = apply(all_rate_matrix, 2, 
        f<-function(rates) quantile_rate_fun(rates, p=0.5)
        )

    # For the output function, if we exceed the bounds on the LHS we can return
    # the most extreme value, and on the RHS we should return zero.
    # We take care of the LHS in approxfun, and the RHS below.
    output_function = approxfun(Mw_seq, final_rates, rule=2)
    upper_function = approxfun(Mw_seq, upper_rate, rule=2)
    lower_function = approxfun(Mw_seq, lower_rate, rule=2)
    median_function = approxfun(Mw_seq, median_rate, rule=2)

    # Function to compute the rate, with a number of useful extra options
    #
    # @param Mw the Mw value at which the rate is desired
    # @param bounds logical. If TRUE, return some summary statistics on the
    # range of rates suggested by the logic tree. More detailed information
    # can be obtained using quantiles or return_all_logic_tree_branches.
    # @param return_all_logic_tree_branches logical. If TRUE, return a list
    # with the all_rate_matrix, Mw_seq, all_par_prob, all_par_prob_prior,
    # and all_par. This option is useful for doing non-standard manipulations
    # of the results, and debugging.
    # @param quantiles numeric vector of probabilities or NULL. If not NULL,
    # then return quantiles of the exceedance rate for Mw (based on the
    # logic tree branches and probabilities)
    # @param return_random_curve TRUE/FALSE. If TRUE, randomly pick a parameter
    # combination from the logic tree (with probability equal to the probability
    # of each branch), and evaluate the exceedance rate of Mw for that
    # parameter combination. This can be useful for some monte-carlo
    # computations.
    # @return Depends on the input arguments, see comments above.
    output_function2<-function(
        Mw, 
        bounds=FALSE, 
        return_all_logic_tree_branches=FALSE, 
        quantiles=NULL,
        return_random_curve=FALSE){

        # Check that at most one 'special option' is taken
        quantiles_not_null = !is.null(quantiles)
        options_sum = bounds + return_all_logic_tree_branches + 
            quantiles_not_null + return_random_curve

        if(options_sum > 1){
            msg = paste0('At most one of bounds, ', 
                'return_all_logic_tree_branches, ',
                '!is.null(quantiles), and return_random_curve can be TRUE')
            stop(msg)
        }


        # Option to just dump the key data. This can be useful for
        # debugging / plotting non-standard outputs, etc 
        if(return_all_logic_tree_branches){
            output = list(all_rate_matrix = all_rate_matrix, 
                          Mw_seq = Mw_seq,
                          all_par_prob = all_par_prob,
                          all_par_prob_prior = all_par_prob_prior,
                          all_par = all_par_combo)
            return(output)
        }

        if(return_random_curve){
            # Grab a curve 'at random' with probabilities based on all_par_prob

            stopifnot(bounds==FALSE & return_all_logic_tree_branches==FALSE)

            par = sample(1:length(all_par_combo[,1]), size=1, 
                prob=all_par_prob)
            output = approx(Mw_seq, all_rate_matrix[par,], xout=Mw, rule=2)$y *
                (Mw <= max_Mw_max)
            return(output)
        }

        # Typical case
        output = (Mw <= (max_Mw_max))*output_function(Mw)

        if(bounds){
            # Compute some summary statistics
            output_up = (Mw <= max_Mw_max)*upper_function(Mw)
            output_low = (Mw <= max_Mw_max)*lower_function(Mw)
            output_med = (Mw <= max_Mw_max)*median_function(Mw)
           
            if(length(Mw) == 1){ 
                output = c(output, output_low, output_up, output_med)
                names(output) = c('rate', 'min', 'max', 'median')
            }else{
                output = cbind(output, output_low, output_up, output_med)
                colnames(output) = c('rate', 'min', 'max', 'median')
            }

        }else if(!is.null(quantiles)){
            # Compute results for a set of quantiles

            if(length(Mw) == 1){
                output = rep(NA, length(quantiles))
                names(output) = quantiles
            }else{
                output = matrix(NA, ncol=length(quantiles), nrow=length(Mw))
                colnames(output) = quantiles
            }
            for(i in 1:length(quantiles)){
                # For each quantile, evaluate the rate curve with interpolation
                quantile_rate = apply(all_rate_matrix, 2, 
                    f<-function(rate) quantile_rate_fun(rate, p=quantiles[i]))
                if(length(Mw) == 1){
                    output[i] = 
                        approx(Mw_seq, quantile_rate, xout=Mw, rule=2)$y * 
                        (Mw <= max_Mw_max)
                }else{
                    output[,i] = 
                        approx(Mw_seq, quantile_rate, xout=Mw, rule=2)$y * 
                        (Mw <= max_Mw_max)
                }
            }
        }

        return(output)
    }

    return(output_function2)

}

#' Evaluate the rate of earthquakes with magnitude > Mw using a doubly
#' truncated Gutenberg Richter distribution.
#'
#' The Gutenberg Richter distribution gives the number of earthquakes with
#' magnitude > Mw as:\cr
#' N_{GR}(x > Mw) = 10^(a-bMw) \cr
#' where a and b are constants. Note 10^(a) is the rate of earthquakes with Mw >
#' 0. \cr
#' By differentiating N_{GR} we can estimate the number of events with 
#' (Mw - dx/2 <= magnitude <= Mw +dx/2) as:\cr
#' N(Mw -dx/2 <= x <= Mw + dx/2) ~= dx * n(Mw) = dx * ([10^(a-bMw)] * bln(10)) \cr
#' where n(Mw) is the negative of the derivative of N(x>Mw). \cr
#' n(Mw) is like a scaled probability-density-function of the number of
#' earthquakes. Unlike a pdf, n(Mw) integrates to the rate of earthquakes,
#' instead of 1.\cr
#' For the truncated Gutenberg Richter distribution, n(Mw) is truncated between
#' lower and upper Mw limits (i.e. set to zero outside these limits). \cr
#' We then have the equivalent of N_{GR} for the truncated distribution as: \cr
#' N_{TGR}(x > Mw) = 10^(-max(Mw, Mw_min)*b + a) - 10^(-Mw_max*b + a) \cr
#' Notice that now, 10^(a) is no-longer the rate of earthquakes with Mw > 0 --
#' instead that rate is :\cr
#' N_{TGR}(x > Mw_min) = 10^(a - b*Mw_min) - 10^(a - b*Mw_max) \cr
#'
#' @param Mw Moment magnitude
#' @param a The a parameter
#' @param b The b parameter
#' @param Mw_min the lower truncated moment magnitude
#' @param Mw_max the upper truncated moment magnitude
#' @return The rate of events with magnitude > Mw
#' @export
Mw_exceedance_rate_truncated_gutenberg_richter<-function(
    Mw, a, b, Mw_min, Mw_max){
   
    N_mag_gt_Mw = 10^(a - b*pmax(Mw, Mw_min)) - 10^(a - b*Mw_max)

    output = N_mag_gt_Mw*(Mw <= Mw_max)
    
    return(output)
}



#' Evaluate the rate of earthquakes with magnitude > Mw using a characteristic
#' variant of the Gutenberg Richter distribution.
#' 
#' The 'pure' Gutenberg Richter distribution gives the number of earthquakes
#' with magnitude > Mw as:\cr
#' N_{GR}(x > Mw) = 10^(a-bMw) \cr
#' where a and b are constants. Note 10^(a) is the rate of earthquakes with Mw
#' > 0. \cr
#' The characteristic formulation used here ajusts this as: \cr
#' N_{GR}(x > Mw) = 10^(a-bMw) if Mw_min <= Mw <= Mw_max \cr
#' ~~~~~~~~~~~~~  = 10^(a-bMw_min) if Mw < Mw_min \cr
#' ~~~~~~~~~~~~~  = 0 if Mw > Mw_max \cr
#' This means that the rate of earthquakes of size EXACTLY Mw_max is finite.
#' @param Mw Moment magnitude
#' @param a The a parameter
#' @param b The b parameter
#' @param Mw_min the lower truncated moment magnitude
#' @param Mw_max the upper truncated moment magnitude
#' @return The rate of events with magnitude > Mw
#' @export
#' 
Mw_exceedance_rate_characteristic_gutenberg_richter<-function(
    Mw, a, b, Mw_min, Mw_max){

    N_mag_gt_Mw = 10^(a - b*pmax(Mw, Mw_min))

    output = N_mag_gt_Mw*(Mw <= Mw_max)
    return(output)
}
