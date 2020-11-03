#
# Multi-gauge scenario selection
#


# Read the get_PTHA_results.R routines in an environment
ptha18 = new.env()
ptha_access_script_location = '../../get_PTHA_results.R'
if(!file.exists(ptha_access_script_location)){
    stop('Error: get_PTHA_results.R not found -- please edit the file path to suit your machine, and try again')
}
source(ptha_access_script_location, local=ptha18, chdir=TRUE)
library(rptha)

#
# Get hazard curves and peak-stages for all events at a specified set of gauges.
# This can later be used for event selection (see get_matching_scenarios).
#
# @param site_gauge_inds indices of the gauges {within the PTHA18 hazard points array}
# @param site_gauges data.frame with lon/lat/elev/gaugeID for the gauges
# corresponding to site_gauge_inds
# @param source_zone name of the source_zone (e.g. 'kermadectonga2')
# @param slip_type either 'stochastic' or 'variable_uniform', denoting the type
# of earthquake slip model to use.
# @return The function environment which contains the stage_exrate_curves at
# all gauges, and the peak_stage values for all events on the source-zone,
# among other things.
#
get_hazard_curves_and_peak_stages<-function(site_gauge_inds, site_gauges, source_zone, slip_type){

    if(length(site_gauge_inds) != nrow(site_gauges)) stop('length(site_gauge_inds) is not equal to nrow(site_gauges)')

    #
    # Read the hazard curve data for each gauge 
    #
    stage_exrate_curves = vector(mode='list', length=length(site_gauge_inds))
    for(i in 1:length(site_gauge_inds)){
        stage_exrate_curves[[i]] = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
            target_index = site_gauge_inds[i],
            non_stochastic_slip_sources=TRUE)
    }

    # Store the peak-stage values for each gauge, for every event on the source-zone
    peak_stage_sourcezone = vector(mode='list', length=length(site_gauge_inds))
    for(i in 1:length(site_gauge_inds)){
        # Get the stochastic slip peak stage values for each gauge
        peak_stage_sourcezone[[i]] = ptha18$get_peak_stage_at_point_for_each_event(
            target_index = site_gauge_inds[i], all_source_names = source_zone,
            slip_type=slip_type, include_earthquake_data=FALSE)
    }

    # Convert into a matrix of peak stage values, with each column corresponding to a gauge,
    # and each row being an event
    peak_stage_sourcezone_matrix = do.call(cbind, 
        lapply(peak_stage_sourcezone, f<-function(x) x[[source_zone]]$max_stage))
    colnames(peak_stage_sourcezone_matrix) = round(site_gauges$gaugeID, 2)

    return(environment())
}

# Choose a default function to compute average slip given Mw.
# For some plotting, we need a function f(Mw) which gives the scaling-relation-mean slip,
# assuming a mean rupture area. While it can be provided, it's a bit inconvenient to create,
# whereas for PTHA18 we know what it should be (depends on the focal mechanism)
get_default_slip_from_Mw_function<-function(source_zone_data){

    if(source_zone_data$unit_source_statistics$rake[1] > 89.9){
        # Thrust faults
        slip_from_Mw_function <-function(x) slip_from_Mw(x)

    }else if(source_zone_data$unit_source_statistics$rake[1] < -89.9){
        # Normal faults -- higher rigidity, difference scaling relation
        slip_from_Mw_function <- function(x){
            slip_from_Mw(x, mu = 6e+10, relation = "Blaser-normal",
                        area_function = function(Mw) {Mw_2_rupture_size(Mw, relation = relation)[1] },
                        constant = 9.05)
            }
    }else{
        stop('Default slip_from_Mw_function only exists for thrust and normal faults')
    }

    return(slip_from_Mw_function)
}

#' Find scenarios matching various criteria using the source_zone_data
#'
#' @param source_zone_data result of function get_hazard_curves_and_peak_stages
#' @param scenario_match a list containing variables which define how we match scenarios.
#' @param slip_from_Mw_function optional function to compute the earthquake
#' slip given Mw. By default use the PTHA18 approach
#'
get_matching_scenarios<-function(source_zone_data, scenario_match, slip_from_Mw_function=NULL){

    if(scenario_match$number_matching_gauges > length(source_zone_data$site_gauge_inds)){
        msg = paste0('You requested scenarios where at least ', 
                     scenario_match$number_matching_gauge, 
                     ' gauges satisfy the target stage(or exrate). This is more than',
                     ' the number of gauges you provided in the source_zone_data (which is ', 
                     length(source_zone_data$site_gauge_inds), '). Please either reduce ', 
                     'scenario_match$number_matching_gauges, or provide more gauges when ', 
                     'creating the source_zone_data')
        stop(msg)
    }

    
    # Copy out key data for the source-zone
    stage_exrate_curves = source_zone_data$stage_exrate_curves
    peak_stage_sourcezone_matrix = source_zone_data$peak_stage_sourcezone_matrix
    slip_type = source_zone_data$slip_type
    source_zone = source_zone_data$source_zone

    # Get the max-stage at the target ARI for each hazard point
    # We have to construct a name that matches the name of a stage_exrate_curve
    if(scenario_match$mu_type == 'variable_mu'){
        # Variable rigidity
        rate_curve_name = paste0('variable_mu_', slip_type, '_slip_', 
                                 scenario_match$hazard_curve_type)
    }else{
        # Constant rigidity
        rate_curve_name = paste0(slip_type, '_slip_', scenario_match$hazard_curve_type)
    }

    # For every gauge, do linear interpolation of the rate-vs-stage curve to
    # get the stage @ scenario_match$exrate.
    max_stage_at_exrate = unlist(lapply(
        stage_exrate_curves, 
        f<-function(x) approx(x[[rate_curve_name]], x$stage, xout=scenario_match$exrate, ties='min')$y))
    # Make a 'max-stage-matrix', with constant values in each column, to
    # facilitate matrix-wise operations below.
    mdim = dim(peak_stage_sourcezone_matrix)
    target_max_stage_matrix = matrix(max_stage_at_exrate, nrow=mdim[1], ncol=mdim[2], byrow=TRUE)

    # Convert the peak-stage matrix to a matrix of integers. 
    # For gauges smaller than the desired exrate stage window, value = 0
    # For gauges within the desired exrate stage window, value = 1
    # For gauges larger than the desired exrate stage window, value = 2
    peak_stage_within_tol = peak_stage_sourcezone_matrix * 0
    peak_stage_within_tol = peak_stage_within_tol + 
        (peak_stage_sourcezone_matrix > 
             target_max_stage_matrix*(1-scenario_match$target_stage_tolerance))
    peak_stage_within_tol = peak_stage_within_tol + 
        (peak_stage_sourcezone_matrix > 
             target_max_stage_matrix*(1+scenario_match$target_stage_tolerance))

    # Count number of gauges within the window
    ngauges_in_window = rowSums(peak_stage_within_tol == 1)
    ngauges_above_window = rowSums(peak_stage_within_tol == 2)

    # Accept scenarios with at least the desired number of gauges in the window,
    # and no gauges that exceed the threshold
    candidate_inds = which((ngauges_in_window >= scenario_match$number_matching_gauges) & 
                           (ngauges_above_window == 0))

    # Get the event data for the candidate events
    candidate_events = ptha18$get_source_zone_events_data(source_zone=source_zone,
        slip_type=slip_type, desired_event_rows=candidate_inds, include_potential_energy=TRUE)

    # We should only consider events that are possible
    k = which(candidate_events$events$rate_annual > 0)

    # Remove the impossible events
    candidate_events$events = candidate_events$events[k,]
    candidate_events$desired_event_rows = candidate_events$desired_event_rows[k]

    # Append the peak stage values
    ii = candidate_inds[k]
    candidate_events$peak_stages = peak_stage_sourcezone_matrix[ii,]

    # Append a matrix that shows whether each gauge is within the tolerance for each event
    candidate_events$peak_stages_within_tol = (
        (peak_stage_sourcezone_matrix[ii,,drop=FALSE] > 
             target_max_stage_matrix[ii,,drop=FALSE]*(1-scenario_match$target_stage_tolerance)) &
        (peak_stage_sourcezone_matrix[ii,,drop=FALSE] < 
             target_max_stage_matrix[ii,,drop=FALSE]*(1+scenario_match$target_stage_tolerance)) ) 

    # Append the peak slip
    candidate_events$peak_slip = sapply(candidate_events$events$event_slip_string, 
        f<-function(x) max(as.numeric(strsplit(x, '_')[[1]])),
        USE.NAMES=FALSE)

    # Useful to also know the ratio of peak slip to the mean slip
    if(is.null(slip_from_Mw_function)){
        # Choose a default approach to computing mean slip from Mw, based on the PTHA18 approach
        slip_from_Mw_function = get_default_slip_from_Mw_function(candidate_events)
    }
    candidate_events$slip_from_Mw_function=slip_from_Mw_function
    # For later plotting it is handy to have the function which computes average slip based on Mw 
    strasser_mean_slip = slip_from_Mw_function(candidate_events$events$Mw)
    candidate_events$peak_slip_ratio = candidate_events$peak_slip/strasser_mean_slip

    # Convenient to know, for each gauge, the max-stage associated with the exceedance-rate.
    candidate_events$max_stage_at_exrate = max_stage_at_exrate

    # Useful to store the criteria used for matching
    candidate_events$scenario_match = scenario_match

    return(candidate_events)

}

# Quick plots to look at the scenarios
summarise_scenarios<-function(candidate_events){

    # What is the distribution of magnitudes for these possible-candidate-events?
    # Better to look at the weighted distribution
    #stem(candidate_events$events$Mw)

    slip_from_Mw_function = candidate_events$slip_from_Mw_function

    #par(mfrow=c(2,2))
    layout(matrix(c(1, 2, 3, 4, 5, 5, 6, 6), byrow=TRUE, ncol=2))
    # Show the event magnitude, and the likeihood they are possible (i.e. chance
    # that Mw-max > Mw)
    plot(candidate_events$events$Mw, candidate_events$events$weight_with_nonzero_rate,
         main='Probability that Mw is possible \n (i.e. Mw < Mw-max )', ylim=c(0, 1),
         xlab="Mw", ylab='Prob',
         cex.main=1.5, cex.lab=1.3, cex.axis=1.3)
    grid(col='orange')
    #
    # Look at the peak-slip for each event
    peak_slip = sapply(candidate_events$events$event_slip_string, 
                       f<-function(x) max(as.numeric(strsplit(x, '_')[[1]])),
                       USE.NAMES=FALSE)
    plot(candidate_events$events$Mw, peak_slip, 
         ylim=c(min(peak_slip)/2, max(peak_slip)),
         main='Mw vs peak-slip', log='y',
         xlab='Mw', ylab='Peak slip', 
         cex.main=1.5, cex.lab=1.3, cex.axis=1.3)
    grid(col='orange')
    strasser_mean_slip = slip_from_Mw_function(candidate_events$events$Mw)
    points(candidate_events$events$Mw, strasser_mean_slip, t='l', col='blue')
    points(candidate_events$events$Mw, strasser_mean_slip*3, t='l', col='orange')
    points(candidate_events$events$Mw, strasser_mean_slip*6, t='l', col='red')
    legend('topleft', c('Scaling-relation mean-slip', ' "" times 3', ' "" times 6'), 
           col=c('blue', 'orange', 'red'), bty='n',
           lty=c(1,1,1))
    add_log_axis_ticks(side=2)

    # Mw vs energy
    plot(candidate_events$events$Mw, 
         candidate_events$events$initial_potential_energy,
         log='y', main="Mw vs initial potential energy \n (Historical event energies are approximate!)", 
         xlab='Mw', ylab='Potential Energy (joules)',
         cex.main=1.5, cex.lab=1.3, cex.axis=1.3)
    grid(col='orange')
    min_mw = min(candidate_events$events$Mw)
    energies = c(1e+15, 3e+15, 4.5e+15, 6.5e+15, 1e+16, 4e+16)
    linecol = c('blue', 'green', 'violet', 'purple', 'red', 'black')
    abline(h=energies, col=linecol)
    labels= c('Chile 2010', 'Tohoku 2011', 'Alaska 1964', 
             'Sumatra2004', 'Chile 1960', '2x Chile 1960 wave size')
    text(min_mw+energies*0, energies, labels, adj=c(0, 0.), col=linecol)
    add_log_axis_ticks(side=2)

    # Weighted Cumulative Mw
    weighted_Mw_distribution = weighted_cdf(candidate_events$events$Mw, 
        candidate_events$events$rate_annual, make_plot=FALSE)
    plot(weighted_Mw_distribution$x, weighted_Mw_distribution$fraction_less_or_equal, t='s',
         xlab='Mw', ylab='Cumulative distribution', ylim=c(0, 1), cex.lab=1.3, cex.axis=1.3)
    points(weighted_Mw_distribution$x, weighted_Mw_distribution$fraction_less_or_equal, pch=19)
    grid(col='orange')
    abline(h=0.5, col='red')
    title(main='Mw cumulative distribution for selected \n scenarios, weighted by annual rate',
          cex.main=1.5)

    # Look at the peak-stage for each event
    boxplot(candidate_events$peak_stages, 
            las=2, cex.axis=1.3)
    title(main=paste0('Peak stages at all gauges for all events \n',
                      'Red points give the target exceedance-rate peak-stage'),
          cex.main=1.7)

    points(candidate_events$max_stage_at_exrate, col='red', pch=15, cex=2)
    grid(col='orange')

    image(candidate_events$peak_stages_within_tol, axes=FALSE, 
          col=c('white', 'red'), zlim=c(0, 1))
    nc = ncol(candidate_events$peak_stages_within_tol)
    nr = nrow(candidate_events$peak_stages_within_tol)
    axis(side=2, at=seq(0, 1, length=nc), 
         labels=colnames(candidate_events$peak_stages), las=1)
    axis(side=1, at=seq(0, 1, length=nr), labels=1:nr)
    title('Are gauges within the target range for each event? Red = yes \n Line gives scenario Mw [RHS axis]', 
          cex.main=1.7)

    # Add magnitudes
    mw_seq = seq(min(candidate_events$events$Mw)-0.1/2, max(candidate_events$events$Mw)+0.1/2, by=0.1)
    mw_line = sapply(mw_seq, f<-function(x) sum(candidate_events$events$Mw <= x))

    min_mw = min(candidate_events$events$Mw)
    max_mw = max(candidate_events$events$Mw)
    points((0:(nr-1))/(nr-1), 
           (candidate_events$events$Mw - min_mw)/(max_mw - min_mw),
           t='o', lty='dashed')
    axis_pts = pretty(c(min_mw, max_mw), n=4)
    axis(side=4, at=(axis_pts-min_mw)/(max_mw-min_mw), labels=axis_pts)
}

# This is also worth looking at -- plot the peak stages for each event/gauge
image_plot_peak_stages<-function(candidate_events){
    library(fields)
    image.plot(candidate_events$peak_stages, axes=FALSE)
    nc = ncol(candidate_events$peak_stages_within_tol)
    nr = nrow(candidate_events$peak_stages_within_tol)
    axis(side=2, at=seq(0, 1, length=nc), 
         labels=colnames(candidate_events$peak_stages), las=1)
    axis(side=1, at=seq(0, 1, length=nr), labels=1:nr)
    title('Peak stage at each gauge for each scenario', cex.main=1.7)
}

weighted_cdf<-function(values, weights, make_plot=TRUE){

    # Sort the values, and store the indices 
    sorted_values = sort(values, index.return=TRUE)
    sorted_weights = weights[sorted_values$ix]
    cumulative_normalised_weights = cumsum(sorted_weights)/sum(weights)

    if(make_plot){
        plot(sorted_values$x, cumulative_normalised_weights, t='s',
             xlab='', ylab='', ylim=c(0, 1))
        points(sorted_values$x, cumulative_normalised_weights, pch=19)
        grid(col='orange')
        mtext(side=1, 'Value', line=2, cex=1.5)
        mtext(side=2, 'Fraction of weight on values less-than-or-equal', 
              line=2, cex=1.5)
    }

    return(data.frame(x=sorted_values$x, fraction_less_or_equal=cumulative_normalised_weights))
}


# This function gives an example of how to use the scripts.
# For software maintanance purposes I include the option to run it
# as a regression test (i.e. check whether the results have changed,
# given the hard-coded parameters).
example_usage<-function(run_as_regression_test=FALSE){
    #
    # Main program below here
    #

    #
    # Define the gauge indices of interest -- you need to make BOTH "site_gauge_inds" and "site_gauges"
    #
    all_gauges = ptha18$get_all_gauges()
    site_gauge_inds = which(all_gauges$lon > 186 & all_gauges$lon < 188.9 &
                            all_gauges$lat > -15 & all_gauges$lat < -12)
    site_gauges = all_gauges[site_gauge_inds,]
    # In this case we have a double-up of one gauge. This is not unusual in
    # PTHA18. For the analysis here we don't want that because it will make 1
    # site count as 2. 
    # Let's remove it by computing a distance matrix, and removing points that
    # have '0' IN THE UPPER TRIANGLE of the distance matrix. Restriction to the
    # upper triangle will mean we don't remove both points
    site_gauges_distm = distm(cbind(site_gauges$lon, site_gauges$lat))
    to_remove = which((upper.tri(site_gauges_distm) & (site_gauges_distm == 0)), arr.ind=TRUE)[,1]
    site_gauge_inds = site_gauge_inds[-to_remove]
    site_gauges = site_gauges[-to_remove,]

    ## For test
    #site_gauges = site_gauges[1,,drop=FALSE]
    #site_gauge_inds = site_gauge_inds[1]

    # 
    # Define the source-zone and slip type. This will be used to download
    # background data for scenario selection.
    #
    # Slip type -- either 'stochastic' or 'variable_uniform'
    slip_type = 'variable_uniform' 
    # Name of source-zone in PTHA18
    source_zone = 'kermadectonga2'

    # Download the hazard curves at gauges, and the peak-stage values for all
    # gauges at the source-zone
    kermadec_source_data = get_hazard_curves_and_peak_stages(
        site_gauge_inds, site_gauges, source_zone, slip_type)


    #
    # Set the criteria for picking scenarios. 
    scenario_match = list()
    # Desired exceedance-rate
    scenario_match$exrate = 1/2500
    # Hazard curve mean or percentile type -- a string identifying the stage-vs-exceedance-rate curve
    # type. Either 'rate' (mean) or 'rate_84pc' (84th percentile) or 'rate_16pc'
    # (16th percentile), 'rate_median' (50th percentile), 'rate_lower_ci' (2.5th
    # percentile), 'rate_upper_ci' (97.5 percentile)
    scenario_match$hazard_curve_type = 'rate_median' # 'rate_84pc' # 'rate_84pc'
    # Rigidity type used for hazard calculations -- either an empty string ''
    # (constant rigidity) or 'variable_mu' (depth varying rigidity). If you use
    # variable_mu, be sure to use the variable_mu_Mw from the event table
    scenario_match$mu_type = '' # 'variable_mu'
    # The 'allowed fractional deviation of the peak-stage from the target value'
    # that is still regarded as a match, e.g. 0.1 = 10%
    scenario_match$target_stage_tolerance = 0.1 #
    # The minimum number of gauges that should satisfy the allowed tolerance, for the
    # scenario to be stored. At all other gauges, the peak-stage must be
    # less than or equal the exrate-derived value
    scenario_match$number_matching_gauges = 3
    #scenario_match$number_matching_gauges = 1

    # Find scenarios which match the above criteria at the previously downloaded gauge+"source-zone"
    kermadectonga_events_2500 = get_matching_scenarios(kermadec_source_data, scenario_match)
    
    # Get summary information
    summarise_scenarios(kermadectonga_events_2500)

    # Have a look at the weighted Mw distribution.
    par(mfrow=c(1,1))
    weighted_Mw_distrbution = weighted_cdf(kermadectonga_events_2500$events$Mw, 
        kermadectonga_events_2500$events$rate_annual)
    title(main='Mw distribution for selected scenarios, weighted by annual rate')

    # This can be used to get some subset of the events
    kermadectonga_events_2500$desired_event_rows

    # For instance if we wanted the earthquake data for 5th entry, in a format
    # easily compatable with our the PTHA18 DETAILED_README, we could do:
    event_id = kermadectonga_events_2500$desired_event_rows[5]

    # Get the event data directly from the PTHA database (easier than hacking it
    # out of the above objects).
    event_of_interest = ptha18$get_source_zone_events_data(
        source_zone, slip_type=slip_type, desired_event_rows=event_id)

    # Get the initial deformation, using the regular approach.
    event_raster = ptha18$get_initial_condition_for_event(
        event_of_interest, event_ID=1)

    par(mfrow=c(1,1))
    plot(event_raster, main='Initial condition for the event of interest')

    if(run_as_regression_test){

        # Let's check that we got the same results as before

        previous_peak_slip = c(18.11, 16.42, 14.25, 16.8, 15.44, 13.8, 13.8,
                               15.28, 21.08, 24.3, 22.29, 22.29, 23.24, 33.5,
                               15.59, 32.99, 15.59, 15.59, 15.59, 15.59, 24.42,
                               42.21, 81.58, 94.9, 56.01, 33.69, 21.55, 22.43,
                               117.6, 66.41, 89.27, 60.43, 73.14)

        if(all((round(kermadectonga_events_2500$peak_slip, 2) - previous_peak_slip) == 0)){
            print('PASS')
        }else{
            print(paste0('FAIL (this just means the selected sceanrios have changed). \n',
                         'This might be because you changed the hard coded parameters in \n', 
                         'this function, or because of other changes to the code or the database.'))
        }

        test_stat = range(as.matrix(event_raster))
        previous_test_stat = c( -1.813732,  4.999213)
        if(all(abs(test_stat - previous_test_stat) < 1.0e-05)){
            print('PASS')
        }else{
            print( paste0('FAIL (the results have changed). This might be because \n',
                          'the scenarios have changed (in which case it should have \n',
                          'been flagged above), or because of changes to the event \n',
                          'raster reconstruction.') )
        }
    }

}

test_select_scenarios<-function(){
    # Convenience function for running the test
    example_usage(run_as_regression_test=TRUE)
}
