#
# Example analysis of rates of earthquake events on the Alaska source-zone.
#
# Assumes the code to produce unit sources for alaska in 
#     ptha/R/source_contours_unit_sources/produce_unit_sources.R
# has been run.
#

library(rptha)

############################################################
#
# INPUT PARAMETERS
#
############################################################

# Name of the source. This must correspond with the name of the corresponding
# input source contours shapefile that was used to make the unit sources, e.g.
# if the latter was 'CONTOURS/manila.shp' then this should be 'manila'
sourcename = 'alaska' 

# Read the discretized source info which was created when making the unit
# sources. This list should have an element with name = sourcename
all_discretized_sources = readRDS(
    '../source_contours_2_unit_sources/OUTPUTS/all_discretized_sources.RDS')

# Define the (m/year) SEISMICALLY COUPLED tectonic convergence on the sourcezone.
# We use same values as global PTHA paper, and account for the interface dip
# deformation_along_interface = (horizontal_deformation/cos(interface_dip_radians)
dip_angle = 9.59 # Mean value -- can extract from unit source summary statistics once they are computed.
slip_rate = (1/1000) * c(40, 45, 50) * 1/cos(dip_angle/180*pi)
slip_rate_prob = c(1/3, 1/3, 1/3)

# Gutenberg Richter b parameter -- same values as global PTHA paper,
# which are based on Berryman et al (2015)
b = c(0.7, 0.95, 1.2)
b_prob = c(1/3, 1/3, 1/3)

# Smallest earthquake that we model. 
# It's not really natural to model this probabilistically as it is more a
# choice of the analyst -- so we just use a single value.
Mw_min = 7.5 
Mw_min_prob = 1

# Maximum possible magnitude on the source zone. Same values as global PTHA
# paper
Mw_max = c(9.3, 9.3, 9.4)
Mw_max_prob = c(0.1, 0.45, 0.45)

# Choice of Mw frequency distributions
Mw_frequency_dists = c("truncated_gutenberg_richter", "characteristic_gutenberg_richter")
Mw_frequency_dists_prob = c(0.5, 0.5)

# Mw increment. This gives the spacing between Mw values in the synthetic event
# set. It should exactly divide max(Mw_max) - min(Mw_min).
dMw = 0.05

# Data to update logic tree weights
# Values are:
#     c("Threshold magnitude", 
#       "Number of events >= threshold magnitude", 
#        "Number of years of observations")
# Use c(NA, NA, NA) to not update logic tree weights with data
#
# (We use mumber of Mw > Mw_min thrust events observed in the CMT catalogue
# (1976-2013))
Mw_count_durations_CMT = c(7.5, 0, 38)

###############################################################
#
# END INPUT PARAMETERS
#
################################################################


# Get a table with unit-source summary statistics
discrete_source_summary_statistics = 
    discretized_source_approximate_summary_statistics(
        all_discretized_sources[[sourcename]])

# Alternative version which uses sub-unit-source points for more accurate
# statistics
# discrete_source_summary_statistics = 
#    discretized_source_summary_statistics(
#        all_discretized_sources[[sourcename]],
#        approx_dx=4000, approx_dy=4000)

# Get all 'uniform slip' earthquake events from Mw 7.5 to Mw 9.6,
# with an Mw increment of 0.1 
#
# Note: We can include 'unrealistically large' events in the table, and later
# assigned them probability zero, following the global PTHA methodology.
#
all_eq_events = get_all_earthquake_events(
    unit_source_statistics = discrete_source_summary_statistics, 
    Mmin = min(Mw_min), 
    Mmax = max(Mw_max), 
    dMw = dMw)

# Source area in km^2
sourcezone_area = sum(discrete_source_summary_statistics$length * 
    discrete_source_summary_statistics$width)

# Conditional probabilities of events for each Mw. Here we ignore variations in
# tectonic rates but let the event rate scale inversely with its slip so that
# all events contribute the same to the convergence. 
# More complex conditional rate models can be implemented by passing a function
# as the 'conditional_probability_model', see:
# ?get_event_probabilities_conditional_on_Mw
event_conditional_prob = get_event_probabilities_conditional_on_Mw(
    all_eq_events, conditional_probability_model = 'inverse_slip')

#
# Get function f(Mw) giving the rate of events >= Mw on the whole sourcezone
#
event_rates_function = rate_of_earthquakes_greater_than_Mw_function(
    slip_rate, slip_rate_prob, 
    b, b_prob, 
    Mw_min, Mw_min_prob, 
    Mw_max, Mw_max_prob, 
    sourcezone_total_area = sourcezone_area,
    event_table = all_eq_events, 
    event_conditional_probabilities = event_conditional_prob,
    Mw_frequency_distribution = Mw_frequency_dists,
    Mw_frequency_distribution_prob = Mw_frequency_dists_prob,
    update_logic_tree_weights_with_data = (!any(is.na(Mw_count_durations_CMT))), 
    Mw_count_duration = Mw_count_durations_CMT, 
    account_for_moment_below_mwmin = TRUE,
    mw_max_posterior_equals_mw_max_prior=FALSE)

# Example: Compute rate of events exceeding Mw 9.0. Should be around 1/540.
# Small variations can arise depending on the unit-source discretization
# and how the unit source summary statistics are computed, but they should
# not be significant.
rate_of_events_gt_9 = event_rates_function(9.0)


#
# Get the (nominal) rates of individual events. 
# Note these are completely determined by the discretization -- we just
# spread the total rate in each Mw bin over all the events with that Mw.
#
all_eq_event_rates = event_conditional_prob * (
        event_rates_function(all_eq_events$Mw - dMw/2) -
        event_rates_function(all_eq_events$Mw + dMw/2))

#
# Make a plot. This code also illustrates how to use the event_rates_function
#
Mws = seq(7.5, 9.4, len=100)
plot(Mws, event_rates_function(Mws), t='l', log='y', xlab='Mw', 
    ylab='Exceedance Rate (events/year)', lwd=5,
    main='Alaska source zone Mw-exceedance-rate curve')
grid(col='brown')

#
# We can use event_rates_function to access curves from all logic trees,
# by passing the argument "return_all_logic_tree_branches=TRUE".
#
# This returns a list with detailed information on each curve
#
all_curves = event_rates_function(return_all_logic_tree_branches=TRUE)
for(i in 1:nrow(all_curves$all_rate_matrix)){
    points(all_curves$Mw_seq, all_curves$all_rate_matrix[i,], col='grey', t='l')
}

legend('topright', 
    c('Weighted mean rate curve (updated prior weights)', 'Logic tree branch'),
    lwd=c(5,1), col=c('black', 'grey'))
