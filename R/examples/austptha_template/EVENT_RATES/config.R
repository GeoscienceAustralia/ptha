#
# Set key user-defined parameters in this file
#


# Set this to FALSE for test/plotting based runs where I don't want to touch
# the main netcdf files
write_to_netcdf = TRUE

#######################################################################################
#
# PARAMETERS CONTROLLING SELECTION OF EVENTS FROM GLOBAL CMT CATALOGUE
#
#######################################################################################
#cmt_catalogue_csv = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/EARTHQUAKE/Earthquake_catalogues/Global_CMT_catalogue/GCMT_1976_2017.csv'
cmt_catalogue_csv = '../../../DATA/EARTHQUAKE/Earthquake_catalogues/Global_CMT_catalogue/GCMT_1976_2017.csv'

# Only use events with depth < depth_threshold
depth_threshold = 71 # km

# Events are treated as inside the source-zone, if either the GCMT hypocentre or GCMT centroid
# is inside the source-zone-polygon, after the polygon is buffered by 'buffer_width' degrees
buffer_width = 0.4

#
# See also 'rake_deviation' and 'strike_deviation' below
#

#isc_gem_catalogue_csv = "/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/EARTHQUAKE/Earthquake_catalogues/isc-gem-catalogue-version4/parsed_gem_catalogue.csv"

#######################################################################################
#
# PARAMETERS CONTROLLING MODELLED RATE OF EARTHQUAKES ON EACH SOURCE-ZONE
#
#######################################################################################

#
# Path to a zipped version of the file derived largely from the
# 'PB2002_steps.dat.txt' from the supplementary material of Peter Bird's 2003
# paper:
#   Bird, P. An updated digital model of plate boundaries. 
#   Geochemistry Geophysics Geosystems, 2003, 4, 1-52
# The file is also derived from similar data based on Jonathan Griffin's work for
# this project (mainly around Indonesia). 
#
# This is used to define spatially variable convergence rates, for source-zones where
# the input parameter file indicates that Bird's model should be used.
#
bird2003_steps_data_zip_file = '../DATA/BIRD_PLATE_BOUNDARIES/sourcezone_traces_table_merged.csv.zip'


#
# When looking at historical earthquake data, interpret 'thrust' (or normal)
# earthquakes as those whose rake deviates from 90 degrees (resp -90) by less
# than this number [e.g. a value of 45 means events within [90-45, 90+45] are
# treated as thrust, and [-90-45, -90+45] are normal]. Likewise,
# treat tectonic convergence in a direction within this much of 'pure thrust' as 
# contributing to the seismic moment which is balanced by earthquakes. 
#
rake_deviation = 50 # Degrees

# As above for strike
strike_deviation = 50 # Degrees

# Should we use a uniform prior for coupling?  (i.e. same at every
# source-zone). Or the input spreadsheet values? Or 50% weight on each? 
coupling_prior_type = 'spreadsheet_and_uniform_50_50' # 'uniform' # 'spreadsheet'
uniform_coupling_prior = c(0.1, 1.3)

# Should we use GCMT to update logic tree weights
update_logic_tree_weights_with_data = TRUE

# Should we apply an edge correction to the conditional probabilities of events?
# This tries to avoid under-estimation of the moment release rate near the edges
# of source-zones, which occurs with methods that do not account for this boundary
edge_correct_event_rates = TRUE

#
# File with parameters for all source-zones, that we use to make source-zone
# specific rate curves.
#
sourcezone_parameter_file = '../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv'

# Increment between Mw values in the earthquake_events table. We will check
# that the table holds the same value
# FIXME: This value must be identical to the value in
# ../SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/config.R
# and correspond to the spacing between alternative Mw values in the earthquake
# event table. It should not have more than 3 figures after the decimal point.
dMw = 0.1

#
# Bounds on the Mw-max parameter used for each source zone.
# For example, Berryman et al (2015) suggest Mw_max should be <= 9.6 everywhere.
# No advice is given for normal faults, however, the largest observed normal
# fault event was Mw 8.5 (Sanriku 1933, see Okal et al., 2016), vs 9.5 for
# thrust (Chile). Doglioni et al. (2015) suggest that all else being equal,
# normal faults generate smaller earthquakes than strike-slip or thrust faults.
# Burbidge et al (2008) used a maximum Mw-max of 9.0 for normal fault sources
# in their analysis.
#
# The minimum_allowed_Mw_max value is not physical, but it should be greater than 
# the smallest earthquakes we simulate on all source-zones, to ensure we have
# at least some earthquake events with non-zero probability!
#
MAXIMUM_ALLOWED_MW_MAX = 9.6
MINIMUM_ALLOWED_MW_MAX = 7.35
MAXIMUM_ALLOWED_MW_MAX_NORMAL = 9.0

#
# Only assign non-zero probability to earthquakes with Mw greater than MW_MIN.
# MW_MIN should be smaller than the largest earthquake in the event set (or they
# will all be assigned a zero rate!). It may also be smaller than the smallest
# Mw event in the event set -- that won't create problems, and can make it
# easier to compare against data for small Mw events. 
MW_MIN = 7.2 - dMw/2

#
# We ensure that (Mw_max >= maximum_observed_mw + mw_observed_perturbation)
# This ensures that no logic-tree curve assigns zero probability to the largest
# observed event
#
mw_observed_perturbation = 0.05

# Weights for both 'truncated' and 'characteristic' Gutenberg Richter model
# These are applied to every source-zone
Mw_frequency_distribution_types = c('truncated_gutenberg_richter', 'characteristic_gutenberg_richter')
Mw_frequency_distribution_weights = c(0.7, 0.3) # Must sum to give 1.0

#
# Interpolate logic-tree parameter variation over this many values, all assumed to have
# equal rate. For example, if we provide source_coupling = c(0.1, 0.2, 0.7), then the
# actual coupling values will be computed as
# >   k = logic_tree_parameter_subsampling_factor*coupling_subsampling_increase
# >   approx(source_coupling, n=k)$y
#
# Beware the computational time scales with a high power of this number, so large values will
# slow the computation greatly.
#
logic_tree_parameter_subsampling_factor = 20
# Use these parameters to increase the number of bins for particular variables
# e.g. a value of 2 would lead to (2*logic_tree_parameter_subsampling_factor) bins
# for the given variable. This can be useful because you may need to use a finer
# discretization for some variables than others.
coupling_subsampling_increase = 1
mwmax_subsampling_increase = 2
b_subsampling_increase = 1

# Inverse quantiles used for lower and upper credible intervals
lower_ci_inv_quantile = 0.025
upper_ci_inv_quantile = 0.975


# Name for a directory where we write log files for each source
sourcezone_log_directory  = 'source_zone_parameter_logs'

#######################################################################################
#
# PARAMETERS CONTROLLING THE STAGE EXCEEDANCE RATE CURVE COMPUTATION AT EVERY STATION
#
#######################################################################################

#
# Apply hazard curve computation / data extraction to chunks of data
# We read in (number_of_events x point_chunk_size) max stage values at once, and then
# compute rate curves for the point_chunk_size gauges before moving to the next
# chunk. 
#
# Small values conserve memory -- large values may be faster (up to a point) but use more
# memory.
#
# Note these values will be scaled by the maximum number of events in any source-zone
# (e.g. southamerica in the 2018 PTHA)
point_chunk_size_uniform = 10000
point_chunk_size_stochastic = 500

# 
# Make 'stage_seq', and ordered set of stage values for which we compute the
# stage exceedance rate at every hazard point. [Later, these values can be used
# to define the curve, e.g. using interpolation]
#
stage_seq_len = 100 # Number of points defining rate curve
stage_seq_max = 20 # Max stage on rate curve, m
stage_seq_min = 0.02 # Min stage on rate curve m
# Stage points have logarithmic spacing
stage_seq = exp(seq(log(stage_seq_min), log(stage_seq_max), len=stage_seq_len))

# Number of cores to use in shared memory parallel
MC_CORES = 16


###########################################################################
#
# Functions to weight individual stochastic and variable_uniform events that
# have the same corresponding uniform event. 
# This allows us to attempt to address biases, such as (hypothetically)
# variable_uniform_events being too low-slip on average, or stochastic_slip
# events being too high-slip on average.
# We do not try to adjust determinstic-size-uniform-slip events
#
# If you do not want to perform adjustment, replace these functions with 
# constant functions (e.g. peak_slip_bias_adjustment_XXXX<-function(x) x=1)
# Currently we do this replacement for normal-fault events (inside compute_rates_all_sources.R)
#
###########################################################################

# Exclude events with {peak_slip > peak_slip_limit_factor*mean-scaling-relation-slip}
# Note we use the 'reference Mw' (i.e. constant shear modulus) for doing this
#peak_slip_limit_factor = 7.5 # Inf
source('config_peak_slip_limit_factor.R', local=TRUE)


peak_slip_bias_adjustment_table = read.csv('peak_slip_quantile_adjustment_factors.csv')
peak_slip_bias_adjustment_stochastic = approxfun(
    peak_slip_bias_adjustment_table$peak_slip_quantile,
    peak_slip_bias_adjustment_table$weight_stochastic)
peak_slip_bias_adjustment_variable_uniform = approxfun(
    peak_slip_bias_adjustment_table$peak_slip_quantile,
    peak_slip_bias_adjustment_table$weight_variable_uniform)
# Variable shear modulus case
peak_slip_bias_adjustment_table_mu_vary = read.csv('peak_slip_quantile_adjustment_factors_varyMu.csv')
peak_slip_bias_adjustment_stochastic_mu_vary = approxfun(
    peak_slip_bias_adjustment_table_mu_vary$peak_slip_quantile,
    peak_slip_bias_adjustment_table_mu_vary$weight_stochastic)
peak_slip_bias_adjustment_variable_uniform_mu_vary = approxfun(
    peak_slip_bias_adjustment_table_mu_vary$peak_slip_quantile,
    peak_slip_bias_adjustment_table_mu_vary$weight_variable_uniform)


#######################################################################################
#
# PARAMETERS POINTING TO PATHS OF FILES ALREADY CREATED IN THE TSUNAMI_EVENTS PHASE
# 
# IF YOU ARE USING THE STANDARD FILE STRUCTURE, THESE PARAMETERS PROBABLY DO NOT NEED
# TO BE CHANGED
#
#######################################################################################

#
# Vector of paths to all unit_source_statistics_SOURCENAME.nc files (one for each source-zone).
# [e.g. unit_source_statistics_puysegur.nc, unit_source_statistics_kermadectonga.nc, ...]. 
#
# These should have been created by running codes in
# ../SOURCE_ZONES/sourcename/TSUNAMI_EVENTS/
#
# 
unit_source_statistics_netcdf_files = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source_statistics*.nc')

#
# Paths to shapefiles containing unit-source grid (one for each source-zone). 
#
unit_source_grid_polygon_shapefiles = Sys.glob(
    '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')
    # '../SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')
    #'/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA/SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')

#
# NetCDF files with uniform slip max_stage for every hazard point, and also a place for rates for every event
#
all_source_uniform_slip_tsunami = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_*.nc')

#
# NetCDF files with stochastic slip max_stage for every hazard point, and also a place rates for every event
#
all_source_stochastic_slip_tsunami = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_*.nc')


#
# NetCDF files with stochastic slip max_stage for every hazard point, and also a place rates for every event
#
all_source_variable_uniform_slip_tsunami = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_variable_uniform_slip_earthquake_events_tsunami_*.nc')


##### END INPUTS




################################################################################
#
# CONSISTENCY CHECKS BELOW HERE
#
################################################################################


#
# Here we check that the parameters exist
# 
stopifnot( all.equal( sum(Mw_frequency_distribution_weights),  1) )
stopifnot(file.exists(sourcezone_parameter_file))
stopifnot(all(file.exists(unit_source_statistics_netcdf_files)))
stopifnot(file.exists(bird2003_steps_data_zip_file))

#
# Check we use enough logic tree branches
#
if(logic_tree_parameter_subsampling_factor < 3){
    stop(
        paste0('The logic_tree_parameter_subsampling_factor should not be less than \n',
        'the number of parameter values for b provided in the \n',
        ' source-zone parameter file. Normally we provide 3 values for each, so this error is \n',
        ' produced if logic_tree_parameter_subsampling_factor < 3') 
    )
}

#
# Here we check for consistency among filenames.
#
source_names_1 = basename(dirname(dirname(unit_source_statistics_netcdf_files)))
source_names_2 = basename(dirname(dirname(dirname(unit_source_grid_polygon_shapefiles))))
source_names_3 = basename(dirname(dirname(all_source_uniform_slip_tsunami)))
source_names_4 = basename(dirname(dirname(all_source_stochastic_slip_tsunami)))
source_names_5 = basename(dirname(dirname(all_source_variable_uniform_slip_tsunami)))

if(!all(source_names_1 == source_names_2)){
    stop('unit_source_statistics_netcdf_files and unit_source_grid_polygon_shapefiles must refer to the same source-zones, in the same order. Ensure the files have a natural pairing. If they do, but their sort order is inconsistent, it might be fixed by setting an environment variable (export LC_ALL=C) before running')
}
if(!all(source_names_1 == source_names_3)){
    stop('unit_source_statistics_netcdf_files and all_source_uniform_slip_tsunami must refer to the same source-zones, in the same order. Ensure the files have a natural pairing. If their sort order is inconsistent, it might be fixed by setting an environment variable (export LC_ALL=C) before running')
}
if(!all(source_names_1 == source_names_4)){
    stop('unit_source_statistics_netcdf_files and all_source_stochastic_slip_tsunami must refer to the same source-zones, in the same order. Ensure the files have a natural pairing. If their sort order is inconsistent, it might be fixed by setting an environment variable (export LC_ALL=C) before running')
}
if(!all(source_names_1 == source_names_5)){
    stop('unit_source_statistics_netcdf_files and all_source_variable_uniform_slip_tsunami must refer to the same source-zones, in the same order. Ensure the files have a natural pairing. If their sort order is inconsistent, it might be fixed by setting an environment variable (export LC_ALL=C) before running')
}

