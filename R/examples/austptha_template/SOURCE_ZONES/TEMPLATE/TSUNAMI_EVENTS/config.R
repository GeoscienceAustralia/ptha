# Key variables required by scripts in this folder

############################################################################
#
# Parameters determining the event magnitudes modelled / number of events, etc.
#
############################################################################

# Make events with magnitudes ranging from Mw_min to Mw_max, separated by dMw
#
# These values do not have to correspond to Mw ranges that are assigned non-zero 
# probability -- but they should fully contain the possible values.
#
# For instance, later on we may decide to assign zero probablity to (e.g.) 
#    Mw >= (Mw_max - 0.8)
# However, we cannot later decide to assign non-zero probability to events
# with Mw >= Mw_max, or events with Mw <= Mw_min, since we will not have
# modelled them!
Mw_min = 7.2
Mw_max = 9.8
dMw = 0.1

# Make at least this many stochastic slip events for each uniform event
number_stochastic_events_for_each_uniform_event = 10
# ... but ensure that there are at least this many stochastic slip events in
# each magnitude category
minimum_number_stochastic_events_for_each_magnitude = 100

##############################################################################
#
# Parameters controlling unit-source summation
#
##############################################################################

# Memory we will allow in unit-sources, in MB. 
# It's a rough estimator -- should be a fraction of 1-nodes memory [e.g. using
# 10 GB on a 32GB machine seemed to work ok]
memory_for_unit_sources = 10 * 1024 

# Number of cores in parallel [shared memory only]
mc_cores = 16

# Pre-tsunami stage initial condition for linear model. Gauges with (elevation > msl)
# are inactive
msl = 0.0 
# Latitude range where gauges should be taken. Used to chop out gauges in the
# North/South boundary condition zone [ 2 cells for our model ]
lat_range = c(-72 + 2/60, 65 - 2/60) 

# Define the wave arrival time as the time at which abs(stage) exceeds this
# threshold in meters. We cannot just use 'zero' because of tiny round-off in
# flow solver outputs. 
stage_threshold_for_arrival_time = 0.001 


#################################################################################
#
# Parameters controlling some output-file cludges we do.
#
#################################################################################

# This is used to represent missing data in some places. It should be a
# 'large-negative' number that is not an integer (to make sure netcdf
# interprets it as floating point real)
null_double = -999.999

# Determine the number of significant figures used when storing stochastic slip
# values. 
#
# We store the stochastic slip events in a table, with one column being a character
# 'event_slip_string' containing the slip values as a string 'slip1_slip2_slip3_slip4_...'.
#
# These slip values which correspond to the slip on unit-sources defined by 
# 'event_index_string' which has a format like 'i1-i2-i3-i4-i5-...'
#
# (The strange character-string format is a cludge to help store unstructured
# data in a table).
#
# By reducing the number of significant figures used to store the slip, we reduce the
# netCDF output file size. However, it should not be too small to represent slip values
# well. 4 seems a good choice [giving O(50m) of slip to 2 decimal places]
stochastic_slip_table_significant_figures = 4

# Temporary directory used when running 'unfinished' netcdf
tmp_RDS_dir = 'R_images_tmp'
