#
# Tide gauge Model vs obs comparison
#

#################################################################################################
## INPUT ARGS 

# Earthquake GMT time
model_start_time = strptime('2010-02-27 06:34:15', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

# event ID used in interface with tidal-gauge data
event_id = '2010-02-27' 

# Pass the model multidomain_dir as a commandline argument
multidomain_dir = commandArgs(trailingOnly=TRUE)[1]

# Maximum elevation for gauges that are extracted.
# Sometimes we have tide gauges at locations that the model considers to be on land
# (due to limited resolution), but in that case we store nearby gauges in wet areas.
# This option can ensure that those gauges are selected
maximum_elevation_of_gauges = -1

## END INPUT 
#################################################################################################
source('_gauges_generic_include.R')

