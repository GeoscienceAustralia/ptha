#
# Tide gauge Model vs obs comparison
#

#################################################################################################
## INPUT ARGS 

# Earthquake GMT time
model_start_time = strptime('2006-07-17 08:19:27', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

# event ID used in interface with tidal-gauge data
event_id = '2006-07-17' 

# Pass the model multidomain_dir as a commandline argument
multidomain_dir = commandArgs(trailingOnly=TRUE)[1]

## END INPUT 
#################################################################################################
source('plot_gauges_generic_include.R')
