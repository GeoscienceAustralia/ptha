#
# Variables that are useful in multiple different scripts
#


# HUNGA-TONGA volcano location
hunga_tonga_volcano = c(360 -175.385, -20.55)

# A very approximate explosion start time
explosion_start_time = strptime('2022-01-15 04:15:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

# Location for output directory, relative to script directory
OUTPUT_DIR = '../post-processed/'

# Location for figures in output directory
OUTPUT_GRAPHICS_DIR = paste0(OUTPUT_DIR, '03_graphical_checks/')

# Location for post-processed tidal time-series, relative to script directory.
OUTPUT_TIDE_DIR = paste0(OUTPUT_DIR, '01_tide_gauges/')

# Location for post-processed MSLP data, relative to script directory
OUTPUT_MSLP_DIR = paste0(OUTPUT_DIR, '02_mslp_sensors/')

# Location of station metadata files, csv format, stored inside OUTPUT_DIR.
MSLP_METADATA_TABLE_FILE = paste0(OUTPUT_DIR, '02_mslp_sensor_locations.csv')
TIDEGAUGE_METADATA_TABLE_FILE = paste0(OUTPUT_DIR, '01_tide_gauge_locations.csv')

# Location of files naming the gauges that were skipped
IGNORED_MSLP_FILE = paste0(OUTPUT_DIR, 'ignored_mslp_sensors.txt')
IGNORED_TIDEGAUGE_FILE = paste0(OUTPUT_DIR, 'ignored_tide_gauges.txt')

# Approximate speed of Lamb wave (used for plotting to estimate theoretical arrivals)
LAMB_WAVE_SPEED = 320 # m/s

# Approximate earth radius (used for plotting to estimate theoretical arrivals)
EARTH_RADIUS = 6378137 # m -- consistent with geosphere::distHaversine

# Only process the gauge data between these start and end times. This was done for speed,
# and the start time also helps remove an artefact in several of the BomPorts stations
START_TIME_LIMIT_GAUGE_DATA = strptime('2022/01/02 00:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC')
END_TIME_LIMIT_GAUGE_DATA = strptime('2022/02/01 00:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC')

# Tsunami travel time from the TTT software, assuming ocean gravity wave propagation from the source.
# The "5.851" hours was provided by Diana Greenslade using TTT.
# GEOWARE (2007), TTT - A tsunami travel-time calculator, http://www.geoware-online.com
CROWDY_HEAD_TTT_ARRIVAL_TIME = explosion_start_time + as.difftime(5.851, units='hours')
