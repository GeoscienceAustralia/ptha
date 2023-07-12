#
# Compute various summary statistics about the data, such as the station sampling frequency,
# but now also many other things.
#

source('global_variables.R')
source('parse_gauge_data.R')

MSLP_metadata = read.csv(MSLP_METADATA_TABLE_FILE)
TG_metadata = read.csv(TIDEGAUGE_METADATA_TABLE_FILE)

# Read all the MSLP and tide-gauge time-series
mslp_data = lapply(MSLP_metadata$postprocessed_file, function(x) read.csv(x, comment.char='#'))
names(mslp_data) = basename(MSLP_metadata$postprocessed_file)
length(mslp_data)

tg_data = lapply(TG_metadata$postprocessed_file, function(x) read.csv(x, comment.char='#'))
names(tg_data) = basename(TG_metadata$postprocessed_file)
length(tg_data)

# Counts by data source
table(MSLP_metadata$Dataset)
table(TG_metadata$Dataset)

# Typical data increment in seconds
mslp_spacing = lapply(mslp_data, function(x) median(diff(x$juliant)*24*3600, na.rm=TRUE))
range(unlist(mslp_spacing))

# Typical data increment in seconds
tg_spacing = lapply(tg_data, function(x) median(diff(x$juliant)*24*3600, na.rm=TRUE))
table(round(unlist(tg_spacing)))

# Gauges with NA values for the tide
tg_fraction_NA = lapply(tg_data, function(x) mean(is.na(x$stage)))
table(unlist(tg_fraction_NA))

# Duration in days
tg_duration = lapply(tg_data, function(x) diff(range(x$juliant, na.rm=TRUE)))
stem(unlist(tg_duration))

sum(round(unlist(tg_duration)) >= 30)
sum(round(unlist(tg_duration)) >= 20)
sum(round(unlist(tg_duration)) >= 10)


mslp_duration = lapply(mslp_data, function(x) diff(range(x$juliant, na.rm=TRUE)))
stem(unlist(mslp_duration))

