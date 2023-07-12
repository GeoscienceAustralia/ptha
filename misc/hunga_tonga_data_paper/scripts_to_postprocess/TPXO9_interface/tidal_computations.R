# Extraction of tides at specific locations

source('predict_tide.R', chdir=TRUE)

site_name = 'Coffs_harbour' # (No spaces in name)
# Site coordinates as decimal degrees (longitude, latitude)
site_coordinates = c(153 + 48/60 + 45.82/(60*60), -(30 + 18/60 + 10.29/(60*60)))

# Timezone MUST be GMT / UTC. Note that -10 means 'ahead 10 hours', so the
# time-zone 'Etc/GMT-10' corresponds to eastern Australia. Confusingly, the
# timezone is often called 'GMT+10'
start_time = strptime('1940-01-01 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-10')
end_time   = strptime('1941-01-01 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-10')
time_interval = '15 min'

coffs_data = get_tidal_prediction(site_name, site_coordinates, 
    start_time, end_time, time_interval)


