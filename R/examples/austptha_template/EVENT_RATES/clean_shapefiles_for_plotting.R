#
# This makes a shapefile with a 'single 1/1000 rate ' peak-stage around
# Australia. It is useful e.g. for standard GA Graphics.
#
# We deliberately remove points 'not near Australia', and also points where the
# tsunami solver is less likely to give good results (e.g. very shallow water).
# Some very shallow points tend to be the most eye-catching and the least valid, because
# our elevation data and hydrodynamic theory are both inadequate for modelling
# in very shallow water, or very near the coast. 
#
# In general we let users decide which points to keep/use. That said, some
# users just want something for a general graphic. Hence this script.
#

elevation_threshold = -30
# Only use points in this polygon
clip_polygon = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/clip_polygon/clip_polygon.shp'
# Shapefile with ARI information
peak_stage = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/tsunami_stages_at_fixed_return_periods'

library(rptha)
clip_pol = readOGR(clip_polygon, layer='clip_polygon')
peak_stg = readOGR(peak_stage,   layer='tsunami_stages_at_fixed_return_periods')

keep = which(c(!is.na(over(peak_stg, clip_pol))) & peak_stg$elev < elevation_threshold)

output = peak_stg[keep, c('lon', 'lat', 'elev', 'ST_1000')]

writeOGR(output, dsn='tsunami_stages_1_in_1000_Australia', 
    layer='tsunami_stages_1_in_1000_Australia', driver='ESRI Shapefile')

system('zip -r tsunami_stages_1_in_1000_Australia.zip tsunami_stages_1_in_1000_Australia')
