# Script to take the multiple outputs produced with 'make_all.sh', and combine
# them (also with the DART buoy locations).

#
# Get DART locations
#
dart_locations = read.csv('INPUTS/DART_LOCATIONS/dart_metadata.csv', 
    colClasses='character', encoding='utf8')
for(i in 1:2) dart_locations[,i] = as.numeric(dart_locations[,i])
dart_xyID = dart_locations[,1:3]

#
# Get other hazard points -- note the files have latitude before longitude,
# which we correct
#
haz_GA250_20 = read.table('OUTPUTS/OUTPUT_GA250_20m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_20[,1:2] = haz_GA250_20[,2:1] 
haz_GA250_100 = read.table('OUTPUTS/OUTPUT_GA250_100m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_100[,1:2] = haz_GA250_100[,2:1] 
haz_GA250_1000 = read.table('OUTPUTS/OUTPUT_GA250_1000m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GA250_1000[,1:2] = haz_GA250_1000[,2:1] 
haz_GEBCO2014_100 = read.table('OUTPUTS/OUTPUT_GEBCO2014_100m/haz_pts_W-180.txt', 
    skip=1, header=FALSE)
haz_GEBCO2014_100[,1:2] = haz_GEBCO2014_100[,2:1] 

haz_gridded = read.csv('OUTPUTS/GRIDDED_POINTS/gridded_points.csv', skip=1, header=FALSE)

#
# Combine the points
#
np = length(dart_xyID[,1]) + length(haz_GA250_100[,1]) + 
    length(haz_GA250_20[,1]) + length(haz_GA250_1000[,1]) + 
    length(haz_GEBCO2014_100[,1]) + length(haz_gridded[,1])
blank = rep(NA, np)
outputs = data.frame(lon=blank, lat=blank, ID=blank, stringsAsFactors=FALSE)
for(i in 1:2){ 
    outputs[,i] = c(haz_GA250_20[,i], haz_GA250_100[,i], haz_GA250_1000[,i], 
        haz_GEBCO2014_100[,i], haz_gridded[,i], dart_xyID[,i])
}

# # Label
# outputs[,3] = c(
#     paste0('GA250_100m_', 1:length(haz_GA250_100[,1])), 
#     paste0('GA250_1000m_', 1:length(haz_GA250_1000[,1])),
#     paste0('GEBCO2014_100m_', 1:length(haz_GEBCO2014_100[,1])), 
#     paste0('dart_', dart_xyID[,3])
#     )
#
# REAL ID -- with a data-source dependent decimal
# This is cruder than using a string, but is more straightforward with netcdf
# output
#
outputs[,3] = c(
    0.0 + 1:length(haz_GA250_20[,1]), # 20m GA250 points end in .0
    0.1 + 1:length(haz_GA250_100[,1]), # 100m ""         end in .1
    0.2 + 1:length(haz_GA250_1000[,1]), # 1000m ""       end in 0.2
    0.3 + 1:length(haz_GEBCO2014_100[,1]), # 100m GEBCO  end in 0.3
    0.5 + 1:length(haz_gridded[,1]), # Gridded GEBCO     end in 0.5
    0.4 + as.numeric(dart_xyID[,3]) # dart               end in 0.4
    )

# We need to translate some longitudes by 360 for everything to work
lower_left_lon = -40
kk = which(outputs$lon < lower_left_lon)
outputs$lon[kk] = outputs$lon[kk]+360

if(length(outputs$ID) == length(unique(outputs$ID))){
    print("Point ID's are unique. Saving to file")
    write.csv(outputs, file='merged_hazard_points.csv', row.names=FALSE, quote=FALSE)
}else{
    stop("ERROR: Non-unique point ID's -- file not saved")
}

