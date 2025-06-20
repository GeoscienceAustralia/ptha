library(rptha)
source('convert_lines_to_lon_lat_maxelev.R')

# Workhorse function here -- for all shapefiles associated with a dem
# (local_raster), make a xyz csv file containing the line (with a spacing of
# the provided distance in m). The z values contain the MINIMUM cell in a
# buffer of size line_spacing around each point.
make_lines_with_max_buffer<-function(all_shp, local_raster, line_spacing, buffer_size){

    all_lines = vector(mode='list', length=length(all_shp))

    # Loop over each shapefile and convert to x/y/z csv file
    for(i in 1:length(all_shp)){
        all_lines[[i]] = convert_lines_to_lon_lat_maxelev(all_shp[i], local_raster, 
                                                          line_spacing_m = line_spacing,
                                                          buffer_width_m = buffer_size,
                                                          aggregation_fun = max)
        colnames(all_lines[[i]][[1]]) = c('lon', 'lat', 'z')
        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, quote=FALSE)
    }
}


params_list <- list(
    list(
        shape = 'shapes/port_alma_rd.shp',
        buffer_size = 5,
        line_spacing = 5,
        local_raster = "/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/clip_neg/Rockhampton_2015_no_water_wgs84.tif"
    ),
    list(
        shape = 'shapes/causeway_lake_bridge.shp',
        buffer_size = 10,
        line_spacing = 5,
        local_raster = "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2737432_56_0001_0001_wgs84.tif"
    )
)


# Make the lines
for (params in params_list){
    print(params)
    make_lines_with_max_buffer(Sys.glob(params$shape), params$local_raster, params$line_spacing, params$buffer_size)
}