#
# Make lines defining "inverts" (i.e. features that we'd like to ensure are
# burned into the elevation grid as a min elevation, irrespective of its resolution).
#

library(rptha)
source('convert_lines_to_lon_lat_maxelev.R')

# Workhorse function here -- for all shapefiles associated with a dem
# (local_raster), make a xyz csv file containing the line (with a spacing of
# the provided distance in m). The z values contain the MINIMUM cell in a
# buffer of size line_spacing around each point.
make_lines_with_min_buffer<-function(all_shp, local_raster, line_spacing, buffer_size){

    all_lines = vector(mode='list', length=length(all_shp))

    # Loop over each shapefile and convert to x/y/z csv file
    for(i in 1:length(all_shp)){
        all_lines[[i]] = convert_lines_to_lon_lat_maxelev(all_shp[i], local_raster, 
                                                          line_spacing_m = line_spacing,
                                                          buffer_width_m = buffer_size,
                                                          aggregation_fun = min)
        colnames(all_lines[[i]][[1]]) = c('lon', 'lat', 'z')
        write.csv(all_lines[[i]][[1]], file=gsub('shp', 'csv', all_shp[i]), row.names=FALSE, quote=FALSE)
    }
}

#
# Site-specific applications below here
#
# Manually use gdal to find the closest tile in the vrt. E.G.:
#   gdallocationinfo -geoloc -wgs84 /g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801_DEM_wgs84.vrt 151.2438771 -23.8571051

# Temporarily Merge two rasters for the Auckland inlet drain
local_raster_1 <- "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3217362_56_0001_0001_wgs84.tif"
local_raster_2 <- "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3217361_56_0001_0001_wgs84.tif"
system("module load gdal")
system(paste("gdalbuildvrt tmp_Gladstone201801_DEM_wgs84.vrt", local_raster_1, local_raster_2))
auckland_inlet_vrt <- "tmp_Gladstone201801_DEM_wgs84.vrt"

# Tempprary merge rasters for barney point drain
barney_raster <- c(
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3227363_56_0001_0001_wgs84.tif",
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3227362_56_0001_0001_wgs84.tif",
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3237362_56_0001_0001_wgs84.tif",
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3237361_56_0001_0001_wgs84.tif"
)
system(paste("gdalbuildvrt tmp_barney_point_drain.vrt", paste(barney_raster, collapse=" ")))
barnery_drain_vrt <- "tmp_barney_point_drain.vrt"

# Kooyong park drain
kooyong_prak_drain <- c(
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3207360_56_0001_0001_wgs84.tif",
    "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3217360_56_0001_0001_wgs84.tif"
)
system(paste("gdalbuildvrt tmp_kooyong_prak_drain.vrt", paste(kooyong_prak_drain, collapse=" ")))
kooyong_prak_drain_vrt <- "tmp_kooyong_prak_drain.vrt"

# park_st_underpass_yeppoon
park_st_underpass_yeppoon <- c(
    "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2687440_56_0001_0001_wgs84.tif",
    "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2687439_56_0001_0001_wgs84.tif"
)
system(paste("gdalbuildvrt tmp_park_st_underpass_yeppoon.vrt", paste(park_st_underpass_yeppoon, collapse=" ")))
park_st_underpass_yeppoon_vrt <- "tmp_park_st_underpass_yeppoon.vrt"


params_list <- list(
    list(
        shape = 'buffer_height_lines/barney_point_drain.shp',
        buffer_size = 15,
        local_raster = barnery_drain_vrt
    ),
    # warning("Check this csv file for NA values.")
    list(
        shape = 'buffer_height_lines/kooyong_park_drain.shp',
        buffer_size = 60,
        local_raster = kooyong_prak_drain_vrt
    ),
    list(
        shape = "buffer_height_lines/arthur_st_drain_boyne_island.shp",
        buffer_size = 40,
        local_raster = "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3327350_56_0001_0001_wgs84.tif"
    ),
    list(
        shape = "buffer_height_lines/taranganba_drain.shp",
        buffer_size = 20,
        local_raster = "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2707438_56_0001_0001_wgs84.tif"
    ),
    list(
        shape = "buffer_height_lines/yeppoon_creek_bridge.shp",
        buffer_size = 20,
        local_raster = "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2687439_56_0001_0001_wgs84.tif"
    ),
    list(
        shape = "buffer_height_lines/apex_park_yeppoon.shp",
        buffer_size = 30,
        local_raster = "/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2687439_56_0001_0001_wgs84.tif"
    ),
    list(
        shape = "buffer_height_lines/auckland_inlet_drain.shp",
        buffer_size = 30,
        local_raster = auckland_inlet_vrt
    ),
    list(
        shape = "buffer_height_lines/lord_st_drain_gladstone.shp",
        buffer_size = 20,
        local_raster = "/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3217362_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/park_st_underpass_yeppoon.shp",
        buffer_size=20,
        local_raster=park_st_underpass_yeppoon_vrt
    ),
    list(
        shape="buffer_height_lines/lakeview_circuit_causeway_lake.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2737433_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/blue_water_boulevard_causeway_lake.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2737433_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/lammermoor_native_gardens.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2717437_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/ross_creek_culvert_1.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2707438_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/ross_creek_culvert_2.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2707438_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/ross_creek_culvert_3.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2707438_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/ross_creek_culvert_4.shp",
        buffer_size=20,
        local_raster="/g/data/w85/tsunami/DATA/ELEVATION/QLD/Gladstone_ELVIS_1m_LIDAR/wgs84_version_clipped_onshore/StLawrence201805_DEM_wgs84/StLawrence201805-DEM-GRID-50cm_2707438_56_0001_0001_wgs84.tif"
    ),
    list(
        shape="buffer_height_lines/willian_miskin_park_ditch.shp",
        buffer_size=3,
        local_raster="/g/data/w85/tsunami/MODELS/inundation/QLD_tsunami_inundation_QFES/Gladstone_2024/elevation/data_links/Gladstone201801_DEM_wgs84/Gladstone201801-DEM-GRID-50cm_3217362_56_0001_0001_wgs84.tif"
    )
)

convert_lines_to_lon_lat <- function(line_shapefile, line_spacing_m) {
    # Read lines denoting approximate levee positions
    lines = readOGR(line_shapefile, gsub('.shp', '', basename(line_shapefile)))
    if(!isLonLat(lines)) stop('line_shapefile must be in lon-lat coordinate system')
    # Interpolate with prescribed spacing
    lines_interp = approxSpatialLines(lines, spacing=5/1000, longlat=TRUE, distinguish_disjoint_line_segments=TRUE)
    # Get the coordinates of these points
    line_coords = coordinates(lines_interp)
    return(line_coords)
}

#' make a csv file with given elevation. 
#' 
#' Containing lat, lon and z, where z is a given scalar elevation value
make_csv_with_elev<-function(shape_file, z, line_spacing=5){
    line_coords = convert_lines_to_lon_lat(shape_file, line_spacing_m=line_spacing)
    colnames(line_coords) = c('lon', 'lat')
    # add z column with fixed values
    z_df <- data.frame(z = rep(z, times=nrow(line_coords)))
    xyz_df <- cbind(line_coords, z_df)
    write.csv(xyz_df, file=gsub('shp', 'csv', shape_file), row.names=FALSE, quote=FALSE)
}

# Make fixed height inverts
make_csv_with_elev("fixed_height_lines/briffney_creek.shp", z=0.0)
make_csv_with_elev("fixed_height_lines/yeppoon_creek.shp", z=0.0)
make_csv_with_elev("fixed_height_lines/joe_joseph_dr_bridge.shp", z=-0.85)
make_csv_with_elev("fixed_height_lines/lake_view_circuit_causeway_lake_1.shp", z=0.8)
make_csv_with_elev("fixed_height_lines/shoreline_cl_rosslyn_bay.shp", z=8.0)

# Make the lines
line_spacing = 3
for (params in params_list){
    make_lines_with_min_buffer(Sys.glob(params$shape), params$local_raster, line_spacing, params$buffer_size)
}

# Clean up the temporary vrt
file.remove(auckland_inlet_vrt, showWarnings=FALSE)
file.remove(barnery_drain_vrt, showWarnings=FALSE)
file.remove(kooyong_prak_drain_vrt, showWarnings=FALSE)
file.remove(park_st_underpass_yeppoon_vrt, showWarnings=FALSE)

# Add the files to a list for use in swals
invert_files <- Sys.glob('*/*.csv')
invert_files_full <- normalizePath(invert_files)
writeLines(invert_files_full, 'swals_invert_files.txt')
