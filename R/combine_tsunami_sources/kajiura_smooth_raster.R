#' Apply kajiura to a raster. 
#'
#' Kajiura filter applied directly to a raster. FIXME Integrate into the main package and test.
#'
#' @param source_raster the raster (object of class RasterLayer). 
#' @param new_origin origin for cartesian coordinate system in which kj filter
#' is applied (if spherical_input)
#' @param elevation_raster_file name of file giving elevation. Negative = below MSL.
#' Depths for kajiura filter are taken from this raster
#' @param kj_filter_grid_dxdy value for grid_dx and grid_dy in kajiura_filter
#' @param kj_filter_def_threshold real (m). Only apply kajiura filter in a
#' region where the deformation exceeds kj_filter_def_threshold. Actually the
#' region is somewhat larger, see kj_cartesian_buffer.
#' @param kj_cartesian_buffer. real (m). Exand the region where kj_filter is
#' applied by this distance to avoid boundary effects.
#' @param minimum_kj_depth. real (m). Minimum possible depth passed to kj
#' filter (irrespective of elevation raster)
#' @param elevation_extraction_longitude_offset Add this number to the 'longitude'
#' or x coordinate of the source_raster, prior to extracting associated elevation
#' values from the elevation_raster. The main use is to set it to +360 or -360 to
#' deal with recentering in lon/lat rasters
#' @param spherical_input TRUE/FALSE. Is input raster in lon/lat coordinates.
#' If NULL, then use isLonLat to determine this
#' @return a raster (source raster with kajiura filter applied)
#'
kajiura_smooth_raster<-function(
    source_raster, 
    new_origin, 
    elevation_raster_file,
    kj_filter_grid_dxdy = 1000, 
    kj_filter_def_threshold=1.0e-03, 
    kj_cartesian_buffer = 10000,
    minimum_kj_depth = 10,
    elevation_extraction_x_offset=0,
    spherical_input = NULL){

    if(is.null(spherical_input)){
        spherical_input = isLonLat(source_raster)
    }

    new_orig = new_origin

    xyz_spherical = rasterToPoints(source_raster)

    xyz_cartesian = xyz_spherical

    if(spherical_input){
        xyz_cartesian[,1:2] = spherical_to_cartesian2d_coordinates(
            xyz_spherical[,1:2], 
            origin = new_orig)
    }

    elevation_raster = raster(elevation_raster_file)
    #xyz_depth = extract(elevation_raster, cbind(xyz_spherical[,1] - 360, xyz_spherical[,2]))
    xyz_depth = extract(elevation_raster, 
        cbind(xyz_spherical[,1] + elevation_extraction_x_offset, xyz_spherical[,2]))
    xyz_depth = pmax(minimum_kj_depth, -xyz_depth)

    # Find indices where we will smooth
    deformation_inds = which(abs(xyz_spherical[,3]) > kj_filter_def_threshold)
    kajiura_xrange = range(xyz_cartesian[deformation_inds,1])
    kajiura_yrange = range(xyz_cartesian[deformation_inds,2])
    
    kajiura_inds = which(
        (xyz_cartesian[,1] >= kajiura_xrange[1] - kj_cartesian_buffer) &
        (xyz_cartesian[,1] <= kajiura_xrange[2] + kj_cartesian_buffer) &
        (xyz_cartesian[,2] >= kajiura_yrange[1] - kj_cartesian_buffer) &
        (xyz_cartesian[,2] <= kajiura_yrange[2] + kj_cartesian_buffer) )

    smoothed_perturbation = kajiura_filter(
        xyz_cartesian[kajiura_inds,], xyz_depth[kajiura_inds],
        grid_dx = kj_filter_grid_dxdy, grid_dy = kj_filter_grid_dxdy,
        verbose=FALSE, fortran_inner_loop=TRUE)

    xyz_spherical[kajiura_inds,3] = smoothed_perturbation[,3]

    smoothed_raster = rptha:::.local_rasterFromXYZ(xyz_spherical,
        res=res(source_raster), crs=CRS(proj4string(source_raster)))

    return(smoothed_raster)
}
