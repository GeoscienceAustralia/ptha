#
# Code for various utility computations in 3D spherical coordinates
#
# Gareth Davies, Geoscience Australia 2015
#

#' Test if 'angle' is within 'dtheta' of 'target_angle'
#'
#' Function to determine whether the difference between 'angle' and 'target
#' angle' is less than or equal to 'dtheta'. ALL ANGLES ARE IN DEGREES.
#' Accounts for circularity of angles. 
#'
#' @param angle vector of angles
#' @param target_angle vector of target angles
#' @param dtheta vector of 'dthetas' or 'angular deviations'
#' @return logical vector with the same length as angle
#' @export
#'
#' @examples
#'    angle = c(340:359, 0:50)
#'    target_angle = 15
#'    dtheta = 25
#'    nearby = angle_within_dtheta_of_target(angle, target_angle, dtheta)
#'
#'    # Should evaluate TRUE for angles in [350, 0] and [0, 40]
#'    expected_result = c(rep(FALSE, 10), rep(TRUE, 51), rep(FALSE, 10))
#'    stopifnot(all(nearby == expected_result))
#'
#'    # Check that more extreme circularity is ok
#'    for(dev in c(-360*2, 360*3)){
#'        angle2 = angle + dev
#'        nearby2 = angle_within_dtheta_of_target(angle2, target_angle, dtheta)
#'        stopifnot(all(nearby == expected_result))
#'    }
#'    
angle_within_dtheta_of_target<-function(angle, target_angle, dtheta){

    ( (target_angle - angle)%%360 <= dtheta ) | 
    ( (angle - target_angle)%%360 <= dtheta )

}


#' Compute distances between 2 points defined by (lon,lat,depth)
#'
#' Assume all depths are relative to a regular datum (MSL -- key point is we
#' ignore topography of earth). \cr
#' The surface projection of the path is assumed to follow a great circle
#' The rate of change of depth with respect to surface distance is assumed constant
#'
#' @param p1 = (lon1, lat1, depth1) and
#' @param p2 = (lon2, lat2, depth2). 
#' @param depth_in_km if TRUE depths are in km (otherwise meters)
#' @param n = number of intermediate points used in distance computation
#' @param r = radius of the earth (m)
#' @return The distance in m
#'
#' @export
distance_down_depth<-function(p1, p2, depth_in_km=TRUE, n = 1e+03, r = 6378137){
    # FIXME: I think there is an analytical solution that can be adapted to
    # this problem (see unit test where one is derived).

    # There will be n+2 points in the path including the start/end
    
    # Get 'surface' coordinates of path
    if(n > 0){
        # gcIntermediate warns about longitudes outside [-180, 180]
        # This is probably because it always returns interpolated longitudes in
        # that range
        # This doesn't affect us here though, because of the cartesian
        # transform a few lines below.
        mypath = suppressWarnings(
            geosphere::gcIntermediate(p1[1:2], p2[1:2], n=n, addStartEnd=TRUE))
    }else{
        stop('n must be large for this approximation')
    }

    if(depth_in_km){
        p1[3] = p1[3]*1000
        p2[3] = p2[3]*1000
    }

    # Assume d(depth)/d(surface_distance) = constant
    depth = seq(p1[3], p2[3], len=n+2)

    deg2rad = pi/180
    mypath = mypath*deg2rad
  
    xs = (r - depth)*sin(mypath[,1])*cos(mypath[,2])
    ys = (r - depth)*cos(mypath[,1])*cos(mypath[,2])
    zs = (r - depth)*sin(mypath[,2])

    total_distance = sum(sqrt(diff(xs)**2 + diff(ys)**2 + diff(zs)**2))

    return(total_distance) 
}

#' Interpolate points along a great circle 
#'
#' Given a path defined by 2 points, interpolate along the great-circle
#' and return a matrix with the points along the great circle (where lon is
#' adjusted for continuity). There will be n+2 points (so n intermediate points)
#'
#' @param surface_path matrix with dim=c(2,3), with columns = (lon, lat, depth)
#' @param n number of points in the interior of the path
#' @return matrix with dim=c(n+2, 3) with columns = (lon,lat,depth), which
#'         interpolate between the 2 points of surface_path
#'
#' @export
interpolate_gc_path<-function(surface_path, n=50){

    # gcIntermediate warns about longitudes outside [-180, 180]
    # This is probably because it always returns interpolated longitudes in
    # that range
    # This doesn't affect us here though, because we correct the longitudes
    # a few lines below
    surface_path_interpolated = suppressWarnings(
        geosphere::gcIntermediate(surface_path[1,], surface_path[2,], n=n, 
            addStartEnd=TRUE))

    # Check whether longitudes have 'wrapped' (jumped) to stay in [-180, 180]
    path_jumps = diff(surface_path_interpolated[,1])

    # Ensure continuity of paths
    if(any(abs(path_jumps) > 180)){
        offset = c(0, path_jumps * (abs(path_jumps) > 180))
        offset = round(offset/360)*360 
        offset = cumsum(offset)

        surface_path_interpolated[,1] = surface_path_interpolated[,1] - offset
    }

    return(surface_path_interpolated) 

}

#' Get points along a path which intersect with depth contours
#'
#' Given a line on the earths surface (defined by 2 points, assumed to follow a
#' great circle path) and some intersecting depth contours, assign depths to
#' path points using the contours. \cr
#' Assumes that the depth is changing at a constant rate with respect to the
#' along-surface distance. 
#'
#' @param surface_path 2 points defining a great-circle path (2x2 matrix with
#' columns lon,lat)
#' @param depth_contours depth_contours as SpatialLinesDataFrame
#' @param n number of points to used to interpolate along the surface path,
#' before we find the contours. Note: The interpolation accounts for the
#' great-circle geometry, while the intersection treats coordinates as
#' cartesian.  So we need n sufficiently dense that the computed intersection is
#' nonetheless reasonable.
#' @param contour_depth_attribute name of the depth_contours attribute giving
#' the depth
#' @param extend_line_fraction Extend contours at the edges of the source zone
#' by (approx) this fraction of their length. This can help ensure
#' intersections at the edges.
#' @return A 3 column matrix with lon/lat/depth giving the intersection with
#' the contours
#'
#' @export
intersect_surface_path_with_depth_contours<-function(
    surface_path, 
    depth_contours, 
    n=200, 
    contour_depth_attribute = 'level', 
    extend_line_fraction=0.0e-03){

    if(length(surface_path[,1]) != 2){
        stop('surface_path can only have 2 points')
    }

    # Optionally extend the lines by a small fraction (to ensure intersections
    # occur)
    if(extend_line_fraction != 0){
        line_dx = surface_path[2,] - surface_path[1,]
        surface_path[1,] = surface_path[1,] - line_dx*extend_line_fraction
        surface_path[2,] = surface_path[2,] + line_dx*extend_line_fraction

        # Do the same to the depth contours
        for(i in 1:length(depth_contours@lines)){
            line = depth_contours@lines[[i]]
            for(j in 1:length(line@Lines)){
                lines = line@Lines[[j]]
                lc = length(lines@coords[,1])
                #dx1 = (lines@coords[lc,] - lines@coords[1,])

                # 17/07/15 Changed for better performance on very curved source zones
                # Compute vector gradient using "endpoints + midpoint" of line
                lch = max(floor(lc/2), 2)
                dx1 = (lines@coords[1,] - lines@coords[lch,])*2
                lch = min(lch, lc - 1)
                dx2 = (lines@coords[lc,] - lines@coords[lch,])*2

                lines@coords  = rbind(lines@coords[1,] + dx1*extend_line_fraction, 
                    lines@coords)
                lines@coords = rbind(lines@coords, 
                    lines@coords[lc+1,] + dx2*extend_line_fraction)

                depth_contours@lines[[i]]@Lines[[j]]@coords = lines@coords
            }
        }

        # Now extend the bounding box of depth_contours to reflect the
        # extension we just did.  The 'trick' to doing that is pass the Lines
        # to the SpatialLines creation function.
        depth_contours = sp::SpatialLinesDataFrame(
            sp::SpatialLines(depth_contours@lines, 
                proj4string = sp::CRS(proj4string(depth_contours))),
            data = depth_contours@data)
    }

    # Get points on the path. This is just done so that gIntersection works
    # well (since the latter is not designed for spherical cases)
    new_path = interpolate_gc_path(surface_path, n=n)
    new_path_sp = sp::SpatialLines(list(sp::Lines(
            list(sp::Line(new_path)), ID="P")), 
        proj4string=depth_contours@proj4string)

    # Compute intersection
    path_contour_intersect = gIntersection(new_path_sp, depth_contours, 
        byid=TRUE)

    pci_indices = match( 
        gsub('P ', '', rownames(path_contour_intersect@coords)),
        rownames(depth_contours@data) )

    path_contour_intersect_spdf = sp::SpatialPointsDataFrame(
        path_contour_intersect@coords,
        data = depth_contours@data[pci_indices,,drop=FALSE], 
        proj4string=depth_contours@proj4string)

    # Make sure the depth attribute is a number
    path_contour_intersect_spdf[[contour_depth_attribute]] = 
        as.numeric(as.character(
            path_contour_intersect_spdf[[contour_depth_attribute]]))

    # Make 3D points
    threeD_points = cbind(path_contour_intersect_spdf@coords, 
        path_contour_intersect_spdf[[contour_depth_attribute]])

    # Check for unique contour levels
    if(length(unique(threeD_points[,3])) != length(threeD_points[,3])){
        msg = 'ERROR in intersect_surface_path_with_depth_contours: Non-unique contour levels in intersection'
        print(msg)
        print(threeD_points)    
        stop()
    }

    # Order by depth
    threeD_points = threeD_points[order(threeD_points[,3]), ]

    return(threeD_points)

}


#' Adjust longitude notation to be near a chosen value
#'
#' Given lon-lat point(s) p0, and reference point(s), change the notation
#' of the longitude of p0 so that it differs by < 180 degrees from the 
#' longitude of the reference point. This can be useful e.g. if you
#' have a raster with longitude ranging from -180 to 180, and need to look up 
#' the values at points where the longitude notation might include values > 180
#'
#' @param p0 numeric matrix with columns lon/lat, or a vector with a single
#' point
#' @param reference_point same as p0.
#' @return matrix with the same shape and coordinates as p0, but with
#' longitudes adjusted as required
#'
#' @export
#'
#' @examples
#' # Adjust a point to have longitude 'close' to 180 degrees
#' adjust_longitude_by_360_deg(c(-90, 10), c(180, 0))
#'
adjust_longitude_by_360_deg<-function(p0, reference_point){

    if(is.null(dim(p0))){
        if(length(p0) != length(reference_point)){
            stop('p0 and reference_point must have the same shape')
        }
        dim(p0) = c(1,2)
        dim(reference_point) = c(1,2)
    }
    
    if(length(p0[,1]) != length(reference_point[,1])){
        stop('p0 and reference_point must have the same shape')
    }


    stopifnot(all(dim(p0) == dim(reference_point)))

    londif = p0[,1] - reference_point[,1]

    p0[,1] = p0[,1] - round(londif/360)*360

    return(p0)
}

#' Interpolate along a 3D path
#'
#' Given a 3D path of points (lon,lat,depth), interpolate a new path with n
#' evenly spaced points (i.e. equal distances between points in 3D)
#'
#' @param threeD_path matrix of lon,lat,depth defining the threeD path
#' @param n desired number of points on the output path
#' @param depth_in_km are depths in km? If FALSE, assume depths are in metres
#' @param ndense number of points used to support the interpolation (just ensure it is high enough)
#' @return a threeD path with the required number of points
#'
#' @export
interpolate_3D_path<-function(threeD_path, n=100, depth_in_km=TRUE, 
    ndense = (n*10)){

    if(length(threeD_path[1,]) != 3){
        stop('threeD_path must have 3 columns')
    }

    if(length(threeD_path[,1]) < 2){
        stop('threeD_path must have at least 2 points')
    }
    
    # Compute the down-dip distance increments along the midline
    path_distance_increments = threeD_path[,1]*0
    for(i in 2:length(threeD_path[,1])){

        path_distance_increments[i] = distance_down_depth(
            threeD_path[i-1,],
            threeD_path[i,],
            depth_in_km=depth_in_km)
    }

    path_distance = cumsum(path_distance_increments)

    # Make a 'dense' path with n points between every pair of points
    fine_path = c(threeD_path[1,], path_distance[1])
    for(i in 2:length(threeD_path[,1])){

        next_2D_path = interpolate_gc_path(
            threeD_path[((i-1):i),1:2], n=ndense)[2:(ndense+1),]

        next_depths = seq(threeD_path[i-1,3], threeD_path[i,3], 
            len=ndense+2)[2:(ndense+1)]

        next_distances = seq(path_distance[i-1], path_distance[i], 
            len=ndense+2)[2:(ndense+1)]

        stopifnot(length(next_2D_path[,1]) == length(next_depths))

        fine_path = rbind(fine_path, 
            cbind(next_2D_path, next_depths, next_distances), 
            c(threeD_path[i,], path_distance[i]))
    }

    # Interpolate the dense path
    out_x = approx(fine_path[,4], fine_path[,1], n=n)$y 
    out_y = approx(fine_path[,4], fine_path[,2], n=n)$y 
    out_z = approx(fine_path[,4], fine_path[,3], n=n)$y

    return(cbind(out_x, out_y, out_z))
}

#' Find the utm zone of a lon/lat point
#'
#' Given a longitude-latitude point in decimal degrees, return its UTM zone
#' with the WGS84 ellipsoid + datum as a proj4string \cr
#' Zone number derived based on discussions on stack-exchange and
#' stack-overflow:
#' -- http://gis.stackexchange.com/questions/13291/computing-utm-zone-from-lat-long-point
#' -- http://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
#' 
#' @param lonlat vector with c(longitude, latitude) in decimal degrees
#' @return proj4string for the UTM Zone 
#'
#' @export
lonlat2utm<-function(lonlat){

    lonlat[1] = lonlat[1]%%360

    zone_number = floor((lonlat[1] + 180)/6)%%60 + 1;

    # Adjustment near Norway
    if( (lonlat[2] >= 56.0) & (lonlat[2] < 64.0) & 
        (lonlat[1] >= 3.0) & (lonlat[1] < 12.0) ){
        zone_number = 32;
    }
    # Adjustments for Svalbard
    if( lonlat[2] >= 72.0 & lonlat[2] < 84.0 ){ 
        if( lonlat[1] >= 0.0  & lonlat[1] <  9.0 ){ 
            zone_number = 31
        }else if( lonlat[1] >= 9.0  & lonlat[1] < 21.0 ){
            zone_number = 33;
        }else if(lonlat[1] >= 21.0 & lonlat[1] < 33.0 ){
            zone_number = 35;
        }else if(lonlat[1] >= 33.0 & lonlat[1] < 42.0 ){ 
            zone_number = 37;
        }
    }

    # +proj=utm +zone=56 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs
    hemisphere = ""
    if(lonlat[2] < 0) hemisphere = '+south '

    proj4_string = paste0('+proj=utm +zone=', zone_number,' ', hemisphere,
        '+ellps=WGS84 +datum=WGS84 +units=m +no_defs', sep="")

    return(proj4_string)
}

#' Convert from spherical coordinates to a local 2D cartesian coordinate system
#'
#' Given lon-lat coordinates for a sphere, and an origin (also lon-lat),
#' compute the coordinates in a 'local' cartesian 2D coordinate system (x,y):\cr
#' dlat = (lat - origin_latitude)*degrees_to_radians \cr
#' dlon = (lon - origin_longitude)*degrees_to_radians \cr
#' x --> r*dlon*cos(origin_latitude * degrees_to_radians) \cr
#' y --> r*dlat \cr
#' This will only have good accuracy for small areas, near to the origin, and
#' not too close to the poles.
#'
#' @param coords_lonlat matrix with 2 columns (lon/lat) in degrees
#' @param origin_lonlat vector of length 2 (or matrix with 2 columns and same number of rows as coords_lonlat) 
#' giving the lon/lat for the origin of
#' the new coordinate system
#' @param r radius of the sphere (earth radius in m as default)
#' @return matrix with x,y in the new coordinate system (units same as r units)
#'
#' @examples
#' # Make some random coordinates near to (20, -19)
#' lonlat = cbind(20 + 0.1*runif(10), -19 + 0.1*runif(10))
#' origin = c(20, -19)
#' # Convert to local xy system
#' new_coords = spherical_to_cartesian2d_coordinates(lonlat, origin)
#'
#' @export
spherical_to_cartesian2d_coordinates<-function(
    coords_lonlat, 
    origin_lonlat = NULL, 
    r = 6378137){

    # Treat case where coords_lonlat is simplified to a vector 
    if(is.null(dim(coords_lonlat))){
        stopifnot(length(coords_lonlat) == 2)
        #coords_lonlat = matrix(coords_lonlat, ncol = 2, nrow = 1)
        dim(coords_lonlat) = c(1, 2)
    }

    if(is.null(origin_lonlat)){
        stop('Must provide origin for cartesian coordinate system in lon/lat)')
    }else{
        if(is.null(dim(origin_lonlat))){
            dim(origin_lonlat) = c(1, length(origin_lonlat))
        }else{
            if(nrow(origin_lonlat) != nrow(coords_lonlat)){
                stop('If origin_lonlat is a matrix, it must have the same number of rows as coords_lonlat')
            }
        }
    }

    deg2rad = pi/180

    dlon = (coords_lonlat[,1] - origin_lonlat[,1])*deg2rad
    dlat = (coords_lonlat[,2] - origin_lonlat[,2])*deg2rad
    
    x = r*dlon*cos(origin_lonlat[,2]*deg2rad)
    y = r*dlat

    return(cbind(x,y))
}


#' Inverse of spherical_to_cartesian2d_coordinates
#'
#' Suppose coords_xy are 'local' cartesian coordinates on a sphere, 
#' defined by a plane tangent to the sphere with origin at a given lon/lat
#' location (see \code{?spherical_to_cartesian2d_coordiantes} for details). 
#' Then this routine converts coords_xy to lon/lat coordinates.
#'
#' @param coords_xy matrix with 2 columns containing x,y coordinates in the
#' local coordinate system
#' @param origin_lonlat vector of length 2 with the lon-lat of the origin of
#' the cartesian coordinate system
#' @param r radius of the sphere with the same units as x and y (earth radius
#' in m as default)
#' @return matrix with lon/lat of coords_xy
#'
#' @export
cartesian2d_to_spherical_coordinates<-function(
    coords_xy, 
    origin_lonlat, 
    r = 6378137){
   
    # Treat case where coords_xy is simplified to a vector 
    if(is.null(dim(coords_xy))){
        stopifnot(length(coords_xy) == 2)
        coords_xy = matrix(coords_xy, ncol = 2, nrow = 1)
    }

    deg2rad = pi/180

    dlat = coords_xy[,2]/(r*deg2rad)
    dlon = coords_xy[,1]/(r*deg2rad*cos(origin_lonlat[2]*deg2rad))

    lon = dlon + origin_lonlat[1]
    lat = dlat + origin_lonlat[2]

    return(cbind(lon,lat))
}


#' Convert spherical lon/lat coordinates to 3D cartesian x,y,z coordinates
#'
#' @param p0 2 column matrix with columns lon,lat
#' @param r radius of sphere (default = earth radius in m)
#' @return matrix with 3 columns x,y,z giving the cartesian coordinates
#'
#' @export
lonlat2xyz<-function(p0, r = 6378137){
    deg2rad = pi/180

    if(is.null(dim(p0))){
        lon = p0[1]*deg2rad
        lat = p0[2]*deg2rad
    }else{
        lon = p0[,1]*deg2rad
        lat = p0[,2]*deg2rad
    }
    
    x = r*cos(lat)*cos(lon)
    y = r*cos(lat)*sin(lon)
    z = r*sin(lat)

    return(cbind(x,y,z))
}

#' Nearest neighbours on the sphere
#' 
#' For each point in 'lonlat1', find the k nearest points in 'lonlat2', assuming
#' spherical (longitude, latitude) coordinates.
#'
#' @param lonlat1 2 column matrix with lon,lat coordinates -- we want to find the neighbour(s)
#' nearest to these points
#' @param lonlat2 2 column matrix with lon,lat coordiantes. We will return the indices of points
#' in this matrix which are nearest points in lonlat1
#' @param k optional integer, return the 'k' nearest neighbours. By default (k=NULL), return a vector 
#' with the nearest neighbours. Otherwise, return a matrix with the k nearest neighbours in k columns. Note that if k=1
#' it returns a 1-column matrix, whereas if k=NULL it returns a vector (for backwards compatibility)
#' @return A vector (or matrix if k is not NULL) giving an integer index into lonlat2 with the nearest neighbours. 
#' @export
#' @examples
#'    # Find points in lonlat2 nearest to points in lonlat1
#'    lonlat1 = matrix(cbind(c(0,0), c(10, 5), c(-10, 5), c(190, 0)), ncol=2, byrow=TRUE)
#'    lonlat2 = matrix(cbind(c(0.5,0.5), c(0,5), c(-180,0)), ncol=2, byrow=TRUE)
#'    index_match = lonlat_nearest_neighbours(lonlat1, lonlat2)
#'    stopifnot(all(index_match == c(1,2,2,3)))
#'
lonlat_nearest_neighbours<-function(lonlat1, lonlat2, k=NULL){

    if(!is.null(k)){
        k_local = k
    }else{
        k_local = 1
    }

    xyz1 = lonlat2xyz(lonlat1)
    xyz2 = lonlat2xyz(lonlat2)

    output = knnx.index(data=xyz2, query=xyz1, k=k_local, algorithm = 'kd_tree')

    # In the past we only supported k=1, and had the function return a vector
    # For backwards compatibility, keep returning a vector if the argument k is not supplied
    if(is.null(k)){
        return(output[,1])
    }else{
        return(output)
    }
}

#' Find the 'straight line' distance of point p0 to p1, i.e. the Euclidean
#' distance between the points in 3D space. Both points are assumed to be on a sphere
#'
#' This will be shorter than the great circle distance (but similar for small
#' differences between the points)
#'
#' @param p0 vector of lon/lat, or matrix with 1st column lon and 2nd column lat
#' @param p1 vector of lon/lat, or matrix with 1st column lon and 2nd column lat
#' @param r radius of the earth
straight_line_distance<-function(p0, p1, r = 6378137){
    xyz0 = lonlat2xyz(p0, r=r)
    xyz1 = lonlat2xyz(p1, r=r)

    distance = rowSums((xyz0 - xyz1)**2)**0.5

    return(distance)
}

#' Compute the mean of a set of angles. 
#'
#' This accounts for the fact that e.g. mean(c(1, 359)) should be 0 for degrees
#' (not 180 as would be true for the standard mean). The 'complex-mean'
#' method seems to be the most widely used. Note that with this method
#' the [mean of c(0,0,90)] != 30. 
#'
#' @param angles numeric vector. Contains angles to compute mean for
#' @param degrees logical. Are input angles in degrees (TRUE) or radians (FALSE)
#' @param method character.  'complex-mean' gives Arg( mean(exp(i * angles))). 
#' @param weights Optionally provide a vector with the same length of 'angles'
#'        giving weights for weighted mean.
#' 
#' @examples
#'   # The counter-intuitive case
#'   m1 = mean_angle(c(0,0,90))
#'   stopifnot(m1 < 30)
#'   # A more typical case
#'   m2 = mean_angle(c(120, 60))
#'   stopifnot(abs(m2 - 90) < sqrt(.Machine$double.eps))
#'
#' @export
mean_angle<-function(angles, degrees=TRUE, method='complex-mean', weights=1){

    deg2rad = pi/180

    if(degrees) angles = angles*deg2rad

    stopifnot(all(weights >= 0))

    # Normalise weights to avoid case with very small weights
    # from catching 'Exact cancellation of complex number' case below
    weights = weights/sum(weights) 

    complex_value = switch(method,
        'complex-mean' = (mean(exp(1i*angles)*weights)),
        # This has problems with e.g. mean_angle(c(-180, 180)) 
        # 'complex-geometric-mean' = ( prod(exp(1i*angles)**(1/length(angles)))),
        stop('method not recognized')
        )

    if(Mod(complex_value) < sqrt(.Machine$double.eps)){
        stop('Exact cancellation of complex number transform')
    }

    mean_angle = Arg(complex_value)

    if(degrees) mean_angle = mean_angle/deg2rad

    return(mean_angle)
}

#' Return a set of points with chosen spacing along a SpatialLines object
#'
#' The spacing is calculated with linear interpolation (no great circles), so it
#' is assumed that the points defining lines in the SpatialLines object SL are
#' close enough together that 'great circle' and 'euclidean' interpolation 
#' are practically identical. However, distances are computed assuming longlat 
#' (and spherical earth) if 'longlat = TRUE'.
#'
#' The units of spacing depend on the projection of SL, and the value of
#' longlat.  If longlat=FALSE, spacing is in units of SL's coordinates (e.g.
#' metres for UTM projections, degrees for longlat projections).  If
#' longlat=TRUE, SL is in longlat degrees, and spacing is in kilometres
#'
#' @param SL a SpatialLines object, along which we desire evenly spaced points.
#' @param spacing numeric. An approximate spacing between points. Alternatively,
#' the user may provide 'n' below.
#' @param n numeric. The number of points in output (if spacing = NULL)
#' @param longlat logical. Is data longlat?
#' @param verbose logical. Print progress information
#' @param distinguish_disjoint_line_segments logical. If FALSE, then the output
#' contains an integer for each point, with the corresponding line-index in the SL.
#' If TRUE, then the latter is replaced with a real-number, defined as
#' "line-index + (segment_index-1)/(num_segments)", which can be used to
#' distinguish between disjoint parts of the same line.
#' @return SpatialPointsDataFrame along SL, with data identifying the
#' corresponding line index in SL
#' @import sp
#' @export
approxSpatialLines<-function(SL, spacing=NULL, n=NULL, longlat=FALSE, 
    verbose=FALSE, distinguish_disjoint_line_segments = FALSE){

    if(is.null(n)) stopifnot(!is.null(spacing))
    if(is.null(spacing)) stopifnot(!is.null(n))

    all_points = c()
    for(j in 1:length(SL@lines)){

        if(verbose) print(paste0('Line ', j))

        # Get the segment distance between points defining each line in SL
        SL_seglength = lapply(SL@lines[[j]]@Lines, 
            myfun<-function(x) sp::LineLength(x,longlat, sum=F))

        # Store x,y points on each segment
        newxlist = list()
        newylist = list()
        newilist = list()
        for(i in 1:length(SL_seglength)){
            seg_coords = SL@lines[[j]]@Lines[[i]]@coords
            seg_lnth = c(0,cumsum(SL_seglength[[i]]))

            if (!is.null(spacing)){
                # Space out_lnth_points evenly around the segment
                out_lnth_num = ceiling(max(seg_lnth)/spacing)

            }else if(!is.null(n)){
                out_lnth_num = n-1
            }

            out_lnth = seq(0, max(seg_lnth), len=out_lnth_num+1) #[1:out_lnth_num]

            # Interpolate the new coordinates
            new_x = approx(seg_lnth, seg_coords[,1], out_lnth, rule=2)$y 
            new_y = approx(seg_lnth, seg_coords[,2], out_lnth, rule=2)$y

            newxlist[[i]] = new_x
            newylist[[i]] = new_y
            # A line-segment indicator
            newilist[[i]] = rep((i-1)/length(SL_seglength), length(new_x))
        }
        
        all_points_local = cbind(unlist(newxlist), unlist(newylist))
        if(distinguish_disjoint_line_segments){
            tmp = unlist(newilist)
            all_points_local = cbind(all_points_local, j+tmp)
        }else{
            all_points_local = cbind(all_points_local, 
                rep(j, length(all_points_local[,1])))
        }
        
        
        all_points = rbind(all_points,all_points_local)
    }

    # Convert to spatial points
    if(length(all_points) > 0){
        out = sp::SpatialPointsDataFrame(all_points[,1:2], 
            data=data.frame(SLID=all_points[,3]),
            proj4string=sp::CRS(proj4string(SL)))
    }else{
        out = NULL
    }
    return(out)
}

#' Point rotation
#'
#' Given a set of cartesian x,y coordinates, create a new rotated coordinate
#' system with a given origin and x-axis direction and return the points in the
#' new coordinate system. With inverse = TRUE you can undo the rotation, see
#' the example.
#'
#' @param points matrix with 2 columns containing point x,y cartesian
#' coordinates
#' @param origin numeric vector of length 2 giving the origin of the new
#' coordinate system 
#' @param x_axis_vector vector of length 2 which defines the positive x
#' @param inverse logical. If TRUE, perform the inverse operation
#' direction of the new coordinate system
#' @return points in the new coordinate system
#'
#' @export
#'
#' @examples
#' library(rptha)
#' p1 = rbind(c(500,300), c(-60, 20), c(300, 12))
#' origin = c(3000,4012)
#' x_axis_vector = c(30,-1)
#' # Get p1 in a new coordinate system
#' p1_rot = rotate_cartesian2d(p1, origin, x_axis_vector)
#' # Back-calculate the original points
#' inverted_points = rotate_cartesian2d(p1_rot, origin, x_axis_vector, inverse=TRUE)
#' # Should be identical to p1
#' stopifnot(isTRUE(all.equal(c(inverted_points), c(p1))))
#'
rotate_cartesian2d<-function(points, origin, x_axis_vector, inverse=FALSE){

    if(!inverse){
        # Change origin
        points = cbind(points[,1] - origin[1], points[,2] - origin[2])
    }else{
        x_axis_vector = c(x_axis_vector[1], -x_axis_vector[2])
    }

    # Ensure x-axis-vector has unit length
    x_axis_vector_norm = sqrt(sum(x_axis_vector**2))

    if(length(x_axis_vector_norm) > 0){
        x_axis_vector = x_axis_vector / x_axis_vector_norm
    }else{
        stop('x_axis_vector cannot have zero length')
    }

    y_axis_vector = c(-x_axis_vector[2], x_axis_vector[1])

    # rotate
    points_rotated_x = points[,1]*x_axis_vector[1] + points[,2]*x_axis_vector[2]
    points_rotated_y = points[,1]*y_axis_vector[1] + points[,2]*y_axis_vector[2]

    if(inverse){
        points_rotated_x = points_rotated_x + origin[1]
        points_rotated_y = points_rotated_y + origin[2]
    }
    

    return(cbind(points_rotated_x, points_rotated_y))

}


#'
#' Determine whether lonlat coordinates are inside a given polygon
#'
#' Takes care of different longitude conventions (i.e. offsets longitude by 360
#' as required to be inside polygon). If buffer_width is provided, then poly is
#' buffered by this amount prior to testing inclusion. Buffer_width should be
#' in degrees, which means the 'buffer distance' is not identical everywhere.
#'
#' @param lonlat 2 column matrix with longitude/latitude of points
#' @param poly SpatialPolygons object
#' @param buffer_width thickness of buffer to apply to poly before testing the point inclusion.
#' @return logical vector with one entry for every row in lonlat.
#' @export
#' @examples
#' # Point p
#'   p = matrix( c(-180, 10, 180, 10), ncol=2, byrow=TRUE) # Same point with different longitude conventions
#'   poly_coords = matrix(c(170, 9, 190, 9, 190, 11, 170, 11), ncol=2, byrow=TRUE)
#'   poly_sp = SpatialPolygons(list(Polygons(list(Polygon(poly_coords)), ID='1')), proj4string=CRS("+init=epsg:4326"))
#'   is_inside = lonlat_in_poly(p, poly_sp)
#'   stopifnot(all(is_inside==TRUE))
#'
lonlat_in_poly<-function(lonlat, poly, buffer_width = 0){

    if(is.null(dim(lonlat))) dim(lonlat) = c(length(lonlat)/2, 2)

    if(buffer_width != 0){
        poly = gBuffer(poly, width=buffer_width, byid=TRUE)
    }

    # Find point on poly, which is used to determine the longitude convention
    # of 'poly'. Note this will not be the exact centroid for lon/lat --
    # doesn't matter
    point_in_poly = as.numeric(coordinates(gCentroid(poly)))

    # Extend to same dimensions as lonlat
    point_in_poly = matrix(point_in_poly, nrow=length(lonlat[,1]), ncol=2, byrow=TRUE)

    # Make sure longitude convention is the same as for the polygon
    lonlat_near_poly = adjust_longitude_by_360_deg(lonlat, point_in_poly)

    lonlat_near_poly = SpatialPoints(lonlat_near_poly, proj4string=CRS(proj4string(poly)))

    inside_poly = over(lonlat_near_poly, as(poly, 'SpatialPolygons'))

    inside_poly = !is.na(inside_poly)

    return(inside_poly)

}

