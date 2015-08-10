#
# Code for various utility computations in 3D spherical coordinates
#
# Gareth Davies, Geoscience Australia 2015
#

suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(rgdal))

## FIXME: At the moment we use a script to fix a bug in geosphere's antipodal
## I have reported the bug, when it it fixed remove this
#library(geosphere)
source('override_antipodal_geosphere.R')

#' Compute distances between 2 points defined by (lon,lat,depth)
#' Assume all depths are relative to a regular datum (MSL -- key point is we
#' ignore topography of earth)
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
        mypath = gcIntermediate(p1[1:2], p2[1:2], n=n, addStartEnd=TRUE)
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

#

.test_distance_down_depth<-function(){

    ## TEST 1 ## Surface distances

    p1 = c(200, 50, 0)
    p2 = c(170, 10, 0)

    dist1 = distance_down_depth(p2, p1, n=1e+04)

    dist_cos = distCosine(p2[1:2], p1[1:2])

    if(!isTRUE(all.equal(dist1, dist_cos))){
        print('FAIL')
    }else{
        print('PASS')
    }

    ## TEST 2 ## Unaffected by + 360 to lat

    p1 = c(200, 50, 0)
    p2 = c(170, 10, 50)

    dist2 = distance_down_depth(p2, p1)

    p1 = c(200 - 360, 50, 0)
    p2 = c(170, 10, 50)

    dist3 = distance_down_depth(p2, p1)

    if(!isTRUE(all.equal(dist2 , dist3))){
        print('FAIL')
    }else{
        print('PASS')
    }

    ## TEST 3 ##
    dist4 = distance_down_depth(p1, p2) # reverse arguments

    if(!isTRUE(all.equal(dist3 , dist4))){
        print('FAIL')
    }else{
        print('PASS')
    }


    ## TEST 4 ## Unaffected by reflecting lat
    p1 = c(200, -50, 0)
    p2 = c(170, -10, 50)

    dist5 = distance_down_depth(p2, p1)
    if(!isTRUE(all.equal(dist5 , dist3))){
        print('FAIL')
    }else{
        print('PASS')
    }
    
   

    ## TEST 5 ##  Analytical solution for a dipping fault in 1D
    ##
    ## Suppose we integrate over a surface S dipping below a land surface (1D
    ## circle with radius R), with the S's depth increasing at a rate ddepth/dtheta = constant
    ## The coordinate can be parameterised by an angular coordinate 'theta'
    ##
    ## At theta=0, depth=0
    ## At any theta, ddepth = dtheta * (ddepth/dtheta) [since ddepth/dtheta = constant]
    ## and the projection of the dipping surface parallel to the land surface is:
    ##               inner_length_element = (R - depth)*dtheta = (R - theta*ddepth/dtheta)*dtheta 
    ## Therefore, 
    ## downdip_length_element**2 = ddepth**2 + inner_length_element**2 
    ##     = dtheta^2 * ( (ddepth/dtheta)^2 + [R^2 - theta*ddepth/dtheta]^2)
    ## so
    ## d(downdip_length)/dtheta = sqrt ( ddepth/dtheta^2 + [R - theta*ddepth/dtheta]^2] )
    ## we can integrate this to get the down-dip length
    ##
    ## (NOTE: I think this can be analytically integrated too, e.g.)
    ## sage: var('R, p')
    ## sage: assume(R>0)
    ## sage: assume(p>0)
    ## sage: assume(R>p)
    ## sage: f = integral( ( (R - x*p)^2 + p^2)^0.5, x)
    ## sage: f.simplify_full()
    ## 1/2*(p^2*arcsinh((p*x - R)/p) + sqrt(p^2*x^2 - 2*R*p*x + R^2 + p^2)*(p*x - R))/p

    # Example:
    # Drop 50km in 1 degree
    p1 = c(0, 0, 0)
    p2 = c(1, 0, 50)

    # Function to integrate (theta in radians)
    f<-function(theta, ddepth_dtheta = 50/(2*pi/360), R = 6378.137){
        sqrt( ddepth_dtheta**2 + (R - theta*ddepth_dtheta)^2 )
    }
    l_analytical = integrate(f, 0, 2*pi/360)

    # Numerically compute
    l_numerical = distance_down_depth(p1, p2)/1000

    if(isTRUE(all.equal(l_analytical$value, l_numerical))){
        print('PASS')
    }else{
        print('FAIL')
    }

}

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

    surface_path_interpolated = gcIntermediate(surface_path[1,], surface_path[2,], 
        n=n, addStartEnd=TRUE)

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


.test_interpolate_gc_path<-function(){

    ## Test 1
    surface_path = matrix(c(200, 5, 170, -5), ncol=2, byrow=T)

    interp_path1 = interpolate_gc_path(surface_path)
   
    if(max(abs(diff(interp_path1[,1]))< 20)){
        print('PASS')
    }else{
        print('FAIL')
    } 


    ## Test 2
    surface_path = surface_path[2:1,]
    interp_path2 = interpolate_gc_path(surface_path)
  
    # Check there is no large spacing 
    if(max(abs(diff(interp_path2[,1]))< 2)){
        print('PASS')
    }else{
        print('FAIL')
    } 

    # Check these paths are the 'same'
    if(any((abs(rev(interp_path2[,1]) - interp_path1[,1]) > 1.0e-03) |
           (abs(rev(interp_path2[,2]) - interp_path1[,2]) > 1.0e-03) )){
        print('FAIL')
    }else{
        print('PASS')
    }


    ## Test 3
    surface_path = matrix(c(-160, 5, 170, -5), ncol=2, byrow=T)

    interp_path1 = interpolate_gc_path(surface_path)
   
    if(max(abs(diff(interp_path1[,1]))< 20)){
        print('PASS')
    }else{
        print('FAIL')
    } 


    ## Test 4
    surface_path = surface_path[2:1,]
    interp_path2 = interpolate_gc_path(surface_path)
  
    # Check there is not any large spacing 
    if(max(abs(diff(interp_path2[,1]))< 2)){
        print('PASS')
    }else{
        print('FAIL')
    } 

    # Check these paths are the 'same'
    # (except now the longitudes are offset by 360
    if(any((abs(rev(interp_path2[,1]) - interp_path1[,1])%%360 > 1.0e-03) |
           (abs(rev(interp_path2[,2]) - interp_path1[,2]) > 1.0e-03) )){
        print('FAIL')
    }else{
        print('PASS')
    }
}

#' Given a line on the earths surface (defined by 2 points, assumed to follow a
#' great circle path) and some intersecting depth contours, assign depths to
#' path points using the contours.
#'
#' Assumes that the depth is changing at a constant rate with respect to the
#' along-surface distance
#'
#' @param surface_path 2 points defining a great-circle path (2x2 matrix with columns lon,lat)
#' @param depth_contours depth_contours as SpatialLinesDataFrame
#' @param n number of points to used to interpolate along the surface path, before we find the contours.
#'        Note: The interpolation accounts for the great-circle geometry, 
#'              while the intersection treats coordinates as cartesian.
#'        So we need n sufficiently dense that the computed intersection 
#'        is nonetheless reasonable.
#' @param contour_depth_attribute name of the depth_contours attribute giving the depth
#' @param extend_line_fraction Extend contours at the edges of the source zone
#'        by (approx) this fraction of their length. This can help ensure intersections at the edges.
#' @return A 3 column matrix with lon/lat/depth giving the intersection with the contours
#'
#' @export
intersect_surface_path_with_depth_contours<-function(surface_path, depth_contours, n=200, 
    contour_depth_attribute = 'level', extend_line_fraction=0.0e-03){

    if(length(surface_path[,1]) != 2){
        stop('surface_path can only have 2 points')
    }

    # Optionally extend the lines by a small fraction (to ensure intersections occur)
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

                lines@coords  = rbind(lines@coords[1,]  + dx1*extend_line_fraction, lines@coords)
                lines@coords = rbind(lines@coords, lines@coords[lc+1,] + dx2*extend_line_fraction)

                depth_contours@lines[[i]]@Lines[[j]]@coords = lines@coords
            }
        }

        # Now extend the bounding box of depth_contours to reflect the extension we just did
        # The 'trick' to doing that is pass the Lines to the SpatialLines creation function
        depth_contours = SpatialLinesDataFrame(
            SpatialLines(depth_contours@lines, proj4string = CRS(proj4string(depth_contours))),
            data = depth_contours@data)
    }

    # Get points on the path. This is just done so that gIntersection works
    # well (since the latter is not designed for spherical cases)
    new_path = interpolate_gc_path(surface_path, n=n)
    new_path_sp = SpatialLines(list(Lines(list(Line(new_path)), ID="P")), 
        proj4string=depth_contours@proj4string)

    # Compute intersection
    path_contour_intersect = gIntersection(new_path_sp, depth_contours, byid=TRUE)

    pci_indices = match( gsub('P ', '', rownames(path_contour_intersect@coords)),
                      rownames(depth_contours@data) )

    path_contour_intersect_spdf = SpatialPointsDataFrame(
        path_contour_intersect@coords,
        data = depth_contours@data[pci_indices,,drop=FALSE], 
        proj4string=depth_contours@proj4string)

    # Make sure the depth attribute is a number
    path_contour_intersect_spdf[[contour_depth_attribute]] = 
        as.numeric(as.character(path_contour_intersect_spdf[[contour_depth_attribute]] ))

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


.test_intersect_surface_path_with_depth_contours<-function(){

    contour1 = readOGR('../MAKE_SOURCE_CONTOURS/OUTPUT_DATA/CONTOURS_FULL/alaska.shp', 
        layer='alaska')

    line1 = matrix(c(210, 57, 208, 59.5), ncol = 2, byrow=TRUE)

    line1_3D = intersect_surface_path_with_depth_contours(line1, contour1)

    # Check that we get nearly the same result with higher 'n'
    line2_3D = intersect_surface_path_with_depth_contours(line1, contour1, n=200)    

    if(all(abs(line1_3D - line2_3D) < 1.0e-03)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Check what happens if we reverse line1_3D
    # Should get the same result
    line3_3D = intersect_surface_path_with_depth_contours(line1[2:1,], contour1)
    if(all(abs(line1_3D - line3_3D) == 0.0)){
        print('PASS')
    }else{
        print('FAIL')
    }

    ## This case was troublesome originally 
    ## The intersection point with the shallowest contour was missed
    contour2 = readOGR('../MAKE_SOURCE_CONTOURS/OUTPUT_DATA/CONTOURS_FULL/sagami.shp', 
        layer='sagami')

    lc = length(contour2) 
    l1 = length(contour2@lines[[1]]@Lines[[1]]@coords[,1] )
    #l2 = length(contour2@lines[[lc]]@Lines[[1]]@coords[,1] )
    l2 = 1 # Note the last line orientation in this file is currently reversed (24/07/2015)
    end_points = rbind(contour2@lines[[1]]@Lines[[1]]@coords[l1,],
                       contour2@lines[[lc]]@Lines[[1]]@coords[l2,])

    line2_3D = intersect_surface_path_with_depth_contours(end_points, contour2, 
        extend_line_fraction=0.1)

    # Originally end_points[1,] was missing from line2_3D    
    if(min(abs(end_points[1,1] - line2_3D[,1]) + 
           abs(end_points[1,2] - line2_3D[,2])) < 1.0e-02){
        print('PASS')
    }else{
        print('FAIL')
    }

}

#' Given lon-lat point(s) p0, and reference point(s), change the notation
#' of the longitude of p0 so that it differs by < 360 degrees from the 
#' longitude of the reference point
#'
#' @param p0 numeric matrix with columns lon/lat, or a vector with a single
#' point
#' @param reference_point same as p0.
#' @return matrix with the same shape and coordinates as p0, but with
#' longitudes adjusted as required
#'
#'@export
adjust_longitude_by_360_deg<-function(p0, reference_point){

    if(is.null(dim(p0))){
        dim(p0) = c(1,2)
        dim(reference_point) = c(1,2)
    }

    stopifnot(all(dim(p0) == dim(reference_point)))

    londif = p0[,1] - reference_point[,1]

    p0[,1] = p0[,1] - round(londif/360)*360

    return(p0)
}

.test_adjust_longitude_by_360_deg<-function(){

    # Test 1 -- good change
    p0 = c(357, 5)
    refpt = c(2, -10)
   
    p1 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p1 == p0 - c(360, 0))){
        print('PASS')
    }else{
        print('FAIL')
    } 
   
    # Test 2 -- good change 
    p0 = c(-154, 5)
    refpt = c(260, -10)
    
    p1 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p1 == p0 + c(360, 0))){
        print('PASS')
    }else{
        print('FAIL')
    }


    # Test 3 -- good change 
    p0 = c(260, -10)
    refpt = c(-154, 5)
    
    p1 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p1 == p0 - c(360, 0))){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Test 4 -- good change 
    p0 = c(260 + 3*360, -10)
    refpt = c(-154, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p11 == p1)){
        print('PASS')
    }else{
        print('FAIL')
    }


    # Test 5 -- good change 
    p0 = c(260 - 6*360, -10)
    refpt = c(-154, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p11 == p1)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Test 6 -- no change
    p0 = c(260 , -10)
    refpt = c(130, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)
    
    if(all(p11 == p0)){
        print('PASS')
    }else{
        print('FAIL')
    }

}


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
interpolate_3D_path<-function(threeD_path, n=100, depth_in_km=TRUE, ndense = (n*10)){

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

        next_2D_path = interpolate_gc_path(threeD_path[((i-1):i),1:2], n=ndense)[2:(ndense+1),]

        next_depths = seq(threeD_path[i-1,3], threeD_path[i,3], len=ndense+2)[2:(ndense+1)]

        next_distances = seq(path_distance[i-1], path_distance[i], len=ndense+2)[2:(ndense+1)]

        stopifnot(length(next_2D_path[,1]) == length(next_depths))

        fine_path = rbind(fine_path, cbind(next_2D_path, next_depths, next_distances), 
            c(threeD_path[i,], path_distance[i]))
    }

    # Interpolate the dense path
    out_x = approx(fine_path[,4], fine_path[,1], n=n)$y 
    out_y = approx(fine_path[,4], fine_path[,2], n=n)$y 
    out_z = approx(fine_path[,4], fine_path[,3], n=n)$y

    return(cbind(out_x, out_y, out_z))
}

.test_interpolate_3D_path<-function(){


    return
}

#' Given a longitude-latitude point in decimal degrees, return its UTM zone
#' with the WGS84 ellipsoid + datum as a proj4string
#'
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
    if( lonlat[2] >= 56.0 & lonlat[2] < 64.0 & lonlat[1] >= 3.0 & lonlat[1] < 12.0 ){
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
#' compute the coordinates in a 'local' cartesian 2D coordinate system (x,y):
#' dlat = lat - origin_latitude
#' dlon = lon - origin_longitude
#' x --> r*dlon*cos(origin_latitude)
#' y --> r*dlat
#' This will only have good accuracy for small areas near to the origin (and not too close to the
#' poles)
#'
#' @param coords_lonlat matrix with 2 columns (lon/lat) in degrees
#' @param origin_lonlat vector of length 2 giving the lon/lat for the origin of
#'        the new coordinate system
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
    coords_lonlat, origin_lonlat = NULL, r = 6378137){

    # Treat case where coords_lonlat is simplified to a vector 
    if(is.null(dim(coords_lonlat))){
        stopifnot(length(coords_lonlat) == 2)
        coords_lonlat = matrix(coords_lonlat, ncol = 2, nrow = 1)
    }

    if(is.null(origin_lonlat)){
        stop('Must provide origin for cartesian coordinate system in lon/lat)')
    }

    deg2rad = pi/180

    dlon = (coords_lonlat[,1] - origin_lonlat[1])*deg2rad
    dlat = (coords_lonlat[,2] - origin_lonlat[2])*deg2rad
    
    x = r*dlon*cos(origin_lonlat[2]*deg2rad)
    y = r*dlat

    return(cbind(x,y))
}


#' Inverse of spherical_to_cartesian2d_coordinates
#'
#' This will only have good accuracy for small areas near to the origin_lonlat
#' (and not too close to the poles)
#'
#' @param coords_xy matrix with 2 columns containing x,y coordinates in the local coordinate system
#' @param origin_lonlat vector of length 2 with the lon-lat of the origin of the cartesian coordinate system
#' @param r radius of the sphere with the same units as x and y (earth radius in m as default)
#' @return matrix with lon/lat of coords_xy
#'
#' @export
cartesian2d_to_spherical_coordinates<-function(
    coords_xy, origin_lonlat, r = 6378137){
   
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


.test_spherical_to_cartesian2d_and_inverse<-function(){

    ## Test 1: Check that the inverse relation holds

    lonlat = cbind(20 + 0.1*runif(10), -19 + 0.1*runif(10))
    origin = c(20, -19)
    # Convert to local xy system
    new_coords = spherical_to_cartesian2d_coordinates(lonlat, 
        origin_lonlat = origin)

    # Back-caculate old coordinates
    old_coords = cartesian2d_to_spherical_coordinates(new_coords, origin)

    if(all(abs(lonlat - old_coords) < 1.0e-06)){
        print('PASS')
    }else{
        print('FAIL')
    }

    ## Test 2: Works for single-point vector input (as well as matrix, tested
    ## above)
    new_coord2 = spherical_to_cartesian2d_coordinates(lonlat[1,1:2], 
        origin_lonlat = origin)
    old_coord2 = cartesian2d_to_spherical_coordinates(new_coord2, origin)
    
    if(all(abs(old_coord2 - old_coords[1,]) == 0)){
        print('PASS')
    }else{
        print('FAIL')
    }

    ## Test 3: Check that distances are similar to UTM projection (cannot be indentical since
    ##         UTM assumes ellipsoidal model of earth)
    ##         Also check distances are similar to spherical distances

    # Range of origins, including determininstic and random points
    origins = list( c(144.96, -37.81), # Melbourne
                    c(144.96, -60),  # Far south
                    c(0, 62), # High latitude
                    c(0, -62), # High latitude south
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59))
                  ) 

    # Relative distance error tolerances for each origin
    # In reality the value can be low for equatorial regions, with larger
    # errors for high latitudes N/S
    tols = c(0.02, rep(0.04, length(origins) - 1))

    for(ind in 1:length(origins)){

        origin = origins[[ind]]
        tol = tols[[ind]]

        # Make a random set of coordinates 'near' origin
        lonlat = cbind(origin[1] + runif(100), origin[2] + runif(100))

        # Convert to our local coordinate system and compute distance matrix
        new_coords = spherical_to_cartesian2d_coordinates(lonlat, 
            origin_lonlat = origin)
        distances0 = as.matrix(dist(new_coords, diag=FALSE, upper=TRUE))

        # Convert to UTM coordinates
        lonlat_wgs84 = SpatialPoints(lonlat, proj4string=CRS("+init=epsg:4326"))
        local_utm_proj4string = lonlat2utm(origin) 
        lonlat_utm = spTransform(lonlat_wgs84, CRS(local_utm_proj4string))
        lonlat_utm_c = coordinates(lonlat_utm)

        # Compute distance matrix for UTM version
        distances1 = as.matrix(dist(lonlat_utm_c, diag=FALSE, upper=TRUE))

        # Check that distance matrix is 'very close'
        if( max( abs(distances0 - distances1)/distances0, na.rm=TRUE ) < tol){
            print('PASS')
        }else{
            print('FAIL')
        }

        # Now compare distances0 against spherical distances
        # Should be similar but not exactly the same (linear vs spherical)
        distances_sphere = distances1*0
        for(i in 1:length(lonlat[,1])){
            for(j in 1:length(lonlat[,1])){
                distances_sphere[i,j] = distHaversine(lonlat[i,], lonlat[j,])
            }
        }

        if( max( abs(distances_sphere - distances0)/distances0, na.rm=TRUE) < tol){
            print('PASS')
        }else{
            print('FAIL')
        }

    }

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
#' @param method. character.  'complex-mean' gives Arg( mean(exp(i * angles))). 
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
        #'complex-geometric-mean' = ( prod(exp(1i*angles)**(1/length(angles)))),
        stop('method not recognized')
        )

    if(Mod(complex_value) < sqrt(.Machine$double.eps)) stop('Exact cancellation of complex number transform')

    mean_angle = Arg(complex_value)

    if(degrees) mean_angle = mean_angle/deg2rad

    return(mean_angle)
}

.test_mean_angle<-function(){

    # By this method the mean of -1,0,90 is < 30 
    m1 = mean_angle(c(0,0,90))

    err = abs(m1 - 26.5650511771)
    if(err > 1.0e-08){
        print('FAIL')
    }else{
        print('PASS')
    }

    # A more typical case
    m2 = mean_angle(c(120, 60))
    err = abs(m2 - 90)

    if(err > sqrt(.Machine$double.eps)){
        print('FAIL')
    }else{
        print('PASS')
    }

    # Check degrees / radians

    a3 = runif(100, -2, 2)*pi
    a4 = a3/pi*180
    m3 = mean_angle(a3, degrees=FALSE)
    m4 = mean_angle(a4, degrees=TRUE)

    err = abs(m3*180/pi - m4)

    if(err > 1.0e-06){
        print('FAIL')
    }else{
        print('PASS')
    }
   
    # Check there is no impact of representing negative angles as positive 
    a5 = a3 + (a3 < 0)*2*pi
    m5 = mean_angle(a5, degrees=FALSE)

    err = abs(m3 - m5)
    if(err > 1.0e-06){
        print('FAIL')
    }else{
        print('PASS')
    }
    
}

#' Convenience plotting
#'
#' @export
scatter3d<-function(x, y, z, colramp = 'cpt-city/ds9/rainbow.cpt', add=FALSE, 
    ...){
    library(rgl)
    library(colorRampPC)

    colfun = colorRampPC(colramp)

    colz = colfun( (z - min(z))/diff(range(z)))
    plot3d(x, y, z, col = colz, add=add, ...)
}


test_all_geometric_util<-function(){
    .test_distance_down_depth()
    .test_intersect_surface_path_with_depth_contours()
    .test_interpolate_gc_path()     
    .test_spherical_to_cartesian2d_and_inverse()
    .test_mean_angle()
    .test_adjust_longitude_by_360_deg()
}



