
#' Normalised distance based interpolation along a line
#
#' Make a function to interpolate along a line defined by sorted x,y
#' coordinates. The user must provide an along-line distance for each point
#' along the line (as the coordinates might not be cartesian, e.g. if x/y define
#' lon/lat coordinates). The output function F:[0,1] --> (x,y) with F(0) being
#' the first point on the line and F(1) being the last point, and e.g.  F(0.5)
#' being halfway along the line.
#'
#' @param coords matrix with at least 2 columns giving the x,y coordinates of the line
#' @param distance numeric vector with (length = nrow(coords)) giving the along
#' line distance of each point
#' @return function F:[0-1] --> (x,y) mapping the unit interval to the line.
#'
#' @examples
#' x = c(0, 1, 2, 1, 0)
#' y = c(0, 1, 2, 3, 2)
#' distance = c(0, cumsum(sqrt(diff(x)**2 + diff(y)**2)))
#' f = make_line_interpolation_fun(cbind(x,y), distance)
#' stopifnot(all(f(0) == c(0,0)))
#' stopifnot(all(f(1) == c(0,2)))
#' stopifnot(all(f(0.5) == c(2,2)))
#'
#' @export
make_line_interpolation_fun<-function(coords, distance){

    maxD = max(distance)
    f1 = approxfun(distance/maxD, coords[,1])
    f2 = approxfun(distance/maxD, coords[,2])

    f3<-function(x){
        cbind(f1(x), f2(x))
    }

    return(f3)
}

#' Functions to interpolate points on source contours
#'
#' Given source contours in lon/lat coordinates, make a list of functions which
#' map [0,1] to points on each contour line using
#' \code{make_line_interpolation_fun}. The list has one such function for each
#' contour line, and is ordered by increasing depth.
#'
#' @param source_contours SpatialLinesDataFrame containing the contour
#' lines (in lon/lat coordinates) and an attribute giving the depth.
#' @param contour_depth_attribute Name of the attribute giving the depth.
#' @return A list of functions (one for each line). See
#' \code{make_line_interpolation_fun} for information on what they do. The list
#' is named based on the corresponding contour depth
#' 
#' @export
make_source_contours_interpolation_function_list<-function(
    source_contours,
    contour_depth_attribute='level'){

    #
    # Establish if x,y point on each contour are properly ordered.
    # (i.e. in along-strike order). We assume they are either in along-strike
    # order, or the reverse. Fix them if required
    #
    num_l = length(source_contours@lines)
    contour_line_order = order(as.numeric(as.character(
        source_contours[[contour_depth_attribute]])))

    top_ind = contour_line_order[1]
    top_line = source_contours@lines[[top_ind]]@Lines[[1]]@coords 
    nc = length(top_line[,1])
    # Get a bearing oriented approximately in the direction of the top line,
    # using a lag of 10 points to make the result less sensitive to 'jitters'
    # in the contour xy data
    suppressWarnings({
        top_bearing = bearing(top_line[1,1:2], top_line[pmin(nc, 10),1:2], 
            f=0)%%360 
    })

    # Find the end-point on the bottom line (BL1) that is nearest to the first
    # point on the top line (TL1). 
    # To decide if the top_line is ordered along strike, we can compare the
    # top_bearing with the bearing of the line joining BL1 to TL1. If top_line
    # is ordered along strike, then the latter line will have (approximately)
    # the bearing of TL1 - 90. But if top_line is reversed, then the latter
    # line will have (approximately) the bearing of TL1 + 90. 

    bot_ind = contour_line_order[num_l]
    bot_line = source_contours@lines[[bot_ind]]@Lines[[1]]@coords
    nb = length(bot_line[,1])
   
    suppressWarnings({ 
        d1 = distHaversine(top_line[1,1:2], bot_line[1,1:2])
        d2 = distHaversine(top_line[1,1:2], bot_line[nb,1:2])
        if(d1 < d2){
            updip_bearing = bearing(bot_line[1,1:2], top_line[1,1:2], f=0)%%360
        }else{
            updip_bearing = bearing(bot_line[nb,1:2], top_line[1,1:2], f=0)%%360
        }
    })

    bearing_diff = (top_bearing - updip_bearing)%%360
    if(bearing_diff < 180){
        # Order is correct
    }else{
        # Needs reversal
        top_line = top_line[nc:1,]
    }

    # Compute interpolation function for top line
    suppressWarnings({
        ds = distHaversine(top_line[1:(nc-1),1:2], top_line[2:nc,1:2])
    })
    coord_distance = c(0, cumsum(ds))
    
    source_contours_interpolator_list = list()
    source_contours_interpolator_list[[1]] =  make_line_interpolation_fun(
        top_line[,1:2], coord_distance)
    counter = 1

    # Compute line interpolation function for all contours
    for(i in contour_line_order[2:num_l]){

        counter = counter+1
        next_line = source_contours@lines[[i]]@Lines[[1]]@coords
        nl = length(next_line[,1])

        # Ensure the line coordinates are sorted correctly (along-strike)
        suppressWarnings({
            d1 = distHaversine(top_line[1,1:2], next_line[1,1:2])
            d2 = distHaversine(top_line[1,1:2], next_line[nl, 1:2])
        })
        if(d2 <= d1){
            # Needs to be reordered
            next_line = next_line[nl:1,]
        }

        # Compute the interpolation function
        suppressWarnings({
            ds = distHaversine(next_line[1:(nl-1),1:2], next_line[2:nl,1:2])
        })
        coord_distance = c(0, cumsum(ds))
        source_contours_interpolator_list[[counter]] = make_line_interpolation_fun(
            next_line[,1:2], coord_distance)

    }
    names(source_contours_interpolator_list) = as.character(
        source_contours[[contour_depth_attribute]])

    return(source_contours_interpolator_list)
}

#
# Given a matrix of normalised distances along each line (s_matrix),
# compute the lon/lat coords by looking up the interpolation function
# This is a convenience function to use alongside
# \code{make_source_contours_interpolation_function_list}
get_xy_coords<-function(s_matrix, source_contours_interpolator_list){

    x_coords = s_matrix*0
    y_coords = s_matrix*0
    num_l = length(s_matrix[,1])
    for(i in 1:num_l){
        if(any(diff(s_matrix[i,]) < 0)) stop('non-increasing s_matrix')
        tmp = source_contours_interpolator_list[[i]](s_matrix[i,])
        x_coords[i,] = tmp[,1]
        y_coords[i,] = tmp[,2]
    }

    return(list(x = x_coords, y=y_coords))
}

#
# Get a measure of 'quality' of the s_matrix. Basically we would like it to
# define a grid that is close to orthogonal, and has nearly the same size
# for each cell side. It is impossible to meet both these constraints
# in general, but here we try to find a 'good' solution by optimization.
#
get_quality_matrix<-function(s_matrix, source_contours_interpolator_list){
    
    # Check s
    tmp = try(
        get_xy_coords(s_matrix, source_contours_interpolator_list), 
        silent=TRUE)
    if(class(tmp) == 'try-error'){
        # In case of error, return values which indicate 'very bad fit'
        q_matrix = s_matrix * 0 + Inf
        return(q_matrix)
    }
    x_coords = tmp$x
    y_coords = tmp$y
    
    np = length(s_matrix[1,])
    num_l = length(s_matrix[,1])
    q_matrix = s_matrix*0

    # Preallocate space for loop
    coords_xy_ip1 = matrix(NA, ncol=2, nrow=np)
    coords_xy_im1 = matrix(NA, ncol=2, nrow=np)
    coords_xy_jm1 = matrix(NA, ncol=2, nrow=np)
    coords_xy_jp1 = matrix(NA, ncol=2, nrow=np)
    coords_xy = matrix(NA, ncol=2, nrow=np)

    # Compute the badness of fit in a vectorized fashion for speed.
    # Vectorization is over 'j' which is 'an indivdual contour line'
    for(i in 1:num_l){
        j = 1:np
        jp1 = pmin(2:(np+1), np)
        jm1 = pmax(0:(np-1), 1)
        ip1 = min(num_l, i+1)
        im1 = max(1, i-1)
    
        # i-plus-one coordinates 
        coords_xy_ip1[,1] = x_coords[ip1, j]
        coords_xy_ip1[,2] = y_coords[ip1, j]
        # i-minus-one coordinates 
        coords_xy_im1[,1] = x_coords[im1, j]
        coords_xy_im1[,2] = y_coords[im1, j]
        # j-minus-one coordinates 
        coords_xy_jm1[,1] = x_coords[i, jm1]
        coords_xy_jm1[,2] = y_coords[i, jm1]
        # j-plus-one coordinates 
        coords_xy_jp1[,1] = x_coords[i, jp1]
        coords_xy_jp1[,2] = y_coords[i, jp1]
        # centred coordinates
        coords_xy[,1] = x_coords[i, j]
        coords_xy[,2] = y_coords[i, j]

        # Convert to local cartesians
        local_xy_coords_ip1 = spherical_to_cartesian2d_coordinates(
            coords_xy_ip1,
            origin_lonlat = coords_xy)
        local_xy_coords_im1 = spherical_to_cartesian2d_coordinates(
            coords_xy_im1,
            origin_lonlat = coords_xy)
        local_xy_coords_jp1 = spherical_to_cartesian2d_coordinates(
            coords_xy_jp1,
            origin_lonlat = coords_xy)
        local_xy_coords_jm1 = spherical_to_cartesian2d_coordinates(
            coords_xy_jm1,
            origin_lonlat = coords_xy)

        # Prepare to take dot product of unit vectors up-dip and along-strike
        v1 = local_xy_coords_ip1 - local_xy_coords_im1
        v2 = local_xy_coords_jp1 - local_xy_coords_jm1
        v1A = local_xy_coords_ip1 # - c(0,0) # Note local_xy_coords is 0,0

        # Local 'badness of fit' measure. It includes a component related to
        # the angle at which the up-dip and along-strike vectors meet
        # (orthogonal == good).
        # It also contains a component that increases if the increment in the
        # s_matrix is not close to 1/(np-1). 
        q_matrix[i,] = 
           # Orthogonal-ness measures
           0.25*abs(rowSums(v1A * v2))/(1.0e-12 + sqrt(rowSums(v1A*v1A)*rowSums(v2*v2)))+ 
           abs(rowSums(v1 * v2))/sqrt(rowSums(v1*v1)*rowSums(v2*v2)) + 
           # Even-distance measures
           0.3*(exp(abs((s_matrix[i,jp1] - s_matrix[i,j]) - 1/(np-1))*(np-1)*(jp1 != j))-1) + 
           0.3*(exp(abs((s_matrix[i,j] - s_matrix[i,jm1]) - 1/(np-1))*(np-1)*(jm1 != j))-1)

    }
    return(q_matrix)
}

#' Downdip lines on source contours
#'
#' Make downdip lines (lon/lat/depth) which can be used to form the boundaries
#' for unit sources on source contours
#'
#' @param source_contours SpatialLinesDataFrame with the source contours in lon/lat coordinates
#' @param desired_unit_source_length Desired top-edge length of the unit sources in km
#' @param contour_depth_attribute Name of attribute in source_contours that contains the depth
#' @param force_even_point_spacing Make the lines consist of points that are evenly spaced along each contour.
#' This will generally lead to low-quality lines, but can be useful for debugging purposes.
#' @param make_plot logical. Plot the contours + downdip lines
#' @return a list with the 3d lines
#' @export
create_downdip_lines_on_source_contours_improved<-function(
    source_contours,
    desired_unit_source_length,
    contour_depth_attribute='level',
    force_even_point_spacing=FALSE,
    make_plot=FALSE
    ){


    source_contour_depths = as.numeric(as.character(
        source_contours[[contour_depth_attribute]]))

    # Get length of top line
    shallow_contour = which.min(source_contour_depths)
    top_line = source_contours[shallow_contour,]
    top_line_length_km = SpatialLinesLengths(top_line, longlat=TRUE)

    desired_num_lines = ceiling(top_line_length_km/desired_unit_source_length + 1)

    # Get functions which can interpolate along
    source_contours_interpolator_list = 
        make_source_contours_interpolation_function_list(
            source_contours, contour_depth_attribute)

    num_l = length(source_contours@lines)

    #
    # We can optionally use a multigrid type method (looping over the 'nps'). 
    # However, this does not seem to offer improvements.
    #
    # nps = 2 + 2**(1:4) 
    nps = desired_num_lines 
    for(i in 1:length(nps)){
        np = nps[i]

        # Define distances of points along each line. 'Distances' are normalised
        # by the line length (i.e. are in [0,1])
        #
        # s_matrix has s[1,] being the 'trench' line, and s[num_l,] being the
        # deepest contour line
        if(i == 1){
            s_matrix = matrix(seq(0,1,len=np), ncol=np, nrow=num_l, byrow=TRUE)
            
            if(force_even_point_spacing){
                # This can be useful for debugging purposes (to look at the initial condition)
                new_s_matrix = s_matrix
                break
            }
        }else{
            # Make finer 's_matrix' by interpolating from the coarser optimum result
            s_matrix = matrix(NA, ncol=np, nrow=num_l, byrow=TRUE)
            for(i in 1:num_l){
                s_matrix[i,] = approx(new_s_matrix[i,], n=np)$y
            }
        }

        # For boundary conditions we need s_matrix to be monotonic along the rows
        # We also fix the first and last column, and might fix the trench (or not?)
        moving_par = c(s_matrix[1:num_l, 2:(np-1)])

        optim_fun<-function(moving_par){
            s_matrix_local = s_matrix
            s_matrix_local[1:num_l, 2:(np-1)] = moving_par
            q_matrix = get_quality_matrix(s_matrix_local, 
                source_contours_interpolator_list)
            return(c(q_matrix))
        }

        #nls_max_iter = 10
        #nls_max_fev = 20 * length(moving_par) # Default is 100 times this

        # Converges when the maximum change in the solution is (relatively)
        # nls_ptol. Set the accuracy to within 1/100 of the typical grid
        # spacing.
        #nls_ptol = 1/(np-1) * 0.01 
        model_fit = nls.lm(moving_par, 
            fn=optim_fun, 
            control=list(
                #maxiter=nls_max_iter, 
                #maxfev=nls_max_fev, 
                #ptol=nls_ptol, 
                ftol = 1.0e-03,
                nprint=0))

        new_s_matrix = s_matrix
        new_s_matrix[1:(num_l), 2:(np-1)] = model_fit$par

    }
    s_matrix = new_s_matrix

    all_x_coords = s_matrix*0
    all_y_coords = s_matrix*0
    for(i in 1:nrow(s_matrix)){
        output_coords = source_contours_interpolator_list[[i]](s_matrix[i,])
        all_x_coords[i,] = output_coords[,1]
        all_y_coords[i,] = output_coords[,2]
    }
    output_3d_lines = list()
    output_depths_sorted = sort(source_contour_depths)
    for(i in 1:ncol(s_matrix)){
    
        output_3d_lines[[i]] = cbind(all_x_coords[,i], all_y_coords[,i], 
            output_depths_sorted)

    }

    if(make_plot){

        plot(source_contours, axes=TRUE)
        for(i in 1:length(output_3d_lines)){
            points(output_3d_lines[[i]][,1:2], t='o', col='red', pch=19, cex=0.5)
        }
    }

    return(output_3d_lines)
}

#' Convert downdip lines to a SpatialLinesDataFrame
#' 
#' Convenience function for making a shapefile from the downdip lines. This
#' is useful for checking the result quality
#' 
#' @param new_xy output from \code{create_downdip_lines_on_source_contours_improved}
#' @return a SpatialLinesDataFrame with the downdip lines
#' @export
downdip_lines_to_SpatialLinesDataFrame<-function(new_xy){

    LinesList = list()
    for(i in 1:length(new_xy)){
        LinesList[[i]] = Lines(list(Line(new_xy[[i]][,1:2])), ID=as.character(i))
    }

    output_sl = SpatialLines(LinesList, proj4string=CRS('+init=epsg:4326'))

    output_sldf = SpatialLinesDataFrame(output_sl, data=data.frame(ID=1:length(new_xy)), match.ID=FALSE)

    return(output_sldf)
}


#' Make '3D' down-dip lines along source contours [deprecated]
#'
#' This is as an intermediate step when defining boundaries of
#' the unit sources. It is not meant to be exported to the main package,
#' and has been replaced by more advanced routines.
#'
#' @param source_contours SpatialLinesDataFrame with the source contours and an
#' attribute 'level' giving the contour depth
#' @param desired_subfault_length Desired length of subfaults along-strike
#' @param contour_depth_attribute character The name of the column in the
#' attribute table giving the contour depth
#' @param extend_line_fraction To ensure that contour lines intersect downdip
#' lines at the left/right edges we extend them by this fraction of the
#' end-to-end source length
#' @param orthogonal_near_trench move unit source points along the trench to enhance
#' orthogonality there. Can reduce numerical artefacts at the trench
#' @param make_plot logical (TRUE/FALSE) Plot the result
#' @return a list with the 3d lines
#'
#' @export
#'
create_downdip_lines_on_source_contours<-function(
    source_contours, 
    desired_subfault_length, 
    contour_depth_attribute, 
    extend_line_fraction,
    orthogonal_near_trench,
    make_plot){

    ## NOTE: This routine is deprecated -- there are better ways to do it ##

    # Previously was a function argument, but now I am considering this the
    # only good choice
    #down_dip_line_type = 'mid'
    down_dip_line_type = 'eq_spacing'

    # Get the deepest/shallowest levels
    contour_levels = as.numeric(as.character(source_contours@data$level))
    shallow_contour = source_contours[which.min(contour_levels),]
    deep_contour = source_contours[which.max(contour_levels),]

    # Get points on the shallow contour with approx desired length spacing
    interp_shallow_line = approxSpatialLines(shallow_contour, longlat=TRUE, 
        spacing=desired_subfault_length)
    interp_shallow_line = coordinates(interp_shallow_line)

    # Get the same number of points on the deep contour
    np = length(interp_shallow_line[,1])
    interp_deep_line = approxSpatialLines(deep_contour, longlat=TRUE,
        spacing=NULL, n = np)
    interp_deep_line = coordinates(interp_deep_line)

    # Order the lines appropriately
    # Note: Here we suppress warnings from geosphere about longitudes in 
    # [-180,180], which are caused by .pointsToMatrix therein
    if(suppressWarnings(
        distCosine(interp_shallow_line[1,], interp_deep_line[1,]) > 
        distCosine(interp_shallow_line[1,], interp_deep_line[np,]) )
        ){
        interp_deep_line = interp_deep_line[np:1,]
    }
    # Now the deep and shallow lines are ordered in the same direction

    # Make sure the line is ordered in the 'along-strike' direction
    # Measure this direction as the angle from the first point to a
    # point at index 'mi' along the line (intended as a rough mid-point)
    #
    # Note: 'bearing' calls 'pointsToMatrix' which warns whenever longitudes
    # are outside [-180, 180], but it is no problem
    mi = max(floor(np/2), 2)
    b1 = suppressWarnings(
        geosphere::bearing(interp_shallow_line[1,], interp_shallow_line[mi,], 
            f=0))
    b2 = suppressWarnings(
        geosphere::bearing(interp_shallow_line[1,], interp_deep_line[1,], 
            f=0))

    # If the interp_shallow_line is ordered along-strike, then b1 + 90 should
    # be similar to b2, accounting for angular addition. This means the
    # difference is close to 0 or 360 or -360, etc
    b1_plus_90 = b1 + 90
    angle_diff = (b2 - b1_plus_90)%%360
    if(!(angle_diff < 90 | angle_diff > 270)){
        print('The top contour does not seem to be oriented in the along-strike direction')
        print('Reordering ... (be careful)')
        interp_shallow_line = interp_shallow_line[np:1,]
        interp_deep_line = interp_deep_line[np:1,]
    }

    # Divide the deep line into finely spaced points
    interp_deep_line_dense = approxSpatialLines(deep_contour, longlat=TRUE,
        spacing=NULL, n=np*50)
    interp_deep_line_dense = coordinates(interp_deep_line_dense)

    # Divide the shallow line into finely spaced points
    interp_shallow_line_dense = approxSpatialLines(shallow_contour, longlat=TRUE,
        spacing=NULL, n=np*50)
    interp_shallow_line_dense = coordinates(interp_shallow_line_dense)

    if(make_plot){

        ## Plot
        plot(source_contours, asp=1, axes=TRUE)
        plot(shallow_contour, col='red', add=T)
        plot(deep_contour, col='red', add=T)
        title(source_shapefile)

    }

    # Find 'nearest' points on shallow/deep lines which can be used to make
    # lines 'cutting' the source along the dip
    dip_cuts = list()
    ll = length(interp_shallow_line[,1])
    for(i in 1:ll){

        dip_cuts[[i]] = list()

        if( (i != 1) & (i != ll)){
           
            # Line 1: Join equally spaced points along the top/bottom contours 
            line_eq_spacing = rbind(interp_shallow_line[i,], 
                interp_deep_line[i,])

            # Line2: Join nearest points (measured along-surface) along the
            # bottom contour to the top
            # Note: Here we suppress warnings from geosphere about longitudes
            # in [-180,180], which are caused by .pointsToMatrix therein
            nearest = suppressWarnings(
                which.min(distCosine(interp_shallow_line[i,], 
                    interp_deep_line_dense)))
            line_nearest_surface = rbind(interp_shallow_line[i,], 
                interp_deep_line_dense[nearest,])

            # Line2B: Join nearest points (measured along-deep-contour) along the
            # top contour to the bottom
            # Note: Here we suppress warnings from geosphere about longitudes
            # in [-180,180], which are caused by .pointsToMatrix therein
            nearestB = suppressWarnings(
                which.min(distCosine(interp_deep_line[i,], 
                    interp_shallow_line_dense)))
            line_nearest_deep = rbind(interp_deep_line[i,], 
                interp_shallow_line_dense[nearestB,])
            

            # Line3: Find a midpoint compromise between the nearest and the
            # equal spacing option
            # Note: Here we suppress warnings from geosphere about longitudes
            # in [-180,180], which are caused by .pointsToMatrix therein
            mid_index_bottom = suppressWarnings(
                round(0.5*(nearest + which.min(
                    distCosine(interp_deep_line[i,], interp_deep_line_dense))))
                )
            mid_index_top = suppressWarnings(
                round(0.5*(nearestB + which.min( 
                    distCosine(interp_shallow_line[i,], 
                    interp_shallow_line_dense)))))

            line_mid = rbind(interp_shallow_line_dense[mid_index_top,], 
                interp_deep_line_dense[mid_index_bottom,])

            # Store results 
            dip_cuts[[i]][['eq_spacing']] = line_eq_spacing
            dip_cuts[[i]][['nearest']] = line_nearest_surface
            dip_cuts[[i]][['mid']] = line_mid

            if(make_plot & FALSE){
                points(interpolate_gc_path(line_eq_spacing), t='l', 
                    col='blue', lwd=2, lty='solid')
                points(interpolate_gc_path(line_nearest_surface), t='l', 
                    col='red', lty='solid')
                points(interpolate_gc_path(line_mid), t='l', col='darkgreen', 
                    lty='solid')
            }

        }else{
            # End points need to match up
            line_eq_spacing = rbind(interp_shallow_line[i,], interp_deep_line[i,])
            
            if(make_plot & FALSE){
                points(interpolate_gc_path(line_eq_spacing), t='l', col='blue',
                    lwd=2, lty='solid')
            }

            # All lines must connect the start/end points
            dip_cuts[[i]][['eq_spacing']] = line_eq_spacing
            dip_cuts[[i]][['nearest']] = line_eq_spacing
            dip_cuts[[i]][['mid']] = line_eq_spacing

        }

    }

    chosen_line = down_dip_line_type

    # Find the 'depth' of points along all dip-cut lines
    mid_line_with_cutpoints = list()
    ll = length(dip_cuts)
    for(i in 1:ll){
        # Split up the dip_cut line as appropriate
        mid_line = dip_cuts[[i]][[chosen_line]]

        # Compute the mid_line as a 3D path
        mid_line_with_cutpoints[[i]] = intersect_surface_path_with_depth_contours(
            mid_line, 
            source_contours, 
            n=1000, 
            contour_depth_attribute=contour_depth_attribute,
            extend_line_fraction=extend_line_fraction) 
    }

    # If desired, adjust the top line to enhance orthogonality
    if(orthogonal_near_trench){
        top_line = matrix(NA, nrow=2, ncol=ll)
        second_line = matrix(NA, nrow=2, ncol=ll)

        for(i in 1:ll){
            top_line[,i] = mid_line_with_cutpoints[[i]][1,1:2]
            second_line[,i] = mid_line_with_cutpoints[[i]][2,1:2]
        }
        new_top_line = orthogonal_near_trench(top_line, second_line)
        for(i in 1:ll){
            # Now make the 'new-top-line' exactly on the contour.
            nearest_index = which.min(
                (interp_shallow_line_dense[,1] - new_top_line[i,1])**2 + 
                (interp_shallow_line_dense[,2] - new_top_line[i,2])**2)
            mid_line_with_cutpoints[[i]][1,1:2] = interp_shallow_line_dense[nearest_index,1:2]
        }
    }

    return(mid_line_with_cutpoints)

}


