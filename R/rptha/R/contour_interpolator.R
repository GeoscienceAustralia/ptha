#' Convert a unit_source_grid to SpatialPolygonsDataFrame
#'
#' Useful to exploit vector operations in sp/rgeos etc
#' @param unit_source_grid 3d array defining the outline of unit sources
#' on some source-zone contours. Typically obtained as the unit_source_grid
#' entry from \code{discretized_source_from_source_contours}.
#' @return SpatialPolygonsDataFrame, with attributes giving the down-dip
#' and along-strike indices of each polygon.
#' @export
#'
unit_source_grid_to_SpatialPolygonsDataFrame<-function(unit_source_grid){

    ndip = length(unit_source_grid[,1,1])-1
    nstrike = length(unit_source_grid[1,1,]) - 1
    np = ndip * nstrike

    poly_list = list()
    poly_data = data.frame(downdip_number = rep(NA, np), alongstrike_number = rep(NA,np))
    counter = 0
    for(i in 1:ndip){
        for(j in 1:nstrike){
            local_polygon = rbind(
                unit_source_grid[i,1:2,j],
                unit_source_grid[i+1,1:2,j], 
                unit_source_grid[i+1,1:2,j+1], 
                unit_source_grid[i,1:2,j+1])
            counter = counter+1
            poly_list[[counter]] = Polygons(list(Polygon(local_polygon)),
                ID=as.character(counter))
            poly_data$downdip_number[counter] = i
            poly_data$alongstrike_number[counter] = j
        }
    }
    
    p1 = SpatialPolygons(poly_list, proj4string=CRS(""))
    p1_df = SpatialPolygonsDataFrame(p1, poly_data)

    return(p1_df)
}

#' Interpolate along downdip lines defining (cartesian) source contours
#'
#' Suppose we have a list of matrices defining a grid of (x,y_depth) points over
#' the source contours, with the along-strike edges aligned with fixed contour
#' levels. This input can be obtained from the 'mid_line_with_cutpoints' entry
#' of the output of \code{discretized_source_from_source_contours}, or from the
#' function \code{create_downdip_lines_on_source_contours_improved}. \cr
#' This function is used to make a function which can assign a depth value to
#' arbitrary x,y points inside the source. Note it can handle conversion from spherical
#' to cartesian coordinates (with convert_to_cartesian=TRUE)
#'
#' @param mid_line_with_cutpoints A list of matrices which is typically obtained
#' from the 'mid_line_with_cutpoints' entry of the output of
#' \code{discretized_source_from_source_contours}, or directly from
#' \code{create_downdip_lines_on_source_contours_improved}. All matrices have
#' the same dimensions (and 3 columns). Each matrix defines a down-dip
#' line along the source contours (x,y,depth). The depth values are increasing,
#' and have the same values in each matrix. The list is ordered along strike. 
#' @param convert_to_cartesian Logical. If TRUE, the mid_line_with_cutpoints
#' is in spherical lon/lat coordinates, and the return function will account
#' for this by transforming x,y data to cartesian prior to fitting, using the
#' provided values of origin_lonlat and r (see
#' \code{spherical_to_cartesian2d_coordinates}).
#' @param origin_lonlat If convert_to_cartesian = TRUE, then this value of
#' origin_lonlat is passed to \code{spherical_to_cartesian2d_coordinates}).
#' @param r radius of the earth. If convert_to_cartesian = TRUE, then this value
#' of origin_lonlat is passed to \code{spherical_to_cartesian2d_coordinates}).
#' @return a function f(xy, allow_outside=FALSE, xy_perturbation_m=NULL) 
#' which can interpolate depth values along the source contours. Non-zero buffer
#' widths can be used to allow the function to be applied to points slightly
#' outside the polygons. xy_perturbation_m can be a matrix of the same shape as xy,
#' in which case it is added to xy before depth computation. This turns out to be convenient
#' for computing dip.
#' @export
#'
make_contour_interpolator<-function(mid_line_with_cutpoints, 
    convert_to_cartesian = FALSE, origin_lonlat = NULL, r = 6378137){

    convert_to_cartesian = convert_to_cartesian
    r = r
    origin_lonlat = origin_lonlat

    if(convert_to_cartesian){
        if(is.null(origin_lonlat)){
            stop('origin_lonlat must be provided if convert_to_cartesian=TRUE')
        }
        for(i in 1:length(mid_line_with_cutpoints)){
            mid_line_with_cutpoints[[i]][,1:2] = 
                spherical_to_cartesian2d_coordinates(
                    mid_line_with_cutpoints[[i]][,1:2],
                    origin_lonlat = origin_lonlat,
                    r = r)
        }
    }

    # Convert mid_line_with_cutpoints to SpatialPolygonsDataFrame, with
    # one polygon for each along-strike region.
    poly_list = list()
    ns = (length(mid_line_with_cutpoints)-1)
    poly_data = data.frame(alongstrike_number = 1:ns)
    for(i in 1:ns){
        nn = length(mid_line_with_cutpoints[[i+1]][,1])
        poly_list[[i]] = Polygons(list(Polygon(
            rbind(mid_line_with_cutpoints[[i]][,1:2], 
                  mid_line_with_cutpoints[[i+1]][nn:1,1:2]))),
            ID = as.character(i))
    }
    p1 = SpatialPolygons(poly_list, proj4string=CRS(""))
    local_polygons = SpatialPolygonsDataFrame(p1, poly_data)

    # Make a function to interpolate from mid_line_with_cutpoints
    #
    # Approach: Find which polygon each 'xy' value (at which we want
    # interpolated depth) is inside
    #
    # NOTE: Sometimes we want to use this function for dip computation. In that
    # case we will compute depth at 2 points (one is down-dip of the other), and
    # it is important that both points are interpolated based on the same mid_line_with_cutpoints
    # polygon. To enforce this, we allow an xy_perturbation_m to be passed to the
    # function directly -- which ensures that 'xy + xy_perturbation_m' will be interpolated
    # using the same polygon as used for 'xy'.
    contour_interpolator<-function(xy, allow_outside=FALSE, xy_perturbation_m=NULL){
        # Ensure xy is a matrix
        if(length(xy) == 2){ 
            dim(xy) = c(1,2)
        }

        if(length(dim(xy)) != 2){
            stop('xy points must be a matrix with 2 columns')
        }

        if(dim(xy)[2] != 2){
            stop('xy points must be a matrix with 2 columns')
        }

        if(!is.null(xy_perturbation_m)){
            if(!(length(xy_perturbation_m) == length(xy))){
                dim(xy_perturbation_m) = dim(xy)
            }
        }

        if(convert_to_cartesian){
            xy = spherical_to_cartesian2d_coordinates(xy, 
                origin_lonlat = origin_lonlat, r=r)
        }


        xy_sp = SpatialPoints(coords=xy)
        local_polygons_sp = as(local_polygons, 'SpatialPolygons')
        poly_id = over(xy_sp, local_polygons_sp)

        if(any(is.na(poly_id))){
            # Some points are outside the polygons
            if(allow_outside){
                # Compute the distance from each point outside the polygons to
                # each polygon, and use the nearest
                kk = which(is.na(poly_id))
                local_dist = gDistance(xy_sp[kk], local_polygons, byid=TRUE)
                local_min_ind = apply(local_dist, 2, which.min)

                # Replace the na polygon id's with the nearest polygon id
                poly_id[kk] = local_min_ind
            }
        }
       
        if(any(is.na(poly_id))){
            stop(paste0(
                ' Some xy points outside the contour mid lines in', 
                ' contour_interpolator. \n',
                'Consider passing allow_outside=TRUE to the function to allow extrapolation'))
        }

        # Do depth interpolation by looping over the polygons which contain 
        # the points
        depths = rep(0, length(poly_id))
        polys_to_use = unique(poly_id)
        for(ptu in polys_to_use){
            # Indices in 'xy' which are inside the "ptu'th" polygon
            #
            # FIXME: Could simplyfy this -- only need polygons defined by pairs
            # of midlines.
            #
            inds = which(poly_id == ptu)
            poly_dip_index = local_polygons$downdip_number[ptu]
            poly_strike_index = local_polygons$alongstrike_number[ptu]

            l0 = mid_line_with_cutpoints[[poly_strike_index]]
            l1 = mid_line_with_cutpoints[[poly_strike_index+1]]

            if(!is.null(xy_perturbation_m)){
                xp = xy[inds,] + xy_perturbation_m[inds,]
            }else{
                xp = xy[inds,]
            }

            interpolated_points = edge_source_interpolator(xp, 
                l0, l1)
            depths[inds] = interpolated_points[,3]
        }

        return(depths)
    }
    return(contour_interpolator)
}

#' Interpolate inside a region defined by 2 downdip lines
#'
#' Interpolate inside a region defined by 2 downdip lines (x,y,z coordinates)
#' based on a given set of xy coordinates
#' 
#' @param xy matrix of xy coordinates ('n' rows, 2 columns)
#' @param edge1 matrix of x,y,z coordinates. 
#' @param edge2 matrix of x,y,z coordinates. 
#' @return xyz coordinates interpolated at xy. Note that x,y should be the same
#' as the input x,y (we return them anyway for debugging/testing purposes)
#' 
edge_source_interpolator<-function(xy, edge1, edge2){
    #interpolating_quad){

    if((is.null(dim(xy))) | (length(dim(xy))!=2) ){
        if(length(xy) == 2){
            dim(xy) = c(1,2)
        }else{
            stop('xy must have 2 columns, or be a vector of length 2')
        }
    }

    # Make function to linearly interpolate along the edge. Allow points
    # outside the edge
    # By using splines, we ensure depth is differentiable, so dip varies
    # smoothly. 
    f_edge<-function(edge){
        s_coord = c(0, cumsum(sqrt(diff(edge[,1])**2 + diff(edge[,2])**2))) 
        s_coord = s_coord/max(s_coord)
        f_x = splinefun(s_coord, edge[,1], method='natural')
        f_y = splinefun(s_coord, edge[,2], method='natural')
        #f_z = splinefun(s_coord, edge[,3], 
        #    method='natural') # Linear extrapolation
        #    #method='monoH.FC') # Monotonic depths
        
        # Use monotonic splines for depth, but enforce linear extrapolation
        # outside [0,1], since monoH splines may not be well behaved with
        # extrapolation
        f_zB = splinefun(s_coord, edge[,3], method='monoH.FC')
        df_dzB_1 = (f_zB(1) - f_zB(1-1.0e-08))/1.0e-08
        df_dzB_0 = (f_zB(1.0e-08) - f_zB(0.0))/1.0e-08
        f_z<-function(alpha){
            (alpha < 0)*(f_zB(0) + (alpha-0)*df_dzB_0) + 
            (alpha > 1)*(f_zB(1) + (alpha-1)*df_dzB_1) + 
            (alpha >= 0)*(alpha <= 1)*f_zB(alpha)
        }

        outfun<-function(alpha){
            out = cbind(f_x(alpha), f_y(alpha), f_z(alpha))
            return(out)
        }
    
        return(outfun)
    }

    edge1_fun = f_edge(edge1)
    edge2_fun = f_edge(edge2)

    # Function to interpolate between the lines (also works outside the lines,
    # which is occasionally useful although not advisable in general)
    f_xyz<-function(alpha_s, xy=NULL, returncoords = FALSE){
        dim(alpha_s) = c(length(alpha_s)/2, 2)
        alpha = alpha_s[,1]
        s = alpha_s[,2]    

        output_coords = s * edge1_fun(alpha) + 
            (1 - s)*(edge2_fun(alpha))

        if(returncoords){
            output = output_coords
        }else{
            # Return a residual
            output = output_coords[,1:2] - xy[,1:2]
            # Coerce to vector since nls.lm works with vectors
            dim(output) = NULL
        }
        return(output)
    }

    # Find 'best' alpha, s values for these coordinates
    alpha_s_best = xy[,1:2,drop=FALSE] * 0
    for(i in 1:length(xy[,1])){
        coordinate_optim = nls.lm(c(0.5, 0.5), lower=c(-Inf,-Inf), upper=c(Inf,Inf), 
            fn=f_xyz, xy=xy[i,1:2, drop=FALSE])
        alpha_s_best[i,] = coordinate_optim$par
    }

    fitted_coordinates = f_xyz(alpha_s_best, xy[,1:2], returncoords=TRUE)
    if(length(fitted_coordinates) == 3) dim(fitted_coordinates) = c(1,3)

    # Check that convergence was ok -- and if not, try to improve the fit
    distance_check = ((fitted_coordinates[,1] - xy[,1])**2 + (fitted_coordinates[,2] - xy[,2])**2)
    refit_threshold = (1.0e-02)**2 # FIXME: Beware hardcoded threshold (but probably ok?)
    refit = which(distance_check > refit_threshold)
    if(length(refit) > 0){
        for(i in 1:length(refit)){
            ri = refit[i]
            coordinate_optim = nls.lm(alpha_s_best[ri,], lower=c(-Inf,-Inf), upper=c(Inf,Inf),
                fn=f_xyz, xy=xy[ri,1:2, drop=FALSE])
            alpha_s_best[ri,] = coordinate_optim$par
        }
        fitted_coordinates[refit,] = f_xyz(alpha_s_best[refit,,drop=FALSE], 
            xy[refit,,drop=FALSE], returncoords=TRUE)
    }

    return(fitted_coordinates)
}

