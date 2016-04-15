# 
# Code to create unit sources from subduction interface contours
#
# Gareth Davies, Geoscience Australia, 2015
#

#suppressPackageStartupMessages(library(rgdal))
#suppressPackageStartupMessages(library(rgeos))
#suppressPackageStartupMessages(library(sp))
#suppressPackageStartupMessages(library(FNN))

#library(geosphere)
## FIXME: At the moment we use a script to fix a bug in geosphere's antipodal
## I have reported the bug, when it it fixed remove this
#source('override_antipodal_geosphere.R')

#cu = new.env()
#gu = new.env()
#source('contour_util.R', local=cu)
#source('geometric_util.R', local=gu)

#' Make discretized_source from subduction interface contours
#'
#' Given a shapefile with contour lines in lon/lat coordinates defining the
#' subduction interface (with an attribute defining their depth), partition into
#' unit sources with chosen approximate length/width. Information on all the
#' unit sources is held in a 'discretized_source' list
#' 
#' @param source_shapefile character name of line shapefile defining the
#' subduction interface contours. It should have an attribute giving the depth
#' @param desired_subfault_length numeric desired length of subfaults (km)
#' @param desired_subfault_width numeric desired width of subfaults (km)
#' @param make_plot logical Make a plot?
#' @param contour_depth_attribute character The name of the column in the
#' attribute table giving the contour depth
#' @param contour_depth_in_km logical Are contour depths given in km? (If False,
#' assume 'm')
#' @param extend_line_fraction To ensure that contour lines intersect downdip
#' lines at the left/right edges we extend them by this fraction of the
#' end-to-end source length
#' @return A list containing: depth_contours The original source contours;
#' unit_source_grid A 3 dimensional array descrbing the unit source vertices;
#' discretized_source_dim A vector of length 2 with number-of-sources-down-dip,
#' number-of-sources-along-strike; fine_downdip_transects A 3 dimensional array
#' containing densly spaced points along the down-dip transects, which might be
#' useful for defining sub-unit-source points for tsunami source integration
#'
#' @export
discretized_source_from_source_contours<-function(
    source_shapefile, 
    desired_subfault_length,
    desired_subfault_width, 
    make_plot=FALSE,
    contour_depth_attribute='level', 
    contour_depth_in_km=TRUE,
    extend_line_fraction=1.0e-01){

    # Previously was a function argument, but now I am considering this the
    # only good choice
    down_dip_line_type = 'eq_spacing'

    # Get the shapefile
    source_contours = rgdal::readOGR(
        dsn = dirname(source_shapefile),
        layer=gsub('.shp', '', basename(source_shapefile)),
        verbose=FALSE)

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
            sphere=TRUE))
    b2 = suppressWarnings(
        geosphere::bearing(interp_shallow_line[1,], interp_deep_line[1,], 
            sphere=TRUE))

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
    strike_cuts = list()   
    fine_strike_cuts = list()
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

    # Find the lengths of the dip cut lines in the down-dip direction
    dip_cut_lengths = unlist(lapply(mid_line_with_cutpoints, 
        f<-function(x) {
            total_distance = 0
            for(i in 1:(length(x[,1])-1)){
                total_distance = total_distance + 
                    distance_down_depth(x[i,], x[i+1,], 
                        depth_in_km=contour_depth_in_km)
            }
            return(total_distance)
        })
        )

    # Decide how many unit-sources there should be (determined by the number of
    # mid-lines cutting the source)
    mean_dip_cut_length = mean(dip_cut_lengths)/1e+03 # km
    number_of_mid_lines = max(
        round(mean_dip_cut_length / desired_subfault_width - 1), 
        0)


    for(i in 1:ll){
        # Now interpolate along the above, with the desired spacing between
        # points
        interpolated_midline = interpolate_3D_path(
            mid_line_with_cutpoints[[i]], 
            n = number_of_mid_lines + 2,
            depth_in_km = contour_depth_in_km)

        fine_interpolated_midline = interpolate_3D_path(
            mid_line_with_cutpoints[[i]], 
            n = 100,
            depth_in_km = contour_depth_in_km)
       
        if(make_plot){
            points(interpolated_midline, col='purple', pch=19)
            points(mid_line_with_cutpoints[[i]][,1:2], col='pink', pch=19, 
                cex=0.2)
            points(mid_line_with_cutpoints[[i]][,1:2], t='l', col='brown', 
                lty = 'dashed') 
        }
        
        strike_cuts[[i]] = interpolated_midline
        fine_strike_cuts[[i]] = fine_interpolated_midline
    }

    # Plot
    for(i in 1:length(interpolated_midline[,1])){
        # Get lines along-strike
        line2plot = matrix( unlist(lapply(strike_cuts, f<-function(x) x[i,])), 
            ncol = 3, byrow = TRUE)
        # Interpolate them along a great circle with the 3D line interpolation
        # routine
        line2plot = cbind(line2plot[,1:2], rep(0, length(line2plot[,1])))
        line2plot = interpolate_3D_path(line2plot)[,1:2]
       
        if(make_plot){ 
            points(line2plot[,1], line2plot[,2], t='l', lty='dashed', col='brown')
        }
    }

    # Convert down-dip lines to 3D arrays (which makes the connectedness
    # properties clearer)
    ll = length(strike_cuts)
    unit_source_grid = array(NA, dim = c(dim(strike_cuts[[1]]), ll))
    fine_unit_source_grid = array(NA, dim = c(dim(fine_strike_cuts[[1]]), ll))
    for(i in 1:ll) unit_source_grid[,,i] = strike_cuts[[i]]
    for(i in 1:ll) fine_unit_source_grid[,,i] = fine_strike_cuts[[i]]

    unit_source_dim = dim(unit_source_grid)[c(1,3)] - 1
    names(unit_source_dim) = c('dip', 'strike')
    
    # Final output data
    unit_source_data = list(
        depth_contours = source_contours, 
        unit_source_grid = unit_source_grid, 
        discretized_source_dim = unit_source_dim,
        fine_downdip_transects = fine_unit_source_grid)

    return(unit_source_data)
}


#' Compute APPROXIMATE summary statistics for all unit sources in a discretized source
#' zone
#'
#' The quantities that we output follow those in the i-invall format in the old
#' URSGA code, although some outputs there (e.g. section number) do not apply
#' with our algorithms
#'
#' @param discretized_source list with an entry 'unit_source_grid' containing a
#'        3 dimensional array defining the vertices of all unit sources for the source
#'        zone (e.g. output of 'discretized_source_from_source_contours')
#' @param default_rake numeric degrees, all unit sources are assigned this rake
#' @param default_slip numeric m, all unit sources are assigned this slip
#' @param make_plot logical, make a plot?
#' @param depth_in_km Are depths in km (TRUE) or meters (FALSE)
#' @return data.frame with key summary statistics
#'
#' @export
discretized_source_approximate_summary_statistics<-function(
    discretized_source,
    default_rake = 90, 
    default_slip = 1, 
    make_plot=FALSE,
    depth_in_km=TRUE){

    #
    # lon, lat, depth, strike, dip, rake, slip, length, width, subfaultnumber, 
    # downdipnumber, alongstrikenumber
    #
    unit_source_grid = discretized_source$unit_source_grid

    deg2rad = pi/180

    along_strike_source_count = dim(unit_source_grid)[3] - 1
    down_dip_source_count = dim(unit_source_grid)[1] - 1
   
    num_sources = along_strike_source_count*down_dip_source_count 

    lon_c = rep(NA, num_sources)
    lat_c = rep(NA, num_sources)
    depth = rep(NA, num_sources)
    strike = rep(NA, num_sources)
    dip = rep(NA, num_sources)
    rake = rep(default_rake, num_sources)
    slip = rep(default_slip, num_sources)
    len = rep(NA, num_sources)
    width = rep(NA, num_sources)
    subfault_number = rep(NA, num_sources)
    downdip_number = rep(NA, num_sources)
    alongstrike_number = rep(NA, num_sources)
    max_depth = rep(NA, num_sources)

    source_coordinates = list() 

    ## BEWARE: The loop ordering here should be the same as for the corresponding
    ## detailed summary statistics routine, so that 'counter' is updated correctly
    counter = 0    
    for(i in 1:along_strike_source_count){
        for(j in 1:down_dip_source_count){

            # Record indices
            counter = counter+1
            subfault_number[counter] = counter
            downdip_number[counter] = j        
            alongstrike_number[counter] = i

            # Get polygon defining unit_source boundary
            source_coords = rbind(unit_source_grid[j, 1:3, i],
                                  unit_source_grid[j+1, 1:3, i],
                                  unit_source_grid[j+1, 1:3, i+1],
                                  unit_source_grid[j, 1:3, i+1])
            # Save for later
            source_coordinates[[counter]] = source_coords

            # Basic depth/lon/lat computation
            lon_c[counter] = mean(source_coords[,1])
            lat_c[counter] = mean(source_coords[,2])
            depth[counter] = mean(source_coords[,3])
            max_depth[counter] = max(source_coords[,3])
           
            # Strike = midpoint bearing along great-circle of top of unit
            # source. 
            # FIXME: Consider averaging with bottom of unit source.
            # NOTE: Here we suppress geosphere warnings about longitudes
            # ranging from [-180, 180], which is caused by their .pointsToMatrix
            midpoint = suppressWarnings(
                midPoint(source_coords[1,1:2], source_coords[4,1:2]))
            st0 = suppressWarnings(
                bearing(midpoint, source_coords[4,1:2], sphere=TRUE)%%360)
            strike[counter] = st0 

            stopifnot( (strike[counter] > 0) & (strike[counter] < 360)) 
  
            # Dip = atan (change in depth / surface distance) 
            # Take average of left/right dipping lines (complex numbers to
            # average an angle)
            #
            # Account for depth being in m or km
            depth_scale = 1
            if(depth_in_km) depth_scale = 1000

            dl = atan2( 
                    (source_coords[2,3] - source_coords[1,3])*depth_scale,
                    suppressWarnings(
                        distHaversine(source_coords[1,1:2], 
                            source_coords[2,1:2])))
            dr = atan2( 
                    (source_coords[3,3] - source_coords[4,3])*depth_scale,
                    suppressWarnings(
                        distHaversine(source_coords[4,1:2], 
                            source_coords[3,1:2])))

            dip[counter] = mean_angle(c(dl, dr), degrees=FALSE)/deg2rad

            stopifnot(dip[counter] > 0)

            # Length = 0.5*( upper length + lower length)
            len0 = 0.5*(
                distance_down_depth(source_coords[1,], source_coords[4,], 
                    depth_in_km=depth_in_km) +
                distance_down_depth(source_coords[2,], source_coords[3,], 
                    depth_in_km=depth_in_km)
            )

            surface_area = areaPolygon(source_coords[,1:2], f = 0)

            # FIXME: Consider that with our definitions there is also a
            # 'horizontal' component to the slope, not considered here
            sloping_area = surface_area*sqrt(1 + tan(dip[counter]*deg2rad)**2)

            # Width := sloping_area / sloping length
            width0 = sloping_area / len0
           
            len[counter] = len0/1000
            width[counter] = width0/1000
 
        }
    }


    output = data.frame(lon_c = lon_c, lat_c = lat_c, depth = depth, strike = strike, 
        dip = dip, rake = rake, slip = slip, length = len, width = width, 
        downdip_number = downdip_number, alongstrike_number = alongstrike_number, 
        subfault_number = subfault_number, max_depth = max_depth)

    # Plotting only below here

    if(make_plot){

        plot_ylim = range( unlist( lapply(source_coordinates, f<-function(x) x[,2])))
        plot_xlim = range( unlist( lapply(source_coordinates, f<-function(x) x[,1])))

        plot(lon_c, lat_c, asp=1, col='purple', pch=19, xlim=plot_xlim, ylim=plot_ylim)

        # Add vertices of unit sources
        for(i in 1:length(source_coordinates)){
            sc = rbind(source_coordinates[[i]], source_coordinates[[i]][1,])
            points(sc[,1:2], t='o', col='red')
        }   

        # Add vector giving direction of slip, using a great-circle
        # Get rid of geosphere warnings about longitudes outside of [-180,180],
        # as we deal with that below
        slip_angle = strike - rake
        end_points = suppressWarnings(
            destPoint(cbind(lon_c, lat_c), slip_angle, (len/3)*1000, f=0))

        # end_points will have lat in [-180, 180] because that's how 'geosphere' works
        # Fix this.
        lap = 1
        while(lap > 0){
            adjust_points = which(abs(lon_c - end_points[,1]) > 180)
            lap = length(adjust_points)
            if(lap > 0){
                end_points[adjust_points, 1] = end_points[adjust_points, 1] + 360
            }
        }

        arrows(lon_c, lat_c, end_points[,1], end_points[,2], col='orange')

        title(main = 'Unit sources with slip direction (strike - rake)')
    }


    return(output)
}


#' Compute summary statistics for all unit sources in a discretized source
#' zone.
#'
#' The quantities that we output follow those in the i-invall format in the old
#' URSGA code, although some outputs there (e.g. section number) do not apply
#' with our algorithms. The calculations are based on filling the unit sources
#' with subgrid points (via unit_source_interior_points_cartesian), so should
#' more accurately account for the interface structure than does the alternative
#' 'approximate summary statistics' routine. The width and length are defined
#' so that (width x length) = (dipping interface area).
#'
#' @param discretized_source list with an entry 'unit_source_grid' containing a
#'        3 dimensional array defining the vertices of all unit sources for the source
#'        zone (e.g. output of 'discretized_source_from_source_contours')
#' @param default_rake numeric degrees, all unit sources are assigned this rake
#' @param default_slip numeric m, all unit sources are assigned this slip
#' @param approx_dx numeric. Value of approx_dx passed to unit_source_interior_points_cartesian
#' @param approx_dy numeric. Value of approx_dy passed to unit_source_interior_points_cartesian
#' @param make_plot logical, make a plot?
#' @param depth_in_km Are depths in km (TRUE) or meters (FALSE)
#' @return data.frame with key summary statistics
#' @export
discretized_source_summary_statistics<-function(
    discretized_source,
    default_rake = 90,
    default_slip = 1, 
    approx_dx=2000,
    approx_dy=2000,
    make_plot=FALSE,
    depth_in_km=TRUE){

    # First get the approximate statistics 
    output = discretized_source_approximate_summary_statistics(
        discretized_source,
        default_rake = default_rake, 
        default_slip = default_slip, 
        make_plot=make_plot,
        depth_in_km=depth_in_km)

    # Correct some of the approximate statistics using the detailed cartesian
    # information
    ndip =    discretized_source$discretized_source_dim[1]
    nstrike = discretized_source$discretized_source_dim[2]
  
    # BEWARE: The loop ordering here should be the same as for the corresponding
    # approximate summary statistics routine, so that the variable 'counter'
    # is updated correctly
    for(ns in 1:nstrike){
        for(nd in 1:ndip){

            counter = which(
                (output$downdip_number == nd) & 
                (output$alongstrike_number == ns))

            # Logical checks
            if(length(counter) != 1){
                stop('Did not find exactly one row corresponding to this subfault in the table')
            }
            if(output$subfault_number[counter] != counter){
                stop('There seems to be a table ordering bug')
            }

            # Compute Cartesian source
            usc = unit_source_interior_points_cartesian(
                discretized_source, 
                c(nd, ns), 
                origin=NULL, 
                approx_dx=approx_dx, 
                approx_dy=approx_dy)

            # Get a good approximation of the rupture interface area
            area = sum(usc$grid_points[,'area_projected']*
                sqrt(1 + tan(usc$grid_points[,'dip']/180*pi)**2))

            # We want length * width = area
            # -- keep length as average of top/bottom lengths
            p0 = usc$unit_source_cartesian[1,1:3]
            p1 = usc$unit_source_cartesian[4,1:3]
            len0 = sqrt(sum((p0 - p1)**2))
            p0 = usc$unit_source_cartesian[2,1:3]
            p1 = usc$unit_source_cartesian[3,1:3]
            len1 = sqrt(sum((p0 - p1)**2))

            # Make 'good' estimates of length and width
            output$length[counter] = (0.5*(len1 + len0))/1000
            output$width[counter] = (area/1e+06)/output$length[counter]

            # Depth in km
            output$depth[counter] = weighted.mean(usc$grid_points[,'depth'], 
                usc$grid_points[,'area_projected'])/1000

            # Strike/dip consistent with subgrid
            output$strike[counter] = mean_angle(usc$grid_points[,'strike'])
            # Ensure strike is > 0
            if(output$strike[counter] < 0){
                output$strike[counter] = output$strike[counter] + 360
            }
            output$dip[counter] = mean_angle(usc$grid_points[,'dip'])
        }
    }

    return(output)
}


#' Get coordinates of a single unit source from the discretized_source list
#'
#' @param discretized_source List holding the discretized_source information (e.g. output
#'        of discretized_source_from_source_contours)
#' @param unit_source_index Vector of length 2 giving the index of the desired
#'        unit source (down-dip, along-strike)
#' @return List containing unit source information
#' 
#' @export
#'
get_unit_source_from_discretized_source<-function(discretized_source, 
    unit_source_index){

    stopifnot(length(unit_source_index) == 2)
    stopifnot(unit_source_index[1] >= 1)
    stopifnot(unit_source_index[1] <= 
        discretized_source$discretized_source_dim[1])
    stopifnot(unit_source_index[2] >= 1)
    stopifnot(unit_source_index[2] <= 
        discretized_source$discretized_source_dim[2])

    d_d = unit_source_index[1] # down-dip
    a_s = unit_source_index[2] # along-strike
    unit_source_coords = rbind(
        discretized_source$unit_source_grid[d_d:(d_d+1), 1:3, a_s], 
        discretized_source$unit_source_grid[(d_d+1):d_d, 1:3, a_s+1])

    fine_downdip_transects = 
        discretized_source$fine_downdip_transects[,,(a_s:(a_s+1))]

    output = list(unit_source_grid = unit_source_coords, 
                  fine_downdip_transects = fine_downdip_transects)
    return(output)
}

#' Extract the lon-lat midpoints of the top edge of all unit sources
#' along the 'shallowest' edge of a discrete source (i.e. near the trench), and
#' compute the strike
#'
#' @param discretized_source A list containing discretized source information
#' (output of discretized_source_from_source_contours or similar)
#' @return matrix with 3 columns giving the lon,lat,strike of the top edge of
#' each unit source
#'
#' @export
get_shallow_unit_source_top_edge_strikes<-function(discretized_source){

    discretized_source_dim = discretized_source$discretized_source_dim
    nstrike = discretized_source_dim[2]
    ndip = discretized_source_dim[1]

    # Make a matrix of all 'top-left' and 'top-right' unit source vertex points
    # (looking up dip)
    vertex_tl = c() # Top left vertex
    vertex_tr = c() # Top right vertex

    # Skip the last down-dip index as it is not the 'top' edge of a unit source
    for(i in 1:nstrike){
         vertex_tl = rbind(vertex_tl, 
            discretized_source$unit_source[1, 1:2, i])
         vertex_tr = rbind(vertex_tr, 
            discretized_source$unit_source[1, 1:2, i+1])
    }

    # Note: Here we suppress warnings from geosphere about longitudes
    # in [-180,180], which are caused by .pointsToMatrix therein
    top_edges = suppressWarnings(midPoint(vertex_tl, vertex_tr))
    strike = suppressWarnings(bearing(top_edges, vertex_tr, sphere=TRUE))

    # Ensure strike in 0-360 [not -180 to 180]
    strike = strike + (strike < 0)*360

    # Correct the top_edges hemisphere
    top_edges = adjust_longitude_by_360_deg(top_edges, vertex_tr)

    return(cbind(top_edges, strike))
}


#' Fill a unit source with points (for tsunami source integration)
#'
#' Make a locally linear coordinate system for a unit source and compute a
#' 'grid' of points inside the source. These can be used later for tsunami
#' source integration
#'
#' @param discretized_source list containing unit sources information (e.g. output of
#' discretized_source_from_source_contours)
#' @param unit_source_index Index (down-dip, along-strike) of the unit source
#' to operate on.
#' @param origin vector with 2 entries (lon,lat) giving the local cartesian coordinate
#' system origin for x,y. The depth origin is at 0. If NULL, the lon,lat of the
#' first coordinate associated with the unit source at unit_source_index is
#' used. 
#' @param r Radius of the earth (default gives radius in m)
#' @param approx_dx Approximate x spacing of grid points inside the unit source
#' (same units as r). If NULL, then half the min depth of the unit source is used.
#' @param approx_dy Approximate y spacing of grid points inside the unit source
#' (same units as r). If NULL, then half the min depth of the unit source is used.
#' @param scale_dxdy approx_dx and approx_dy are scaled by this constant.
#' Most useful in the case that the latter are set to NULL since it allows
#' rescaling of defaults.
#' @param depths_in_km logical. If TRUE, input depths are assumed in km, and
#' are transformed to m (so output depths are in m). Otherwise input depths are
#' assumed in m
#' @param plot_source logical. make a simple plot?
#' @return list containing: 'unit_source': matrix with unit_source_coords in the
#' cartesian coordinate system + an additional depth column (m); 'origin' the
#' lon/lat origin of the local coordinate system; r (input argument); 'dx, dy'
#' the actual dx, dy of the grid points (not exactly the same as the approximate
#' values); 'grid_points' grid points inside the source in the cartesian
#' coordinate system; 'grid_point_polygons' a polygon associated with the grid points
#'
#'@export
unit_source_interior_points_cartesian<-function(
    discretized_source,
    unit_source_index,
    origin = NULL,
    r = 6378137,
    approx_dx = NULL,
    approx_dy = NULL,
    scale_dxdy = 1,
    depths_in_km = TRUE,
    plot_source=FALSE){

    ## Get the unit source coordinates, and down-dip transects we can use for
    ## interpolation
    unit_source_info = get_unit_source_from_discretized_source(
        discretized_source, 
        unit_source_index)

    # We will change these variables to cartesian coordinates
    unit_source_coords = unit_source_info$unit_source_grid
    fine_downdip_transects_cartesian = unit_source_info$fine_downdip_transects

    if(is.null(origin)){
        origin = unit_source_coords[1,1:2]
    }

    ## Convert the unit source lon/lat coordinates to a local 2D cartesain system
    unit_source_cartesian = unit_source_coords
    unit_source_cartesian[,1:2] = spherical_to_cartesian2d_coordinates(
        unit_source_coords[,1:2], origin_lonlat = origin, r = r)

    if(depths_in_km) unit_source_cartesian[,3] = 1000*unit_source_cartesian[,3]

    ## Adjust approx_dx, approx_dy
    if(is.null(approx_dx)){
        approx_dx = min(unit_source_cartesian[,3])*0.5
    }
    if(is.null(approx_dy)){
        approx_dy = min(unit_source_cartesian[,3])*0.5
    }
    approx_dx = approx_dx*scale_dxdy
    approx_dy = approx_dx*scale_dxdy

    ## Replace lon/lat spherical coordinates in 'fine-transects' with local cartesians
    fine_downdip_transects_cartesian[,1:2,1] = spherical_to_cartesian2d_coordinates(
        fine_downdip_transects_cartesian[,1:2,1], origin_lonlat = origin, r=r)
    fine_downdip_transects_cartesian[,1:2,2] = spherical_to_cartesian2d_coordinates(
        fine_downdip_transects_cartesian[,1:2,2], origin_lonlat = origin, r=r)

    if(depths_in_km){
        fine_downdip_transects_cartesian[,3,1:2] = 
            1000*fine_downdip_transects_cartesian[,3,1:2]
    }

    ## Get a 'grid' of points inside the unit source (later used for
    ## integration)
    grid_point_data = compute_grid_point_areas_in_polygon(
        unit_source_cartesian[,1:2], approx_dx = approx_dx, 
        approx_dy = approx_dy)

    grid_points = grid_point_data$grid_points

    grid_point_areas = grid_point_data$area

    grid_point_polygon = grid_point_data$grid_point_polygon

    # Test for consistency between the grid point areas and the original unit
    # source
    a0 = sum(grid_point_areas)
    a1 = areaPolygon(unit_source_coords[,1:2], f=0)
    rel_err = abs(a0-a1)/((a0+a1)*0.5)
    if( rel_err > 0.01 ){
        msg = paste0('Projected and spherical unit source area differ by a fraction ',
            rel_err, sep=" ")
        warning(msg)
    }

    
    ## Get strike at the grid points

    # Convert lon-lat centroids of ALL unit sources in the discrete source to
    # the local coordinate system and use them to make a continuous function of
    # strike

    ds1_lonlatstrike = get_shallow_unit_source_top_edge_strikes(discretized_source) 
    ds1_stats_points_cartesian = spherical_to_cartesian2d_coordinates(
        ds1_lonlatstrike[,1:2], origin_lonlat=origin, r=r)

    # Ideally we interpolate strike using a weighted nearest-neighbours
    # However, the code used here requires that we have > 1 point
    np = length(ds1_lonlatstrike[,1])
    k = min(4, np)
    # Compute an inverse distance weighted angular mean of k nearest neighbours
    inds = knnx.index(data = ds1_stats_points_cartesian[,1:2, drop=FALSE], 
        query = grid_points[,1:2, drop=FALSE], k=k)
    dists = knnx.dist(data = ds1_stats_points_cartesian[,1:2, drop=FALSE],
        query = grid_points[,1:2, drop=FALSE], k=k)
    mean_strike = dists[,1]*NA 
    for(ii in 1:length(mean_strike)){
        mean_strike[ii] = mean_angle(ds1_lonlatstrike[inds[ii, ], 3], 
            weights = 1/(dists[ii, ])**2)
    }
    strike = mean_strike

    # Ensure > 0
    strike = strike + (strike < 0)*360

    #print('GD: Deliberately breaking strike')
    #strike = strike*0 + mean_angle(strike)

    ## Get dip/depth etc at the grid points
    grid_points = get_depth_dip_at_unit_source_interior_points(
        unit_source_cartesian, grid_points,
        fine_downdip_transects_cartesian, strike)

    grid_points = cbind(grid_points, grid_point_areas, strike)

    # Fix row/column names for later ease of access
    colnames(grid_points) = c('x', 'y', 'depth', 'dip', 'alpha', 's', 'area_projected', 'strike')

    rownames(grid_points) = NULL

    # Make a list to hold the unit source interior points information
    output = list(unit_source_cartesian=unit_source_cartesian, 
        origin_lonlat=origin, r=r, 
        dx=grid_point_data$dx, dy=grid_point_data$dy, 
        grid_points=grid_points,
        grid_point_polygon=grid_point_polygon,
        fine_downdip_transects_cartesian = fine_downdip_transects_cartesian)

    if(plot_source){
        plot_unit_source_interior_points_cartesian(output)
    }

    return(output)
}
#'
#' Fill a polygon with grid points (typically used for numerical integration)
#'
#' Given a polygon, we fill it with grid points with the provided approximate
#' spacing, compute an area associated with each grid point such that the total
#' area = polygon area. The area is less for points near the boundary, and the
#' spacing here will generally deviate from a pure regular grid
#' 
#' @param polygon matrix with x,y coordinates defining a polygon
#' @param approx_dx approximate x spacing of points
#' @param approx_dy approximate y spacing of points
#' @param rotation_origin Origin about which the polygon is rotated before
#' filling with x-y aligned grid points. If NULL, the final point in the polygon
#'is used.
#' @param rotation_x_axis_vector Vector defining the rotated x-axis that the
#' polygon is rotated to before filling with x-y aligned grid points. If NULL, the
#' vector joining the final and first points in the polygon is used.
#' @return A list containing grid points in the polygon and other useful information.
#'
#' @export
compute_grid_point_areas_in_polygon<-function(polygon, approx_dx, approx_dy,
    rotation_origin = NULL, rotation_x_axis_vector = NULL){

    lp = length(polygon[,1])

    # Set default rotation origin
    if(is.null(rotation_origin)){
        new_origin = polygon[lp,]
    }else{
        new_origin = rotation_origin
    }

    # Default rotation x-axis vector
    if(is.null(rotation_x_axis_vector)){
        x_axis_vector = polygon[lp,] - polygon[1,]
    }else{
        x_axis_vector = rotation_x_axis_vector
    }
    
    # Ensure that the first/last points are not the same
    stopifnot(!isTRUE(all.equal(x_axis_vector, c(0,0))))

    # Rotate polygon. This helps us ensure the grid points are aligned with the
    # top of the unit source
    rotated_polygon = rotate_cartesian2d(points = polygon, origin = new_origin, 
        x_axis_vector = x_axis_vector)

    # Get extent of the unit source
    grid_extent = rbind(range(rotated_polygon[,1]), 
        range(rotated_polygon[,2]))

    # Make a 'grid' of points inside the unit source
    grid_xs = seq(grid_extent[1,1] - approx_dx/2, 
        grid_extent[1,2] + approx_dx/2, 
        len=max(3, ceiling(diff(grid_extent[1,])/approx_dx)+2))
    grid_ys = seq(grid_extent[2,1] - approx_dy/2, 
        grid_extent[2,2] + approx_dy/2, 
        len=max(3, ceiling(diff(grid_extent[2,])/approx_dy)+2))

    # The grid spacing may not be exactly dx,dy
    dx_local = grid_xs[2] - grid_xs[1]
    dy_local = grid_ys[2] - grid_ys[1]

    grid_points = expand.grid(grid_xs, grid_ys)
    
    # Convert polygon to SpatialPolygons, so we can use rgeos::gIntersection
    p0 = SpatialPolygons(list(Polygons(list(Polygon(polygon)), ID='P')),
        proj4string=CRS(""))

    # Create polygon around each grid point
    grid_pol = list()
    for(i in 1:length(grid_points[,1])){
        gp = grid_points[i,]

        local_polygon_rotated = rbind( 
                   gp + c(-dx_local, -dy_local)/2,
                   gp + c(-dx_local, dy_local)/2,
                   gp + c( dx_local, dy_local)/2,
                   gp + c( dx_local, -dy_local)/2)

        # Get the polygon in 'unrotated' coordiantes
        local_polygon = rotate_cartesian2d(local_polygon_rotated, 
            origin = new_origin, x_axis_vector = x_axis_vector, inverse=TRUE)

        grid_pol[[i]] = Polygons(list(Polygon(local_polygon)),
            ID=as.character(i))
    }

    p1 = SpatialPolygons(grid_pol, proj4string=CRS(""))

    # Alternative methods for computing the intersection of p1 and p0
    exact_intersection = TRUE
    if(exact_intersection){
        p_intersect = gIntersection(p1, p0, byid=TRUE, drop_lower_td = TRUE) 
    }else{
        # Approximate method which however keeps the centroids = rotated grid_points
        pc = coordinates(p1)
        pip = point.in.polygon(pc[,1], pc[,2], polygon[,1], polygon[,2]) 
        keepers = which(pip == 1)
        p_intersect = p1[keepers,]
    }



    # If there are 'just touching' relations and similar, then
    # p_intersect may not by SpatialPolygons. Throw an error for now
    stopifnot(class(p_intersect)=='SpatialPolygons')

    areas = unlist(lapply(p_intersect@polygons, 
        f<-function(x) x@Polygons[[1]]@area))
    indices = unlist(lapply(p_intersect@polygons, f<-function(x) x@ID))
    indices = as.numeric(gsub('P', '', indices))

    centroids = coordinates(p_intersect)

    stopifnot(length(areas) == length(indices))

    return(list(grid_points = centroids, index = indices, area = areas, dx = dx_local, dy = dy_local,
        grid_point_polygon=p_intersect))
}

#' Utility to assist with interpolating depth at points inside a unit source
#'
#' Given one down-dip edge of a unit source, and an overlapping
#' fine_downdip_transect , make a function f(alpha): [0-1] --> (x,y,depth) which
#' gives the cartesian coordinates and depth at fraction 'alpha' of the distance
#' down the unit source edge.
#' 
#' @param unit_source_edge matrix with 2 rows and x,y,depth columns, giving the
#' coordiantes of an edge of the unit source. x,y,depth MUST be local cartesian coordinates (m)
#' @param fine_downdip_transect matrix with many rows and x,y,depth columns
#' which covers the unit source edge, and gives detail to the depth variation.
#' Must be in local cartesian coordinates, with x, y in m, and depth in km
#' @return a function f(alpha) --> (x, y, depth), where alpha is a coordinate
#' in [0,1] giving the fraction distance along that part of the
#' fine_downdip_transect that covers the unit source edge (might be slightly
#' larger than the unit source edge)
#'
#' @export
make_edge_approxfun<-function(unit_source_edge, fine_downdip_transect){
   
    # Extract that part of the 'fine_downdip_transect' that just covers
    # the unit source edge 
    n0 = max(1, 
             min(which(fine_downdip_transect[,3] > unit_source_edge[1,3])) - 1)
    n1 = min(length(fine_downdip_transect[,1]), 
             max(which(fine_downdip_transect[,3] < unit_source_edge[2,3])) + 1)
    t1_subset = fine_downdip_transect[n0:n1,1:3]

    # Compute the distance, normalised to 0-1, along the transect
    dlength = ( diff(t1_subset[,1])**2 + diff(t1_subset[,2])**2 + 
        (diff(t1_subset[,3]))**2 )**0.5
    t1_length = c(0, cumsum(dlength))
    t1_length = t1_length/max(t1_length)

    # Make lookup function to return
    fx = approxfun(t1_length, t1_subset[,1])
    fy = approxfun(t1_length, t1_subset[,2])
    fdepth = approxfun(t1_length, t1_subset[,3])

    fout<-function(alpha){
        return(cbind(fx(alpha), fy(alpha), fdepth(alpha)))
    }

    return(fout)
}


#' Interpolate depth and dip at interior points in the unit source
#'
#' @param unit_source matrix with 4 rows and 3 columns (x,y,depth) defining the
#' unit source boundary. The 4 points should be ordered so that unit_source[1,]
#' and unit_source[4,] are the top edge (shallower), and unit_source[2,] is on
#' the same side of the source as unit_source[1,] (so plotted as a line it looks
#' sensible). The unit source x,y,depth should be in m (local cartesian coords)
#' @param grid_points data.frame with x,y coordinates of the points inside the
#' unit source where we want to find depth/dip/etc
#' @param fine_downdip_transects array of dimension (N, 3, 2).
#' fine_downdip_transects[,,1] contains a line (x,y, depth) which covers the
#' edge unit_source[1:2,], while fine_downdip_transects[,,2] contains a similar
#' line which covers the edge unit_source[4:3,]. The x,y,depth coordinates should be in
#' m in local cartesian coordinates
#' @param grid_points_strike The strike for each grid point
#' @return matrix with x, y, depth, dip,.. etc at points inside the unit source
#' 
#' @export
get_depth_dip_at_unit_source_interior_points<-function(
    unit_source, grid_points, fine_downdip_transects, 
    grid_points_strike){

    # Make functions which take alpha [0,1] and return
    # x,y,depth along the unit source edge. We will use these to interpolate
    # the depth
    f1 = make_edge_approxfun(unit_source[1:2,], fine_downdip_transects[,,1])
    f2 = make_edge_approxfun(unit_source[4:3,], fine_downdip_transects[,,2])

    # Make interpolation function across unit source
    # Input (alpha, s), both in [0,1]
    line_across_source<-function(alpha, s){

        p1 = f1(alpha)
        p2 = f2(alpha)
    
        output = p1 + s*(p2 - p1)

        return(output)
    }

    #library(minpack.lm)
    # Define a function that gives alpha/s associated with a given grid point
    find_alpha_s<-function(initial_guess, grid_point){

        # Function giving the distance**2 between the desired x,y point
        # and the line at alpha,s. We minimise this to find the best alpha/s
        minme<-function(alpha_s){
            #sum((line_across_source(alpha_s[1], alpha_s[2])[1:2] - grid_point[1:2])**2)
            (line_across_source(alpha_s[1], alpha_s[2])[1:2] - grid_point[1:2])
        }

        # Find (alpha, s) at grid_points[i,1:2]
        #fit_alpha_s = optim(initial_guess, minme, method='BFGS', 
        #    control = list(reltol=1e-8, abstol=1.0e-8))
        fit_alpha_s = nls.lm(initial_guess, lower=c(-0.1, -0.1), 
            upper=c(1.1, 1.1), minme)

        return(fit_alpha_s)
    }

    # Find the 'alpha' and 's' values corresponding to each grid point inside
    # the unit source (and also a perturbation of that grid point in the down-dip direction)
    alpha_s = as.matrix(grid_points[,1:2])*0 + 0.5 # Initial values -- reasonable for search
    alpha_s_perturbation = alpha_s

    deg2rad = pi/180
    strike_perp = (grid_points_strike + 90)*deg2rad
    perturbation_scale = 100
    grid_points_perturb = grid_points[,1:2] + cbind(sin(strike_perp), cos(strike_perp))*perturbation_scale

    for(i in 1:length(alpha_s[,1])){
        fit_alpha_s = find_alpha_s(alpha_s[i,], grid_points[i,1:2])
        alpha_s[i,] = fit_alpha_s$par

        alpha_s_perturbation[i,] = 
            find_alpha_s(fit_alpha_s$par, grid_points_perturb[i,1:2])$par
    }

    # Check for alpha/s which are outside [0,1]. Very small excursions are ok
    # (due to the great-circle fine_downdip_transects having slight disagreement
    # with the straight-line along the unit source edge -- the latter are used
    # for point-in-polygon checking).
    #
    # It's ok to keep points affected by that, since they will be missing from the
    # neighbouring unit source (which arguably should have them), and they will
    # have a depth/dip similar to the edge transect depth/dip
    adjusters = which(alpha_s < 0)
    if(length(adjusters) > 0){
        if(min(alpha_s[adjusters]) > -5.0e-02){
            #alpha_s[adjusters] = 0
        }else{
            msg = paste0('bad interpolation a:', min(alpha_s[adjusters]))
            stop(msg)
        }
    }
    adjusters = which(alpha_s > 1)
    if(length(adjusters) > 0){
        if(max(alpha_s[adjusters]) < (1.0 + 5.0e-02) ){
            #alpha_s[adjusters] = 1
        }else{
            msg = paste0('bad interpolation b:', max(alpha_s[adjusters]))
            stop(msg)
        }
    }

    # Get the x,y,depth points and perturbed points
    output_points = line_across_source(alpha_s[,1], alpha_s[,2])
    output_points_perturb = line_across_source(alpha_s_perturbation[,1], 
        alpha_s_perturbation[,2])

    # Now compute the dip numerically [ tan(dip) = change_in_depth/distance ]
    dip = atan((output_points_perturb[,3] - output_points[,3])/perturbation_scale)/deg2rad

    #print('GD: Deliberately breaking dip')
    #dip = dip*0 + mean_angle(dip)

    output_points = cbind(output_points, dip, alpha_s)

    return(output_points)

}

#' Quick plot of a unit source interior points
#'
#' @param us List with unit source interior points information: Output of
#' 'unit_source_interior_points_cartesian' or similar
#' @return nothing but make a plot
#'
#' @export
plot_unit_source_interior_points_cartesian<-function(us){

    unit_source = us$unit_source_cartesian
    output_points = us$grid_points

    col_value = (max(unit_source[,3]) - unit_source[,3])/diff(range(unit_source[,3])) 
    plot(unit_source[,1:2], t='o', 
        col = rgb(col_value**0.5, (1-col_value)**2, col_value**2), 
        pch=1, asp=1, 
        main='Unit source integration points: color ~ depth; \n Orange angle ~ dip ; black line ~ (strike - rake)')

    col_value = (max(unit_source[,3]) - output_points[,3])/diff(range(unit_source[,3])) 
    points(output_points[,1:2], pch=19, 
        col = rgb(col_value**0.5, (1-col_value)**2, col_value**2),
        cex = 1)


    # Add lines in the up-dip direction
    dx = us$dx    
    dy = us$dy
    ds = 0.5 * sqrt(dx**2 + dy**2)

    perturb = ds * cbind(sin(us$grid_points[,'strike']/180*pi - pi/2),
        cos(us$grid_points[,'strike']/180*pi - pi/2))
    arrows(output_points[,1], output_points[,2], output_points[,1] + perturb[,1],
        output_points[,2] + perturb[,2], length=0)

    # Add angles which show the dip
    for(i in 1:length(output_points[,1])){
        
        l1 = cbind( c(output_points[i,1], output_points[i,1] + 2*ds/3),
                    c(output_points[i,2], output_points[i,2]))
        points(l1, t='l', col='orange', lwd=0.5)

        sindip = sin(output_points[i,'dip']/180*pi)
        l2 = cbind(c(output_points[i,1], output_points[i,1] + 2*ds/3),
                   c(output_points[i,2], output_points[i,2] - 2*ds/3*sindip))
        points(l2, t='l', col='orange', lwd=0.5)
    }

}

#'
#' 3D rgl plot of the unit source interior points
#'
#' @param us List with unit source interior points information: Output of
#' 'unit_source_interior_points_cartesian' or similar
#' @param aspect aspect ratio of plot (passed to rgl::plot3d)
#' @param add logical Add to an existing plot?
#' @param add_zero_plane Draw a plane at depth = 0
#' @param ... further arguments to plot3d (or points3d if add=TRUE)
#' @return Nothing, but make a plot
#' 
#' @export
plot3d_unit_source_interior_points_cartesian<-function(us, aspect='iso', 
    add=FALSE, add_zero_plane=TRUE, ...){

    if(!add){
        rgl::plot3d(us$grid_points[,1], us$grid_points[,2], -us$grid_points[,3], 
            aspect=aspect, add=add, ...)
    }else{
        rgl::points3d(us$grid_points[,1], us$grid_points[,2], 
            -us$grid_points[,3], ...)
    }

    rgl::polygon3d(us$unit_source_cartesian[,1], us$unit_source_cartesian[,2], 
            -us$unit_source_cartesian[,3], fill=FALSE)

    m = us$fine_downdip_transects_cartesian
    rgl::lines3d(m[,1,1], m[,2,1], -m[,3,1], col='red')
    rgl::lines3d(m[,1,2], m[,2,2], -m[,3,2], col='red')

    if(add_zero_plane){
        rgl::polygon3d(x = us$unit_source_cartesian[,1], 
            y = us$unit_source_cartesian[,2],
            z = rep(0, 4), fill=TRUE, col='blue', alpha=0.2)
    }

    # Get the pan3d function -- suppress printing of the example where it is
    # defined
    # To run the example we still need 'require'
    require(rgl)
    tmp = capture.output(example(rgl.setMouseCallbacks))
    tmp = capture.output(pan3d(3))

}

