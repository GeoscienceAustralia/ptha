# 
# Code to create unit sources from subduction interface contours
#
# Gareth Davies, Geoscience Australia, 2015
#

#' Enhance near-trench orthogonality on discretized source [deprecated]
#'
#' Suppose top-line is a set of lon/lat points defining the top of a grid of
#' unit sources, and second line is a set of lon/lat points defining the bottom
#' of the unit sources. If top-line is at the trench, then numerically it can be
#' desirable for the unit sources to be highly orthogonal where they intersect
#' the trench. One way to do this is by moving points along the top line. That's
#' what this function does. Note the approach has been superceded by an improved
#' algorithm to generate the unit-sources, which is invoked by passing 
#' improved_downdip_lines=TRUE to \code{discretized_source_from_source_contours}
#'
#' @param top_line matrix with 2 rows and n columns. First row is longitude,
#' second row is latitude
#' @param second_line as top_line, but for the second line
#' @return an alternative top_line which should have better orthogonality
orthogonal_near_trench<-function(top_line, second_line){

    top_line_bearing = top_line[1,]*0
    n = length(top_line_bearing)

    top_line_bearing[2:(n-1)] = bearing(t(top_line[1:2,1:(n-2)]), t(top_line[1:2, 3:n]), f=0)
    top_line_bearing[1] = bearing(top_line[1:2,1], top_line[1:2,2], f=0)
    top_line_bearing[n] = bearing(top_line[1:2,n-1], top_line[1:2,n], f=0)

    # Make function to interpolate along the line, and also give the tangent
    # bearing
    make_approx_top_line<-function(){
        approx_top_line_xfun = approxfun(1:n, top_line[1,])
        approx_top_line_yfun = approxfun(1:n, top_line[2,])
        approx_bearing = approxfun(1:n, top_line_bearing)
        
        approx_top_line<-function(n){
            return(cbind(approx_top_line_xfun(n), approx_top_line_yfun(n), approx_bearing(n)))
        }

        return(approx_top_line)
    }

    approx_top_line = make_approx_top_line()

    # Find the new top of grid point
    best_n = rep(NA, n)
    for(i in 2:(n-1)){

        bearing_end_point<-function(n){
            p1 = approx_top_line(n)
            p2 = destPoint(p1[1:2], (p1[3] + 90), d = 1000*500)
            point_distance = dist2Line(second_line[1:2,i], rbind(p1[1:2], p2))[1]
            return(point_distance)
        }

        best_n[i] = optimize(bearing_end_point, interval = c(i-0.4, i+0.4))$minimum
    }
    best_n[1] = 1
    best_n[n] = n

    return(approx_top_line(best_n)[,1:2])
}


#' Make discretized_source from subduction interface contours
#'
#' Given a shapefile with contour lines in lon/lat coordinates defining the
#' source-zone subduction interface depth, partition the source-zone into
#' unit sources with chosen approximate length/width. Information on all the
#' unit sources is held in a 'discretized_source' list. \cr 
#' Optionally, the user can provide a line shapefile defining the along-strike
#' boundaries of the unit sources. Otherwise, they are automatically created.
#' In any case care may be required to ensure good quality unit-source
#' boundaries are created. A good strategy is to make them automatically as
#' a first iteration, then save the corresponding down-dip lines as a line shapefile.
#' This can be achieved by passing the mid_lines_with_cutpoints output of the
#' current routine to  \code{downdip_lines_to_SpatialLinesDataFrame}. Then you
#' may optionally edit the latter shapefile in GIS, and subsequently pass it
#' to this routine directly.
#' 
#' @param source_shapefile character filename of the line shapefile defining the
#' subduction interface contours. It should have an attribute giving the depth.
#' Alternatively, source_shapefile can be a \code{SpatialLinesDataFrame} obtained
#' by reading the file with \code{rgdal::readOGR}. 
#' @param desired_subfault_length numeric desired length of subfaults (km)
#' @param desired_subfault_width numeric desired width of subfaults (km)
#' @param make_plot logical Make a plot?
#' @param contour_depth_attribute character The name of the column in the
#' source_shapefile attribute table giving the contour depth
#' @param contour_depth_in_km logical Are contour depths given in km? (If False,
#' assume 'm')
#' @param extend_line_fraction To ensure that contour lines intersect downdip
#' lines at the left/right edges we extend (or buffer) them by this fraction of the
#' end-to-end source length. Not required if 'improved_downdip_lines=TRUE' 
#' @param orthogonal_near_trench move unit source points along the trench to enhance
#' orthogonality there. Can reduce numerical artefacts at the trench. 
#' Not used if 'improved_downdip_lines=TRUE' or if downdip_lines is provided (we
#' encourage the latter approaches)
#' @param improved_downdip_lines If TRUE, use an iterative algorithm for
#' computing the downdip lines (which define unit source boundaries). This
#' should lead to unit-sources that are more orthogonal, which is normally a 'good thing'. 
#' The iterative algorithm can take a few minutes if there are hundreds of unit-sources in
#' each along-trench direction. In general we do not suggest setting this to
#' FALSE (except for compatibility with early versions of the code).
#' @param downdip_lines Either NULL, or the name of a Line Shapefile, or a 
#' SpatialLinesDataFrame derived by reading the latter. It should contain lines
#' which cross the source contours in a down-dip direction, with a single attribute 
#' giving the line order in the along strike direction. If not NULL, then these
#' lines define the along-strike boundaries of the unit-sources. If NULL and 
#' 'improved_downdip_lines=TRUE', then the code uses
#' \code{create_downdip_lines_on_source_contours_improved} to auto-generate downdip
#' lines. The latter routine can be called separately to make downdip
#' lines, which can be manually edited before passing to this function. 
#' @return A list containing: depth_contours The original source contours;
#' unit_source_grid A 3 dimensional array descrbing the unit source vertices;
#' discretized_source_dim A vector of length 2 with number-of-sources-down-dip,
#' number-of-sources-along-strike; fine_downdip_transects A 3 dimensional array
#' containing densly spaced points along the down-dip transects, which might be
#' useful for defining sub-unit-source points for tsunami source integration
#'
#' @export
#'
#' @examples
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip,
#'    make_plot=TRUE)
#'
discretized_source_from_source_contours<-function(
    source_shapefile, 
    desired_subfault_length,
    desired_subfault_width, 
    make_plot=FALSE,
    contour_depth_attribute='level', 
    contour_depth_in_km=TRUE,
    extend_line_fraction=1.0e-06,
    orthogonal_near_trench = FALSE,
    improved_downdip_lines = TRUE,
    downdip_lines=NULL){

    # Get the shapefile
    if(class(source_shapefile) == 'SpatialLinesDataFrame'){
        source_contours = source_shapefile
    }else{
        if(!file.exists(source_shapefile)) stop('cannot find shapefile')
        source_contours = rgdal::readOGR(
            dsn = dirname(source_shapefile),
            layer=gsub('.shp', '', basename(source_shapefile)),
            verbose=FALSE)
    }

    if(is.null(downdip_lines)){
        # Make the mid_line_with_cutpoints automatically

        if(!improved_downdip_lines){
            mid_line_with_cutpoints = create_downdip_lines_on_source_contours(
                source_contours,
                desired_subfault_length, 
                contour_depth_attribute, 
                extend_line_fraction,
                orthogonal_near_trench,
                make_plot)
        }else{
            mid_line_with_cutpoints = create_downdip_lines_on_source_contours_improved(
                source_contours,
                desired_subfault_length, 
                contour_depth_attribute,
                make_plot=make_plot)
        }

    }else{
        # Use the provided downdip_lines to make the mid_line_with_cutpoints

        if(is.character(downdip_lines)){
            if(!file.exists(downdip_lines)){
                stop(paste0('Could not find file defined by downdip_lines ', downdip_lines))
            }
            downdip_lines = readOGR(dsn=downdip_lines, 
                layer=gsub('.shp', '', basename(downdip_lines)), verbose=FALSE)
        }

        mid_line_with_cutpoints = mid_line_with_cutpoints_from_downdip_sldf_and_source_contours(
            source_contours, 
            downdip_lines, 
            contour_depth_attribute=contour_depth_attribute, 
            buffer_width=extend_line_fraction,
            make_plot=make_plot)

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

    # Decide how many unit-sources there should be in the down-dip direction
    # (the number along-strike is already determined from mid_line_with_cutpoints)
    mean_dip_cut_length = mean(dip_cut_lengths)/1e+03 # km
    number_of_mid_lines = max(
        round(mean_dip_cut_length / desired_subfault_width - 1), 
        0)

    strike_cuts = list()   
    ll = length(mid_line_with_cutpoints)
    for(i in 1:ll){
        # Now interpolate along the above, with the desired spacing between
        # points
        interpolated_midline = interpolate_3D_path(
            mid_line_with_cutpoints[[i]], 
            n = number_of_mid_lines + 2,
            depth_in_km = contour_depth_in_km)
        strike_cuts[[i]] = interpolated_midline
    }

    fine_strike_cuts = list()
    for(i in 1:ll){

        fine_interpolated_midline = interpolate_3D_path(
            mid_line_with_cutpoints[[i]], 
            n = 100,
            depth_in_km = contour_depth_in_km)
       
        if(make_plot){
            points(strike_cuts[[i]], col='purple', pch=19)
            points(mid_line_with_cutpoints[[i]][,1:2], col='pink', pch=19, 
                cex=0.2)
            points(mid_line_with_cutpoints[[i]][,1:2], t='l', col='brown', 
                lty = 'dashed') 
        }
        
        fine_strike_cuts[[i]] = fine_interpolated_midline
    }

    # Correct the 'fine' interpolated midline to match the strike_cuts
    if(orthogonal_near_trench){
        for(i in 1:ll){
            fine_strike_cuts[[i]][1,1:2] = strike_cuts[[i]][1,1:2]        
        }
    }

    # Plot
    if(make_plot){
        for(i in 1:length(interpolated_midline[,1])){
            # Get lines along-strike
            line2plot = matrix( unlist(lapply(strike_cuts, f<-function(x) x[i,])), 
                ncol = 3, byrow = TRUE)
            # Interpolate them along a great circle with the 3D line interpolation
            # routine
            line2plot = cbind(line2plot[,1:2], rep(0, length(line2plot[,1])))
            line2plot = interpolate_3D_path(line2plot)[,1:2]
           
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
        fine_downdip_transects = fine_unit_source_grid,
        mid_line_with_cutpoints = mid_line_with_cutpoints)

    return(unit_source_data)
}


#' Compute APPROXIMATE summary statistics for all unit sources in a discretized source
#' zone
#'
#' The quantities that we output follow those in the i-invall format in the old
#' URSGA code, although some outputs there (e.g. section number) do not apply
#' with our algorithms. \cr
#' The approach used in this routine is faster than discretized_source_summary_statistics,
#' at the expense of accuracy [the latter draws on the sub-unit-source grid points].
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
#' @examples
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip)
#'
#' # Make summary statistics
#' puysegur_summary_statistics = discretized_source_approximate_summary_statistics(
#'     puysegur_discretized_source)
#' print(head(puysegur_summary_statistics))
#'
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
                bearing(midpoint, source_coords[4,1:2], f=0)%%360)
            strike[counter] = st0 
            stopifnot( (strike[counter] >= 0) & (strike[counter] < 360)) 
  
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

        arrows(lon_c, lat_c, end_points[,1], end_points[,2], code=3, length=0, 
            col='orange')

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
#' so that (width x length) = (dipping interface area). \cr
#' Note that the statistics will not account for edge tapering that may be
#' optionally applied to the slip (i.e. the area is entirely based on the interior
#' of the earthquake unit source, so does not include any extra regions
#' affected by smoothing). 
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
#' @examples
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip)
#'
#' puysegur_stats = discretized_source_summary_statistics(puysegur_discretized_source, 
#'    approx_dx = 5000, approx_dy = 5000)
#' head(puysegur_stats)
#'
#' # Compare with 'approximate' summary statistics, which don't use sub-grid points
#' # Most obvious difference occurs in depth for up-dip unit-sources, which have a fairly
#' # non-linear (concave-down parabolic-ish) profile. This causes the
#' # 'approximate' mean depth to be larger than the mean depth based on sub-grid points.
#' puysegur_stats_approx = discretized_source_approximate_summary_statistics(
#'     puysegur_discretized_source)
#' summary(puysegur_stats_approx/puysegur_stats)
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

#' Coordinates of a single unit source
#'
#' Convenience function to extract coordinates of a single unit source from the
#' discretized_source list
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
                  fine_downdip_transects = fine_downdip_transects
                  )
    return(output)
}

#' Get strike along top edge of shallowest unit sources
#'
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
    strike = suppressWarnings(bearing(top_edges, vertex_tr, f=0))

    # Ensure strike in 0-360 [not -180 to 180]
    strike = strike + (strike < 0)*360

    # Correct the top_edges hemisphere
    top_edges = adjust_longitude_by_360_deg(top_edges, vertex_tr)

    return(cbind(top_edges, strike))
}

#' Outline of the discretized source
#'
#' Convenience function to get the outline of the discretized source from its
#' unit_source_grid
#'
#' @param discretized_source list containing unit sources information (e.g.
#' output of discretized_source_from_source_contours)
#' @return matrix defining the 'outline' polygon of the discretized source
#'
#' @export 
get_discretized_source_outline<-function(discretized_source){
    
    unit_source_grid = discretized_source$unit_source_grid
    ndip = discretized_source$discretized_source_dim['dip']
    nstrike = discretized_source$discretized_source_dim['strike']

    top_line = t(unit_source_grid[1,,])
    bottom_line = t(unit_source_grid[ndip+1,,])
    left_line = unit_source_grid[,,1]
    right_line = unit_source_grid[,, nstrike+1]

    # Drop the 'first' point from the second matrix in rbind to avoid repeats 
    output_grid = rbind(top_line, right_line[-1,, drop=FALSE])
    # (note bottom_line has (nstrike+1) rows, so we are still dropping a point here)
    output_grid = rbind(output_grid, bottom_line[nstrike:1,, drop=FALSE])
    # Need to drop 2 points from left_line
    if(ndip > 1){
        # (note left_line has (ndip+1) rows, so we are dropping 2 points here)
        output_grid = rbind(output_grid, left_line[ndip:2,, drop=FALSE])
    }
    
    return(output_grid) 

}


#' Fill a unit source with "sub-unit-source" points (for tsunami source integration)
#'
#' Make a locally linear coordinate system for a unit source and compute a
#' 'grid' of points inside the source. These can be used later for tsunami
#' source integration
#'
#' @param discretized_source list containing unit sources information (e.g. output of
#' discretized_source_from_source_contours). Coordinates must be in lon/lat
#' @param unit_source_index vector of length 2 with integer index (down-dip,
#' along-strike) of the unit source to operate on.
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
#' @param edge_taper_width m. Use tapering of the slip boundaries. This extends
#' the sub-unit-source grid points beyond the unit-source 'box' by a distance
#' edge_taper_width. The slip reduces linearly from 1 to 0 for points near the boundary
#' (ranging from -edge_taper_width to +edge_taper_width from the unit-source edges)
#' @return list containing: 'unit_source': matrix with unit_source_coords in the
#' cartesian coordinate system + an additional depth column (m); 'origin' the
#' lon/lat origin of the local coordinate system; r (input argument); 'dx, dy'
#' the actual dx, dy of the grid points (not exactly the same as the approximate
#' values); 'grid_points' grid points inside the source in the cartesian
#' coordinate system; 'grid_point_polygons' a polygon associated with the grid points
#'
#' @export
unit_source_interior_points_cartesian<-function(
    discretized_source,
    unit_source_index,
    origin = NULL,
    r = 6378137,
    approx_dx = NULL,
    approx_dy = NULL,
    scale_dxdy = 1,
    depths_in_km = TRUE,
    plot_source=FALSE,
    edge_taper_width = 0){

    # Get the unit source coordinates, and down-dip transects we
    # can use for interpolation
    unit_source_info = get_unit_source_from_discretized_source(
        discretized_source, 
        unit_source_index)

    # We will change these variables from spherical to cartesian coordinates
    unit_source_coords = unit_source_info$unit_source_grid

    if(is.null(origin)){
        origin = unit_source_coords[1,1:2]
    }

    # Convert the unit source lon/lat coordinates to a local 2D
    # cartesain system
    unit_source_cartesian = unit_source_coords
    unit_source_cartesian[,1:2] = spherical_to_cartesian2d_coordinates(
        unit_source_coords[,1:2], origin_lonlat = origin, r = r)

    # We need the 'full' unit source grid if we use slip tapering, since
    # we need to be able to find sub-grid points inside neighbouring polygons
    full_unit_source_grid_cartesian = discretized_source$unit_source_grid
    for(i in 1:dim(full_unit_source_grid_cartesian)[3]){
        full_unit_source_grid_cartesian[,1:2,i] = 
            spherical_to_cartesian2d_coordinates(
                full_unit_source_grid_cartesian[,1:2,i], origin_lonlat = origin, 
                r = r)
    }

    if(depths_in_km){
        unit_source_cartesian[,3] = 1000*unit_source_cartesian[,3]
        full_unit_source_grid_cartesian[,3,] = 1000 * 
            full_unit_source_grid_cartesian[,3,]
    }

    # Adjust approx_dx, approx_dy
    if(is.null(approx_dx)){
        approx_dx = min(unit_source_cartesian[,3])*0.5
    }
    if(is.null(approx_dy)){
        approx_dy = min(unit_source_cartesian[,3])*0.5
    }

    approx_dx = approx_dx*scale_dxdy
    approx_dy = approx_dy*scale_dxdy

    # Replace lon/lat spherical coordinates in 'fine-transects'
    # with local cartesians
    #
    # Currently only used for output (consider removing)
    fine_downdip_transects_cartesian = unit_source_info$fine_downdip_transects
    fine_downdip_transects_cartesian[,1:2,1] = spherical_to_cartesian2d_coordinates(
        fine_downdip_transects_cartesian[,1:2,1], origin_lonlat = origin, r=r)
    fine_downdip_transects_cartesian[,1:2,2] = spherical_to_cartesian2d_coordinates(
        fine_downdip_transects_cartesian[,1:2,2], origin_lonlat = origin, r=r)
    if(depths_in_km){
        fine_downdip_transects_cartesian[,3,1:2] = 
            1000*fine_downdip_transects_cartesian[,3,1:2]
    }
    

    # Get a 'grid' of points inside the unit source (later used
    # for integration)
    grid_point_data = compute_grid_point_areas_in_polygon(
        unit_source_cartesian[,1:2], 
        approx_dx = approx_dx, 
        approx_dy = approx_dy,
        edge_taper_width=edge_taper_width,
        full_unit_source_grid = full_unit_source_grid_cartesian)

    # NOTE: grid points which might fall outside the unit source boundaries,
    # because of edge tapering (if edge_taper_width > 0). This is deliberate.
    grid_points = grid_point_data$grid_points_buffer

    # Test for consistency between the grid point areas and the original unit
    # source. Do NOT use the buffered areas
    a0 = sum(grid_point_data$area)
    a1 = areaPolygon(unit_source_coords[,1:2], f=0)
    rel_err = abs(a0-a1)/((a0+a1)*0.5)
    if( rel_err > 0.01 ){
        msg = paste0('Cartesian and spherical unit source area ',
            "(surface only, not 3d) differ by a fraction ",
            signif(rel_err,3), sep=" ")
        warning(msg)
    }

    
    # Get strike at the grid points
    # 

    ## APPROACH 1: Strike is discontinuous between neighbouring unit-sources
    #top_edge_vector = unit_source_cartesian[4,1:2] - unit_source_cartesian[1,1:2]
    ## Since the sub-unit-sources are created to be aligned with the top edge vector,
    ## we should use that strike. Degrees from north
    #strike_top_edge = atan2(top_edge_vector[1], top_edge_vector[2])/pi*180
    #strike = rep(strike_top_edge, length(grid_points[,1]))
    ## Ensure strike is > 0
    #strike = strike + (strike < 0)*360

    #
    # Convert lon-lat centroids of ALL unit sources in the discrete source to
    # the local coordinate system and use them to make a continuous function of
    # strike

    ds1_lonlatstrike = get_shallow_unit_source_top_edge_strikes(discretized_source) 
    ds1_stats_points_cartesian = spherical_to_cartesian2d_coordinates(
        ds1_lonlatstrike[,1:2], origin_lonlat=origin, r=r)

    # Ideally we interpolate strike using a weighted nearest-neighbours
    # However, the code used here requires that we have > 1 point
    np = length(ds1_lonlatstrike[,1])
    k = min(11, np)
    # Compute a distance**(negative_power) weighted angular mean of k nearest neighbours
    inds = knnx.index(data = ds1_stats_points_cartesian[,1:2, drop=FALSE], 
        query = grid_points[,1:2, drop=FALSE], k=k)
    dists = knnx.dist(data = ds1_stats_points_cartesian[,1:2, drop=FALSE],
        query = grid_points[,1:2, drop=FALSE], k=k)
    mean_strike = dists[,1]*NA 
    for(ii in 1:length(mean_strike)){
        mean_strike[ii] = mean_angle(ds1_lonlatstrike[inds[ii,], 3], 
            weights = 1/(dists[ii,])**3)
    }
    strike = mean_strike

    # Ensure strike is > 0
    strike = strike + (strike < 0)*360

    #
    # Get dip/depth etc at the grid points
    # Use the contour_function
    # Easiest if we convert the mid_line_with_cutpoints to cartesian
    # coordinates before starting 
    #
    mid_line_with_cutpoints_cartesian = discretized_source$mid_line_with_cutpoints
    for(i in 1:length(mid_line_with_cutpoints_cartesian)){
        mid_line_with_cutpoints_cartesian[[i]][,1:2] = 
            spherical_to_cartesian2d_coordinates(
                mid_line_with_cutpoints_cartesian[[i]][,1:2],
                origin_lonlat = origin, 
                r=r)
        if(depths_in_km){
            mid_line_with_cutpoints_cartesian[[i]][,3] = 1000 * 
                mid_line_with_cutpoints_cartesian[[i]][,3]
        }
    }

    contour_fun = make_contour_interpolator(mid_line_with_cutpoints_cartesian)

    # NOTE: It is possible for points to be outside the grid defined by
    # mid_line_with_cutpoints_cartesian, but inside the grid defined by the
    # unit_source edges. This is because the former can be jagged (particularly
    # at lateral boundaries). So better use 'allow_outside=TRUE'
    new_depths = contour_fun(grid_points, allow_outside=TRUE)

    # Compute dip -- first make some 'down-dip' points, then find their depth,
    # and do the computation
    dx_numerical = approx_dx/2
    deg2rad = pi/180
    grid_points_perturb = grid_points * 0
    grid_points_perturb[,1] = dx_numerical*cos(strike*deg2rad)
    grid_points_perturb[,2] = -dx_numerical*sin(strike*deg2rad)
    depth_perturb = contour_fun(
        grid_points, 
        allow_outside=TRUE,
        xy_perturbation_m = grid_points_perturb)
    if(any(depth_perturb < new_depths)) stop('Negative dip')
    if(any(new_depths < 0) | any(depth_perturb < 0)) stop('Negative depth')
    dip = atan((depth_perturb - new_depths)/dx_numerical)/deg2rad

    # Now that we know dip, we can properly normalise the unit_slip_scale to
    # preserve seismic moment. Since we record dip in the 'buffer' points, we
    # indirectly compute area inside the original unit source by multiplying
    # the buffer areas by their fraction inside the unit source, before summing.
    grid_point_unit_slip_scale = grid_point_data$unit_slip_scale
    a0 = sum(grid_point_data$area_buffer * 
        grid_point_data$area_buffer_fraction_inside_unit_source * 
        sqrt(1 + tan(dip*deg2rad)**2))
    a1 = sum(grid_point_unit_slip_scale * 
        grid_point_data$area_buffer * 
        sqrt(1 + tan(dip*deg2rad)**2))
    grid_point_unit_slip_scale = grid_point_unit_slip_scale * a0/a1

    # Make output with x,y,depth,dip, area_projected, strike, unit_slip_scale
    grid_point_areas = grid_point_data$area_buffer
    grid_point_polygon = grid_point_data$grid_point_polygon_buffer
    fraction_area_inside_unit_source = grid_point_data$area_buffer_fraction_inside_unit_source

    grid_points = cbind(grid_points, new_depths, dip)
    grid_points = cbind(grid_points, grid_point_areas, strike, grid_point_unit_slip_scale, 
        fraction_area_inside_unit_source)
    # Fix row/column names for later ease of access
    colnames(grid_points) = c('x', 'y', 'depth', 'dip', 'area_projected', 'strike', 
        'unit_slip_scale', 'fraction_area_inside_unit_source')

    rownames(grid_points) = NULL

    # Make a list to hold the unit source interior points information
    output = list(
        unit_source_cartesian=unit_source_cartesian, 
        origin_lonlat=origin, 
        r=r, 
        dx=grid_point_data$dx, 
        dy=grid_point_data$dy, 
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
#' area = polygon area. If edge_tapering is used, the points may extend outside
#' the polygon, and they are assigned a unit_slip_scale giving the slip corresonding
#' to an overall 1m slip (not exactly 1m around edges so that result is smooth).
#' 
#' @param polygon matrix with x,y coordinates defining a polygon
#' @param approx_dx approximate x spacing of points in same units as x,y
#' @param approx_dy approximate y spacing of points in same units as x,y
#' @param rotation_origin Origin about which the polygon is rotated before
#' filling with x-y aligned grid points. If NULL, the final point in the polygon
#' is used.
#' @param rotation_x_axis_vector Vector defining the rotated x-axis that the
#' polygon is rotated to before filling with x-y aligned grid points. If NULL, the
#' vector joining the final and first points in the polygon is used.
#' @param edge_taper_width Length scale over which slip should be smoothed at 
#' the edges of the unit-source. This will expand slip outside the original polygon,
#' and reduce slip near the polygon edges.
#' @param bounding_polygon If not NULL, then a 2 column matrix of coordinates 
#' defining a polygon in which the output grid_points are forced to be inside. 
#' Should be used in conjunction with edge_taper_width > 0 to ensure points near 
#' boundaries are weighted correctly (note: This could just be recomputed with
#' full_unit_source_grid)
#' @param full_unit_source_grid Polygons for the entire sourcezone, stored in an 
#' array as a unit_source_grid. This is needed if edge_taper_width>0, since we
#' need to find sub-grid points inside neighbouring polygons. 
#' @return A list containing grid points in the polygon and other useful information.
#'
#' @export
compute_grid_point_areas_in_polygon<-function(polygon, approx_dx, approx_dy,
    rotation_origin = NULL, rotation_x_axis_vector = NULL,
    edge_taper_width=0, bounding_polygon=NULL, full_unit_source_grid=NULL){

    p0 = SpatialPolygons(list(Polygons(list(Polygon(polygon)), ID='P')),
        proj4string=CRS(""))

    if(!is.null(full_unit_source_grid)){
        bp = as(
            unit_source_grid_to_SpatialPolygonsDataFrame(full_unit_source_grid), 
            'SpatialPolygons')
        bp = gUnaryUnion(bp)
    }

    lp = length(polygon[,1])

    if(edge_taper_width == 0){
        # Simple case -- just fill 'polygon' with points, all having unit_slip_scale=1
        p1 = fill_polygon_with_grid_points(polygon, rotation_origin, 
            rotation_x_axis_vector, approx_dx, approx_dy)
    
        unit_slip_scale = rep(1, length(p1))
    }else{
        # Difficult case.
        #
        # Loop over all polygons in the full unit source grid, and find if they
        # intersect with ours. If they do, fill them with points, and do some
        # work to figure out unit_slip_scale weightings
        #
        if(is.null(full_unit_source_grid)){
            stop('Must provide full unit source grid with edge_taper_width>0')
        }
        p0_buf = gBuffer(p0, width=edge_taper_width)
        ndip = dim(full_unit_source_grid)[1] - 1
        nstrike = dim(full_unit_source_grid)[3] - 1
        p1_list = list()
        counter = 0
        for(i in 1:ndip){
            for(j in 1:nstrike){
                local_pol = rbind(
                        full_unit_source_grid[i:(i+1),1:2,j],
                        full_unit_source_grid[(i+1):(i),1:2,j+1])

                local_pol_sp = SpatialPolygons(list(Polygons(list(Polygon(
                    local_pol)), ID='1')))

                if(gIntersects(p0_buf, local_pol_sp)){
                    # Need to include points in this polygon
                    # FIXME: Might need approx_dx, approx_dy from the other
                    # polygon for accuracy (which might not be the same as the
                    # current values). Check.
                    local_p1 = fill_polygon_with_grid_points(
                        local_pol, rotation_origin, rotation_x_axis_vector, 
                        approx_dx, approx_dy)
                    counter = counter + 1
                    p1_list[[counter]] = local_p1
                }
            }
        }

        # Merge all the points into one SpatialPolygons object
        if(counter == 0){
            stop('Mismatch between unit-source-polygon and unit-source-grid')
        }else{
            p1 = p1_list[[1]]
            if(counter > 1){
                for(i in 2:counter){
                    p1 = rbind(p1, p1_list[[i]], makeUniqueIDs=TRUE)
                }
            }
        }

        # Now for every point, take a circular buffer with radius = edge_taper_width, 
        # and find the fraction of its area inside the original p0 polygon.
        # This gives the unit-slip-scale. Note this is basically a circular
        # moving-average filter, applied to a set of points with '1' inside p0, 
        # and 0 elsewhere
        unit_slip_scale = rep(NA, length(p1))
        point_bufs = gBuffer(gCentroid(p1, byid=TRUE), width=edge_taper_width, quadsegs=5, byid=TRUE)
        point_bufs_intersect = gIntersection(point_bufs, bp, byid=TRUE)
        if(length(point_bufs) != length(point_bufs_intersect)){
            stop('BUG: Some sub-unit-source points not inside bounding polygon')
        }

        for(i in 1:length(p1)){
            #point_buf = gBuffer(gCentroid(p1[i]), width=edge_taper_width, 
            #    quadsegs=10)

            # Need this to prevent points near the trench from having
            # slip suppressed
            #point_buf = gIntersection(point_buf, bp)
            #if(is.null(point_buf)) stop('BUG: Logically, point_buf should be inside bounding_polygon.')

            point_buf_intersect_i = gIntersection(p0, point_bufs_intersect[i])
            if(is.null(point_buf_intersect_i)){
                unit_slip_scale[i] = 0
            }else{
                unit_slip_scale[i] = gArea(point_buf_intersect_i)/gArea(point_bufs_intersect[i])
    
                #if(unit_slip_scale[i] > 1.001){
                #    # I have seen a case where order reversal seems to lead to
                #    # a bug in gIntersection, but could not reproduce it in an example.
                #    point_buf_intersect_i = gIntersection(point_bufs_intersect[i], p0)
                #    unit_slip_scale[i] = gArea(point_buf_intersect_i)/gArea(point_bufs_intersect[i])
                #}
            }
        }

        if(any(unit_slip_scale < 0)){
            stop('Error: Unit_slip_scale < 0')
        }
        
        # Remove points with no slip
        kk = which(unit_slip_scale > 0)
        if(length(kk) == 0) stop('Error in unit_slip_scale computation')
        p1 = p1[kk]
        unit_slip_scale = unit_slip_scale[kk]
    
        if(any(unit_slip_scale > 1)){
            if(any(unit_slip_scale > 1.0001)){
                print(max(unit_slip_scale))
                stop('Error: Unit_slip_scale > 1')
            }else{
                unit_slip_scale = pmin(unit_slip_scale, 1)
            }
        }
    }

    p1_points = SpatialPoints(coords = coordinates(p1))
    p_intersect = gIntersection(p1, p0, byid=TRUE, drop_lower_td = TRUE) 
    p_intersect_buf = p1 

    # If there are 'just touching' relations and similar, then
    # p_intersect may not be SpatialPolygons. Throw an error for now
    stopifnot(class(p_intersect)=='SpatialPolygons')

    areas = unlist(lapply(p_intersect@polygons, 
        f<-function(x) x@Polygons[[1]]@area))
    indices = unlist(lapply(p_intersect@polygons, f<-function(x) x@ID))
    indices = as.numeric(gsub('P', '', indices))

    centroids = coordinates(p_intersect)

    stopifnot(length(areas) == length(indices))

    areas_buffer = unlist(lapply(p_intersect_buf@polygons, 
        f<-function(x) x@Polygons[[1]]@area))

    # Find fraction of each 'buffer' point area that is inside the original
    # unit source region
    # gBuffer calls are used to try to work around topology exceptions from rgeos
    #unit_source_region = gUnaryUnion(gBuffer(gBuffer(p_intersect,
    #    width=1.0e-06), width=-1e-06))
    unit_source_region = p0
    area_buffer_fraction_inside_unit_source = areas_buffer * 0
    for(i in 1:length(p_intersect_buf)){
        ptmp = gIntersection(p0, p_intersect_buf[i])
        if(!is.null(ptmp)){
            area_buffer_fraction_inside_unit_source[i] = gArea(ptmp) /
                gArea(p_intersect_buf[i])
        }
    }

    return(list(
        grid_points = centroids, 
        index = indices, 
        area = areas, 
        dx = NA, dy = NA,
        grid_point_polygon=p_intersect,
        grid_point_polygon_buffer=p_intersect_buf,
        grid_points_buffer = coordinates(p_intersect_buf),
        area_buffer = areas_buffer,
        area_buffer_fraction_inside_unit_source = area_buffer_fraction_inside_unit_source,
        unit_slip_scale = unit_slip_scale
        ))
}

#
# Workhorse function used inside compute_grid_point_areas_in_polygon
#
# @param polygon 
# @param new_origin
# @param x_axis_vector
# @param approx_dx
# @param approx_dy
fill_polygon_with_grid_points<-function(polygon, new_origin, x_axis_vector,
        approx_dx, approx_dy){

    lp = length(polygon[,1])

    # Default rotation origin
    if(is.null(new_origin)){
        new_origin = polygon[lp,]
    }

    # Default rotation x-axis vector
    if(is.null(x_axis_vector)){
        x_axis_vector = polygon[lp,] - polygon[1,]
    }

    # Ensure that the first/last points are not the same
    stopifnot(!isTRUE(all.equal(x_axis_vector, c(0,0))))


    # Rotate polygon. This helps us ensure the grid points are aligned with the
    # top of the unit source. 
    rotated_polygon = rotate_cartesian2d(points = polygon, origin = new_origin, 
        x_axis_vector = x_axis_vector)
    
    f_edge<-function(edge){
        s_coord = c(0, 1) 
        f_x = approxfun(s_coord, edge[,1], rule=2)
        f_y = approxfun(s_coord, edge[,2], rule=2)
        outfun<-function(s_coord){
            out = cbind(f_x(s_coord), f_y(s_coord))
            return(out)
        }
        return(outfun)
    }

    # Functions to interpolate from [0-1] to xy on edges. Oriented
    # left-to-right or top-to-bottom
    f_top = f_edge(rotated_polygon[c(4,1),1:2])
    f_bot = f_edge(rotated_polygon[3:2,1:2])

    f_all<-function(alpha, s){
        out = (1-alpha) * f_top(s) + alpha*f_bot(s)
        return(out)
    }

    # Choose grid point spacing
    xv = rotated_polygon[4,1:2] - rotated_polygon[1,1:2]
    xn = rotated_polygon[4,1:2] - rotated_polygon[3,1:2]
    npx = ceiling(sqrt(sum(xv^2))/approx_dx)
    npy = ceiling(sqrt(sum(xn^2))/approx_dy)
    grid_normalised = expand.grid(
        seq(0+0.5/npy,1-0.5/npy,len=npy), 
        seq(0+0.5/npx,1-0.5/npx,len=npx))

    grid_points = f_all(grid_normalised[,1], grid_normalised[,2]) 
    
    grid_pol = list()
    for(i in 1:length(grid_points[,1])){
        gp = grid_points[i,]
        gp_as = c(grid_normalised[i,1], grid_normalised[i,2])

        local_polygon_rotated = rbind( 
            f_all(gp_as[1]-0.5/npy, gp_as[2]-0.5/npx),
            f_all(gp_as[1]-0.5/npy, gp_as[2]+0.5/npx),
            f_all(gp_as[1]+0.5/npy, gp_as[2]+0.5/npx),
            f_all(gp_as[1]+0.5/npy, gp_as[2]-0.5/npx))

        # Get the polygon in 'unrotated' coordiantes
        local_polygon = rotate_cartesian2d(local_polygon_rotated, 
            origin = new_origin, x_axis_vector = x_axis_vector, inverse=TRUE)

        grid_pol[[i]] = Polygons(list(Polygon(local_polygon)),
            ID=as.character(i))
    }

    p1 = SpatialPolygons(grid_pol, proj4string=CRS(""))

    return(p1)

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


#' Interpolate depth and dip at interior points in the unit source [deprecated]
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

    # Define a function that gives alpha/s associated with a given grid point
    find_alpha_s<-function(initial_guess, grid_point){

        # Function giving the (dx,dy) between the desired x,y point and the
        # line at alpha,s. We minimise the norm of this with the lm-algorithm
        # to find the best alpha/s
        minme<-function(alpha_s){
            (line_across_source(alpha_s[1], alpha_s[2])[1:2] - grid_point[1:2])
        }

        # Find (alpha, s) at grid_points[i,1:2]
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
    perturbation_scale = 100 # m
    grid_points_perturb = grid_points[,1:2] + 
        cbind(sin(strike_perp), cos(strike_perp))*perturbation_scale

    for(i in 1:length(grid_points[,1])){
        fit_alpha_s = find_alpha_s(alpha_s[i,], grid_points[i,1:2])
        alpha_s[i,] = fit_alpha_s$par

        alpha_s_perturbation[i,] = 
            find_alpha_s(fit_alpha_s$par, grid_points_perturb[i,1:2])$par
    }

    # Get the x,y,depth points and perturbed points
    output_points = line_across_source(alpha_s[,1], alpha_s[,2])
    output_points_perturb = line_across_source(alpha_s_perturbation[,1], 
        alpha_s_perturbation[,2])

    # Now compute the dip numerically [ tan(dip) = change_in_depth/distance ]
    dip = atan((output_points_perturb[,3] - output_points[,3])/perturbation_scale)/deg2rad

    output_points = cbind(output_points, dip, alpha_s)

    return(output_points)

}

#' Quick plot of a unit source interior points cartesian
#'
#' @param us List with unit source interior points information: Output of
#' 'unit_source_interior_points_cartesian' or similar. 
#' @return nothing but make a plot
#'
#' @export
plot_unit_source_interior_points_cartesian<-function(us){

    unit_source = us$unit_source_cartesian
    output_points = us$grid_points

    if(is.null(unit_source)){
        msg = paste0('"us" is not a CARTESIAN unit source. \n',
            'It should be the output of "unit_source_interior_points_cartesian" or similar')
        stop(msg)
    }


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

    deg2rad = pi/180
    perturb = ds * cbind(sin(us$grid_points[,'strike']*deg2rad - pi/2),
        cos(us$grid_points[,'strike']*deg2rad - pi/2))
    arrows(output_points[,1], output_points[,2], output_points[,1] + perturb[,1],
        output_points[,2] + perturb[,2], length=0)

    # Add angles which show the dip
    for(i in 1:length(output_points[,1])){
        
        l1 = cbind( c(output_points[i,1], output_points[i,1] + 2*ds/3),
                    c(output_points[i,2], output_points[i,2]))
        points(l1, t='l', col='orange', lwd=0.5)

        sindip = sin(output_points[i,'dip']*deg2rad)
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

    if(is.null(us$unit_source_cartesian)){
        msg = paste0('us is not a CARTESIAN unit source. \n',
            'It should be the output of "unit_source_interior_points_cartesian" or similar')
        stop(msg)
    }

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


#' Find the unit source containing a given xy coordinate
#'
#' @param point_xy numeric vector of length 2 giving c(lon, lat)
#' @param unit_source_geometry SpatialPolygonsDataFrame with unit source
#' geometry, e.g. output of \code{unit_source_grid_to_SpatialPolygonsDataFrame}
#' @param unit_source_statistics data.frame with unit source statistics, e.g.
#' output of \code{discretized_source_summary_statistics}
#' @return integer giving the row index in unit_source_statistics of the unit
#' source containing point_xy
#' @export
#' @examples
#'
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip)
#' 
#' # Get geometry
#' puysegur_geometry = unit_source_grid_to_SpatialPolygonsDataFrame(
#'     puysegur_discretized_source$unit_source_grid)
#' 
#' # Get summary statistics
#' puysegur_summary_statistics_approx = discretized_source_approximate_summary_statistics(
#'    puysegur_discretized_source)
#' 
#' # Find the unit source containing this lon/lat point
#' pp = c(166.5, -45.9)
#' 
#' # Plot the situation
#' plot(puysegur_geometry, axes=TRUE, asp=1)
#' points(pp[1], pp[2], col='red', pch=19)
#' # Find the index of the unit source containing pp [i.e. index in the summary statistics]
#' unit_source_index = find_unit_source_index_containing_point(
#'     pp, puysegur_geometry, puysegur_summary_statistics_approx)
#' # Here are the summary statistics
#' target_unit_source = puysegur_summary_statistics_approx[unit_source_index,]
#' 
#' # Plotting
#' points(target_unit_source$lon_c, target_unit_source$lat_c, col='green', pch=19, cex=0.3)
#' geo_ind = which(target_unit_source$downdip_number == puysegur_geometry$downdip_number &
#'     target_unit_source$alongstrike_number == puysegur_geometry$alongstrike_number)
#' plot(puysegur_geometry[geo_ind,], border='green', add=TRUE)
#' 
find_unit_source_index_containing_point<-function(point_xy, 
    unit_source_geometry, 
    unit_source_statistics){

    # Do the spatial operation
    xx = SpatialPoints(coords = matrix(point_xy, ncol=2), 
        proj4string=CRS(proj4string(unit_source_geometry)))
    dd_as = over(xx, unit_source_geometry)

    if(nrow(dd_as) != 1){
        stop('Count not find unique unit source containing point')
    }

    # Find the corresponding data.frame row 
    output_ind = which( 
        (unit_source_statistics$downdip_number == dd_as[1,1]) & 
        (unit_source_statistics$alongstrike_number == dd_as[1,2]))

    return(unit_source_statistics$subfault_number[output_ind])
}
