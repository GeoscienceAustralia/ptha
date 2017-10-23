# 
# Routines to extract data from the CMT catalogue within polygons that define
# our source-zones
#
library(rptha)
config = new.env()
source('config.R', local=config)

# Key inputs
cmt_catalogue_csv = config$cmt_catalogue_csv 
unit_source_grid_polygon_shapefiles = config$unit_source_grid_polygon_shapefiles 


# Properties of earthquake events we extract
mw_threshold = config$MW_MIN 
depth_threshold = config$depth_threshold #70 # depth < depth_threshold
buffer_width = config$buffer_width # Events are inside polygon, after polygon is buffered by 'buffer_width' degrees
allowed_rake_deviation = config$rake_deviation

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


#
# Read all unit-source grid shapefiles into a list, with names corresponding to
# source-zones
#
unit_source_grid_poly = lapply(
    as.list(unit_source_grid_polygon_shapefiles), 
    f<-function(x) readOGR(x, layer=gsub('.shp', '', basename(x))))

names(unit_source_grid_poly) = basename(dirname(dirname(dirname(
    unit_source_grid_polygon_shapefiles))))

#
# Read earthquake observational datasets
#
gcmt = read.csv(cmt_catalogue_csv)
days_in_year = 365.25
cmt_start_time_julian = julian(
    strptime('1976-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC'),
    origin=strptime('1900-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S',
        tz='Etc/UTC'))
cmt_start_time_julianDay1900 = as.numeric(cmt_start_time_julian, units='days')
cmt_end_time_julian = julian(
    strptime('2017-03-01 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC'),
    origin=strptime('1900-01-01 00:00:00', format='%Y-%m-%d %H:%M:%S',
        tz='Etc/UTC'))
cmt_end_time_julianDay1900 = as.numeric(cmt_end_time_julian, units='days')
cmt_duration_years = (cmt_end_time_julianDay1900 - cmt_start_time_julianDay1900)/days_in_year 
# Check the times are accurate by comparing with time between first/last
# earthquake
if( (cmt_duration_years < diff(range(gcmt$julianDay1900))/days_in_year) |
    (cmt_duration_years > 1.0005* diff(range(gcmt$julianDay1900))/days_in_year) ){
    stop('GCMT date/time inconsistencies')
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
#' @param buffer width thickness of buffer to apply to poly before testing the point inclusion.
#' @return logical vector with one entry for every row in lonlat.
#'
lonlat_in_poly<-function(lonlat, poly, buffer_width = 0){

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

#'
#' Extract earthquake events > threshold for a source-zone or segment
#'
#' Spatial inclusion is judged by testing both hypocentres and centroids [if
#' either is judged as 'inside', then the event is treated as inside.
#'
#' @param source_name name of a source-zone in unit_source_grid_poly
#' @param alongstrike_index_min NULL, or an integer alongstrike index where the
#'   segment of interest begins
#' @param alongstrike_index_min NULL, or an integer alongstrike index where the
#'   segment of interest ends
#' @param local_mw_threshold keep earthquakes with magnitude >= this
#' @param local_depth_threshold keep earthquakes with depth <= this
#' @param target_rake_value Keep events within 'allowed_rake_deviation' of 'target_rake_value'
#' @param allowed_rake_deviation See 'target_rake_value' above
#' @param local_buffer_width buffer polygon by this many degrees before running
#'   point-in-polygon test.
#' @param use_both_gcmt_solutions. Use either rake1 or rake2 when testing for
#'   event inclusion.
#' @param filter_by_strike Logical. If TRUE, also consider the strike when
#'   accepting/rejecting gcmt data
#' @param allowed_strike_deviation. Accepted events must have strike within
#'   this many degrees of the nearest unit source strike. Only matters if filter_by_strike==TRUE
#' @param unit_source_table Unit source table for this source zone (only required if filter_by_strike=TRUE)
#' @return subset of gCMT catalogue
#'
get_gcmt_events_in_poly<-function(source_name, 
    alongstrike_index_min = NULL, 
    alongstrike_index_max = NULL,
    local_mw_threshold = mw_threshold, 
    local_depth_threshold = depth_threshold, 
    target_rake_value = 90,
    allowed_rake_deviation = config$rake_deviation, 
    local_buffer_width=buffer_width,
    use_both_gcmt_solutions=FALSE,
    filter_by_strike = FALSE,
    allowed_strike_deviation = NA,
    unit_source_table = NULL){

    target_rake_value = target_rake_value

    if(!(source_name %in% names(unit_source_grid_poly))){
        stop(paste0('No matching source_name for ', source_name))
    }

    # Get the polygon
    poly = unit_source_grid_poly[[source_name]]

    # If only looking at a subset of the polygon, we will use this
    # variable to keep track of other regions [to avoid point double-counting]
    poly_neighbour_segments = NULL

    # Get a subset of 'poly' with the right alongstrike indices
    if(!is.null(alongstrike_index_min) | !is.null(alongstrike_index_max)){

        # Check input args
        if(is.null(alongstrike_index_min) | is.null(alongstrike_index_max)){
            stop('Must provide BOTH alongstrike_index_min and alongstrike_index_max, or neither')
        }
        # Ensure name-mangling in shapefile has not changed
        if( !('alngst_' %in% names(poly)) ){
            stop('Polygon does not have"alngst_" attribute')
        }

        # Get unit-source table
        if(!is.null(unit_source_table)){
            keep = which(unit_source_table$alongstrike_number >= alongstrike_index_min &
                unit_source_table$alongstrike_number <= alongstrike_index_max)
            unit_source_table = unit_source_table[keep,]
        }

        poly_alongstrike_index = as.numeric(as.character(poly[['alngst_']]))

        keep = which((poly_alongstrike_index >= alongstrike_index_min) & 
            (poly_alongstrike_index <= alongstrike_index_max))

        if(length(keep) == 0) stop('No indices to keep')

        # Get the subset of interest 
        poly_orig = poly
        poly = poly_orig[keep,]

        # Keep the 'remainder' of the polygon, to check whether events
        # would be included in multiple segments
        if(length(poly_alongstrike_index) > length(keep)){
            poly_neighbour_segments = poly_orig[-keep,]
        }

    }

    # Count as inside if EITHER the hypocentre, or the centroid, is inside 
    inside_events_hypo = lonlat_in_poly(
        gcmt[,c('hypo_lon', 'hypo_lat')], poly, buffer_width=local_buffer_width)
    inside_events_centroid = lonlat_in_poly(
        gcmt[,c('cent_lon', 'cent_lat')], poly, buffer_width=local_buffer_width)
   
    inside_events = (inside_events_hypo | inside_events_centroid)

    # Other criteria we have (might not be used, see below)
    mw_depth_ok = (gcmt$Mw >= local_mw_threshold) & (gcmt$depth <= local_depth_threshold)
    # Note we need to be careful about circularity of angles when checking for rakes
    rake1_ok = angle_within_dtheta_of_target(gcmt$rake1, target_rake_value, allowed_rake_deviation)
    rake2_ok = angle_within_dtheta_of_target(gcmt$rake2, target_rake_value, allowed_rake_deviation)

    if(filter_by_strike){
        # Determine whether strike1 and strike2 are within some threshold of
        # the nearest unit-source strike value
        
        if(is.null(unit_source_table)){
            stop('Must provide unit_source_table if filter_by_strike==TRUE')
        }

        if(is.na(allowed_strike_deviation)){
            stop('Must provide allowed_strike_deviation if filter_by_strike==TRUE')
        }

        #
        # Compute the strike of the nearest unit source for each point that
        # might be inside
        #       
        # Variable to hold the strike value of the nearest unit source
        nearest_strike = rep(NA, length(inside_events))
        # For every earthquake in the 'inside-events' group, find the strike of
        # the nearest unit-source
        for(i in 1:length(inside_events)){
            if(inside_events[i] & mw_depth_ok[i]){
                # Unit source coordinates in this segment
                p2 = cbind(unit_source_table$lon_c, unit_source_table$lat_c)
                # Earthquake coordinates (centroid based)
                p1 = p2*0
                p1[,1] = gcmt$cent_lon[i]
                p1[,2] = gcmt$cent_lat[i]
                # Nearest unit source strike
                nearest_uss = which.min(distHaversine(p1, p2))
                nearest_strike[i] = unit_source_table$strike[nearest_uss]%%360
            }
        }
        
        # Is the strike close to the desired value?
        strike1_ok = angle_within_dtheta_of_target(gcmt$strk1, nearest_strike, 
            allowed_strike_deviation)
        strike1_ok[which(is.na(strike1_ok))] = FALSE
        strike2_ok = angle_within_dtheta_of_target(gcmt$strk2, nearest_strike, 
            allowed_strike_deviation)
        strike2_ok[which(is.na(strike2_ok))] = FALSE

    }else{
        # Do not limit the selected events based on strike
        strike1_ok = rep(TRUE, length(inside_events))
        strike2_ok = rep(TRUE, length(inside_events))
    }

    # If we don't use both gcmt solutions, then ensure 'strike2/rake2' the same
    # as strike1/rake1
    if(!use_both_gcmt_solutions){
        strike2_ok = strike1_ok
        rake2_ok = rake1_ok
    }
    
    # Criterion for point selection - note we will keep the event if either
    # rake1 or rake2 is within a given range of pure thrust
    inside_keep = which( inside_events & mw_depth_ok &
            ((rake1_ok & strike1_ok) | (rake2_ok & strike2_ok) )
        )


    if(length(inside_keep) > 0){

        output_gcmt = gcmt[inside_keep,]

        # Keep track of double-counted points
        output_gcmt$double_counted = rep(FALSE, length=length(inside_keep))
        if(!is.null(poly_neighbour_segments)){
            # Identify points that are also inside neighbour segments. 
            inside_events_hypo = lonlat_in_poly(
                output_gcmt[,c('hypo_lon', 'hypo_lat')], poly_neighbour_segments, 
                buffer_width=local_buffer_width)
            inside_events_centroid = lonlat_in_poly(
                output_gcmt[,c('cent_lon', 'cent_lat')], poly_neighbour_segments, 
                buffer_width=local_buffer_width)
            inside_events = (inside_events_hypo | inside_events_centroid)
            output_gcmt$double_counted = inside_events
            
        }
    }else{
        # If there is no data, return an empty data.frame, still with the
        # correct names
        output_gcmt = as.data.frame(matrix(0, nrow=0, ncol=ncol(gcmt)+1))
        names(output_gcmt) = c(names(gcmt), 'double_counted')
    }

    return(output_gcmt)
}


source_zone_events_plot<-function(source_name, eq_events, focsize_t1 = 6.5, focsize_t2 = 6.0){

    library(RFOC)

    plot(unit_source_grid_poly[[source_name]], main=source_name, border='grey', axes=TRUE)
    points(eq_events$cent_lon, eq_events$cent_lat, cex=(eq_events$Mw-focsize_t1)**2, col='red')
    # Make sure we don't miss points due to longitude convention
    points(eq_events$cent_lon-360, eq_events$cent_lat, cex=(eq_events$Mw-focsize_t1)**2, col='red')
    points(eq_events$cent_lon+360, eq_events$cent_lat, cex=(eq_events$Mw-focsize_t1)**2, col='red')

    #
    # Same plot as above, with beachballs
    #

    plot(unit_source_grid_poly[[source_name]], main=source_name, border='grey', axes=TRUE)

    pol_centroid = as.numeric(coordinates(gCentroid(unit_source_grid_poly[[source_name]])))

    for(j in 1:length(eq_events$cent_lon)){

        lon_lat = c(eq_events$cent_lon[j], eq_events$cent_lat[j])
        lon_lat = adjust_longitude_by_360_deg(lon_lat, pol_centroid)
    
        # Beach-ball details
        mec = SDRfoc(eq_events$strk1[j], eq_events$dip1[j], eq_events$rake1[j], PLOT=FALSE)

        fcol =  c( rgb(1,0,0, alpha=1), 'lightblue', 'green', 'orange', 'yellow', 'purple', 'black')[foc.icolor(mec$rake1)]

        par(lwd=0.2) # Try to avoid strong lines on beachballs
        try(
            justfocXY(mec, lon_lat[1], lon_lat[2],
                focsiz=4*((eq_events$Mw[j]-focsize_t2)/10)**2, xpd=TRUE, fcol=fcol,
                fcolback='white'),
            silent=TRUE
            )
        par(lwd=1)

    }

}
