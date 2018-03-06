##
## Code to support manipulation of the source zone fault traces and slab1.0
## data
##

library(rgdal)
library(sp)
library(geosphere)
library(raster)

#' Read all shapefiles in 'source_zones_geom', clean up for later processing,
#' and combine into a single shapefile
#'
#' @return SpatialLinesDataFrame with all traces of source-zones
#'
collate_source_traces<-function(){
    sourcefiles = Sys.glob('./source_zones_geom/*.shp')

    for(i in 1:length(sourcefiles)){
        # Read in the shapefile
        ith_layer = strsplit(basename(sourcefiles[i]), '\\.')[[1]][1]
        ithgeo = readOGR(dsn=sourcefiles[i], layer=ith_layer)

        # We need unique ids to do the merge
        # Get the present id
        ithgeoID = unlist(lapply(ithgeo@lines, myfun<-function(x) x@ID))
        # Add the id to the data frame
        ithgeo@data = cbind(ithgeo@data, data.frame(ID1=ithgeoID))
        # Make a new id
        newids = paste0(ith_layer, '-', ithgeoID)
        ithgeo@data = cbind(ithgeo@data,data.frame(newID=newids))
        ithgeo = spChFIDs(ithgeo,newids)

        # Some faults have different attrubute columns. Make sure we get only the
        # ones we need
        needed_data =data.frame(Segnum=ithgeo@data$Segnum, Dip=ithgeo@data$Dip, 
            Rake=ithgeo@data$Rake, ID1=ithgeo@data$ID1, newID=ithgeo@data$newID)
        rownames(needed_data) = newids
        ithgeo@data = needed_data
        

        # Collate
        if(i == 1){
            all_traces = ithgeo
        }else{
            all_traces = rbind(all_traces,ithgeo)
        }
        
    }

    writeOGR(all_traces, dsn='./OUTPUT_DATA/allTraces', 
             layer='allTraces', driver="ESRI Shapefile",
             overwrite=TRUE)
    return(all_traces)
}

###############################################################################

#
# Given a fault trace known to be dipping to the right of the line
# orientation, with dip provided in the attribute table, and a known
# maximum depth, extend 'normals' from the fault trace down the dip
# to the max depth. We can then use this to make a polygon / interpolate
# from the trace
#
# INPUTS: fault_trace = SpatialLines of the fault trace, with appropriate
#                       attribute table. 
#         maxdepth/mindepth == max/min depths of the slab, in KM
#
# OUTPUT: List with 1) Polygon of slab region; 2) Raster of slab depths; 3)
#         Contour from slab depths
# @param fault_trace SpatialLinesDataFrame with the fault trace, and an
# appropriate attribute table
# @param maxdepth The maximum depth to which the contours should go (constant)
# @param mindepth The minimum depth to which the contours should go (constant)
# @param maxdepth_offset Do intermediate computations with 
#   maxdepth = maxdepth+maxdepth_offset. This is useful e.g. to ensure
#   smoothness of deepest contours (constant)
# @param in_0_360 logical. If TRUE, assume longitudes are in [0-360] degrees. 
#   Otherwise assume longitudes are in [-180, 180] degrees.
# @param downdip_profile character with 'linear' or 'parabolic', defining the
#   type of down-dip profile to use
# @param dip_depth (useful if downdip_profile = 'parabolic') depth in km at
#   which to make trace dip = dip, when fitting the parabola. This is not 
#   required for downdip_profile = 'linear', because the dip is constant. If
#   provided, it should be a vector with length equal to the length of fault_trace.
# @param mindepth_dip assumed dip of the interface at the 'trench' or other 
#   min depth point. Only required for downdip_profile == 'parabolic'. If
#   provided, should have length equal to the length of fault_trace
#
extend_trace_to_depth_contours<-function(fault_trace, maxdepth, mindepth, 
    maxdepth_offset=0, in_0_360=TRUE, downdip_profile='linear', 
    dip_depth=NULL, mindepth_dip=NULL){

    #
    # Check inputs
    #
    lft = length(fault_trace)
    stopifnot(length(maxdepth) == 1)
    stopifnot(length(mindepth) == 1)
    stopifnot(length(maxdepth_offset) == 1)
    if(downdip_profile != 'linear'){
        stopifnot(length(dip_depth) == lft)
        stopifnot(length(mindepth_dip) == lft)
    }
    

    # To make the contours well behaved we need to extend the max depth. 
    maxdepth_extended = maxdepth + maxdepth_offset

    if(maxdepth_extended > 250){
        stop('maxdepth_extended > 250, but I have assumed maxdepth is in km. Suggests an input error')
    }
    if(maxdepth_extended<mindepth){
        stop('maxdepth_extended < mindepth ! Input error!')
    }

    # Get points along the trace, and remove repeated points (caused by
    # snapping line segments of the input shapefile)
    fault_trace_points = as(fault_trace,'SpatialPointsDataFrame')
    keep_points = c( seq(1, length(fault_trace_points), 2), 
        length(fault_trace_points))
    fault_trace_points = coordinates(fault_trace_points)[keep_points,]
    
    # Compute bearing of down dip direction, and find point at the maximal
    # depth along this bearing
    store_seg_destpoints = c()
    for(i in 1:length(fault_trace)){

        segment_pts = as(fault_trace[i,], 'SpatialPointsDataFrame')

        if(length(segment_pts) != 2){
            stop('Fault segment is not defined by a single line. Problem with input data')
        }

        fault_bearing = bearingRhumb(segment_pts[1,], segment_pts[2,])

        # Make normal to trace, varies from [0-360]
        fault_downdip_bearing = (fault_bearing+90)%%360

        # Compute the along-earth-surface distance that we extend the geometric
        # subduction interface in the down-dip direction. This assumes a linear
        # dip, so the distance will be long enough so long as the assumed depth
        # profile is linear or super-linear. In sub-linear cases (concave-up),
        # the distance will be enough if maxdepth_offset is sufficiently large.
        # FIXME: Consider a more robust approach
        max_distance = (maxdepth_extended - mindepth) * 1000 /
            tan(segment_pts[1,]@data$Dip/180*pi) # units meters

        # Find the 'end' points of a line with direction
        # 'fault-downdip-bearing' and length 'max_distance' from each segment
        # point
        seg_destpoint1 = destPoint(segment_pts[1,], fault_downdip_bearing,
            d=max_distance)
        if(in_0_360 & (seg_destpoint1[1]<0)){
            seg_destpoint1[1] = 360+seg_destpoint1[1]
        }
        seg_destpoint2 = destPoint(segment_pts[2,], fault_downdip_bearing,
            d=max_distance)
        if(in_0_360 & (seg_destpoint2[1]<0)){
            seg_destpoint2[1] = 360+seg_destpoint2[1]
        }

        store_seg_destpoints = rbind(store_seg_destpoints, seg_destpoint1)
        store_seg_destpoints = rbind(store_seg_destpoints, seg_destpoint2)

    }

    # We now have 2 seg_destpoints for each fault trace point, except at the
    # end points. Let's use the arithmetic average to remove the repeated
    # points, and thus determine the fault boundary
    ll = length(fault_trace)+1
    thinned_seg_destpoints = matrix(NA,nrow=ll,ncol=2)
    for(i in 1:(ll)){
        if(i == 1){
            # Special case -- first point 
            thinned_seg_destpoints[i,] = store_seg_destpoints[1,]
        }else if(i == ll){
            # Special case -- last point
            # Get the last seg-destpoint. There are 2*(length(fault_trace))
            # in total
            thinned_seg_destpoints[i,] = store_seg_destpoints[2*i-2,]
        }else{
            thinned_seg_destpoints[i,] = midPoint(store_seg_destpoints[2*i-2,],
                                                  store_seg_destpoints[2*i-1,])
            if(in_0_360 & (thinned_seg_destpoints[i,1]<0)){
                thinned_seg_destpoints[i,1] = 360+thinned_seg_destpoints[i,1]
            }
        }
    }
    
    # Make a polygon surrounding the fault in areas where we have depths
    fault_seg_poly = rbind(fault_trace_points,
                         thinned_seg_destpoints[ll:1,])
    fault_pt_depths = c(rep(mindepth,ll), rep(maxdepth_extended,ll))
    fault_seg_poly_closed = rbind(fault_seg_poly, fault_seg_poly[1,]) # Close the ring
    fault_seg_poly_sp =
        SpatialPolygons(list(Polygons(
                             list(Polygon(fault_seg_poly_closed, hole=FALSE))
                                , ID='0')),
                        proj4string=CRS('+init=epsg:4326'))

    # Now, we want to interpolate these depths onto a raster, then contour
    # 
    # Inefficient good-enough approach: Let's generate loads of long,lat,depth
    # points, then just do nearest neighbour or something
    plot(fault_seg_poly_sp, col='green', axes=T, asp=1)
    points(fault_seg_poly, col=as.factor(fault_pt_depths))
    store_pts = c()
    box_horiz_n = 100
    box_down_n = 100
    lft = length(fault_trace)
    for(i in 1:lft){

        segment_pts = as(fault_trace[i,], 'SpatialPointsDataFrame')

        ddepth_dx_reference = tan(segment_pts[1,]@data$Dip/180*pi)

        box_start = fault_seg_poly[i:(i+1),]
        box_end = fault_seg_poly[(2*ll-(i-1)):(2*ll-1-(i-1)),]

        # Plot it
        mybox=rbind(box_start,box_end[2:1,])
        points(mybox,t='l')
        
        for(j in 1:box_horiz_n){
            bp1 = box_start[1,] + (box_start[2,]- box_start[1,])*(j-1)/(box_horiz_n-1)
            bp2 = box_end[1,] + (box_end[2,]- box_end[1,])*(j-1)/(box_horiz_n-1)
            local_distance = distHaversine(bp1, bp2)/1000

            if(downdip_profile == 'linear'){

                local_pts = cbind(
                    gcIntermediate(bp1, bp2, n=box_down_n-2, addStartEnd=TRUE),
                    seq(mindepth, maxdepth_extended, len=box_down_n))

            }else if(downdip_profile == 'parabolic'){

                if(is.null(dip_depth[i])){
                    stop('Must provide dip_depth for parabolic interpolation')
                }

                ## Parabolic variation:
                ##
                # depth = mindepth + beta * distance + alpha * (distance^2)
                ##
                ## where: 
                ##   distance := 0; at depth=mindepth
                ##   slope = mindepth_interface_slope; at depth=mindepth
                ##   slope = ddepth_dx_reference; at depth=dip_depth
                ##
                ## Denote depth1, d1, s1 as the depth, distance, and slope at the mindepth,
                ## and    depth2, d2, s2 as the depth, distance, and slope at the reference_depth.
                ## We know that:
                ##   depth1 = mindepth
                ##   d1 = 0
                ##   s1 = mindepth_interface_slope
                ##   depth2 = dip_depth
                ##   s2 = ddepth_dx_reference
                ## We need to compute 'beta, alpha' in the parabolic equation above.
                ## Differentiating the parabola, we know:
                ##     ddepth/dx = beta + 2*alpha*distance
                ## and so from the mindepth location, we know
                ##     s1 = beta   
                ## Further, using the gradient at dip_depth, we know:
                ##     s2 = beta + 2 * alpha * d2
                ## or
                ##     d2 = (s2 - s1)/(2*alpha)
                ## Also, at dip_depth:
                ##     depth2 = mindepth + beta * d2 + alpha * d2^2
                ##            = mindepth + s1 * (s2-s1)/(2*alpha) + alpha * (s2-s1)^2 / (4 * alpha^2)
                ## Multiplying by alpha and simplifying, we get:
                ##     0 = (mindepth-depth2)*alpha + s1*(s2-s1)/2 + (s2-s1)^2/4
                ## So:
                ##    alpha = (s1*(s2-s1)/2 + (s2-s1)^2 / 4) / (depth2 - mindepth)
                ## We now know 'alpha' and 'beta', so know the shape of the parabola.
                ##
                mindepth_interface_slope = tan(mindepth_dip[i] * pi/180)
                s1 = mindepth_interface_slope
                s2 = ddepth_dx_reference
                beta = s1
                alpha = ( s1 * (s2 - s1)/2 + 1/4 * (s2-s1)^2 ) / (dip_depth[i] - mindepth)
                dists = seq(0, local_distance, length=box_down_n)
                depth_vals = mindepth + beta*dists + alpha * dists*dists

                if(max(depth_vals) < maxdepth ){
                    xv = -beta/(2*alpha)
                    yv = mindepth + beta * xv + alpha * xv * xv
                    # Check if there is a turning point
                    if(alpha < 0 & yv < maxdepth){
                        msg1 = paste0('Based on your inputs, the parabolic profile has a turning ',
                            ' \n point near to depth=', yv, ' which is < max_depth=', maxdepth, 
                            ' \n This suggests a data input error, and is not allowed. You should either',
                            ' \n reduce the dip at the trench, or increase the deeper dip value, or',
                            ' \n decrease the maxdepth. The parameters were:',
                            ' \n mindepth = ', mindepth, 
                            ' \n dip@mindepth = ', mindepth_dip[i], 
                            ' \n dip_depth = ', dip_depth[i], 
                            ' \n dip@dip_depth = ', segment_pts[1,]@data$Dip,
                            ' \n Leading to a parabola like " depth = mindepth + beta * distance + alpha * distance^2 "',
                            ' \n with parameters alpha=', alpha, ' beta=', beta,
                            ' \n where distance=0 at depth=mindepth' 
                            )
                        stop(msg1)
                    }else{
                        msg2 = paste0('Based on your inputs, the depth never reaches maxdepth=', maxdepth, 
                            '\n This might mean that the maxdepth_offset is too small, try increasing it.',
                            '\n If that doesn"t work, debugging is required')
                        stop(msg2)
                    }
                }

                # For concave-down profiles, the depth can decrease after the initial increase. Do not allow that
                if(alpha < 0 & min(diff(depth_vals)) < 0){
                    xv = -beta/(2*alpha)
                    yv = mindepth + beta * xv + alpha * xv * xv
                    k = which(dists > xv & depth_vals < yv)
                    # Set the depth vals after the turning point to a 'deep' value, that will be cut off
                    depth_vals[k] = yv + 1000
                }


                #depth_vals = parabola_const * seq(0, local_distance, length=box_down_n)**2 + mindepth

                # Higher values above maxdepth (which are typically later truncated)
                local_pts = cbind(
                    gcIntermediate(bp1, bp2, n=box_down_n-2, addStartEnd=TRUE),
                    depth_vals)
            }else{

                stop(
                    paste0('downdip_profile value "', downdip_profile, 
                        '" not recognized'))

            }

            if(in_0_360 & any(local_pts[,1]<0)){
                local_pts[,1] = local_pts[,1]*(local_pts[,1]>=0) + 
                    (360+local_pts[,1])*(local_pts[,1]<0)
            }
            store_pts = rbind(store_pts,local_pts)
        }
    }

    # Push onto raster, in a crude way
    myrast = raster(extent(fault_seg_poly_sp)) #,nrow=50,ncol=50)
    res(myrast) = c(1,1)/10 # 1/10 of a degree
    proj4string(myrast) = CRS('+init=epsg:4326')
    myrast = setValues(myrast,NA)

    store_pts_coords = SpatialPoints(store_pts[,1:2], 
        proj4string=CRS('+init=epsg:4326'))
    store_pts_indices = extract(myrast,store_pts_coords,cellnumbers=T)[,1]

    rast_vals = aggregate(store_pts[,3],list(store_pts_indices),FUN=mean,na.rm=T)

    myrast[rast_vals[,1]] = rast_vals[,2]
    myrast = extend(myrast, 10) # Add a buffer so that contouring works ok

    ## Now, let's ensure the trace is burned in with the mindepth
    trace_pt_coords = extract(myrast, 
        spsample(fault_trace,n=3000,type='regular'), 
        small=T, cellnumbers=T)
    trace_pt_coords = trace_pt_coords[,1]
    myrast[trace_pt_coords] = mindepth
    
    plot(myrast, col=rev(terrain.colors(255)), alpha=0.85,add=T)

    # Make contour with 5km spacing
    mycontours = rasterToContour(myrast,
        levels=c(mindepth+0.001, 
            seq(mindepth+5, maxdepth_extended-5,by=5), 
            maxdepth_extended-0.001))

    output = list(fault_seg_poly_sp, myrast, mycontours)

    return(output)
}

############################################################################

# For all source zones, get the min/max depth from the source-trace-metata
# and gem subduction data as appropriate.
# Method of selection based on discussions with Jonathan Griffin, October 2013
# Put in a function to keep the namespace clean
get_min_max_depth<-function(source_traces, source_trace_metadata, 
    gem_subduction_data){

    layer_min_max_store = data.frame(
        layername = rep(NA,length(source_traces)),
        mindepth = rep(NA,length(source_traces)),
        maxdepth = rep(NA,length(source_traces)),
        stringsAsFactors=FALSE)

    for(i in 1:length(source_traces)){

        mylayer = strsplit(basename(source_traces[i]), ".shp")[[1]][1]

        mindepth=0

        # Find reference to the trace in source_trace_metadata
        meta_inds=grep(mylayer, source_trace_metadata[,1])

        # Find reference to the trace in the gem data.
      
        if(length(meta_inds)>=1){
            # We may have had a partial match -- make sure it is an exact match
            meta_inds_refine = which(source_trace_metadata[meta_inds,1] %in% mylayer)
            meta_inds = meta_inds[meta_inds_refine]
        
            if(length(meta_inds)==0) stop('No match with mylayer')
    
            if(!(any(source_trace_metadata[meta_inds,2]=="NULL"))){
                # If the second column of source_trace_metadata is not null
                # then look up the GEM table for the max-depth data
                GEM_match=match(source_trace_metadata[meta_inds,2], 
                    gem_subduction_data[,1])
                # If the layer name matches multiple rows in
                # source-trace-metadata, then take the maximum depth 
                maxdepth=max(gem_subduction_data[GEM_match,]$Down.dip.depth.max)
            }else{
                # If the second column is 'NULL', then the 6th column should
                # contain the max depth
                if( length(meta_inds) > 1){
                    stop('meta_inds must be unique for NULL gem match')
                }
                maxdepth = as.numeric(as.character(
                    source_trace_metadata[meta_inds,7]))
            }
        }else{
            stop(paste0('No match for ', mylayer))
        }

        # NOTE: gem_subduction_data$Updip_min is giving the trench depth below sea
        #       level. So it is zero anyway in our notation (and identical to the
        #       trench-depth column)
        #
        layer_min_max_store[i,] = data.frame(
            layername = mylayer,
            mindepth = mindepth,
            maxdepth = maxdepth,
            stringsAsFactors=FALSE)
    }

    write.table(layer_min_max_store, file='layer_mindepth_maxdepth.csv',
        sep=",", row.names=F)

    return(layer_min_max_store)
}


##############################################################################################

# Parse the SLAB 1.0 data, and generate polygons around its extent
#
get_slab<-function(slab_depth_files){

    s_outline_poly_store = list()
    SLAB_rasters = list()
    wgs84= CRS("+init=epsg:4326")

    for(i in 1:length(slab_depth_files)){

        print(i)

        # Read in rasters
        s_Depth = raster(slab_depth_files[i])

        # Make sure the projection info is there
        proj4string(s_Depth) = wgs84

        # Make polygon of raster non-na extent
        # First buffer the raster, so that the contour lines don't break at edges
        s_Depth = extend(s_Depth, 10) # This seems enough
        s_Depth_nona = !is.na(s_Depth) # 0 with nodata or 1 where there is data
        s_outline = rasterToContour(s_Depth_nona, level=0.5) # Outline data

        # There can be multiple loops in the contour -- deal with these by
        # finding the longest one
        s_outline_coords = coordinates(s_outline)
        poly_index = which.max(unlist(lapply(s_outline_coords[[1]], length)))

        # Add last point to coordinates to 'close the polygon'
        s_outline_coords = s_outline_coords[[1]][[poly_index]]
        s_outline_coords = rbind(s_outline_coords,s_outline_coords[1,])
        # Finally convert to polygon
        s_outline_poly = SpatialPolygons(
            list(Polygons(list(Polygon(s_outline_coords, hole=FALSE)),
                ID=slab_depth_files[i])), 
            proj4string=CRS(proj4string(s_Depth)))

        # Save key info into lists
        s_outline_poly_store[[i]] = s_outline_poly
        SLAB_rasters[[i]] = s_Depth
    }

    # Stick all the outline polygons into a single shapefile
    outline_polys = s_outline_poly_store[[1]]
    for(i in 2:length(slab_depth_files)){
        outline_polys = rbind(outline_polys, s_outline_poly_store[[i]])
    }

    SLAB_polys = SpatialPolygonsDataFrame(
        outline_polys,
        data=data.frame(ID=slab_depth_files),
        match.ID=FALSE)

    return(list(SLAB_polys, SLAB_rasters))

}

#################################################################

# Reverse the order of a SpatialLinesDataFrame
#
# Reverses the line segment order as well as line segments themselves
# @param kuril  A SpatialLinesDataFrame to be reversed
# @return the reversed SLDF
reverse_SLDF_line_order<-function(kuril){

    # Flip coordinate order
    kuril2 = kuril[length(kuril):1,]
    # Flip internal line order
    kuril2_lines = list()

    for(i in 1:length(kuril2)){

        tmp = kuril2[i,]@lines[[1]]@Lines[[1]]@coords

        ll = length(tmp[,1])

        kuril2_lines[[i]] = Lines(list(Line(tmp[ll:1,])), 
            ID=kuril2[i,]@lines[[1]]@ID)

    }
    kuril3 = SpatialLines(kuril2_lines, proj4string=CRS(proj4string(kuril)))
    kuril3SLDF = SpatialLinesDataFrame(kuril3, data=kuril2@data)
}

##################################################################

# Make a raster suitable for 'smooth' contouing
#
# Given a raster newrast, fit a smoother to the xyz points. Then, makes a new
# raster with a resolution finer than the input raster by a factor of
# resShrink, and predicts the smoother at these points
#
# @param newrast: A raster
# @param resShrink: Ratio of the input raster resolution to the new raster
# resolution
# @param minPred lower bound on the predicted value
# @param smoother string, either 'loess' or 'tps'
# @param myspan span for loess smooth
# @return A raster with resolution res(newrast)/resShrink
#
smooth_raster_for_contouring<-function(newrast, resShrink=1, minPred=-Inf, 
    smoother='loess', myspan=0.05,...){

    newrast_xyz = as.data.frame(rasterToPoints(newrast))

    if(smoother == 'loess'){

        # Apply a 'loess' local regression to the raster
        newrast_smoother = loess(layer~x*y, data=newrast_xyz, span=myspan, order=1, 
            family='symmetric', control=loess.control(iterations=20))

    }else if(smoother=='tps'){

        # Apply a gam with thin-plate smoothing, degrees of freedom
        # automatically guessed
        # Experimented with thin plate (s) + tensor product (te) + modified
        # tensor product (t2) -- the first looked to make more continuous
        # contours and deform the original raster less, with default settings

        library(mgcv)
        newrast_smoother = gam(layer ~ s(x,y, bs='tp',k=20), data=newrast_xyz)
        #newrast_smoother=gam(layer~t2(x,y),data=newrast_xyz)

    }else{
        stop('Smoother not recognized in smooth_raster_for_contouring')
    }

    # Make a (possibly finer) raster to predict on
    newrast_fine = raster(newrast) 
    res(newrast_fine) = res(newrast)/resShrink
    newrast_fine = resample(newrast,newrast_fine)

    # Get points from the fine raster, and predict at them with the smoother
    newrast_xyz2 = as.data.frame(rasterToPoints(newrast_fine))
    newrast_xyz2[,3] = pmax(minPred, 
        predict(newrast_smoother, newdata=newrast_xyz2[,1:2]))

    # Get our smoother raster -- working around an old bug in 'raster'
    # (probably fixed now)
    library(rptha)
    newrast_smoothed = rptha:::.local_rasterFromXYZ(cbind(newrast_xyz2[,1:2], 
        pmax(newrast_xyz2[,3], minPred)), digits=3)

    proj4string(newrast_smoothed) = proj4string(newrast)

    return(newrast_smoothed)

}

#####################################################################

# Code to help repair broken contour lines
#
# Sometimes we have 'broken' contourlines, because of data-merging artefacts
# Our methods can only have one contour. Experience suggests we just want the
# first one. Do that here.
#
# @param mycontours SpatialLinesDataFrame made from contouring a raster
# @return new version of mycontours, with repeats removed
strip_small_repeat_contours<-function(mycontours){

    myLines = mycontours@lines

    for(i in 1:length(myLines)){
        # Find 'distances^2' of start/end point on each contour line
        all_lengths = unlist(lapply(myLines[[i]]@Lines, 
            myfun<-function(x){
                sum((x@coords[1,]-x@coords[length(x@coords[,1]),])**2)
                }))
        j = which.max(all_lengths)
        myLines[[i]]@Lines = list(myLines[[i]]@Lines[[j]])
    }
   
    newSL = SpatialLines(myLines, proj4string=CRS(proj4string(mycontours))) 

    newSLDF = SpatialLinesDataFrame(newSL, data=mycontours@data)

}


######################################################################
#' 
#' Append the source trace to contours defining the subduction interface
#' The trace is assigned a user-provided depth
#'
#' @param contours SpatialLinesDataFrame with source trace contours
#' @param trace SpatialLinesDataFrame with source trace.
#' @param trace_depth the depth of the trace
#' @return Contours which include the trace with the required depth
append_trace_to_contours<-function(contours, trace, trace_depth){

    # Make the trace a single line    
    trace_line_data = lapply(coordinates(trace), f<-function(x) x[[1]])

    trace_line_matrix = matrix(trace_line_data[[1]][1,], ncol=2)

    lt = length(trace_line_data)
    if(lt > 1){
        for(i in 2:lt){
            # Check lines are organised continuously
            stopifnot(isTRUE(all.equal(trace_line_data[[i-1]][2,], 
                                       trace_line_data[[i]][1,])))
            trace_line_matrix = rbind(trace_line_matrix, trace_line_data[[i]][1,])
        }
    }
    trace_line_matrix = rbind(trace_line_matrix, trace_line_data[[lt]][2,])

    new_trace = SpatialLinesDataFrame(
        SpatialLines(list(Lines(list(Line(trace_line_matrix)), ID='C_0')), 
            proj4string = CRS(proj4string(contours)) ),
        data = data.frame(level = trace_depth), 
        match.ID=FALSE )

    rownames(new_trace@data) = 'C_0'

    final_contours = rbind(new_trace, contours)

    return(final_contours)
}
