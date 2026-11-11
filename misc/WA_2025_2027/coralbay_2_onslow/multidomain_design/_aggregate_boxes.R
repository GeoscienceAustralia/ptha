
# Returns TRUE if the right-edge of box b1 is
# equal to the left edge of box b2
matching_edge_EW<-function(b1, b2){
    return( # Bottom right b1 == Bottom left b2
     (b1$urx == b2$llx & b1$lly == b2$lly) &
      # Top right b1 == Top left b2
     (b1$urx == b2$llx & b1$ury == b2$ury))
}

# Symmetric version of matching_edge_EW
matching_side_edge<-function(b1, b2){
    return(matching_edge_EW(b1, b2) | matching_edge_EW(b2, b1))
}

# Returns TRUE if the top-edge of box b1 is
# equal to the bottom edge of box b2. 
matching_edge_NS<-function(b1, b2){
    return(# Top left b1 == bottom left b2
     (b1$llx == b2$llx & b1$ury == b2$lly) &
     # Top right b1 == Bottom right b2
     (b1$urx == b2$urx & b1$ury == b2$lly) )
}

# Symmetric version of matching_edge_NS
matching_top_edge<-function(b1, b2){
    return((matching_edge_NS(b1, b2) | matching_edge_NS(b2, b1)))
}

merge_boxes<-function(b1, b2){
    b3 = b1
    b3$llx = min(b1$llx, b2$llx)                
    b3$lly = min(b1$lly, b2$lly)
    b3$urx = max(b1$urx, b2$urx)
    b3$ury = max(b1$ury, b2$ury)
    return(b3)
}

expected_box_substepping<-function(max_depth, depth_lower_bound=40){

    grid_size = 1 # This factors out but is useful conceptually
    # Gravity-wave time-step on the box 
    expected_box_ts = grid_size/sqrt(9.8 * pmax(max_depth, depth_lower_bound))    
    # Gravity-wave time-step on the box if its depth were equal to the lower bound 
    max_box_ts = grid_size/sqrt(9.8 * depth_lower_bound)
    return(ceiling(max_box_ts/expected_box_ts))

}

#' Box merging
#'
#' Combine boxes if they have a 'similar' min elevation relative to our
#' depth lower bound
#'
#' @param boxes a data.frame with columns llx, lly, urx, ury (in order) giving
#' the lower left and upper right box coordinates
#' @param min_elev the minimum elevation in each box
#' @param depth_lower_bound parameter influencing the merging of boxes with similar depth.
#' @return a list with new boxes made by merging the old ones, and associated min_elev values. 
aggregate_boxes<-function(boxes, min_elev, depth_lower_bound=40){

    original_boxes = boxes

    # Estimate of how many sub-time-steps a domain would take, if the global
    # time step corresponds to the gravity wave speed at depth_lower_bound
    values = expected_box_substepping(-min_elev, depth_lower_bound)
    
    box_area = rep(1, length=nrow(boxes))

    counter = 0
    while(TRUE){
        counter = counter+1

        have_deleted_boxes = 0

        # Merge side edges that are identical
        nb = nrow(boxes)
        delete_box = rep(FALSE, nb)
        touched_box = rep(FALSE, nb)

        i_seq = seq(1,(nb-1), by=1)
        for(i in i_seq){
            if(touched_box[i]) next
            b1 = boxes[i,]

            j_seq = seq(i+1, nb, by=1)

            for(j in j_seq){
                if(touched_box[j]) next
                b2 = boxes[j,]
                if(matching_side_edge(b1, b2) & values[i] == values[j] &
                   (box_area[i]*values[i] + box_area[j]*values[j]) <= 9){
                    b3 = merge_boxes(b1, b2)
                    boxes[i,] = b3
                    boxes[j,] = b3
                    delete_box[i] = TRUE
                    min_elev[i] = min(min_elev[i], min_elev[j])
                    min_elev[j] = min_elev[i]
                    touched_box[i] = TRUE
                    touched_box[j] = TRUE
                    box_area[i] = box_area[i] + box_area[j]
                    box_area[j] = box_area[i]
                }
            }
        }

        # Remove 'double-up' boxes
        k = which(!delete_box)
        boxes = boxes[k,]
        values = values[k]
        min_elev = min_elev[k]
        box_area = box_area[k]

        have_deleted_boxes = have_deleted_boxes + sum(delete_box)

        #print(mean(delete_box))

        # Merge top edges that are identical
        nb = nrow(boxes)
        delete_box = rep(FALSE, nb)
        touched_box = rep(FALSE, nb)

        i_seq = seq(1,(nb-1), by=1)

        for(i in i_seq){
            if(touched_box[i]) next
            b1 = boxes[i,]

            j_seq = seq(i+1, nb, by=1)

            for(j in j_seq){
                if(touched_box[j]) next
                b2 = boxes[j,]
                if(matching_top_edge(b1, b2) & values[i] == values[j] &
                   (box_area[i]*values[i] + box_area[j]*values[j]) <= 9){
                    b3 = merge_boxes(b1, b2)
                    boxes[i,] = b3
                    boxes[j,] = b3
                    delete_box[i] = TRUE
                    min_elev[i] = min(min_elev[i], min_elev[j])
                    min_elev[j] = min_elev[i]
                    touched_box[i] = TRUE
                    touched_box[j] = TRUE
                    box_area[i] = box_area[i] + box_area[j]
                    box_area[j] = box_area[i]
                }
            }
        }

        k = which(!delete_box)
        boxes = boxes[k,]
        values = values[k]
        min_elev = min_elev[k]
        box_area = box_area[k]

        #print(mean(delete_box))

        have_deleted_boxes = have_deleted_boxes + sum(delete_box)

        # If we didn't change anything, there is no more deletion to do.
        if(have_deleted_boxes == 0) break
    }

    # Double check that all of the 'original box centres are only in one new box.
    # If not, there is a logic error in the algorithm. 
    # If yes, then the new boxes should by OK for the model.
    xy = cbind(0.5*(original_boxes[,1] + original_boxes[,3]),
               0.5*(original_boxes[,2] + original_boxes[,4]))
    for(i in 1:nrow(xy)){
        in_poly = (xy[i,1] > boxes[,1] & xy[i,2] > boxes[,2] &
                   xy[i,1] < boxes[,3] & xy[i,2] < boxes[,4])
        sip = sum(in_poly) 
        if(sip != 1){
            print(paste0('problem in box deaggregation at point: ', xy[i,1], ',', xy[i,2]))
            print(paste0('sip = ', sip))
            print(paste0('A zero value indicates "no box contains point", a value > 1 indicates "multiple boxes contain it"'))
            stop()
        }
    }

    return(list(boxes=boxes, min_elev=min_elev))
}

#' Make a spatial polygon from the domain_boundaries data.frame
#'
#' @param domain_boundaries a 4-column data.frame with lower-left-x, lower-left-y, upper-right-x, upper-right-y
#' for each rectangle
#' @param ID_flag Character to appear in the ID
#' @return Spatial Polygons
domain_boundaries_to_spatialpolygons<-function(domain_boundaries, ID_flag=''){
    srl = vector(mode='list', length=nrow(domain_boundaries))
    for(i in 1:length(srl)){
        tmp = as.numeric(domain_boundaries[i,])
        coords = matrix(c(tmp[c(1,2)], 
                          tmp[c(3,2)], 
                          tmp[c(3,4)], 
                          tmp[c(1,4)]),
                        ncol=2, byrow=TRUE)
        srl[[i]] = Polygons(list(Polygon(coords)), ID=paste0(ID_flag, i))
    }
    fine_sp = SpatialPolygons(srl)
    return(fine_sp)
}

# Convert the domain_boundaries data.frame to a shapefile
domain_boundaries_to_shapefile<-function(domain_boundaries, file_dsn, 
    file_layer, ID_flag='', proj4string=CRS("epsg:4326")){

    spdf = SpatialPolygonsDataFrame(
        domain_boundaries_to_spatialpolygons(domain_boundaries, ID_flag=ID_flag),
        data=cbind(domain_boundaries, 
                   data.frame(ID=paste0(ID_flag, 1:nrow(domain_boundaries)))),
        match.ID=FALSE)
    proj4string(spdf) = proj4string
    writeOGR(spdf, dsn=file_dsn, layer=file_layer, 
        driver='ESRI Shapefile', overwrite=TRUE)

    return(spdf)
}

#' Fill a polygon with a set of boxes (domains) with a given size.
#'
#' Convert a polygonal region to a set of boxes (rectangles) with size domain_size[1],
#' domain_size[2], which collectively cover the polygon. These boxes are aligned so that
#' "logically" if we extended the grid infinitely, the "grid_alignment_point" would be
#' a corner of a box.
#'
#' The idea is that we manually define a polygon that includes regions where we
#' need high-resolution domains, and use this function to define the resulting
#' domain boundaries. By using a consistent grid_alignment_point, we can make multiple
#' levels of grids nest nicely with each other.
#'
cover_polygon_with_rectangular_domains<-function(region_poly, domain_size, grid_alignment_point=c(0,0), make_plot=TRUE){

    # Define a grid that will cover the polygon
    x_range = range(region_poly[,1])
    y_range = range(region_poly[,2])

    # Ensure that the rectangles lie on a grid that "in-principle" could
    # include a domain with corner touching 'grid_alignment_point'.
    x_range[1] = grid_alignment_point[1] + floor( (x_range[1] - grid_alignment_point[1])/domain_size[1] )*domain_size[1]
    y_range[1] = grid_alignment_point[2] + floor( (y_range[1] - grid_alignment_point[2])/domain_size[2] )*domain_size[2]

    # Convert to sequence
    x_grid = seq(x_range[1], x_range[1] + ceiling(diff(x_range)/domain_size[1])*domain_size[1], by=domain_size[1])
    y_grid = seq(y_range[1], y_range[1] + ceiling(diff(y_range)/domain_size[2])*domain_size[2], by=domain_size[2])

    # Set of lower-left coordinates of 'possible' domains.
    # Not all will be inside the region_poly
    domain_ll = expand.grid(x_grid, y_grid)
    
    lower_left_upper_right = data.frame(
        llx = domain_ll[,1], lly = domain_ll[,2], 
        urx = domain_ll[,1] + domain_size[1], ury = domain_ll[,2] + domain_size[2])

    # Use rgeos for the intersection computation, to find domains inside or touching the region_poly
    # Thus we have to convert the data to SpatialPolygons
    pol_sp = domain_boundaries_to_spatialpolygons(lower_left_upper_right)
    pol_region = SpatialPolygons(list(Polygons(list(Polygon(coords=region_poly)), ID='a')))

    is_required = gIntersects(pol_sp, pol_region, byid=TRUE)
    keep = which(is_required)

    if(make_plot){
        plot(region_poly, t='l', asp=1)
        abline(v=x_grid, col='red')
        abline(h=y_grid, col='red')
        points(domain_ll[,1] + domain_size[1]/2, domain_ll[,2] + domain_size[2]/2, 
               col=c('white', 'red')[is_required+1], pch=19, cex=0.1)
    }

    
    return(lower_left_upper_right[keep,])
}

