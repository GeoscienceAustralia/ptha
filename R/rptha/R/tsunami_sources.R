#' unit_source water surface deformation 
#'
#' Function to make a tsunami for the ij unit source in a discrete source. 
#' 
#' @param i integer. The down-dip index of the unit source
#' @param j integer. The along-strike index of the unit source
#' @param discrete_source List with discrete source information (e.g. output of
#' \code{discretized_source_from_source_contours})
#' @param rake Rake (degrees) of the slip vector. 0 is along-strike slip, 90 is pure thrust,..
#' @param tsunami_surface_points_lonlat matrix with lon/lat coordinates of
#' points where tsunami initial condition will be evaluated.
#' @param approx_dx Approximate x spacing of integration points inside unit
#' source. See unit_source_interior_points_cartesian.
#' @param approx_dy Approximate x spacing of integration points inside unit
#' source. See unit_source_interior_points_cartesian.
#' @param scale_dxdy Adjustment factor for approx_dx and approx_dy
#' @param depths_in_km logical. Are input depths in the discrete source in km
#' (TRUE) or m (FALSE).
#' @param kajiura_smooth logical. Should we apply kajiura smoothing to the
#' z-displacement computed by tsunami_function?
#' @param surface_point_ocean_depths numeric vector giving ocean depths (in m)
#' at the tsunami surface points. Used only for Kajiura filtering
#' @param kajiura_grid_spacing Value used for grid_dx and grid_dy in kajiura filter.
#' If NULL, max(surface_point_ocean_depths)/2 is used
#' @param kajiura_volume_change_error_threshold Value of volume_change_error_threshold
#' passed to kajiura_filter. 
#' @param minimal_output Logical. If TRUE, set the unit source and
#' tsunami_surface_points_lonlat to NA in the outputs. These are memory heavy so
#' in some settings they are best removed.
#' @param kajiura_where_deformation_exceeds_threshold Numeric. The Kajiura
#' filter is applied to a rectangular region including all points where the
#' absolute value of the initial deformation (computed with
#' \code{tsunami_function}) exceeds the given threshold. Use of non-zero values
#' much smaller than the significant deformation (e.g. 1.0e-04) can lead to a
#' large reduction in the memory and time required for the Kajiura filter
#' computation, although if the threshold is too high then the filter might
#' not be applied to important areas.
#' @param tsunami_function Function used to make the tunami source from
#' a cartesian unit source with interior points, the rake, and the surface
#' points.  Must return a list containing zdsp (a vector of all z displacements
#' in m), and perhaps also edsp, ndsp. 
#' @param verbose logical. Print some information on progress
#' @param edge_taper_width Distance over which to taper slip on unit sources. 
#' Values > 0 can lead to less peaked slip on sources with a shallow top edge
#' depth. Note values>0 also imply some slip occurring outside the unit source
#' location (due to smoothing), although we ensure seismic moment is conserved
#' @param ... further arguments to tsunami function
#' @return a list with the cartesian unit source, tsunami source, i,j indices
#' and tsunami surface points, smoothed deformation, and rake
#'
#' @export
make_tsunami_unit_source<-function(
    i, 
    j, 
    discrete_source, 
    rake,
    tsunami_surface_points_lonlat, 
    approx_dx = NULL, 
    approx_dy = NULL, 
    scale_dxdy = 1, 
    depths_in_km = TRUE, 
    kajiura_smooth=FALSE, 
    surface_point_ocean_depths=NULL, 
    kajiura_grid_spacing=NULL,
    kajiura_volume_change_error_threshold = 0.1, 
    kajiura_where_deformation_exceeds_threshold = 0.0,
    minimal_output=FALSE,
    tsunami_function = unit_source_cartesian_to_okada_tsunami_source,
    verbose=FALSE,
    edge_taper_width=0,
    ...){

    if(verbose){
        print('Making sub-unit-source points...')
        t0 = Sys.time()
    }
    
    ## Get unit source in local cartesian coordinates + unit source statistics
    # By default the first us coordinate will be the origin
    us = unit_source_interior_points_cartesian(
        discrete_source, 
        unit_source_index = c(i,j),
        approx_dx=approx_dx, 
        approx_dy=approx_dy, 
        scale_dxdy=scale_dxdy,
        depths_in_km=depths_in_km,
        edge_taper_width=edge_taper_width)

    if(verbose){
        print(Sys.time() - t0)
        print('')
    }

    # Convert tsunami surface points from lonlat to local cartesian
    # coordinates
    tsunami_surface_points_cartesian = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat=us$origin_lonlat, r=us$r)

    if(verbose){
        print('Computing co-seismic deformation..')
        t0 = Sys.time()
    }

    # Compute the tsunami deformation (without Kajiura smoothing)
    ts = tsunami_function(us, rake, tsunami_surface_points_cartesian, ...)

    if(verbose){
        print(Sys.time() - t0)
        print('')
    }

    # Optionally apply Kajiura smoothing
    if(kajiura_smooth){

        if(verbose){
            print('Applying Kajiura filter')
            t0 = Sys.time()
        }

        stopifnot(length(surface_point_ocean_depths) == length(tsunami_surface_points_lonlat[,1]))

        if(is.null(kajiura_grid_spacing)){
            grid_dx = max(surface_point_ocean_depths)/2
            grid_dy = grid_dx
        }else{
            grid_dx = kajiura_grid_spacing
            grid_dy = kajiura_grid_spacing
        }

        # Rotate the surface points so the strike direction is aligned with
        # pixels, to get better behaviour near discontinuities in Kajiura.
        # Artefacts can be caused by the re-interpolation to a regular grid in that routine
        strike_vec = us$unit_source_cartesian[4,1:2] - us$unit_source_cartesian[1,1:2]
        strike_vec = strike_vec/sqrt(sum(strike_vec**2))

        new_xy = rotate_cartesian2d(tsunami_surface_points_cartesian[,1:2], 
            origin = c(0,0), x_axis_vector = strike_vec)

        # For computational efficiency we only apply Kajiura to a region covering
        # the extent over which the absolute value of the initial deformation is > threshold
        # This should include the important parts but ignore large zero areas
        tmp = which(abs(ts$zdsp) > kajiura_where_deformation_exceeds_threshold)
        k_bbox = rbind(range(new_xy[tmp,1]), range(new_xy[tmp,2]))
        kajiura_inds = which(
            (new_xy[,1] >= k_bbox[1,1]) & (new_xy[,1] <= k_bbox[1,2]) &
            (new_xy[,2] >= k_bbox[2,1]) & (new_xy[,2] <= k_bbox[2,2]))
        rm(tmp, k_bbox)
        

        # Call Kajiura, interpolating separately over (rotated) y > 0 and y <=0
        # This can help reduce artefacts if we have a discontinuity in the deformation
        # at y=0 (as occurs for ruptures to the trench, e.g. Goda 2015 BSSA)
        kajiura_source = kajiura_filter(
            cbind(new_xy[kajiura_inds,], ts$zdsp[kajiura_inds]),
            surface_point_ocean_depths[kajiura_inds], 
            grid_dx = grid_dx, grid_dy=grid_dy,
            volume_change_error_threshold = kajiura_volume_change_error_threshold,
            interpolator='linear', interpolator_categories=function(xy){xy[,2]>0})

        smooth_tsunami_displacement = ts$zdsp 
        smooth_tsunami_displacement[kajiura_inds] = kajiura_source[,3]

        # Careful with memory usage (in parallel, auto garbage collection not so good)
        rm(kajiura_source, new_xy, kajiura_inds); gc()
        if(verbose){
            print(Sys.time() - t0)
            print('')
        }
    }else{
        # In this case just repeat the tsunami source displacement
        smooth_tsunami_displacement = ts$zdsp
    }

    # Option to save memory
    if(minimal_output){
        us = NA
        tsunami_surface_points_lonlat = NA
        ts = NA
    }

    tsunami_unit_source = list(unit_source_interior_points = us, 
        smooth_tsunami_displacement = smooth_tsunami_displacement,
        tsunami_source = ts, i=i, j = j, 
        tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
        rake = rake)

    # Force garbage collection since in parallel, this might not detect
    # high overall memory usage
    rm(us, smooth_tsunami_displacement, ts, tsunami_surface_points_lonlat, rake)
    gc()

    return(tsunami_unit_source)
}

#' I introduced this to work-around an old bug in library(raster),
#' but that was reported and fixed some time ago. 
#' @import raster
#' @import stats
.local_rasterFromXYZ<-function(xyz, crs = NA, res=c(NA,NA), ...){

    xs = sort(unique(xyz[,1]))
    ys = sort(unique(xyz[,2]))

    dxs = mean(diff(xs))
    dys = mean(diff(ys))

    # If res is provided, use it
    if(!is.na(res[1])){
        stopifnot(abs(dxs - res[1])/res[1] < 1.0e-05)
        dxs = res[1]
    }
    if(!is.na(res[2])){
        stopifnot(abs(dys - res[2])/res[2] < 1.0e-05)
        dys = res[2]
    }

    # Check that x and y are evenly spaced
    sd_x = sd(diff(xs))
    sd_y = sd(diff(ys))
    stopifnot(sd_x/dxs < 1.0e-05)
    stopifnot(sd_y/dys < 1.0e-05)

    # Find raster extents
    xmn = xs[1] - dxs*0.5
    xmx = xs[length(xs)] + dxs*0.5
    ymn = ys[1] - dys*0.5
    ymx = ys[length(ys)] + dys*0.5

    output = raster(xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs, res=c(dxs, dys))
    cells = cellFromXY(output, xyz[,1:2])

    output[cells] = xyz[,3]
    return(output)
}


#' Convert the tsunami unit source to a z-displacement raster
#' 
#' @param tsunami_unit_source output of make_tsunami_unit_source or similar.
#' @param filename Name for output raster file. If NULL, no file is saved.
#' @param saveonly logical. If FALSE, return the raster object as well as
#' saving. If TRUE but filename is not NULL, then save the file, but return NULL.
#' @param tsunami_surface_points_lonlat matrix with 2 columns containing surface points
#' at which tsunami_unit_source$smooth_tsunami_displacement occurs. If NULL, look for
#' this in tsunami_unit_source -- however, to save memory, the latter may be set to NA. In
#' which case this argument must be provided.
#' @param res optional argument 'res' to pass to \code{rasterFromXYZ}
#' @return Either a raster, or NULL. Can save the raster to a file as a side effect.
#' @export
tsunami_unit_source_2_raster<-function(tsunami_unit_source, filename=NULL, 
    saveonly=FALSE, tsunami_surface_points_lonlat=NULL, res=c(NA,NA)){

    #library(raster)
    #library(rgdal)

    if(is.null(tsunami_surface_points_lonlat)){

        stopifnot(length(tsunami_unit_source$tsunami_surface_points_lonlat[,1]) == 
                  length(tsunami_unit_source$smooth_tsunami_displacement))

        xyz = cbind(tsunami_unit_source$tsunami_surface_points_lonlat, 
            tsunami_unit_source$smooth_tsunami_displacement)
    }else{

        stopifnot(length(tsunami_surface_points_lonlat[,1]) == 
                  length(tsunami_unit_source$smooth_tsunami_displacement))

        xyz = cbind(tsunami_surface_points_lonlat, 
            tsunami_unit_source$smooth_tsunami_displacement)
    }

    # Currently use local version of rasterFromXYZ to deal with some bugs in it
    # (which have been reported)
    outrast = .local_rasterFromXYZ(xyz, crs=CRS('+init=epsg:4326'), res=res)
    #outrast = rasterFromXYZ(xyz, crs=CRS('+init=epsg:4326'), res=res)

    # Free up memory
    rm(xyz); gc()

    if(!is.null(filename)){
        writeRaster(outrast, filename, driver='GTiff', 
            options='COMPRESS=DEFLATE', overwrite=TRUE)
    }

    if(saveonly){
        return()
    }else{
        return(outrast)
    }
}

#' Convert a cartesian unit source with interior points to an okada tsunami source
#'
#' @param us A unit source with interior points in cartesian coordinates.
#' @param rake The rake of the slip in degrees
#' @param tsunami_surface_points_cartesian Points at which to compute the
#' tsunami deformation in cartesian coordinates
#' @param dstmx For each sub-source, only compute the deformation to a
#' horizontal distance of dstmx*source depth. Using a low value of dstmax can
#' speed up the algorithm if the displacement for each sub-source does not have
#' to be computed for all tsunami_surface_points_cartesian. e.g. I have seen a
#' value of 50 work well in one instance.
#' @param upper_depth_limit Limit for the top-depth of any unit source (km)
#' @param cell_integration_scale vector of length 2, giving a dx/dy value. 
#' If c(0,0), then do nothing. Otherwise we compute a 16 point smooth within each cell where the
#' initially computed abs(deformation) exceeds 10% of its maximum. The smooth
#' is the average okada deformation, computed over the tsunami_surface_points_cartesian 
#' perturbed by dx and +-dy multiplied by c(-1, -1/3, 1/3, 1) [i.e. 4x4 = 16 points] 
#' @return List with edsp, ndsp, zdsp giving the displacements at the
#' tsunami_surface_points_cartesian.
#' @export
unit_source_cartesian_to_okada_tsunami_source<-function(us, rake,
    tsunami_surface_points_cartesian, dstmx = 9e+20,
    upper_depth_limit = 0.0e-03, cell_integration_scale=c(0,0)){

    deg2rad = pi/180

    src = us$grid_points
    nsrc = length(src[,'x'])
    strike = src[, 'strike']
    dip = src[, 'dip']
    depth = src[,'depth']/1000 # depth in km
    # Unit slip in m, accounting for tapering (variable unit_slip_scale)
    thrust_slip = src[,'unit_slip_scale'] * sin(rake*deg2rad) 
    strike_slip = src[,'unit_slip_scale'] * cos(rake*deg2rad)

    dest = tsunami_surface_points_cartesian


    # Choose equivalent rectangular length and width for each grid point.
    #
    # This is non-trivial near boundaries of the unit source, where the area
    # associated with each grid point may be cut so it does not exceed the unit
    # source boundaries.
    #
    # A careful approach to the integration is required to avoid oscillation in
    # the tsunami source, unless an extremely fine grid is used. The approach
    # here is:
    # Generate a 'surface width' scale for the grid point area, by taking the
    # dot product of the vectors making the sides of each grid-points polygon
    # with a unit vector pointing up-dip, and summing their absolute values. To
    # get an 'alongstrike-length' scale, do the same for a unit vector pointing
    # along strike. The ratio of these gives the chosen ratio of the surface
    # width/length for the equivalent rectangular element.

    cs = cos(strike*deg2rad)
    ss = sin(strike*deg2rad)
    updip_vec = cbind(-cs, ss)
    alongstrike_vec = cbind(ss, cs) 

    updip_scale = updip_vec[,1]*NA
    alongstrike_scale = alongstrike_vec[,1]*NA
    for(i in 1:length(updip_scale)){
        poly_coords = us$grid_point_polygon@polygons[[i]]@Polygons[[1]]@coords

        # Matrix with rows giving edge vectors
        region_boundary_vectors = diff(poly_coords)

        # Compute scale for 'up-dip' dimension
        updip_scale[i] = sum(abs(region_boundary_vectors%*%updip_vec[i,]))
        # Compute scale for 'along-strike' dimension
        alongstrike_scale[i] = sum(abs(region_boundary_vectors%*%alongstrike_vec[i,]))

    }

    area_scale = updip_scale*alongstrike_scale

    # Make a 'length' for the point source in km
    src_len = (sqrt(src[, 'area_projected']/area_scale)*alongstrike_scale) * 1/1000
    # Make a 'width' for the point source in km, adjusted for down-dip distance
    area_downslope_adjust = sqrt(1 + tan(dip*deg2rad)**2)
    src_wdt = (sqrt(src[, 'area_projected']/area_scale)*updip_scale)*area_downslope_adjust * 1/1000

    stopifnot(isTRUE(all.equal(
        sum(src_len*src_wdt), sum(src[,'area_projected']*area_downslope_adjust)/1e+06
        )))

    # If sources are protruding from the earth, adjust their width 
    width_limit = 2*(depth-upper_depth_limit)/sin(dip*deg2rad) 
    too_shallow = which(src_wdt > width_limit)
    if(length(too_shallow) > 0){
        warning('Reducing source widths to prevent negative depths')
        #src_len[too_shallow] = src_len[too_shallow]*src_wdt[too_shallow]/width_limit[too_shallow]
        src_wdt[too_shallow] = width_limit[too_shallow]
    }

    # Our Okada function integrates over each rectangular source with constant
    # depth along-strike
    ts = okada_tsunami(
        elon = src[,'x'], 
        elat = src[,'y'], 
        edep = depth,
        strk = strike, 
        dip = dip,
        lnth = src_len, 
        wdt = src_wdt,
        disl1 = strike_slip, 
        disl2 = thrust_slip,
        rlon = dest[,1], 
        rlat = dest[,2], 
        dstmx = dstmx,
        verbose=FALSE)

    # Optionally apply sub-grid sampling to points with significant deformation
    # This is used to avoid issues with 'spikes' etc around the trench
    if(!all(cell_integration_scale == 0)){


        # Find points in a 'box' where deformation is significant
        # FIXME: Hardcoded to 10% of maximum
        largedef = which(abs(ts$zdsp) > 0.1*max(abs(ts$zdsp)))
        lonrange = range(dest[largedef,1])
        latrange = range(dest[largedef,2])

        sp = which(dest[,1] >= lonrange[1] & dest[,1] <= lonrange[2] &
            dest[,2] >= latrange[1] & dest[,2] <= latrange[2])


        counter = 0
        di = cell_integration_scale[1] * c(-1, -1/3, 1/3, 1)
        dj = cell_integration_scale[2] * c(-1, -1/3, 1/3, 1)
        for(i in di){
            for(j in dj){

                # Compute the deformation at coordiantes perturbed by i/j
                ts_local = okada_tsunami(
                    elon = src[,'x'], 
                    elat = src[,'y'], 
                    edep = depth,
                    strk = strike, 
                    dip = dip,
                    lnth = src_len, 
                    wdt = src_wdt,
                    disl1 = strike_slip, 
                    disl2 = thrust_slip,
                    rlon = dest[sp,1] + i, 
                    rlat = dest[sp,2] + j, 
                    dstmx = dstmx,
                    verbose=FALSE)

                if(counter == 0){
                    ts_sp = ts_local
                }else{
                    for(ii in 1:length(ts)){
                        # Summation [divide later to get average]
                        ts_sp[[ii]] = ts_sp[[ii]] + ts_local[[ii]]

                        # Store result for median
                        #ts_sp[[ii]] = cbind(ts_sp[[ii]], ts_local[[ii]])
                    }
                }
                counter = counter+1
            }
        }
        rm(ts_local); gc()
        # Divide to get the average
        for(i in 1:length(ts)) ts_sp[[i]] = ts_sp[[i]]/counter

        # row-wise median 
        #for(i in 1:length(ts)) ts_sp[[i]] = apply(ts_sp[[i]], 1, median)

        # Update original ts
        for(i in 1:length(ts)) ts[[i]][sp] = ts_sp[[i]]
    }

    return(ts)
}


#' Make plot of all tsunami unit sources and their sum
#' 
#' @param sourcename name for source zone
#' @param all_tsunami list with the tsunami unit sources for all source zones
#' @param all_tsunami_rast list of rasters corresponding to all_tsunami
#' @param discrete_source discrete source information for sourcename
#' @return Nothing, but make a plot
#'
#' @export
plot_all_tsunami_unit_sources<-function(sourcename, all_tsunami, 
    all_tsunami_rast, discrete_source){

    ds1 = discrete_source

    # Make raster with 'sum of all unit source tsunami'
    tsunami_sum = all_tsunami_rast[[1]]*0
    for(i in 1:length(all_tsunami_rast)) tsunami_sum = tsunami_sum + all_tsunami_rast[[i]]

    # Plot all
    pdf(paste0('Unit_source_data/', sourcename, '_plots.pdf'), width=10,
        height=10)
    plot(tsunami_sum, asp=1, main='Sum all unit sources')
    plot(tsunami_sum, asp=1, main='Sum all unit sources')

    # quick function to add source info to the plot
    add_discrete_source2plot<-function(){

        # Plot along-strike lines
        for(i in 1:dim(ds1$unit_source_grid)[1]){
            points(ds1$unit_source_grid[i,1,], ds1$unit_source_grid[i,2,],
                t='l', col='brown', lty='dotted')
        }
        # Add the shallowest one in a different colour
        points(ds1$unit_source_grid[1,1,], ds1$unit_source_grid[1,2,],
            t='l', col='red', lty='dashed')

        # Plot down-dip lines
        for(i in 1:dim(ds1$unit_source_grid)[3]){
            points(ds1$unit_source_grid[,1:2,i], t='l', col='brown', 
                lty='dotted')
        }

        plot(ds1$depth_contours, add=TRUE, col = 'grey', alpha=0.3)

    }
    add_discrete_source2plot()

    # Plot ever unit source
    for(i in 1:length(all_tsunami_rast)){
        titlei = paste0('Dip index: ', all_tsunami[[i]]$i, ' ; Strike index: ',
            all_tsunami[[i]]$j)

        if(!is.na(all_tsunami[[i]]$unit_source_interior_points[1])){
            # This information might not be stored (if minimal_output=TRUE in
            # make_tsunami_unit_source)
            plot_unit_source_interior_points_cartesian(
                all_tsunami[[i]]$unit_source_interior_points)
            title(main = titlei, line = -2)
        }

        plot(all_tsunami_rast[[i]], asp=1, main = titlei)
        add_discrete_source2plot()
    }

    dev.off()
}

