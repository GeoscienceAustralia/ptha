#source('unit_sources.R')

#' Function to make a tsunami for the ij unit source in a discrete source. 
#' 
#' @param i integer. The down-dip index of the unit source
#' @param j integer. The along-strike index of the unit source
#' @param discrete_source List with discrete source information
#' @param tsunami_surface_points_lonlat matrix with lon/lat coordinates of
#' points where tsunami initial condition will be evaluated.
#' @param approx_dx Approximate x spacing of integration points inside unit
#' source. See unit_source_interior_points_cartesian.
#' @param approx_dy Approximate x spacing of integration points inside unit
#' source. See unit_source_interior_points_cartesian.
#' @param scale_dxdy Adjustment factor for approx_dx and approx_dy
#' @param depths_in_km logical. Are input depths in the discrete source in km
#' (TRUE) or m (FALSE).
#' @param tsunami_function Function used to make the tunami source from
#' a cartesian unit source with interior points
#' @return a list with the cartesian unit source, tsunami source, i,j indices
#' and tsunami surface points
#'
#' @export
make_tsunami_unit_source<-function(i,j, discrete_source, 
    tsunami_surface_points_lonlat, approx_dx = NULL, approx_dy = NULL, 
    scale_dxdy = 1, depths_in_km = TRUE, 
    tsunami_function = unit_source_cartesian_to_okada_tsunami_source){

    ## Get unit source in local cartesian coordinates + unit source statistics
    # By default the first us coordinate will be the origin
    us = unit_source_interior_points_cartesian(discrete_source, 
        unit_source_index = c(i,j),
        approx_dx=approx_dx, approx_dy=approx_dy, scale_dxdy=scale_dxdy,
        depths_in_km=depths_in_km)

    ## Convert tsunami surface points from lonlat to local cartesian
    ## coordinates
    tsunami_surface_points_cartesian = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat=us$origin_lonlat, r=us$r)

    ts = tsunami_function(us, tsunami_surface_points_cartesian)

    tsunami_unit_source = list(unit_source_interior_points = us, 
        tsunami_source = ts, i=i, j = j, 
        tsunami_surface_points_lonlat = tsunami_surface_points_lonlat)
    return(tsunami_unit_source)
}

#' Convert the tsunami unit source to a z-displacement raster
#' 
#' @param tsunami_unit_source output of make_tsunami_unit_source or similar
#' @param filename Name for output raster file. If NULL, no file is saved.
#' @param saveonly logical. If FALSE, return the raster object as well as
#' saving. If TRUE but filename is not NULL, then save the file, but return NULL
#' @return Either a raster, or NULL. Can save the raster to a file as a side effect
#' @export
tsunami_unit_source_2_raster<-function(tsunami_unit_source, filename=NULL, 
    saveonly=FALSE){

    #library(raster)
    #library(rgdal)

    xyz = cbind(tsunami_unit_source$tsunami_surface_points_lonlat, 
        tsunami_unit_source$tsunami_source$zdsp)

    outrast = rasterFromXYZ(xyz, crs=CRS('+init=epsg:4326'))

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
#' @param tsunami_surface_points_cartesian Points at which to compute the
#' tsunami deformation in cartesian coordinates
#' @param point_scale Multiply the length/width of area sources by point scale
#' prior to convolution, then divide by this afterwoulds. Allows moving between
#' Okada's area source representation at each grid point, and a point source
#' representation. Set to 1 for standard rectangular area source representation
#' @return List with edsp, ndsp, zdsp giving the displacements at the
#' tsunami_surface_points_cartesian.
#' @export
unit_source_cartesian_to_okada_tsunami_source<-function(us,
    tsunami_surface_points_cartesian, point_scale=1){

    deg2rad = pi/180

    src = us$grid_points
    nsrc = length(src[,'x'])
    strike = src[, 'strike']
    dip = src[, 'dip']
    depth = src[,'depth']/1000 # depth in km
    thrust_slip = depth*0 + 1 # Slip in m

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


    # If sources are protruding from the earth, adjust their width and length
    # to preserve area
    width_limit = 2*depth/sin(dip*deg2rad) - thrust_slip*1/1000
    too_shallow = which(src_wdt > width_limit)
    if(length(too_shallow) > 0){
        warning('Reducing source widths to prevent negative depths')
        #src_len[too_shallow] = src_len[too_shallow]*src_wdt[too_shallow]/width_limit[too_shallow]
        src_wdt[too_shallow] = width_limit[too_shallow]
    }

    # Our Okada function is for a rectangular source with constant
    # depth along-strike.
    # We can rescale length/width to be like a point source
    # Note okada_tsunami uses depth in km. 
    ts = okada_tsunami(
        elon = src[,'x'], elat = src[,'y'], edep = depth,
        strk = strike, dip = dip,
        lnth = src_len*point_scale, wdt = src_wdt*point_scale,
        disl1 = rep(0, len=nsrc), disl2 = thrust_slip,
        rlon = dest[,1], rlat = dest[,2],
        verbose=FALSE)

    # Rescale deformation to give correct area
    for(nn in names(ts)){
        ts[[nn]] = ts[[nn]]/(point_scale**2)
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
    tsunami_sum = all_tsunami[[1]]$tsunami_source$zdsp*0
    for(i in 1:length(all_tsunami)){ 
        tsunami_sum = tsunami_sum + all_tsunami[[i]]$tsunami_source$zdsp
    }
    #library(raster)
    tsunami_sum = rasterFromXYZ(
        cbind(all_tsunami[[1]]$tsunami_surface_points_lonlat, tsunami_sum),
        crs=CRS('+init=epsg:4326'))

    # Plot all
    pdf(paste0('Unit_source_data/', sourcename, '_plots.pdf'), width=10,
        height=10)
    plot(tsunami_sum, asp=1, main='Sum all unit sources')
    plot(tsunami_sum, asp=1, main='Sum all unit sources')

    # quick function to add source info to the plot
    add_discrete_source2plot<-function(){

        for(i in 1:dim(ds1$unit_source_grid)[1]){
            points(ds1$unit_source_grid[i,1,], ds1$unit_source_grid[i,2,],
                t='l', col='brown', lty='dotted')
        }

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
        plot_unit_source_interior_points_cartesian(
            all_tsunami[[i]]$unit_source_interior_points)
        title(main = titlei, line = -2)
        plot(all_tsunami_rast[[i]], asp=1, main = titlei)
        add_discrete_source2plot()
    }

    dev.off()
}

