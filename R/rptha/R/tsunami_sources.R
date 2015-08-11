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
#' @return List with edsp, ndsp, zdsp giving the displacements at the
#' tsunami_surface_points_cartesian.
#' @export
unit_source_cartesian_to_okada_tsunami_source<-function(us,
    tsunami_surface_points_cartesian){
    ## Call the tsunami source function

    #library(EqSim)

    src = us$grid_points
    nsrc = length(src[,'x'])
    # Make a 'length' for the point source in km
    src_len = rep( sqrt(us$dx*us$dy)/1000, len=nsrc) 
    # Make a 'width' for the point source in km
    area_projector = sqrt(1 + tan(src[,'dip']/180*pi)**2)
    src_wdt = (src[,'area']/1e+06)/src_len*area_projector

    dest = tsunami_surface_points_cartesian

    # Strike/Dip were already computed
    strike = src[, 'strike']
    dip = src[, 'dip']

    # Our Okada function is for a rectangular source with constant
    # depth along-strike.
    # We will rescale length/width to be like a point source
    point_scale = 0.01
    # Note okada_tsunami uses depth in km. 
    ts = okada_tsunami(
        elon = src[,'x'], elat = src[,'y'], edep = src[,'depth']/1000,
        strk = strike, dip = dip,
        lnth = src_len*point_scale, wdt = src_wdt*point_scale,
        disl1 = rep(0, len=nsrc), disl2 = rep(1, len=nsrc),
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

