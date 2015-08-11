
context('test_unit_source_interior_points_cartesian')

test_that('test_unit_source_interior_points_cartesian', {
#test_unit_source_interior_points_cartesian<-function(){

    # Consider a source with width 20km, length 40km, top-edge-depth
    # of 6km, bottom-edge-depth of 10km
    
    d0 = 6000
    d1 = 10000
    width = 20000
    len = 40000
    dip = atan((d1 - d0)/width)*180/pi
    strike = 270    

    unit_source_coords_cartesian=rbind(c(0, 0, d0),
                                       c(0, width, d1),
                                       c(-len, width, d1),
                                       c(-len, 0, d0))

    # Back-calculate lon/lat
    unit_source_coords_lonlat = cartesian2d_to_spherical_coordinates(
        unit_source_coords_cartesian[,1:2], origin=c(0,0))

    # Make discrete source consisting of a single unit source
    ds = list()
    ds$unit_source_grid = array(dim=c(2,3,2))
    ds$unit_source_grid[,,1] = cbind(unit_source_coords_lonlat[1:2,1:2], 
        unit_source_coords_cartesian[1:2,3]/1000)
    ds$unit_source_grid[,,2] = cbind(unit_source_coords_lonlat[4:3,1:2], 
        unit_source_coords_cartesian[4:3,3]/1000)
   
    ds$discretized_source_dim = c(1,1) 

    ds$fine_downdip_transects = ds$unit_source_grid    

    us = unit_source_interior_points_cartesian(ds, unit_source_index = c(1,1),
        approx_dx = NULL, approx_dy = NULL, depths_in_km=TRUE)

    # Check strike
    expect_that(all(us$grid_points[,'strike'] == strike), is_true())

    # Check dip (small inaccuracies allowed due to numerical optimization)
    expect_that(all(abs(us$grid_points[,'dip'] - dip)<0.02), is_true())

    # Check depth (small inaccuracies allowed due to numerical optimization)
    pred_depth = (us$grid_points[,'y']/width)*(d1 - d0) + d0

    expect_that(all(abs(us$grid_points[,'depth'] - pred_depth) < 1), is_true())

    # Check area in m^2
    expect_that(abs(sum(us$grid_points[,'area']) - len*width) < 0.1, is_true())

    # Check x/y
    expect_that( (min(us$grid_points[,'x']) > -len) &
                 (max(us$grid_points[,'x']) < 0) &
                 (min(us$grid_points[,'y']) > 0) &
                 (max(us$grid_points[,'y']) < width), is_true() )
}

