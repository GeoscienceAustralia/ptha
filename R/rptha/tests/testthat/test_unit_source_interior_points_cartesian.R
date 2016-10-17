
context('test_unit_source_interior_points_cartesian')

test_that('test_unit_source_interior_points_cartesian', {

    # Consider a source with width 20km, length 40km, top-edge-depth
    # of 6km, bottom-edge-depth of 10km
    
    d0 = 6000
    d1 = 10000
    width = 20000
    len = 40000
    dip = atan((d1 - d0)/width)*180/pi
    strike = 270    

    unit_source_coords_cartesian = rbind(
        c(0, 0, d0),
        c(0, width, d1),
        c(-len, width, d1),
        c(-len, 0, d0))

    # Back-calculate lon/lat
    unit_source_coords_lonlat = cartesian2d_to_spherical_coordinates(
        unit_source_coords_cartesian[,1:2], origin_lonlat=c(0,0))

    # Make discrete source consisting of a single unit source
    ds = list()
    ds$unit_source_grid = array(dim=c(2,3,2))
    ds$unit_source_grid[,,1] = cbind(unit_source_coords_lonlat[1:2,1:2], 
        unit_source_coords_cartesian[1:2,3]/1000)
    ds$unit_source_grid[,,2] = cbind(unit_source_coords_lonlat[4:3,1:2], 
        unit_source_coords_cartesian[4:3,3]/1000)

    mid_line_with_cutpoints = list()
    mid_line_with_cutpoints[[1]] = cbind(unit_source_coords_lonlat[1:2,], c(d0, d1)/1000)
    mid_line_with_cutpoints[[2]] = cbind(unit_source_coords_lonlat[4:3,], c(d0, d1)/1000)
    ds$mid_line_with_cutpoints = mid_line_with_cutpoints
    
    ds$discretized_source_dim = c(1,1) 
    names(ds$discretized_source_dim) = c('dip', 'strike')

    ds$fine_downdip_transects = ds$unit_source_grid    

    us = unit_source_interior_points_cartesian(ds, unit_source_index = c(1,1),
        approx_dx = NULL, approx_dy = NULL, depths_in_km=TRUE)

    # Check strike
    expect_true(all(us$grid_points[,'strike'] == strike))

    # Check dip (small inaccuracies allowed due to numerical optimization)
    expect_true(all(abs(us$grid_points[,'dip'] - dip)<0.02))

    expect_true(all(us$grid_points[,'unit_slip_scale'] == 1))
    expect_true(all(us$grid_points[,'fraction_area_inside_unit_source'] == 1))

    # Check depth (small inaccuracies allowed due to numerical optimization)
    pred_depth = (us$grid_points[,'y']/width)*(d1 - d0) + d0

    expect_true(all(abs(us$grid_points[,'depth'] - pred_depth) < 1))

    # Check area in m^2
    expect_true(abs(sum(us$grid_points[,'area_projected']) - len*width) < 0.1)

    # Check x/y
    expect_true( (min(us$grid_points[,'x']) > -len) &
                 (max(us$grid_points[,'x']) < 0) &
                 (min(us$grid_points[,'y']) > 0) &
                 (max(us$grid_points[,'y']) < width))

    us2 = unit_source_interior_points_cartesian(ds, unit_source_index = c(1,1),
        approx_dx = NULL, approx_dy = NULL, depths_in_km=TRUE, 
        edge_taper_width=3000,
        allow_points_outside_discrete_source_outline=TRUE)


    # Check that the moment-normalization has worked
    a0 = sum(us$grid_points[,'area_projected'] * 
        sqrt(1 + atan(us$grid_points[,'dip']/180 * pi)**2))
    a1 = sum(us2$grid_points[,'unit_slip_scale'] * 
        us2$grid_points[,'area_projected'] *
        sqrt(1 + atan(us2$grid_points[,'dip']/180 * pi)**2))
    expect_true(abs(a0 - a1) < 1.0e-12 * a1)

    # Check that the area normalization has worked
    a0 = sum(us$grid_points[,'area_projected'])
    a1 = sum(us2$grid_points[,'area_projected'] * 
        us2$grid_points[,'fraction_area_inside_unit_source'])
    expect_true(abs(a0 - a1) < 1.0e-12 * a1)
       
})

