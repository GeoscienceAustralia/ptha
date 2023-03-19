
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

    l0 = rbind(c(0, 0, d0),
               c(0, width, d1),
               c(0, 2*width, d1+(d1-d0)))

    l1 = l0
    l1[,1] = -len

    l2 = l1
    l2[,1] = -2*len

    unit_source_coords_cartesian = array(NA, dim=c(dim(l0), 3))
    unit_source_coords_cartesian[,,1] = l0
    unit_source_coords_cartesian[,,2] = l1
    unit_source_coords_cartesian[,,3] = l2
    

    # Back-calculate lon/lat
    origin = c(0,0)
    unit_source_coords_lonlat = unit_source_coords_cartesian
    for(i in 1:3){
        unit_source_coords_lonlat[,1:2,i] = cartesian2d_to_spherical_coordinates(
            unit_source_coords_cartesian[,1:2,i], origin=origin)
    }

    # Make discrete source consisting of a 2x2 set of unit sources
    ds = list()
    ds$unit_source_grid = unit_source_coords_lonlat
    ds$unit_source_grid[,3,] = ds$unit_source_grid[,3,]/1000
   
    ds$discretized_source_dim = c(2,2) 
    names(ds$discretized_source_dim) = c('dip','strike')

    mid_line_with_cutpoints = list()
    for(i in 1:3) mid_line_with_cutpoints[[i]] = ds$unit_source_grid[,,i]
    ds$mid_line_with_cutpoints = mid_line_with_cutpoints

    ds$fine_downdip_transects = ds$unit_source_grid

    us = unit_source_interior_points_cartesian(ds, unit_source_index = c(1,1),
        approx_dx = NULL, approx_dy = NULL, depths_in_km=TRUE)

    # Check strike
    expect_true(all(us$grid_points[,'strike'] == strike))

    # Check dip (small inaccuracies allowed due to numerical optimization)
    expect_true(all(abs(us$grid_points[,'dip'] - dip)<0.02))

    expect_true(all(us$grid_points[,'unit_slip_scale'] == 1))
    #expect_true(all(us$grid_points[,'fraction_area_inside_unit_source'] == 1))
    # Following updated rgeos functions
    expect_true(all(abs(us$grid_points[,'fraction_area_inside_unit_source'] - 1) < 1e-10))

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

    ##
    ## Case with edge tapering
    ##

    us2 = unit_source_interior_points_cartesian(ds, unit_source_index = c(1,1),
        approx_dx = 1000, approx_dy = 1000, depths_in_km=TRUE, 
        edge_taper_width=3000)

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

