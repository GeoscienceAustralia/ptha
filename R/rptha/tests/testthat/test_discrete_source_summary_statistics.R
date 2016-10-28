
context('test_discrete_sources')

test_that('test_discrete_sources_and_summary_statistics', {

    interface_shapefile = 'testshp/sagami.shp'

    desired_subfault_length = 100
    desired_subfault_width = 50
    
    # Create sagami source using 'old' discretization method 
    discrete_source1 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=FALSE)

    # Test that our outline of the discrete source is ok
    # Do this by taking the outline, buffering slightly, and checking that all
    # unit_source_grid points are inside it
    discrete_source1_outline = get_discretized_source_outline(discrete_source1)

    p0 = SpatialPolygons(list(
            Polygons(list(Polygon(discrete_source1_outline[,1:2])), ID='P')),
        proj4string=CRS(""))

    p0_buf = gBuffer(p0, width=1)

    grid_pts = discrete_source1$unit_source_grid[,1:2,1]
    for(i in 2:(dim(discrete_source1$unit_source_grid)[3])){
        grid_pts = rbind(grid_pts, discrete_source1$unit_source_grid[,1:2,i])
    }
    grid_pts_sp = SpatialPoints(coords = grid_pts)

    test_result = gCovers(p0_buf, grid_pts_sp)
    expect_true(test_result)

    # Compute summary stats in the approximate way

    output1 = discretized_source_approximate_summary_statistics(discrete_source1)
    output2 = discretized_source_summary_statistics(discrete_source1, approx_dx=5000, approx_dy=5000)

    # Check that the results don't differ too much
    output_diff = abs(output1/output2 - 1)
    expect_true(all(abs(output_diff) < 0.2))


    # Do it again with improved orthogonality

    discrete_source2 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=TRUE)
    output2.1 = discretized_source_approximate_summary_statistics(discrete_source2)
    output2.2 = discretized_source_summary_statistics(discrete_source2, approx_dx=5000, approx_dy=5000)
    
    output_diff = abs(output2.1/output2.2 - 1)
    expect_true(all(abs(output_diff) < 0.2))

    # Changes should be small
    expect_true(all(abs(output1/output2.1 - 1) < 0.1 ))
    expect_true(all(abs(output2/output2.2 - 1)  < 0.1))

    #
    # Do a similar test for a more complex case, to confirm that
    # results with/without improved downdip lines are reasonably consistent
    # 

    interface_shapefile2 = 'testshp/alaska.shp'
    desired_subfault_length = 50
    desired_subfault_width = 50

    alaska_source1 = discretized_source_from_source_contours(
        interface_shapefile2, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=FALSE)
    output_a1.1 = discretized_source_approximate_summary_statistics(alaska_source1)
    alaska_source2 = discretized_source_from_source_contours(
        interface_shapefile2, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=TRUE)
    output_a2.1 = discretized_source_approximate_summary_statistics(alaska_source2)

    # Top length can vary quite a bit if we try to allow for orthogonality
    expect_true(all(abs(output_a1.1$length/output_a2.1$length - 1) < 0.25))
    expect_true(all(abs(output_a1.1$width/output_a2.1$width - 1) < 0.10))


    # Check that we can make a discrete_source while providing downdip lines
    new_downdip_lines = downdip_lines_to_SpatialLinesDataFrame(
        discrete_source2$mid_line_with_cutpoints)
    discrete_source3 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, downdip_lines = new_downdip_lines)
    # Check that mid_lines_with_cutpoints is hardly affected
    for(i in 1:length(new_downdip_lines)){
        expect_true(max(abs(
            discrete_source3$mid_line_with_cutpoints[[i]] - 
            discrete_source2$mid_line_with_cutpoints[[i]])) < 1.0e-05)
    }


})

test_that('test_sub_unit_source_grid_point_creation', {

    # Test compute_grid_point_areas_in_polygon

    ############################################################
    # Setup
    #
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

    # The tested routine uses cartesian info
    p1 = get_unit_source_from_discretized_source(ds, c(1,1))
    polygon = p1$unit_source_grid[,1:2]
    polygon = spherical_to_cartesian2d_coordinates(polygon, origin_lonlat=c(0,0))
    full_unit_source_grid = unit_source_coords_cartesian

    # 
    # End preliminary setup
    #########################################################

    #polygon = rbind(c(0, 0), c(0, 10000), c(10000, 10000), c(10000, 0))

    xx = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000)

    c1 = coordinates(xx$grid_point_polygon)
    c2 = coordinates(xx$grid_point_polygon_buffer)

    # Tapered and un-tapered outputs should be the same, but since some trivial
    # geometric operations were involved, we accept round-off error
    expect_true(max(abs(c1-c2)) < 1.0e-10)
    expect_true(all(abs(xx$area - xx$area_buffer) < (xx$area * 1.0e-12)))

    ##################################################################
    #
    # Examine results with edge_taper
    #

    xx2 = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000, edge_taper_width = 1000,
        full_unit_source_grid = full_unit_source_grid)

    # Check that this has not affected the 'untapered' results
    c3 = coordinates(xx2$grid_point_polygon)
    test_result = all(c1 == c3)
    expect_true(test_result)
    
    # Check there are more points in the tapered case    
    c4 = coordinates(xx2$grid_point_polygon_buffer)
    test_result = all(dim(c4)[1] > dim(c3)[1])
    expect_true(test_result)

    # Actually we should have a layer of 1 cell around 2 edges
    expect_true(dim(c4)[1] - dim(c3)[1] == len/1000 + width/1000 + 1)

    # General sanity check
    expect_true(sum(is.nan(xx2$unit_slip_scale)) == 0)

    # The peak unit_source_slip_scale might be slightly different to 1 due to
    # normalisation
    err = 1 - max(xx2$unit_slip_scale)
    expect_true( abs(err) < 5.0e-10)

    # 'Unit Moment' conservation -- sum(slip * down-dip-area) is same in 
    # buffered and unbuffered cases.
    # Since the dip is constant, we expect the following relation to hold
    # (in general we would have to adjust the 'surface areas' to reflect the
    # 'down-dip areas' by multiplying by sqrt(1+tan(dip)**2 -- but dip is not
    # computed by this routine)
    expect_true(
        abs(sum(xx2$unit_slip_scale * xx2$area_buffer) - sum(xx$area)) < 
            (1.0e-12*sum(xx$area)))

    # Ensure that all points are in the 'full' bounding polygon
    bounding_polygon = get_discretized_source_outline(ds)
    bounding_polygon[,1:2] = spherical_to_cartesian2d_coordinates(bounding_polygon[,1:2],
        origin_lonlat=c(0,0))
    
    m1 = point.in.polygon(xx2$grid_points_buffer[,1], xx2$grid_points_buffer[,2],
        bounding_polygon[,1], bounding_polygon[,2])
    expect_true(all(m1 == 1))

    
    ###################################################################
    # Larger edge taper

    xx2 = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000, edge_taper_width = 3000,
        full_unit_source_grid = full_unit_source_grid)

    # Check that this has not affected the 'untaperd' results
    c3 = coordinates(xx2$grid_point_polygon)
    test_result = all(c1 == c3)
    expect_true(test_result)

    # Check there are more points in the tapered case    
    c4 = coordinates(xx2$grid_point_polygon_buffer)
    test_result = all(dim(c4)[1] > dim(c3)[1])
    expect_true(test_result)

    # Actually we should have a layer of 3 cells around 2 edges, except
    # the top corner is missing (doesn't make it inside the buffer)
    expect_true(dim(c4)[1] - dim(c3)[1] == (3*(len/1000 + width/1000 + 3) -1))

    expect_true(sum(is.nan(xx2$unit_slip_scale)) == 0)

    # The peak unit_source_slip_scale might be slightly different to 1 due to
    # normalisation
    err = 1 - max(xx2$unit_slip_scale)
    expect_true( abs(err) < 5.0e-10)

    # 'Unit Moment' conservation -- sum(slip * down-dip-area) is same in 
    # buffered and unbuffered cases.
    # Since the dip is constant, we expect the following relation to hold
    # (in general we would have to adjust the 'surface areas' to reflect the
    # 'down-dip areas' by multiplying by sqrt(1+tan(dip)**2 -- but dip is not
    # computed by this routine)
    expect_true(
        abs(sum(xx2$unit_slip_scale * xx2$area_buffer) - sum(xx$area)) < 
            (1.0e-12*sum(xx$area)))

    # Ensure that all points are in the 'full' bounding polygon
    bounding_polygon = get_discretized_source_outline(ds)
    bounding_polygon[,1:2] = spherical_to_cartesian2d_coordinates(bounding_polygon[,1:2],
        origin_lonlat=c(0,0))
    
    m1 = point.in.polygon(xx2$grid_points_buffer[,1], xx2$grid_points_buffer[,2],
        bounding_polygon[,1], bounding_polygon[,2])
    expect_true(all(m1 == 1))


})
