
context('test_discrete_source_summary_statistics')

test_that('test_discrete_source_summary_statistics', {

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
    expect_true(all(abs(output_a1.1$length/output_a2.1$length - 1) < 0.5))
    expect_true(all(abs(output_a1.1$width/output_a2.1$width - 1) < 0.10))

})

test_that('test_sub_unit_source_grid_point_creation', {

    # Test compute_grid_point_areas_in_polygon

    polygon = rbind(c(0, 0), c(0, 10000), c(10000, 10000), c(10000, 0))

    xx = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000)

    c1 = coordinates(xx$grid_point_polygon)
    c2 = coordinates(xx$grid_point_polygon_buffer)
    test_result = all(c1 == c2)
    expect_true(test_result)

    l1 = length(xx$unit_slip_scale)
    l2 = length(c2[,1])
    expect_true(l1 == l2)

    expect_true(all(xx$area == xx$area_buffer))

    ##################################################################
    #
    # Examine results with edge_taper
    #

    xx2 = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000, edge_taper_width = 1000)

    # Check that this has not affected the 'untaperd' results
    c3 = coordinates(xx2$grid_point_polygon)
    test_result = all(c1 == c3)
    expect_true(test_result)

    expect_true(sum(is.nan(xx2$unit_slip_scale)) == 0)

    # The peak unit_source_slip_scale might be slightly different to 1 due to
    # normalisation
    err = 1 - max(xx2$unit_slip_scale)
    expect_true( abs(err) < 5.0e-03)

    ##################################################################
    #
    # Examine results with edge_taper and bounding_polygon
    #
    bounding_polygon = rbind(c(0, 0), c(0, 10000), c(0, 20000), 
        c(20000, 20000), c(20000, 0))
    xx3 = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000, edge_taper_width = 1000,
        bounding_polygon = bounding_polygon)

    # Ensure that the bounding polygon did clip as expected
    m1 = point.in.polygon(xx3$grid_points_buffer[,1], xx3$grid_points_buffer[,2],
        bounding_polygon[,1], bounding_polygon[,2])
    expect_true(all(m1 == 1))

    
    ###################################################################
    # Larger edge taper

    xx2 = compute_grid_point_areas_in_polygon(polygon, 
        approx_dx=1000, approx_dy=1000, edge_taper_width = 3000)

    # Check that this has not affected the 'untaperd' results
    c3 = coordinates(xx2$grid_point_polygon)
    test_result = all(c1 == c3)
    expect_true(test_result)

    expect_true(sum(is.nan(xx2$unit_slip_scale)) == 0)



})
