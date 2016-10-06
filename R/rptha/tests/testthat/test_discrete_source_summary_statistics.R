
context('test_discrete_source_summary_statistics')

test_that('test_discrete_source_summary_statistics', {

    interface_shapefile = 'testshp/sagami.shp'

    desired_subfault_length = 100
    desired_subfault_width = 50
    
    # Create sagami source using 'old' discretization method 
    discrete_source1 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=FALSE)

    # Compute summary stats in the approximate way

    output1 = discretized_source_approximate_summary_statistics(discrete_source1)
    output2 = discretized_source_summary_statistics(discrete_source1)

    # Check that the results don't differ too much
    output_diff = abs(output1/output2 - 1)
    expect_that(all(abs(output_diff) < 0.2), is_true())


    # Do it again with improved orthogonality

    discrete_source2 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE, improved_downdip_lines=TRUE)
    output2.1 = discretized_source_approximate_summary_statistics(discrete_source2)
    output2.2 = discretized_source_summary_statistics(discrete_source2)
    
    output_diff = abs(output2.1/output2.2 - 1)
    expect_that(all(abs(output_diff) < 0.2), is_true())

    # Changes should be small
    expect_that(all(abs(output1/output2.1 - 1) < 0.1 ), is_true())
    expect_that(all(abs(output2/output2.2 - 1)  < 0.1), is_true())

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
    expect_that(all(abs(output_a1.1$length/output_a2.1$length - 1) < 0.3), is_true())
    expect_that(all(abs(output_a1.1$width/output_a2.1$width - 1) < 0.07), is_true())

})
