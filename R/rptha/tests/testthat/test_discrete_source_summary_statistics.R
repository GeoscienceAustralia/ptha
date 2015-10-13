
context('test_discrete_source_summary_statistics')

test_that('test_discrete_source_summary_statistics', {

    interface_shapefile = 'testshp/sagami.shp'

    desired_subfault_length = 100
    desired_subfault_width = 50
     
    discrete_source1 = discretized_source_from_source_contours(
        interface_shapefile, desired_subfault_length, desired_subfault_width,
        make_plot=FALSE)

    # Compute summary stats in the approximate way

    output1 = discretized_source_approximate_summary_statistics(discrete_source1)
    output2 = discretized_source_summary_statistics(discrete_source1)

    # Check that the results don't differ too much
    output_diff = abs(output1/output2 - 1)
    expect_that(all(abs(output_diff) < 0.2), is_true())

})
