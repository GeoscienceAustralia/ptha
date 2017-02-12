context('test_intersect_surface_path_with_depth_contours')

test_that('test_intersect_surface_path_with_depth_contours', {
#test_intersect_surface_path_with_depth_contours<-function(){

    contour1 = readOGR('testshp/alaska.shp', layer='alaska')

    line1 = matrix(c(210, 57, 208, 59.5), ncol = 2, byrow=TRUE)

    line1_3D = intersect_surface_path_with_depth_contours(line1, contour1)

    # Check that we get nearly the same result with higher 'n'
    line2_3D = intersect_surface_path_with_depth_contours(line1, contour1, n=200)    

    expect_that(all(abs(line1_3D - line2_3D) < 1.0e-03), is_true())

    # Check what happens if we reverse line1_3D
    # Should get the same result
    line3_3D = intersect_surface_path_with_depth_contours(line1[2:1,], contour1)

    expect_that(all(abs(line1_3D - line3_3D) == 0.0), is_true())

    ## This case was troublesome originally 
    ## The intersection point with the shallowest contour was missed
    contour2 = readOGR('testshp/sagami.shp', layer='sagami')

    lc = length(contour2) 
    l1 = length(contour2@lines[[1]]@Lines[[1]]@coords[,1] )
    l2 = 1 # Note the last line orientation in this file is currently reversed (24/07/2015)
    end_points = rbind(contour2@lines[[1]]@Lines[[1]]@coords[l1,],
                       contour2@lines[[lc]]@Lines[[1]]@coords[l2,])

    line2_3D = intersect_surface_path_with_depth_contours(end_points, contour2, 
        extend_line_fraction=0.1)

    # Originally end_points[1,] was missing from line2_3D
    expect_that(
        (min(abs(end_points[1,1] - line2_3D[,1]) + 
             abs(end_points[1,2] - line2_3D[,2])) < 1.0e-02),
        is_true())

})
