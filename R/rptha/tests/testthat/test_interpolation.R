context('test_interpolation')

test_that('test_interpolation', {

    # Make a single triangle in the plane z=x+y, and interpolate from it
    xy = matrix(c(0,0,0,1,1,1),ncol=2,byrow=T)
    vals = c(0, 1, 2) # z=x+y
    newPts = matrix(c(0.5, 0.5, 0.3, 0.3), ncol=2, byrow=T)

    out = triangular_interpolation(xy, vals, newPts)

    expect_that(all(isTRUE(all.equal(out, c(1.0,0.6)))), is_true())

    # Re-order triangle
    xy = xy[3:1,]
    vals = vals[3:1]
    out = triangular_interpolation(xy, vals, newPts)
    
    expect_that(all(isTRUE(all.equal(out,c(1.0,0.6)))), is_true())

    #another one, with formula z=0.5*x+0.2*y+7
    xy = matrix(c(-1, -1, 1, -0.5, 0.5, 1), ncol=2,byrow=2)
    vals = 0.5*xy[,1]+0.2*xy[,2]+7
    newPts = matrix(c(0,0, 0.5, 0.3),ncol=2,byrow=T)
    expectedVals = 0.5*newPts[,1]+0.2*newPts[,2]+7
    out=triangular_interpolation(xy,vals,newPts)

    expect_that(all(isTRUE(all.equal(out,expectedVals))), is_true())

    # A point outside the triangle 
    newPts = matrix(c(-1,0, -1, 1),ncol=2,byrow=T)
    out = triangular_interpolation(xy,vals,newPts)

    expect_that((all(is.na(out))| formals(triangular_interpolation)$useNearestNeighbour), is_true())

    # A single point
    newPts = matrix(c(0,0),ncol=2)
    out = triangular_interpolation(xy,vals,newPts)

    expect_that(out == 7, is_true())

    # Points on the triangle
    newPts = xy
    out = triangular_interpolation(xy,vals,newPts)
    
    expect_that(all(out == vals), is_true())

    # Point on an edge
    newPts = matrix(0.5*(xy[1,]+xy[2,]),ncol=2)
    out = triangular_interpolation(xy,vals,newPts)

    expect_that(all(out==0.5*(vals[1]+vals[2])), is_true())


    # Test  discontinuous interpolation
    for(interpolator in list(triangular_interpolation, nearest_neighbour_interpolation)){

        xy = matrix(runif(200), ncol=2)
        vals = (xy[,1] > 0.5) + (xy[,2] > 0.5)

        newPts = matrix(runif(3000), ncol=2)
        # Make 4 categories
        category_function<-function(xy){
            (xy[,1] > 0.5) + 10*(xy[,2] > 0.5)
        }

        newPts_vals = interpolation_discontinuous(xy, vals, newPts, category_function, 
            interpolator=interpolator)

        theoretical_answer = (newPts[,1] > 0.5) + (newPts[,2] > 0.5)

        expect_that(all(newPts_vals == theoretical_answer), is_true())


        ##  Test with a matrix of input arguments
        vals = cbind(vals, 1-vals)
        newPts_vals = interpolation_discontinuous(xy, vals, newPts, category_function, 
            interpolator=interpolator)

        theoretical_answer = cbind(theoretical_answer, 1 - theoretical_answer)
        expect_that(all(newPts_vals == theoretical_answer), is_true())

    }

})

test_that('test_contour_interpolation', {

    contours = readOGR('testshp/alaska.shp', layer='alaska')
    dsc = discretized_source_from_source_contours('testshp/alaska.shp', 50, 50, 
        improved_downdip_lines=TRUE)

    # Test some 'general' points with previously confirmed reasonable values
    p1 = rbind(c(210, 58), c(212, 60), c(208, 57))
    expected_p1_depths = c(14.669, 27.898, 12.471) + 6 # These contours are offset by 6km from trench=0
    contour_fun = make_contour_interpolator(dsc$mid_line_with_cutpoints, 
        convert_to_cartesian=TRUE, origin_lonlat = c(208,58))
    p1_depths = contour_fun(p1)
    expect_true(max(abs(p1_depths - expected_p1_depths)) < 1.0e-02)

    ## Test points exactly on the contours
    mid_line_first = dsc$mid_line_with_cutpoints[[1]]
    expected_depths = mid_line_first[,3]
    interpolated_z = contour_fun(mid_line_first[,1:2])
    expect_true(max(abs(expected_depths - interpolated_z)) < 1.0e-05)

})
