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

})

