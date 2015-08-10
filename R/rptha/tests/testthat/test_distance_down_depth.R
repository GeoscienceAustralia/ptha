test_distance_down_depth<-function(){

    ## TEST 1 ## Surface distances

    p1 = c(200, 50, 0)
    p2 = c(170, 10, 0)

    dist1 = distance_down_depth(p2, p1, n=1e+04)

    dist_cos = distCosine(p2[1:2], p1[1:2])

    expect_that(all.equal(dist1, dist_cos), is_true())

    ## TEST 2 ## Unaffected by + 360 to lat

    p1 = c(200, 50, 0)
    p2 = c(170, 10, 50)

    dist2 = distance_down_depth(p2, p1)

    p1 = c(200 - 360, 50, 0)
    p2 = c(170, 10, 50)

    dist3 = distance_down_depth(p2, p1)

    expect_that(isTRUE(all.equal(dist2 , dist3)), is_true()) 
    
    ## TEST 3 ##
    dist4 = distance_down_depth(p1, p2) # reverse arguments

    expect_that(isTRUE(all.equal(dist3 , dist4)), is_true())

    ## TEST 4 ## Unaffected by reflecting lat
    p1 = c(200, -50, 0)
    p2 = c(170, -10, 50)

    dist5 = distance_down_depth(p2, p1)

    expect_that(isTRUE(all.equal(dist5 , dist3)), is_true())

    ## TEST 5 ##  Analytical solution for a dipping fault in 1D
    ##
    ## Suppose we integrate over a surface S dipping below a land surface (1D
    ## circle with radius R), with the S's depth increasing at a rate ddepth/dtheta = constant
    ## The coordinate can be parameterised by an angular coordinate 'theta'
    ##
    ## At theta=0, depth=0
    ## At any theta, ddepth = dtheta * (ddepth/dtheta) [since ddepth/dtheta = constant]
    ## and the projection of the dipping surface parallel to the land surface is:
    ##               inner_length_element = (R - depth)*dtheta = (R - theta*ddepth/dtheta)*dtheta 
    ## Therefore, 
    ## downdip_length_element**2 = ddepth**2 + inner_length_element**2 
    ##     = dtheta^2 * ( (ddepth/dtheta)^2 + [R^2 - theta*ddepth/dtheta]^2)
    ## so
    ## d(downdip_length)/dtheta = sqrt ( ddepth/dtheta^2 + [R - theta*ddepth/dtheta]^2] )
    ## we can integrate this to get the down-dip length
    ##
    ## (NOTE: I think this can be analytically integrated too, e.g.)
    ## sage: var('R, p')
    ## sage: assume(R>0)
    ## sage: assume(p>0)
    ## sage: assume(R>p)
    ## sage: f = integral( ( (R - x*p)^2 + p^2)^0.5, x)
    ## sage: f.simplify_full()
    ## 1/2*(p^2*arcsinh((p*x - R)/p) + sqrt(p^2*x^2 - 2*R*p*x + R^2 + p^2)*(p*x - R))/p

    # Example:
    # Drop 50km in 1 degree
    p1 = c(0, 0, 0)
    p2 = c(1, 0, 50)

    # Function to integrate (theta in radians)
    f<-function(theta, ddepth_dtheta = 50/(2*pi/360), R = 6378.137){
        sqrt( ddepth_dtheta**2 + (R - theta*ddepth_dtheta)^2 )
    }
    l_analytical = integrate(f, 0, 2*pi/360)

    # Numerically compute
    l_numerical = distance_down_depth(p1, p2)/1000

    expect_that(isTRUE(all.equal(l_analytical$value, l_numerical)), is_true())

}

