test_interpolate_gc_path<-function(){

    ## Test 1
    surface_path = matrix(c(200, 5, 170, -5), ncol=2, byrow=T)

    interp_path1 = interpolate_gc_path(surface_path)

    expect_that(max(abs(diff(interp_path1[,1]))< 20), is_true())

    ## Test 2
    surface_path = surface_path[2:1,]
    interp_path2 = interpolate_gc_path(surface_path)
  
    # Check there is no large spacing 
    expect_that(max(abs(diff(interp_path2[,1])) < 2), is_true())

    expect_that(
        any((abs(rev(interp_path2[,1]) - interp_path1[,1]) > 1.0e-03) |
            (abs(rev(interp_path2[,2]) - interp_path1[,2]) > 1.0e-03) ),
        is_false())

    ## Test 3
    surface_path = matrix(c(-160, 5, 170, -5), ncol=2, byrow=T)

    interp_path1 = interpolate_gc_path(surface_path)

    expect_that(max(abs(diff(interp_path1[,1])) < 20), is_true())


    ## Test 4
    surface_path = surface_path[2:1,]
    interp_path2 = interpolate_gc_path(surface_path)
  
    # Check there is not any large spacing 
    expect_that(max(abs(diff(interp_path2[,1]))< 2), is_true())

    # Check these paths are the 'same'
    # (except now the longitudes are offset by 360
    expect_that(
        (any((abs(rev(interp_path2[,1]) - interp_path1[,1])%%360 > 1.0e-03) |
             (abs(rev(interp_path2[,2]) - interp_path1[,2]) > 1.0e-03) )),
        is_false())
}

