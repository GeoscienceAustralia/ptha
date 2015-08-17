context('test_adjust_longitude_by_360_deg')

test_that("test_adjust_longitude_by_360_deg", {
#test_adjust_longitude_by_360_deg<-function(){

    # Test 1 -- good change
    p0 = c(357, 5)
    refpt = c(2, -10)
   
    p1 = adjust_longitude_by_360_deg(p0, refpt)

    expect_that(all(p1 == p0 - c(360, 0)), is_true())
   
    # Test 2 -- good change 
    p0 = c(-154, 5)
    refpt = c(260, -10)
    
    p1 = adjust_longitude_by_360_deg(p0, refpt)
  
    expect_that(all(p1 == p0 + c(360, 0)), is_true())

    # Test 3 -- good change 
    p0 = c(260, -10)
    refpt = c(-154, 5)
    
    p1 = adjust_longitude_by_360_deg(p0, refpt)

    expect_that(all(!is.null(p1)), is_true())

    expect_that(all(p1 == p0 - c(360, 0)), is_true())

    # Test 4 -- good change 
    p0 = c(260 + 3*360, -10)
    refpt = c(-154, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)

    expect_that(all(!is.null(p11)), is_true())

    expect_that(all(p11 == p1), is_true())

    # Test 5 -- good change 
    p0 = c(260 - 6*360, -10)
    refpt = c(-154, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)

    expect_that(all(!is.null(p11)), is_true())

    expect_that(all(p11 == p1), is_true())

    # Test 6 -- no change
    p0 = c(260 , -10)
    refpt = c(130, 5)
    
    p11 = adjust_longitude_by_360_deg(p0, refpt)

    expect_that(all(!is.null(p11)), is_true())
    expect_that(all(p11 == p0), is_true())
#}
})
