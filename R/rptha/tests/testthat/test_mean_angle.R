
context('test_mean_angle')

test_that('test_mean_angle', {
#test_mean_angle<-function(){

    # By this method the mean of -1,0,90 is < 30 
    m1 = mean_angle(c(0,0,90))

    err = abs(m1 - 26.5650511771)
    
    expect_that(err < 1.0e-08, is_true())

    # A more typical case
    m2 = mean_angle(c(120, 60))
    err = abs(m2 - 90)

    expect_that(err < sqrt(.Machine$double.eps), is_true())

    # Check degrees / radians

    a3 = runif(100, -2, 2)*pi
    a4 = a3/pi*180
    m3 = mean_angle(a3, degrees=FALSE)
    m4 = mean_angle(a4, degrees=TRUE)

    err = abs(m3*180/pi - m4)

    expect_that(err < 1.0e-06, is_true())
   
    # Check there is no impact of representing negative angles as positive 
    a5 = a3 + (a3 < 0)*2*pi
    m5 = mean_angle(a5, degrees=FALSE)

    err = abs(m3 - m5)
    
    expect_that(err < 1.0e-06, is_true())
#}
})
