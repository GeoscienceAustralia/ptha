
context('test_okada_tsunami')

test_that("test_okada_tsunami", {
    # Test against Okada's (1985) table

    library(rptha)


    ## Case 2 in Okada's table
    ## Note that Okada's coordinate system differs from ours
    # x axis --> negative okada y axis; y axis --> x axis

    okada_x_origin = 0.0
    okada_y_origin = 0.0
    okada_d = 4
    okada_x = -3 # 2 in okada's frame
    okada_y = 2  # 3 in okada's frame 

    strike=0.
    dip = 70.
    L = 3.
    W = 2.
    slip=1.0

    # Compute centroid of rupture in our reference frame (not the same as
    # Okada's)
    deg2rad = pi/180
    centroid_y = L/2 * 1000
    centroid_x = -(W/2)*cos(dip*deg2rad) * 1000
    centroid_d = okada_d - (W/2)*sin(dip*deg2rad)

    ## Dip slip case ##
    ans = okada_tsunami(centroid_x, centroid_y, centroid_d, strike, dip, L, W,
        0, slip, okada_x*1000, okada_y*1000)
    # Note that x axis --> negative okada y axis, y axis --> x axis
    okada_ans = c(3.527e-02, -4.682e-03, -3.564e-02)

    expect_that(all(abs(unlist(ans) - okada_ans) < 1.0e-05), is_true())

    ## Strike slip case ##
    ans = okada_tsunami(centroid_x, centroid_y, centroid_d, strike, dip, L, W,
        slip, 0, okada_x*1000, okada_y*1000)
    # Note that x axis --> negative okada y axis, y axis --> x axis
    okada_ans = c(4.298e-03, -8.689e-03, -2.747e-03)

    expect_that(all(abs(unlist(ans) - okada_ans) < 1.0e-05), is_true())
})
