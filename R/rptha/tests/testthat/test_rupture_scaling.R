
context('test_rupture_scaling')

test_that('test_rupture_scaling', {
    
    # Example quoted in Bird et al (2009) pg 3102, 2nd last paragraph
    M01 = M0_2_Mw(5.66, inverse=TRUE)

    expect_that(signif(M01, 3) == 3.47e+17, is_true())

    # Check inverse relationship holds
    Mw1 = M0_2_Mw(M01)
    expect_that(signif(Mw1, 3) == 5.66, is_true())

    # Check independetly that Mw can be computed
    Mw1b = M0_2_Mw(3.47e+17)
    expect_that(signif(Mw1b,3) == 5.66, is_true())


    # Test Mw to rupture size 
    output_simple = Mw_2_rupture_size(9.0, relation='Strasser')
    output_complex = Mw_2_rupture_size(9.0, relation='Strasser', detailed=TRUE)

    expect_that(all(output_simple == output_complex$values), is_true())

    expect_that(all(output_complex$log10_sigmas == c(0.304, 0.173, 0.18)), is_true()) 

    expect_that(all(round(output_simple,0) == c(123595, 189, 614)), is_true())


    # Test slip from Mw area mu

    slip0 = slip_from_Mw_area_mu(9.0, output_simple[1], mu=3e+10)
    expect_that(round(slip0,0) == 10, is_true())

    slip1 = slip_from_Mw_area_mu(9.0, output_simple[1], mu=4e+10)

    expect_that( isTRUE(all.equal(slip1, slip0*3/4)), is_true())

    # Test slip_from_Mw
    slip2 = slip_from_Mw(9.0)
    # remove name from slip0 for equality test
    names(slip0) = NULL
    expect_that(isTRUE(all.equal(slip0, slip2)), is_true())

})
