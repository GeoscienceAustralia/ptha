context('test_kajiura_filter')

test_that('test_kajiura_g', {

    # Test that kajiura_g is correctly coded 
    testr = 4 
    # (1/2pi)* (The integral of this) is another representation of the G function 
    f<-function(m) m*besselJ(m*testr,nu=0)/cosh(m) 
    fInt = integrate(f, 0, Inf) 
    fInt = fInt$value/(2*pi) 
    ref_val = kajiura_g(testr) 
    out = abs(fInt-ref_val) 
    #print('Compare theoretical methods 1 and 2')
    #print('Difference is: ')
    #print(out)
    #print('Differences < 1.0e-03 are ok')

    expect_that( out < 1.0e-03, is_true())

    #print('Compare theoretical and empirical')

    kg_E = kajiura_g_empirical()
    # Add small perturbation to x point so that we don't evaluate at spline knot
    e_ref = kg_E(testr+0.001)
    ref_val = kajiura_g(testr+0.001)
    out = abs(ref_val-e_ref) 
    #print('Difference is: ')
    #print(out)
    #print('Differences < 1.0e-06 are ok')
    expect_that(out < 1.0e-06, is_true())

    x = seq(0,8,len=50)
    kgV = Vectorize(kajiura_g)
    theory = kgV(x)
    splineApprox = kg_E(x)
    par(mfrow=c(2,2))
    plot(x,theory,t='o',log='y',main='Kajuira function (log y-scale)')
    points(x,splineApprox,t='l',col=2)
    legend('topright',c('Theory', 'Splines'), lty=c(1,1),col=c(1,2)) 

    plot(theory-splineApprox, main='ABS Difference',t='l')
    plot((theory-splineApprox)/theory, main='REL Difference',t='l')

})



test_that('test_kajiura_filter_qualitative', {

    #
    # Smoothing of step function
    # Qualitative check only
    #

    lScale = 50000
    n = 100
    pts = as.matrix(expand.grid(seq(0,1,len=n), seq(0,2,len=n)))*lScale
    pts_z = pts[,1] > (lScale/2.)
    xyDef = cbind(pts,pts_z)
    
    # Make the depth constant at the top / bottom, and otherwise
    # linearly varying
    depth = 200 + pmin( pmax(xyDef[,2]-0.2*lScale, 0), 1.6*lScale)/100

    new_xyDef = kajiura_filter(xyDef, depth, grid_dx=400,grid_dy=300)
    m1 = matrix(new_xyDef[,3], ncol=n)

    new_xyDef = kajiura_filter(xyDef, depth, grid_dx=200,grid_dy=200)
    m2 = matrix(new_xyDef[,3], ncol=n)
    
    depthMat = matrix(depth,ncol=n)

    # Put a check to see if they are the 'same' away from edges
    testStat=max(abs(m1[,20:80]-m2[,20:80]))

    expect_that(testStat <= 0.1, is_true())

})


test_that('test_kajiura_filter_1', {
    
    # Check that the smoothing of a 'point' source with constant depth
    # leads to a bump with the shape of the Kajiura G function.
    # Use an uneven point cloud.
    #
    # NOTE: The volume of the 'point' source is not well defined.
    #       Since the input xyDef data is interpreted as unstructured points,
    #        the volume of a 'point' source depends on how you interpolate the data.
    #       So, we account for this in the test.
    #
    # Regardless, the repeated re-gridding in the routine will still introduce some errors.

    # Define initial xyDef + depth
    xrange = seq(-1,1,len=201)*50000
    yrange = seq(-1,1,len=401)*50000

    xyGrid = as.matrix(expand.grid(xrange,yrange))

    mm = which((xyGrid[,1]==0) & (xyGrid[,2]==0))
    if(length(mm)!=1) stop('Did not find origin')

    # Make deformation mostly zero, except for a 'spike' of 1 at 0,0
    depths = xyGrid[,1]*0+2000
    def = xyGrid[,1]*0
    def[mm] = 1.0

    xyDef = cbind(as.matrix(xyGrid), def)

    # Interpolate / smooth
    # NOTE: We need to ensure grid_dx, grid_dy are small enough to 'see'
    # the input data correctly. 
    new_xyDef = kajiura_filter(xyDef, depths, grid_dx=500,grid_dy=250)

    # Unsmoothed
    m1 = matrix(xyDef[,3],ncol=201,byrow=T)
    # Smoothed
    m2 = matrix(new_xyDef[,3],ncol=201,byrow=T)
    #
    d1 = matrix(depths, ncol=201,byrow=T)

    # Find the 'r' value to apply to kajiura G
    rMax = matrix(sqrt(xyGrid[,1]^2+xyGrid[,2]^2), ncol=201,byrow=T)/d1

    # Make kajiura function
    kgE = kajiura_g_empirical()

    # Limit the bounds of the kajiura function inputs
    rMax = pmin(rMax,8)

    k2 = matrix(kgE(rMax)*(rMax<8),ncol=201)
    k2 = k2/sum(k2)*sum(m2) # Normalise to sum, to get around the 'volume' interpretation issue for a point source

    testStat = max(abs(m2-k2))

    expect_that(testStat < 1.0e-05, is_true())
})


test_that('test_kajiura_filter_3', {

    # Check that the smoothing of a 'point' source    
    # is as desired, with an uneven point cloud.
    # Uses variable depths. 
    # Recall that in this situation, the implementation of the kajiura filter
    # is not really well defined (since the underlying theory applies to
    # constant depths). However, for slowly varying depths, application of the
    # filter with locally varying filter width should be quite accurate.
    #
    # NOTE: The volume of the 'point' source is not well defined.
    #       Since the input xyDef data is interpreted as unstructured points,
    #        the volume of a 'point' source depends on how you interpolate the data
    #       So, we account for this in the test, by ensuring sum(k2)=sum(m2)
    #
    # However, the repeated re-gridding in the routine will still introduce some errors

    # Define initial xyDef + depth
    xrange = seq(-1,1,len=201)*50000
    yrange = seq(-1,1,len=401)*50000

    xyGrid=as.matrix(expand.grid(xrange,yrange))

    mm = which((xyGrid[,1]==0) & (xyGrid[,2]==0))
    if(length(mm) != 1) stop('Did not find origin')

    # Make deformation mostly zero, except for a 'spike' of 1 at 0,0
    def = xyGrid[,1]*0
    def[mm] = 1.0
   
    # Make 'varying' depths 
    # If they vary too rapidly, the volume error becomes too large (since our kajiura filter)
    # is not volume conservative with spatially varying depths.
    depths = 2000.+100.*(sin(xyGrid[,1]/3000)+cos(xyGrid[,2]/3000*0.8))
    xyDef = cbind(xyGrid, def)

    # Interpolate / smooth
    # NOTE: We need to ensure grid_dx, grid_dy are small enough to 'see'
    # the input data correctly. 
    new_xyDef = kajiura_filter(xyDef, depths, grid_dx=500,grid_dy=250)

    # Unsmoothed
    m1 = matrix(xyDef[,3],ncol=201,byrow=T)
    # Smoothed
    m2 = matrix(new_xyDef[,3],ncol=201,byrow=T)
    #
    #
    # Find the 'r' value to apply to kajiura G
    rMax = sqrt(xyGrid[,1]^2+xyGrid[,2]^2)/depths
    # Limit the bounds of the kajiura function inputs
    rMax = pmin(rMax,8)

    # Make kajiura function
    kgE = kajiura_g_empirical()

    k2 = matrix(kgE(rMax)*(rMax<8), ncol=201,byrow=T)
    k2 = k2/max(k2)*max(m2) # Normalise to sum, to get around the 'volume' interpretation issue for a point source
    

    # Test = difference in surfaces
    testStat = max(abs(m2-k2))

    # Because both our 'analytical' approach in this test and numerical
    # approach can both only be considered approximate with varying depths, we
    # don't expect perfect agreement here.
    expect_that(testStat < 1.0e-04, is_true())

    #
    # Now, let's do it using the raster kajiura code, and check that results
    # are 'the same'
    #

    elevation_rast = rptha:::.local_rasterFromXYZ(cbind(xyDef[,1:2], -depths))
    deformation_rast = rptha:::.local_rasterFromXYZ(xyDef)

    deformation_kj = kajiura_smooth_raster(
        deformation_rast,
        new_origin=c(0,0),
        elevation_raster=elevation_rast,
        kj_filter_grid_dxdy = c(500, 250),
        kj_filter_def_threshold = 0,
        spherical_input=FALSE)

    # Should be 'the same' as m2
    deformation_kj_mat = as.matrix(deformation_kj)
    testStat = max(abs(deformation_kj_mat - m2))
    expect_that( testStat < 1.0e-15, is_true())

})
