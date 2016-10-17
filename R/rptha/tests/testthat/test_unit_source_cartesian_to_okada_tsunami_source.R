context('test_unit_source_cartesian_to_okada_tsunami_source')

test_that('test_unit_source_cartesian_to_okada_tsunami_source', {
#.test_unit_source_cartesian_to_okada_tsunami_source<-function(){

    # Here we copy code from a test of unit source interior points cartesian to
    # set up the problem
    d0 = 6000
    d1 = 10000
    width = 20000
    len = 40000
    dip = atan((d1 - d0)/width)*180/pi
    strike = 270    

    unit_source_coords_cartesian=rbind(c(0, 0, d0),
                                       c(0, width, d1),
                                       c(-len, width, d1),
                                       c(-len, 0, d0))

    # Back-calculate lon/lat
    origin = c(0,0)
    unit_source_coords_lonlat = cartesian2d_to_spherical_coordinates(
        unit_source_coords_cartesian[,1:2], origin=origin)

    # Make discrete source consisting of a single unit source
    ds = list()
    ds$unit_source_grid = array(dim=c(2,3,2))
    ds$unit_source_grid[,,1] = cbind(unit_source_coords_lonlat[1:2,1:2], 
        unit_source_coords_cartesian[1:2,3]/1000)
    ds$unit_source_grid[,,2] = cbind(unit_source_coords_lonlat[4:3,1:2], 
        unit_source_coords_cartesian[4:3,3]/1000)
   
    ds$discretized_source_dim = c(1,1) 
    names(ds$discretized_source_dim) = c('dip','strike')

    mid_line_with_cutpoints = list()
    mid_line_with_cutpoints[[1]] = cbind(unit_source_coords_lonlat[1:2,], c(d0, d1)/1000)
    mid_line_with_cutpoints[[2]] = cbind(unit_source_coords_lonlat[4:3,], c(d0, d1)/1000)
    ds$mid_line_with_cutpoints = mid_line_with_cutpoints

    ds$fine_downdip_transects = ds$unit_source_grid    

    # Get tsunami surface points
    tsunami_surface_points_lonlat = expand.grid(seq(-1,1,len=200), 
        seq(-1,1,len=200)) 

    # Make tsunami with our routines
    tsunami1 = make_tsunami_unit_source(1, 1, ds, rake=90,
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy=3000,
        tsunami_function = unit_source_cartesian_to_okada_tsunami_source)

    # Make tsunami directly with Okada routine
    tsunami_surface_points_cartesian = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat = origin)

    # [Remember to correct width to be down-dip]
    wc = sqrt(1 + atan(dip/180*pi)**2) # width-correction for down-dip
    tsunami2 = okada_tsunami(elon = -len/2, elat = width/2, 
        edep = 0.5*(d0+d1)/1000,
        strk = strike, dip = dip, lnth = len/1000, wdt = width/1000*wc,
        disl1 = 0, disl2 = 1, rlon = tsunami_surface_points_cartesian[,1],
        rlat = tsunami_surface_points_cartesian[,2])

    # Check both displacements are sufficiently close (can get closer with more
    # interior points)
    expect_true(max(abs(range(tsunami1$tsunami_source$zdsp - tsunami2$zdsp))) < 1.0e-02)

    # Check that edge_taper_width > 0 runs
    tsunami3 = make_tsunami_unit_source(1, 1, ds, rake=90,
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy=3000,
        tsunami_function = unit_source_cartesian_to_okada_tsunami_source,
        edge_taper_width = 6000,
        allow_points_outside_discrete_source_outline=TRUE)

    # Simple check -- edge_tapering decreases the extremes of the displacement
    r1 = diff(range(tsunami1$smooth_tsunami_displacement))
    r3 = diff(range(tsunami3$smooth_tsunami_displacement))
    
    expect_true(r1 > r3)
    
    tsunami4 = make_tsunami_unit_source(1, 1, ds, rake=90,
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy=3000,
        tsunami_function = unit_source_cartesian_to_okada_tsunami_source,
        edge_taper_width = 9995, # 10000 will induce an error
        allow_points_outside_discrete_source_outline=TRUE)
    
    r4 = diff(range(tsunami4$smooth_tsunami_displacement))
    expect_true(r3 > r4)
    
    #
    # Check that edge tapering only effects the edges of sums of unit sources
    #
    sagami = readOGR('testshp/sagami.shp', 'sagami')
    sagami_source = discretized_source_from_source_contours('testshp/sagami.shp',
        desired_subfault_length = 50, desired_subfault_width=50, make_plot=FALSE)
    tsunami_surface_points_lonlat = expand.grid(seq(138,144,len=200), 
        seq(33,37,len=200)) 
    
    # Case 1 -- no edge tapering
    tsunami11 = make_tsunami_unit_source(1,1,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000)
    tsunami21 = make_tsunami_unit_source(2,1,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000)
    tsunami12 = make_tsunami_unit_source(1,2,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000)
    tsunami22 = make_tsunami_unit_source(2,2,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000)

    m11 = tsunami_unit_source_2_raster(tsunami11)
    m12 = tsunami_unit_source_2_raster(tsunami12)
    m21 = tsunami_unit_source_2_raster(tsunami21)
    m22 = tsunami_unit_source_2_raster(tsunami22)

    sum1 = m11 + m12
    sum_2 = m11 + m21
    sum_all = m11 + m12 + m21 + m22
  
    # Case 2 -- edge tapering 
    tsunami11 = make_tsunami_unit_source(1,1,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000,
        edge_taper_width=10000, allow_points_outside_discrete_source_outline=TRUE)
    tsunami21 = make_tsunami_unit_source(2,1,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000,
        edge_taper_width=10000, allow_points_outside_discrete_source_outline=TRUE)
    tsunami12 = make_tsunami_unit_source(1,2,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000,
        edge_taper_width=10000, allow_points_outside_discrete_source_outline=TRUE)
    tsunami22 = make_tsunami_unit_source(2,2,sagami_source, rake=90, 
        tsunami_surface_points_lonlat, approx_dx = 3000, approx_dy = 3000,
        edge_taper_width=10000, allow_points_outside_discrete_source_outline=TRUE)

    m11B = tsunami_unit_source_2_raster(tsunami11)
    m12B = tsunami_unit_source_2_raster(tsunami12)
    m21B = tsunami_unit_source_2_raster(tsunami21)
    m22B = tsunami_unit_source_2_raster(tsunami22)

    sum1B = m11B + m12B
    sum_2B = m11B + m21B
    sum_allB = m11B + m12B + m21B + m22B

    # Note: This 'test' is really better illustrated with plotting -- we see a
    # strong ridge in the m11 raster (and similar), and it is gone in the m12 raster.

    r1 = diff(range(as.matrix(sum_all)))
    r2 = diff(range(as.matrix(sum_allB)))
    expect_true(r1*0.9 > r2)
    
    r1 = diff(range(as.matrix(m11)))
    r2 = diff(range(as.matrix(m11B)))
    expect_true(r1*0.8 > r2)
    
    r1 = diff(range(as.matrix(m12)))
    r2 = diff(range(as.matrix(m12B)))
    expect_true(r1*0.8 > r2)

    # png('Slip_tapering_effects.png', width=11,height=10,units='in',res=300)
    # par(mfrow=c(2,2))
    # nc = c(123, 137)
    # plot(sum_2, zlim=c(-0.2, 0.5), col=rainbow(255))
    # abline(v=xFromCol(sum_2,nc), lty='longdash')
    # title('No tapering')
    # plot(sum_2B, zlim=c(-0.2, 0.5), col=rainbow(255))
    # abline(v=xFromCol(sum_2,nc), lty='longdash')
    # title('Tapering of slip')
    # plot(sum_2[,nc[1],],t='o')
    # points(sum_2B[,nc[1],],t='l',col='red', lwd=2)
    # title('Left transect')
    # plot(sum_2[,nc[2],],t='o')
    # points(sum_2B[,nc[2],],t='l',col='red', lwd=2)
    # title('Right transect')
    # dev.off()

})
