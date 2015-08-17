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

    ds$fine_downdip_transects = ds$unit_source_grid    

    # Get tsunami surface points
    tsunami_surface_points_lonlat = expand.grid(seq(-1,1,len=200), 
        seq(-1,1,len=200)) 

    # Make tsunami with our routines
    tsunami1 = make_tsunami_unit_source(1, 1, ds, 
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
    expect_that(max(abs(range(tsunami1$tsunami_source$zdsp - tsunami2$zdsp))) < 1.0e-02, is_true())
    
})
