context('test_spherical_to_cartesian2d_and_inverse')

test_that('test_spherical_to_cartesian2d_and_inverse', {
#test_spherical_to_cartesian2d_and_inverse<-function(){

    ## Test 1: Check that the inverse relation holds

    lonlat = cbind(20 + 0.1*runif(10), -19 + 0.1*runif(10))
    origin = c(20, -19)
    # Convert to local xy system
    new_coords = spherical_to_cartesian2d_coordinates(lonlat, 
        origin_lonlat = origin)

    # Back-caculate old coordinates
    old_coords = cartesian2d_to_spherical_coordinates(new_coords, 
        origin_lonlat = origin)

    expect_that(all(abs(lonlat - old_coords) < 1.0e-06), is_true())

    ## Test 2: Works for single-point vector input (as well as matrix, tested
    ## above)
    new_coord2 = spherical_to_cartesian2d_coordinates(lonlat[1,1:2], 
        origin_lonlat = origin)
    old_coord2 = cartesian2d_to_spherical_coordinates(new_coord2,
        origin_lonlat = origin)

    expect_that(all(abs(old_coord2 - old_coords[1,]) == 0), is_true())

    ## Test 3: Check that distances are similar to UTM projection (cannot be indentical since
    ##         UTM assumes ellipsoidal model of earth)
    ##         Also check distances are similar to spherical distances

    # Range of origins, including determininstic and random points
    origins = list( c(144.96, -37.81), # Melbourne
                    c(144.96, -60),  # Far south
                    c(0, 62), # High latitude
                    c(0, -62), # High latitude south
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59)),
                    c(runif(1, -180, 359), runif(1, -60, 59))
                  ) 

    # Relative distance error tolerances for each origin
    # In reality the value can be low for equatorial regions, with larger
    # errors for high latitudes N/S
    tols = c(0.02, rep(0.04, length(origins) - 1))

    for(ind in 1:length(origins)){

        origin = origins[[ind]]
        tol = tols[[ind]]

        # Make a random set of coordinates 'near' origin
        lonlat = cbind(origin[1] + runif(100), origin[2] + runif(100))

        # Convert to our local coordinate system and compute distance matrix
        new_coords = spherical_to_cartesian2d_coordinates(lonlat, 
            origin_lonlat = origin)
        distances0 = as.matrix(dist(new_coords, diag=FALSE, upper=TRUE))

        # Convert to UTM coordinates
        lonlat_wgs84 = SpatialPoints(lonlat, proj4string=CRS("+init=epsg:4326"))
        local_utm_proj4string = lonlat2utm(origin) 
        lonlat_utm = spTransform(lonlat_wgs84, CRS(local_utm_proj4string))
        lonlat_utm_c = coordinates(lonlat_utm)

        # Compute distance matrix for UTM version
        distances1 = as.matrix(dist(lonlat_utm_c, diag=FALSE, upper=TRUE))

        # Check that distance matrix is 'very close'
        expect_that(
            (max( abs(distances0 - distances1)/distances0, na.rm=TRUE ) < tol),
            is_true())

        # Now compare distances0 against spherical distances
        # Should be similar but not exactly the same (linear vs spherical)
        distances_sphere = distances1*0
        for(i in 1:length(lonlat[,1])){
            for(j in 1:length(lonlat[,1])){
                distances_sphere[i,j] = distHaversine(lonlat[i,], lonlat[j,])
            }
        }

        expect_that(
            (max( abs(distances_sphere - distances0)/distances0, na.rm=TRUE) < tol),
            is_true())

    }

})

