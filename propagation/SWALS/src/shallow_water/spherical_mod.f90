
module spherical_mod
    !!
    !! Useful routines for computation on sphere
    !!
    use global_mod, only: dp, ip, radius_earth, force_double, pi
    implicit none

    real(dp), parameter :: DEG2RAD = pi /180.0_dp
    !! Convert degrees to radians
#ifdef LEGACY_CORIOLIS_PARAMETER
    real(dp), parameter :: EARTH_ANGULAR_FREQ = 2.0_dp * pi / (3600.0_dp * 24.0_dp)
    !! Radians/second of earth rotation -- this assumes one rotation per day,
    !! which is very slightly incorrect, but matches JAGURS and COMCOT.
#else
    real(dp), parameter :: EARTH_ANGULAR_FREQ = 7.292115e-05_dp
    !! Radians/second of earth rotation, with slightly more than one rotation per day
    !! (which is correct).
#endif

    contains

    elemental function area_lonlat_rectangle(lon1, lat1, dlon, dlat, flat) result (area_sp)
        !! Compute the area on a sphere of a 'rectangle' bounded by lines of
        !! constant latitude and longitude, defined by lon1, lon1+dlon, lat1, lat1+dlat.
        !! Here lon1, lat1 must be the lower-left corner.
        !!
        !! The input parameters are in degrees, with latitude ranging from -90 to 90.
        !!
        !! Solution is explained here
        !! http://mathforum.org/library/drmath/view/63767.html
        !!
        !! Note similarity with area = R^2 cos(lat) * dlat * dlon.
        !! The equation here is similar if we note [sin(lat+dlat)-sin(lat)]/dlat = cos(lat)
        !! as dlat --> 0.
        !! If flat=.true. then we use the 'cos' formula, but if flat=.FALSE. (default) we do not.
        real(dp), intent(in):: lon1, lat1, dlon, dlat
        logical, optional, intent(in):: flat
        real(dp) :: area_sp
        logical:: flat_flag
        real(force_double):: sin1, sin2

        if(present(flat)) then
            flat_flag = flat
        else
            flat_flag = .FALSE.
        end if

        if(flat_flag) then
            area_sp = radius_earth * radius_earth * &
                cos((lat1 + 0.5_dp*dlat)*DEG2RAD) * abs(dlat) * DEG2RAD *&
                mod(abs(dlon), 360.0_dp)*DEG2RAD
        else
            ! Here with single precision, we make an effort to avoid
            ! cancellation in the diff(sin) term. Otherwise in the unit-tests,
            ! we would need different error tol's for single and double dp
            sin1 = sin(real(lat1, force_double)*real(DEG2RAD, force_double))
            sin2 = sin((real(lat1, force_double) + real(dlat, force_double)) * real(DEG2RAD, force_double))
            area_sp = radius_earth * radius_earth * &
                abs( sin1 - sin2 ) * &
                mod(abs(dlon), 360.0_dp)*DEG2RAD
        end if

    end function

    pure elemental function distance_haversine(p1_lon, p1_lat, p2_lon, p2_lat, radius) result(dist)
        !! Haversine distance formula giving the distance (meters) between points p1, p2
        !! on a sphere with the specified radius (by default = radius_earth in meters)
        real(dp), intent(in) :: p1_lon, p1_lat
            !! Lon/lat coordinates of input point
        real(dp), intent(in) :: p2_lon, p2_lat
            !! Lon/lat coordinates 
        real(dp), optional, intent(in) :: radius
            !! Sphere radius. The output distance will have the same units as radius. By default use radius_earth
        real(dp) :: dist

        ! Try to avoid round-off error in single precision case by using double-precision here
        real(force_double) :: a, dlat, dlon, r_

        r_ = radius_earth * 1.0_force_double
        if(present(radius)) r_ = radius * 1.0_force_double

        dlat = (p2_lat - p1_lat) * 1.0_force_double * DEG2RAD
        dlon = (p2_lon - p1_lon) * 1.0_force_double * DEG2RAD
        a = sin(dlat*0.5_force_double)**2 + cos(p1_lat*1.0_force_double*DEG2RAD)*&
            cos(p2_lat*1.0_force_double*DEG2RAD)*sin(dlon*0.5_force_double)**2
        a = min(a, 1.0_force_double)
        dist = 2.0_dp * atan2(sqrt(a), sqrt(1.0_dp - a)) * r_

    end function

    subroutine test_spherical_mod()
        !! Unit tests
        real(dp):: lon1, lon2, lat1, lat2, area
        real(dp):: ans, val
        ! the test cases assume this radius
        real(dp), parameter:: test_radius = 6378137.0_dp
        real(dp), parameter:: radius_adjustment_factor = (test_radius/radius_earth)**2

        !! Test 1
        ans = 12391399902.0_dp
        lon1 = 0.0_dp
        lon2 = 1.0_dp
        lat1 = 0.0_dp
        lat2 = 1.0_dp

        area = area_lonlat_rectangle(lon1, lat1, 1.0_dp, 1.0_dp)

        if(abs(area*radius_adjustment_factor - ans) < 1.0e-06_dp * ans ) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Test 1.5 -- 'flat' case
        ans = radius_earth * radius_earth * cos(0.5_dp * DEG2RAD) * 1.0_dp * DEG2RAD * 1.0_dp * DEG2RAD
        area = area_lonlat_rectangle(lon1, lat1, 1.0_dp, 1.0_dp, .TRUE.)
        if(abs(area - ans) < 1.0e-06_dp * ans) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !! Test 2
        ans = 3379082.32316_dp
        lon1 = 108.0_dp
        lon2 = 108.0_dp + 1.0_dp/60.0_dp
        lat1 = -11.0_dp
        lat2 = -11.0_dp + 1.0_dp/60.0_dp

        area = area_lonlat_rectangle(lon1, lat1, 1.0_dp/60.0_dp, 1.0_dp/60.0_dp)

        if(abs(area*radius_adjustment_factor - ans) < 1.0e-06_dp * ans) then
            print*, 'PASS'
        else
            print*, 'FAIL', area * radius_adjustment_factor, ans
        end if

        !! Test 3
        ans = 1689517.29722_dp
        lon1 = 108.0_dp
        lon2 = 108.0_dp + 1.0_dp/60.0_dp
        lat1 = -11.0_dp
        lat2 = -11.0_dp + 1.0_dp/120.0_dp

        area = area_lonlat_rectangle(lon1, lat1, 1.0_dp/60.0_dp, 1.0_dp/120.0_dp)

        if(abs(area*radius_adjustment_factor - ans) < 1.0e-06_dp * ans) then
            print*, 'PASS'
        else
            print*, 'FAIL', area * radius_adjustment_factor, ans
        end if

        !! Test 1 of distance_haversine
        lon1 = 0.0_dp
        lon2 = 1.0_dp
        lat1 = 0.0_dp
        lat2 = 0.0_dp
        ans = 111319.490793_dp ! In R, " distHaversine(c(0, 0), c(1, 0)) "
        val = distance_haversine(lon1, lat1, lon2, lat2, test_radius)
        if(abs(val - ans) < 1.0e-5_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL distance_haversine 1', ans, val
        end if

        !! Test 2 of distance_haversine
        lon1 = 0.0_dp
        lon2 = 1.0_dp
        lat1 = 0.0_dp
        lat2 = 1.0_dp
        ans = 157425.537108_dp ! In R, " distHaversine(c(0, 0), c(1, 1)) "
        val = distance_haversine(lon1, lat1, lon2, lat2, test_radius)
        if(abs(val - ans) < 1.0e-5_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL distance_haversine 2', ans, val
        end if

        !! Test 3 of distance_haversine
        lon1 = -179.0_dp
        lon2 = 12.0_dp
        lat1 = 56.0_dp
        lat2 = -72.0_dp
        ans = 18184349.8625_dp  ! In R, " distHaversine(c(-179, 56), c(12, -72)) "
        val = distance_haversine(lon1, lat1, lon2, lat2, test_radius)
        if(abs(val - ans) < 1.0e-3_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL distance_haversine 3', ans, val
        end if

        !! Test 4 of distance_haversine -- using longitudes well outside [-180, 360]
        lon1 = -179.0_dp - 360.0_dp * 3
        lon2 = 12.0_dp + 360.0_dp * 2
        lat1 = 56.0_dp
        lat2 = -72.0_dp
        ans = 18184349.8625_dp  ! In R, " distHaversine(c(-179, 56), c(12, -72)) "
        val = distance_haversine(lon1, lat1, lon2, lat2, test_radius)
        if(abs(val - ans) < 1.0e-3_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL distance_haversine 4', ans, val
        end if

    end subroutine

end module
