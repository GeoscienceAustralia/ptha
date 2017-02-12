PROGRAM unit_tests
    ! Test code for read_raster_mod
    !USE ISO_C_BINDING
    !USE global_mod, only: charlen, ip, dp
    USE read_raster_mod, only: test_read_raster1 ! Important to include the module here as well
    USE spherical_mod, only: test_spherical_mod
    USE points_in_poly_mod, only: test_points_in_poly_mod
    USE which_mod, only: test_which
    USE point_gauge_mod, only: test_point_gauge_mod
    USE linear_interpolator_mod, only: test_linear_interpolator_mod
    !USE coarray_utilities_mod, only: test_coarray_utilities_mod
    !USE nested_grid_comms_mod, only: test_nested_grid_comms_mod
    IMPLICIT NONE

    print*, 'Testing read raster'
    call test_read_raster1()

    print*, 'Testing spherical mod'
    call test_spherical_mod()

    print*, 'Testing points_in_poly'
    call test_points_in_poly_mod()

    print*, 'Testing which'
    call test_which()
    
    print*, 'Testing point_gauge_mod'
    call test_point_gauge_mod()

    print*, 'Testing linear_interpolator_mod'
    call test_linear_interpolator_mod()

    !print*, 'Testing coarray_utilities_mod'
    !call test_coarray_utilities_mod()

    !print*, 'Testing nested_grid_comms_mod'
    !call test_nested_grid_comms_mod()
END PROGRAM
    
