PROGRAM unit_tests
    USE coarray_utilities_mod, only: test_coarray_utilities_mod
    USE coarray_point2point_comms_mod, only: test_coarray_point2point_comms_mod
    USE nested_grid_comms_mod, only: test_nested_grid_comms_mod
    IMPLICIT NONE

    integer, parameter :: stdout = 6

    sync all
    write(stdout, *) 'Testing coarray_point2point_coms_mod'
    sync all

    call test_coarray_point2point_comms_mod()

    sync all
    write(stdout,*) 'Testing coarray_utilities_mod'
    sync all 

    call test_coarray_utilities_mod()

    sync all
    write(stdout,*) 'Testing nested_grid_comms_mod'
    sync all 

    call test_nested_grid_comms_mod()



END PROGRAM
