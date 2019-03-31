program unit_tests

    use coarray_utilities_mod, only: test_coarray_utilities_mod
    use coarray_point2point_comms_mod, only: test_coarray_point2point_comms_mod
    use nested_grid_comms_mod, only: test_nested_grid_comms_mod
    use multidomain_mod, only: test_multidomain_mod
    use iso_fortran_env, only: output_unit
    implicit none

    integer, parameter :: stdout = output_unit

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

    sync all
    write(stdout,*) 'Testing multidomain_mod'
    sync all 

    call test_multidomain_mod()


end program
