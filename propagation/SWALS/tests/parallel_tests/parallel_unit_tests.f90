program parallel_unit_tests
    !! Run unit-tests which employ coarrays and/or MPI (i.e. distributed-memory parallel)

    use coarray_utilities_mod, only: test_coarray_utilities_mod
    use coarray_point2point_comms_mod, only: test_coarray_point2point_comms_mod
    use coarray_intrinsic_alternatives, only: sync_all_generic, swals_mpi_init, swals_mpi_finalize
    use nested_grid_comms_mod, only: test_nested_grid_comms_mod
    use multidomain_mod, only: test_multidomain_mod
    use iso_fortran_env, only: output_unit
    implicit none

    integer, parameter :: stdout = output_unit
    integer :: ierr

    call swals_mpi_init

    call sync_all_generic
    write(stdout, *) 'Testing coarray_point2point_coms_mod'
    call sync_all_generic

    call test_coarray_point2point_comms_mod()

    call sync_all_generic
    write(stdout,*) 'Testing coarray_utilities_mod'
    call sync_all_generic 

    call test_coarray_utilities_mod()

    call sync_all_generic
    write(stdout,*) 'Testing nested_grid_comms_mod'
    call sync_all_generic 

    call test_nested_grid_comms_mod()

    call sync_all_generic
    write(stdout,*) 'Testing multidomain_mod'
    call sync_all_generic 

    call test_multidomain_mod()

    call swals_mpi_finalize

end program
