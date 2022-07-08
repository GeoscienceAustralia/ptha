module local_routines 
    !!
    !! Setup for the dam-break problem
    !!
    use global_mod, only: dp, ip, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    contains 

    subroutine set_initial_conditions_dam(domain, h_upstream, h_downstream)
        !!
        !! Initial conditions for the dam-break problem
        !!
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: h_upstream, h_downstream

        integer(ip):: i,j
        real(dp):: x, y, cx, cy, initial_stage_1, initial_stage_2, radius

        initial_stage_1 = h_downstream !0.0001_dp 
        initial_stage_2 = h_upstream !1.0_dp

        ! Stage
        domain%U(:,:,STG) = initial_stage_1 
        do j = 1, domain%nx(2)
            where(domain%x > 0.0_dp) domain%U(:,j,STG) = initial_stage_2
        end do

        ! Elevation
        domain%U(:,:,ELV) = 0._dp

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Frictionless
        domain%manning_squared = 0.0_dp

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program dam_break
    !!
    !! 1D Dam-break problem with a flat bed.
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    use boundary_mod, only: transmissive_boundary
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.2_dp
    real(dp), parameter :: final_time = 30.0_dp * 1

    ! length/width
    real(dp), parameter :: global_lw(2) = [200._dp, 20._dp] 
    ! lower-left corner coordinate
    real(dp), parameter :: global_ll(2) = -global_lw/2.0_dp
    ! grid size (number of x/y cells)
    integer(ip), parameter :: global_nx(2) = [200, 20] * 4 ! [400, 400] 

    ! analytical solution
    real(dp), allocatable :: analytical_solution(:,:)
    character(len=charlen):: analytical_solution_file, input_char
    real(dp) :: h_upstream, h_downstream, timestep

    ! Get the upstream/downstream initial depth from the command line
    call get_command_argument(1, input_char)
    read(input_char, *) h_upstream
    call get_command_argument(2, input_char)
    read(input_char, *) h_downstream

    ! Set up the domain
    allocate(md%domains(1))
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method
    ! rk2 still works well with this CFL (with a variable timstep), 
    ! which would usually be seen as voilating stability constraints.
    !domain%cfl = 1.49_dp
    ! Deliberately use a non-TVD limiter (> 2.0) to stress-test that approach. 
    ! It turns out to work well.
    md%domains(1)%theta = 4.0_dp 
    call md%setup

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_dam(md%domains(1), h_upstream, h_downstream)

    timestep = md%stationary_timestep_max() * 0.5_dp

    md%domains(1)%boundary_subroutine => transmissive_boundary

    ! Evolve the code
    do while (.TRUE.)

        ! Avoid storing grids often
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(timestep)

    end do

    call md%finalise_and_print_timers

end program
