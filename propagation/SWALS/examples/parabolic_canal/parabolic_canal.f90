module local_routines 
    !
    ! Setup for the parabolic canal problem
    !
    use global_mod, only: dp, ip, wall_elevation, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use logging_mod, only: log_output_unit
    implicit none

    real(dp), parameter :: problem_a = 3000.0_dp
    real(dp), parameter :: problem_h0 = 10.0_dp
    real(dp), parameter :: problem_g = 10.0_dp

    real(dp) :: problem_ohm ! parameter that appears in the solution

    contains 

    subroutine set_initial_conditions(domain)
        !
        ! Initial conditions
        !
        class(domain_type), target, intent(inout):: domain

        integer(ip):: i,j
        real(dp):: x, y, cx, cy, initial_stage_1, initial_stage_2, radius
        real(dp) :: imid, jmid, coriolis_f

        
        !We assume coriolis is constant over these small spatial scales
        write(log_output_unit,*) 'Coriolis range: ', minval(domain%coriolis), maxval(domain%coriolis)
        coriolis_f = domain%coriolis(domain%nx(2)/2)
        write(log_output_unit,*) 'Forcing coriolis = ', coriolis_f
        domain%coriolis = coriolis_f

        ! This parameter appears in the solution repeatedly. 
        problem_ohm = sqrt(coriolis_f**2 + 2 * gravity * problem_h0/problem_a**2)

        ! Define parabolic elevation
        imid = domain%nx(1)/2.0_dp
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                ! A 'cartesian' x coordinate
                x = (i-imid)*domain%distance_bottom_edge(j)
                ! Elevation
                domain%U(i,j,ELV) = - problem_h0 * (1.0_dp - (x/problem_a)**2)
            end do
        end do

        !@! Walls to the north and south and east and west
        !domain%U(:,1:2, ELV) = 100.0_dp
        !domain%U(:, (domain%nx(2)-1):domain%nx(2), ELV) = 100.0_dp
        !domain%U(1:2,:, ELV) = 100.0_dp
        !domain%U((domain%nx(1)-1):domain%nx(1),:, ELV) = 100.0_dp

        ! Compute the initial stage, specify the initial VH (nonzero)
        ! and the initial UH (zero).
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                ! A 'cartesian' x coordinate
                x = (i-imid)*domain%distance_bottom_edge(j)
                ! Stage
                domain%U(i,j,STG) = -problem_g**2 * problem_h0 / problem_a**2 + &
                    2 * problem_G * problem_h0 / problem_a**2 * x
                ! Ensure stage >= elevation
                domain%U(i,j,STG) = max(domain%U(i,j,STG), domain%U(i,j,ELV))

                ! VH
                domain%U(i,j,VH) = - problem_G * coriolis_f * (domain%U(i,j,STG) - domain%U(i,j,ELV))

                ! UH begins as zero
                domain%U(i,j,UH) = 0.0_dp

            end do
        end do

        domain%manning_squared = 0.0_dp


    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program parabolic_canal
    !!
    !! Flow in a parabolic canal with coriolis. This problem has an analytical solution
    !! by Thacker, which is reviewed in Sampson et al. (2006) Moving boundary shallow water flow above
    !! parabolic bottom topography, ANZIAM J. 47 (EMAC2005) pp.C373–C387
    !!
    use global_mod, only: ip, dp, charlen, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use file_io_mod, only: read_csv_into_array
    use local_routines
    use boundary_mod, only: transmissive_boundary
    implicit none

    integer(ip):: i, nsteps, yi
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 30.0_dp
    real(dp), parameter :: final_time = 7200.0_dp * 1

    real(dp), parameter :: mesh_refine = 3.0_dp

    real(dp), parameter, dimension(2):: global_ll = [0.0_dp, 30.0_dp] ! Latitude so we have Coriolis
        ! lower-left corner coordinate

    ! A 'wide' setup in the NS direction. This performs poorly for the VH component -- perhaps because of
    ! the gradients in grid-size in the NS direction (due to lon/lat coordinates) introduces genuinely 
    ! non-NS-periodic effects (and the VH flow is very slight compared to the UH flow, so easily corrupted 
    ! by any numerical issues)
    !real(dp), parameter, dimension(2):: global_lw = [0.10372884337_dp, 0.0044915764_dp*20] 
    !    ! length/width (10 km x 10km @ 30 degrees latitude)
    !integer(ip), parameter, dimension(2):: global_nx = nint([200, 200]*mesh_refine) !
    !    ! grid size (number of x/y cells)

    ! A 'narrow' setup in the NS direction.
    real(dp), parameter, dimension(2):: global_lw = [0.10372884337_dp, 0.0044915764_dp]
        ! length/width (10 km x 500 m @ 30 degrees latitude)
    integer(ip), parameter, dimension(2):: global_nx = nint([200, 5]*mesh_refine) !
        ! grid size (number of x/y cells)

    real(dp) :: non_adaptive_ts 
        ! Fixed timestep
    character(len=charlen) :: timestepping_method

#ifndef SPHERICAL
    stop 'Program assumes spherical coordinates - you must compile with -DSPHERICAL'
#endif

    call get_command_argument(1, timestepping_method)

    ! Only one domain
    allocate(md%domains(1))

    md%domains(1)%timestepping_method = timestepping_method !default_nonlinear_timestepping_method
    md%domains(1)%lw = global_lw
    md%domains(1)%nx = global_nx
    md%domains(1)%lower_left = global_ll

    ! Make the domain NS periodic
    md%periodic_ys = [global_ll(2), global_ll(2) + global_lw(2)]

    ! Setup 
    call md%setup()

    ! Set initial conditions
    CALL set_initial_conditions(md%domains(1))

    !do i = 1, size(md%domains)
    !    ! Enforce a constant coriolis -- noting multiple parameters affect the coriolis treatment
    !    yi = size(md%domains(i)%y)/2 + 1
    !    if(allocated(md%domains(i)%coriolis)) md%domains(i)%coriolis = md%domains(i)%coriolis(yi)
    !    if(allocated(md%domains(i)%coriolis_bottom_edge)) md%domains(i)%coriolis = md%domains(i)%coriolis_bottom_edge(yi)
    !    if(allocated(md%domains(i)%tanlat_on_radius_earth)) md%domains(i)%coriolis = md%domains(i)%coriolis_bottom_edge(yi)
    !end do

    ! Choose the fixed timestep
    if(md%domains(1)%timestepping_method == 'cliffs') then
        ! Currently cliffs does not perform well on this problem for the VH component.
        ! While there are some wet-dry artefacts, possibly it is also more sensitive to the
        ! NS-periodic approximation than the other solvers?
        non_adaptive_ts = md%domains(1)%stationary_timestep_max() * 0.7_dp
        md%domains(1)%cliffs_minimum_allowed_depth = 0.03_dp
    else
        non_adaptive_ts = md%domains(1)%stationary_timestep_max() * 0.9_dp
    end if

    call md%make_initial_conditions_consistent()


    ! Evolve the code
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency = approximate_writeout_frequency)

        if(md%domains(1)%time > final_time) then
            exit 
        end if

        call md%evolve_one_step(non_adaptive_ts)

    end do

    call md%finalise_and_print_timers

end program
