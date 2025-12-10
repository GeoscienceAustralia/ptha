module local_routines 
    !!
    !! Stationary nesting
    !!
    use global_mod, only: dp, ip, charlen, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    implicit none

    real(dp), parameter :: bed_elev = -4000.0_dp

    contains 

    subroutine set_initial_conditions_nesting_stability(domain)
        !
        ! This function sets the initial conditions in a domain
        !
        class(domain_type), target, intent(inout):: domain

        type(multi_raster_type):: stage_data
        integer(ip) :: i, j
        real(dp), allocatable :: x(:), y(:)
        real(dp), parameter :: wall = 10.0_dp

        ! Initial fluxes are zero
        domain%U(:,:,UH:VH) = 0.0_dp

        allocate(x(domain%nx(1)), y(domain%nx(1)))

        ! Stationary with uneven topography. Used to study nesting instabilities 

        ! Flat stage
        domain%U(:,:,STG) = domain%msl_linear

        ! Make a bump in the topography, including over the nesting boundary.
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            domain%U(:,j,ELV) = bed_elev - &
                (0.8_dp*bed_elev)*abs(sin(x*2*pi/60000.0_dp)*cos(y*2*pi/90000.0_dp))*exp(-(x/50000.0_dp)**2 - (y/50000.0_dp)**2)
        end do

        ! Boundary walls
        domain%U(1:2,:,ELV) = wall
        domain%U(:,1:2,ELV) = wall
        domain%U(domain%nx(1)-1:domain%nx(1),:,ELV) = wall
        domain%U(:, domain%nx(2)-1:domain%nx(2),ELV) = wall

        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV)+1.0e-07_dp)

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program nesting_stability
    !!
    !! Stationary flow with nesting, used to detect unstable nesting
    !!
    use global_mod, only: ip, dp, charlen, default_linear_timestepping_method, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use domain_mod, only: UH, VH
    use logging_mod, only: log_output_unit
    use local_routines
    use coarray_intrinsic_alternatives, only: swals_mpi_init, swals_mpi_finalize
    implicit none

    integer(ip):: i, nsteps
    real(dp):: last_write_time
    type(multidomain_type) :: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 5000.0_dp
    real(dp), parameter :: final_time = 5000.0

    real(dp), parameter :: mesh_refine = 4.0_dp

    ! length/width
    real(dp), parameter, dimension(2):: global_lw = 200000.0_dp * [1.0_dp, 1.0_dp]*2
    ! lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = -global_lw/2.0_dp
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = nint([100, 50]*mesh_refine + 1) ! Deliberately uneven

    integer(ip), parameter :: timestepping_refinement_factor = 1_ip ! Overall timestep scale
    real(dp) :: global_dt = 1.2_dp * 401.0_dp/(global_nx(1) * timestepping_refinement_factor) !* (1.0_dp/3.0_dp)
    integer(ip), parameter :: nested_timestepping_refinement_factor = 3_ip * timestepping_refinement_factor

    !integer(ip), parameter :: mnlt = 40 ! Min nesting layer thickness:
    real(dp), parameter :: nesting_layer_thick_on_depth = 2.5_dp

    ! Misc
    integer(ip) :: j, nd
    character(len=charlen) :: buffer, outer_grid_ts_method

    call swals_mpi_init

    ! Read load balance file as first commandline argument -- can be empty ""
    ! in which case defaults will apply.
    call get_command_argument(1, outer_grid_ts_method)
  
    nd = 2 
    allocate(md%domains(2))

    !
    ! Set the domain properties
    !
    md%domains(1)%timestepping_method = outer_grid_ts_method
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%msl_linear = 0.0_dp
    md%domains(1)%minimum_nesting_layer_thickness = nint(nesting_layer_thick_on_depth * (-bed_elev)/maxval(global_lw/global_nx)) !mnlt
    md%domains(1)%theta = 4.0_dp
    md%domains(1)%timestepping_refinement_factor = timestepping_refinement_factor
    md%domains(1)%use_dispersion = .true.

    ! Nested domain
    call md%domains(2)%match_geometry_to_parent(&
        parent_domain=md%domains(1), &
        lower_left = global_ll + global_lw/4.0_dp, &
        upper_right = global_ll + global_lw*(13.0_dp/24.0_dp), & ! Deliberately not symmetric
        dx_refinement_factor = 3_ip, &
        timestepping_refinement_factor = nested_timestepping_refinement_factor)
    md%domains(2)%timestepping_method = 'midpoint' !
    md%domains(2)%theta = 4.0_dp
    md%domains(2)%msl_linear = 0.0_dp
    md%domains(2)%minimum_nesting_layer_thickness = nint(&
        nesting_layer_thick_on_depth * (-bed_elev)/maxval(global_lw/global_nx)) !mnlt
    md%domains(2)%use_dispersion = .true.

    ! Output variables to store
    do j = 1, size(md%domains)
        md%domains(j)%time_grids_to_store = [character(len=charlen):: 'stage', 'uh', 'vh']
        md%domains(j)%nontemporal_grids_to_store = [character(len=charlen):: 'max_stage', 'max_speed', 'elevation0']
    end do

    call md%setup()

    ! Set initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions_nesting_stability(md%domains(j))
    end do
    call md%make_initial_conditions_consistent()

    write(log_output_unit, *) 'Max timestep: ', md%stationary_timestep_max()

    ! Evolve the solution
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

    end do

    call md%finalise_and_print_timers
    call swals_mpi_finalize

end program
