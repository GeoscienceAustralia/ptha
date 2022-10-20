
module local_routines 
    !!
    !! NTHMP tsunami-currents benchmark problem, Conical island on a triangular shelf
    !!
    use global_mod, only: dp, ip, charlen, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type

    implicit none

    real(dp), parameter :: mesh_refine = 1.0_dp ! Factor to increase resolution 
    
    ! Grid geometry, using extended domain following Macais et al. 2020 so that initial conditions can
    ! be imposed (rather than boundary conditions)
    real(dp), parameter    :: global_ll(2) = [-9.0_dp, -13.0_dp] 
    real(dp), parameter    :: global_lw(2) = [44.6_dp,  13.0_dp] - global_ll
    integer(ip), parameter :: global_nx(2) = nint(global_lw*10*mesh_refine)

    real(dp), parameter :: initial_sea_level = 0.78_dp

    contains 

    ! Main geometry setup routine
    subroutine set_initial_conditions(domain)
        class(domain_type), target, intent(inout):: domain

        type(multi_raster_type):: elevation_data
        real(dp), allocatable:: x(:), y(:)
        integer(ip):: j
        ! Wave parameters, Macais et al. 2020
        real(dp), parameter :: wa=0.39_dp, wx0=-3.3_dp, wH=0.78_dp

        real(dp), parameter :: wall_t = 0.2_dp, wall_h = 10.0_dp
        
        !
        ! Setup elevation
        !
        call elevation_data%initialise([character(len=charlen) :: './bathy/bathy_with_cone.tif'])
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x,y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip)

            ! Add walls on 3 boundaries -- 20cm thick
            if(domain%y(j) < global_ll(2) + wall_t)                domain%U(:,j,ELV) = wall_h
            if(domain%y(j) > global_ll(2) + global_lw(2) - wall_t) domain%U(:,j,ELV) = wall_h
            where(x > global_ll(1) + global_lw(1) - wall_t)        domain%U(:,j,ELV) = wall_h 
        end do

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), &
            maxval(domain%U(:,:,ELV))

        !
        ! Setup wave forcing
        !
        domain%msl_linear = initial_sea_level

        do j = 1, domain%nx(2)
            y = domain%y(j)
            ! Solitary wave
            domain%U(:,j,STG) = initial_sea_level + wa / cosh( (x - wx0)*sqrt(3.0_dp*wa/(4.0_dp*wH**3)))
            domain%U(:,j,UH) = &
                max(domain%U(:,j,STG) - domain%U(:,j,ELV), 0.0_dp) * & ! depth
                max(domain%U(:,j,STG) - initial_sea_level, 0.0_dp) * & ! wave perturbation
                sqrt(gravity/wH) ! sqrt(g/H)
        end do

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        ! ( Manning coefficient )^2
        domain%manning_squared = 0.005_dp**2

        deallocate(x,y)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Seaside_OSU
    !!
    !! NTHMP tsunami-currents benchmark problem, Seaside OSU Experiment
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type, setup_multidomain
    use boundary_mod, only: flather_boundary
    use timer_mod, only: timer_type
    use logging_mod, only: log_output_unit
    use local_routines
    implicit none

    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter ::  global_dt = 0.0125_dp / mesh_refine 


    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.05_dp
    real(dp), parameter :: final_time = 40._dp
    integer(ip), parameter :: print_every_nth_step = 4_ip, &
        write_grids_every_nth_step = 4_ip

    ! Useful misc variables
    integer(ip):: j, i, nd, forcing_case

    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 1
    allocate(md%domains(nd))

    !
    ! Setup basic metadata
    !

    ! Main domain
    md%domains(1)%lower_left =global_ll
    md%domains(1)%lw = global_lw
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method
    !md%domains(1)%use_dispersion = .true.
    !if(md%domains(1)%use_dispersion) md%domains(1)%minimum_nesting_layer_thickness = 25_ip

    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do

    ! FIXME: Add point gauges
    call md%set_point_gauges_from_csv('point_gauge_locations.csv', skip_header=1_ip)

    ! Boundary conditions
    md%domains(1)%boundary_subroutine => flather_boundary

    call md%make_initial_conditions_consistent()
    call md%set_null_regions_to_dry()

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        print*, 'domain: ', j, 'ts: ',  md%domains(j)%stationary_timestep_max()
    end do

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)
        ! Write gauges every time-step, but print and write grids less often
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency,&
            write_grids_less_often=write_grids_every_nth_step, &
            print_less_often = print_every_nth_step,&
            write_gauges_less_often= 1_ip)
        call program_timer%timer_end('IO')

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    print*, ''
    call program_timer%print(output_file_unit=log_output_unit)

end program
