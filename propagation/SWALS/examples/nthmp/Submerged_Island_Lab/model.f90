
module local_routines 
    !!
    !! NTHMP velocities Submerged Island benchmark
    !!
    use global_mod, only: dp, ip, charlen, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV

    implicit none

    !
    ! Parameters defining the domain geometry.
    !

    ! Length/width -- FULL
    real(dp), parameter :: global_lw(2) = [9.75_dp, 1.152_dp]
    real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp]

    ! Conical island details
    real(dp), parameter :: island_centre(2) = [5.0_dp, global_lw(2)/2.0_dp]
    real(dp), parameter :: island_radius_base = 0.75_dp / 2.0_dp
    real(dp), parameter :: island_radius_top = 0.05_dp / 2.0_dp
    real(dp), parameter :: island_max_height = 0.049_dp
    real(dp), parameter :: island_slope = island_max_height/(island_radius_base - island_radius_top)

    ! Manning friction
    real(dp), parameter :: friction_manning = 0.025_dp
    real(dp), parameter :: wall_friction_manning = 0.05_dp ! Option to change friction along the walls.

    ! The velocity (in the x direction) as specified in the experiment.
    real(dp), parameter :: target_velocity = 0.115_dp

    ! To get the right depth, we should have the downstream boundary stage lower than 0.0_dp
    ! This can be implemented by specifying the boundary_function for the flather boundary, or 
    ! setting domain%msl_linear (which is used in the default boundary_function). The latter 
    ! is simpler
    real(dp), parameter :: msl_linear = -0.011_dp ! with low-friction (0.01), can use -0.007_dp
    ! "Background Depth", so that MSL = 0. 
    real(dp), parameter :: depth_ocean = 0.054_dp 
   
    contains 

    ! Main setup routine
    subroutine set_initial_conditions(domain)
        type(domain_type), intent(inout):: domain

        integer(ip):: i, j
        real(dp), allocatable:: x(:), y(:)
        real(dp) :: gauges_1_to_4_x_coord, dd, gauge_xy(3,2)

        ! Stage
        domain%U(:,:,STG) = 0.0e-0_dp
        
        ! Elevation (background)
        domain%U(:,:,ELV) = -depth_ocean

        !
        ! Setup conical island
        !
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            do i = 1, domain%nx(1)

                ! Add the island
                dd = sqrt( (x(i) - island_centre(1))**2 + (y(i) - island_centre(2))**2 )
                if(dd < island_radius_base) then
                    domain%U(i,j,ELV) = -depth_ocean + &
                        min( island_slope * (island_radius_base - dd), island_max_height)
                end if

            end do
        end do
        deallocate(x,y)

        ! Add a wall on 3 sides -- the right side will be an (outflow) boundary condition
        domain%U(1:2,   :, ELV) = 1.0_dp
        domain%U(:  , 1:2, ELV) = 1.0_dp
        domain%U(:  , domain%nx(2)-1:domain%nx(2), ELV) = 1.0_dp

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))


        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        ! Manning friction
        domain%manning_squared = friction_manning*friction_manning
        ! with artificial "wall friction"
        domain%manning_squared(:, 3) = wall_friction_manning*wall_friction_manning
        domain%manning_squared(:, domain%nx(2) - 2) = wall_friction_manning*wall_friction_manning

        ! Gauges have fixed locations
        gauge_xy(1:3, 1) = [ island_centre(1) + island_radius_base + 1.0_dp, global_lw(2)/2.0_dp          , 1.0_dp]
        gauge_xy(1:3, 2) = [ island_centre(1) + island_radius_base + 1.0_dp, global_lw(2)/2.0_dp + 0.27_dp, 2.0_dp]
        call domain%setup_point_gauges(xy_coords = gauge_xy(1:2,:), gauge_ids=gauge_xy(3,:))

    end subroutine

    !
    ! This is the discharge source term
    !
    subroutine apply_discharge_forcing(domain, dt)
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt
        integer(ip) :: j

        ! Volume_in_per_unit_time * dt = dStage_in_a_timestep * cell_area = Discharge_in * dt
        !    and
        ! Discharge_in = target_velocity * target_depth * distance_left_edge
        !    so 
        ! dStage = target_velocity * target_depth * distance_left_edge * dt / cell_area

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            ! Add a volume so that the x-directed velocity = target_velocity
            domain%U(3, j, STG) = domain%U(3, j,STG) + &
                target_velocity * depth_ocean * domain%distance_left_edge(j) * dt / &
                    (domain%distance_bottom_edge(3)*domain%distance_left_edge(j))
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Submerged_Island
    !!
    !! NTHMP velocities Submerged Island benchmark
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: flather_boundary
    use timer_mod
    use logging_mod, only: log_output_unit
    use local_routines
    implicit none

    ! Useful misc variables
    integer(ip):: j, i, nd, forcing_case

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 8.0_dp ! Increase resolution by this amount
    
    real(dp), parameter ::  global_dt = 0.06_dp / mesh_refine ! * 0.5_dp

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 5.0_dp ! 0.5_dp 
    real(dp), parameter :: final_time = 1500._dp

    ! Use this to read command line arguments
    character(len=20) :: tempchar

    !
    ! Key geometric parameters are defined in the 'local routines' module
    !

    ! Grid size (number of x/y cells) in outer domain
    integer(ip), parameter:: global_nx(2) = int(global_lw*10, ip) * mesh_refine

    integer(ip), parameter :: write_grids_and_print_every_nth_step = ceiling(approximate_writeout_frequency/global_dt)

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
    md%domains(1)%cliffs_minimum_allowed_depth = 0.01_dp

    !md%domains(1)%eddy_visc_constants = [0.000_dp, 0.0_dp]
    !md%domains(1)%eddy_visc_constants = [0.000_dp, 0.5_dp]
    !md%domains(1)%use_eddy_viscosity = .true.
    md%domains(1)%msl_linear = msl_linear

    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do

    ! Build boundary conditions
    md%domains(1)%boundary_subroutine => flather_boundary 
    md%domains(1)%forcing_subroutine => apply_discharge_forcing

    call md%make_initial_conditions_consistent()
    
    ! For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        print*, 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)

        ! Write gauges every time-step, but print and write grids less often
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=0.0_dp,&
            write_grids_less_often=write_grids_and_print_every_nth_step, &
            print_less_often = write_grids_and_print_every_nth_step,&
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
