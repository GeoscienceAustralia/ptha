!
! NTHMP Benchmark problem 6
!

module local_routines 
    use global_mod, only: dp, ip, charlen, wall_elevation, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: gdal_raster_dataset_type
    use file_io_mod, only: count_file_lines, read_csv_into_array
    use linear_interpolator_mod, only: linear_interpolator_type
    use boundary_mod, only: flather_boundary, transmissive_boundary

    implicit none

    !
    ! Parameters defining the domain geometry.
    !

    ! Length/width -- FULL
    real(dp), parameter :: global_lw(2) = [29.3_dp, 30.0_dp]
    ! Length/width, neglecting the absorbing boundary region
    real(dp), parameter :: non_absorbing_lw(2) = [25.0_dp, 28.2_dp]
    ! The top/bottom of the wavemaker is offset from the model absorbing boundary by this amount
    real(dp), parameter :: wavemaker_edge_offset = (global_lw(2) - non_absorbing_lw(2))/2.0_dp + 0.38_dp
    ! Lower-left corner coordinate. Account for buffers, and the 2m x-offset between the wavemaker and the buffer
    real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp] - (global_lw - non_absorbing_lw)/2.0_dp - [2.0_dp, 0.0_dp]
    ! "Background Depth", so that MSL = 0
    real(dp), parameter :: depth_ocean = 0.32_dp

    ! Centre of conical island
    real(dp), parameter :: island_centre(2) = [12.96_dp, 13.80_dp]
    real(dp), parameter :: island_radius_base = 7.2_dp / 2.0_dp
    real(dp), parameter :: island_radius_top = 2.2_dp / 2.0_dp
    real(dp), parameter :: island_slope = 0.25_dp
    real(dp), parameter :: island_max_height = 0.625_dp

    ! Wavemaker x position
    real(dp), parameter :: wavemaker_x_start = 0.0_dp

    ! Use the wavemaker forcing, or a gauge-based estimate?
    logical :: use_wavemaker_forcing = .false.
    ! Both approaches force the velocity near the wavemaker, and also adjust the stage based
    ! on a linear plane wave approximation.
    ! The wavemaker_forcing approach estimates the velocity from the paddle displacement time-series.
    ! See "convert_obs_to_wavemaker.R" for details of the gauge-based forcing - basically we roughly estimate
    ! the wave profile 'propagated backward in time', and use that at the wavemaker forcing.
    ! Both approaches give fairly similar results, but the gauge-based estimate can
    ! better match the gauge-1-4 initial waves (likely it is correcting for some dispersion
    ! that is not in the model itself). 

    ! Manning friction
    real(dp), parameter :: low_friction_manning  = 0.01_dp
    real(dp), parameter :: high_friction_manning = 0.01_dp
   
    !
    ! Hold some data for the boundary condition
    !
    type :: boundary_information_type

        ! Wavemaker data
        character(len=charlen):: wavemaker_file = &
            '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/fdbk2abc.txt'
        real(dp), allocatable:: wavemaker_data(:,:)
        type(linear_interpolator_type):: wavemaker_position
        real(dp):: wavemaker_t0 = 0.0_dp
        ! Information on the forcing case
        integer(ip) :: forcing_case

        ! Gauge data. This has been processed to give an estimate of the velocity that should be used to
        ! force the model, like with the wavemaker. See "convert_obs_to_wavemaker.R"
        character(len=charlen) :: gauge_vel_file(3) = [ &
              './gauge_forcing_nonlinear_case1.csv',&
              './gauge_forcing_nonlinear_case2.csv',&
              './gauge_forcing_nonlinear_case3.csv']
        real(dp), allocatable:: gauge_velocity_data(:,:)
        type(linear_interpolator_type) :: gauge_velocity(3)
        real(dp):: gauge_velocity_t0 = 0.0_dp
    end type

    ! This will hold the information -- is seen by other parts of the module
    type(boundary_information_type):: boundary_information

    contains 

    !
    ! Read files with wavemaker info, and make a subroutine which applies wavemaker forcing
    ! in the boundary condition
    !
    subroutine setup_boundary_information(forcing_case)
        integer(ip), intent(in) :: forcing_case

        integer(ip):: bc_unit, nr, nc, skip, i, extra, j
        ! Checks
        real(dp) :: test_times(4) = [0.0, 1.0, 2.0, 3.0], test_vels(4)

        boundary_information%forcing_case = forcing_case

        !
        ! Wavemaker paddle information
        !
        ! Read from file
        ! First column is time, next three columns contain the paddle
        ! positions for case A, B, C.
        open(newunit=bc_unit, file=boundary_information%wavemaker_file)
        nr = count_file_lines(bc_unit)
        nc = 4
        skip = 1
        extra = 1 ! One more data point to avoid exceeding time 
        allocate(boundary_information%wavemaker_data(nr - skip + extra, nc))
        do i = 1, nr
            if(i > skip) then
                read(bc_unit, *) boundary_information%wavemaker_data(i - skip,:)
            else
                read(bc_unit, *) 
            end if 
        end do
        close(bc_unit)
        ! Convert from centimeters to meters
        boundary_information%wavemaker_data(:,2:4) = boundary_information%wavemaker_data(:,2:4)/100.0_dp
        ! Extend the time-series with a constant value, so that time does not
        ! exceed model run time
        boundary_information%wavemaker_data(nr - skip + extra,1) = 1.0e+06_dp + &
            boundary_information%wavemaker_data(nr - skip + extra - 1, 1)
        boundary_information%wavemaker_data(nr - skip + extra,2) = &
            boundary_information%wavemaker_data(nr - skip + extra - 1, 2)
        !
        ! Make the linear interpolation function for the correct forcing case
        call boundary_information%wavemaker_position%initialise(&
                boundary_information%wavemaker_data(:,1), &
                boundary_information%wavemaker_data(:,forcing_case+1))
        boundary_information%wavemaker_t0 = boundary_information%wavemaker_data(1,1)

        !
        ! Gauge information
        !
        ! Loop over each case
        do j = 1, 3 

            ! Get the data
            call read_csv_into_array(boundary_information%gauge_velocity_data, &
                boundary_information%gauge_vel_file(j), skip_header=1)

            ! Make an interpolation function
            call boundary_information%gauge_velocity(j)%initialise(&
                boundary_information%gauge_velocity_data(1,:), &
                boundary_information%gauge_velocity_data(2,:))

            ! Test
            !call boundary_information%gauge_velocity(j)%eval(test_times, test_vels)
            !do i = 1, 4
            !    print*, 'case ', j, test_times(i), test_vels(i)
            !end do

            ! Cleanup (only really required on the last j ...)
            deallocate(boundary_information%gauge_velocity_data)
        end do

    end subroutine
   
    ! 
    ! Pass this as a "boundary condition", also affecting the interior (wavemaker)
    !
    subroutine forcing_subroutine(domain)
        type(domain_type), intent(inout) :: domain

        ! Get the wavemaker position in the past/future, and differentiate to get the
        ! velocity which is imposed at the right location
        real(dp) :: time, time_future, time_past, time1(1), pos_x(1), pos_x_future(1), pos_x_past(1)
        real(dp) :: forcing_vel, forcing_vel1(1)
        real(dp), parameter :: time_lag = 0.2_dp
        real(dp) :: inflation_factor

        real(dp) :: time_limit, dstage, stage_forced
        integer(ip) :: j, i, fc, ii

        !
        ! Firstly, we do the regular boundary update
        !
        ! FIXME: Really should prevent mass flowing out.
        call flather_boundary(domain)
        !call transmissive_boundary(domain)

        !
        ! Secondly, we force the model velocities internally to simulate a wavemaker.
        ! This is imperfect, but can be done based on the wavemaker paddle, OR, based
        ! on a velocity inferred from the gauge observations.
        !
        fc = boundary_information%forcing_case

        if(use_wavemaker_forcing) then
            ! Here we differentiate the wavemaker paddle velocity to estimate the forcing velocity

            ! The forcing is imperfect -- with these parameters we try to 'calibrate' the result.
            ! See the alternative gauge_based approach.
            if(fc == 1) then
                inflation_factor = 0.9_dp
            else
                inflation_factor = 1.0_dp
            end if

            time = domain%time
            time_future = time + time_lag
            time_past = max(time - time_lag, 0.0_dp)

            ! Wavemaker offset position. 
            call boundary_information%wavemaker_position%eval([time + boundary_information%wavemaker_t0], pos_x)
            call boundary_information%wavemaker_position%eval([time_future + boundary_information%wavemaker_t0], pos_x_future)
            call boundary_information%wavemaker_position%eval([time_past   + boundary_information%wavemaker_t0], pos_x_past)
          
            ! Numerical estimate of the velocity. Time-lag needs to give enough smoothing!
            forcing_vel = (pos_x_future(1) - pos_x_past(1)) / (2.0_dp * time_lag) * inflation_factor
            time_limit = HUGE(1.0_dp) 

        else
            !
            ! Here we use a gauge-based estimate of the forcing
            ! Preprocessing is already done.
            !
            time = domain%time
            ! We do not want to force for too long (because otherwise the gauges become affected by reflections etc)
            time_limit = 1.0_dp + 6.0_dp/fc

            if(time < time_limit) then
                ! Get the forcing velocity
                call boundary_information%gauge_velocity(fc)%eval( [time], forcing_vel1)
                forcing_vel = forcing_vel1(1)
            else
                forcing_vel = 0.0_dp
            end if
            pos_x = 0.0_dp ! Apply boundary at regular wavemaker position
        end if

        ! Now impose this velocity at the wavemaker x-locations
        ! Note the gauges suggest imperfections in the experiment. For instance
        ! the leading wave at gauges 1-4 should be very similar, but there are differences.
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, forcing_vel, pos_x) 
        !$OMP DO
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                if( &
                    ! X coordinate within some distance of wavemaker position
                    (abs(domain%x(i) - wavemaker_x_start - pos_x(1)) < 0.5*domain%dx(1)) .and. &
                    ! Y coordinate is not too close to edge, as specified by geometry
                    (domain%y(j) - global_ll(2) > wavemaker_edge_offset) .and. & 
                    (domain%y(j) - global_ll(2) < global_lw(2) - wavemaker_edge_offset) &
                    ) then

                    ! Velocity = paddle velocity
                    domain%U(i,j,UH) = forcing_vel * depth_ocean
                    ! Stage as for a linear plane wave -- beware mass conservation violation
                    stage_forced = max( domain%U(i,j,UH) / sqrt(gravity * depth_ocean), domain%U(i,j,ELV) )
                    dstage = stage_forced - domain%U(i,j,STG)
                    domain%U(i,j,STG) = stage_forced
                    ! Enforce mass conservation -- spread it over a 3 cells (or we can force cells dry!)
                    do ii = 1, 3
                        domain%U(i-ii,j,STG) = domain%U(i-ii,j,STG) - dstage/3
                    end do

                end if
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine 

    ! Main setup routine
    subroutine set_initial_conditions(domain, forcing_case)
        class(domain_type), target, intent(inout):: domain
        integer(ip), intent(in) :: forcing_case

        integer(ip):: i, j
        character(len=charlen):: input_elevation, input_stage
        real(dp), allocatable:: x(:), y(:)
        real(dp) :: wall, dd, gauges_1_to_4_x_coord
        real(dp) :: gauge_xy(3,8), leftmost_x(3)
        real(dp) :: absorbing_boundary_width(2)

        ! Thickness of absorbing boundary region
        absorbing_boundary_width = (global_lw - non_absorbing_lw)/2.0_dp

        ! Stage
        domain%U(:,:,STG) = 0.0e-0_dp
        
        ! Elevation (background)
        domain%U(:,:,ELV) = -depth_ocean

        !
        ! Setup detailed elevation / manning 
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

                ! Friction -- try modelling absorbing boundaries with high friction. 
                if(domain%timestepping_method /= 'linear') then
                    if( (x(i) - global_ll(1) < absorbing_boundary_width(1)) .or. &
                        (x(i) - global_ll(1) > non_absorbing_lw(1) + absorbing_boundary_width(1)) .or. &
                        (y(i) - global_ll(2) < absorbing_boundary_width(2)) .or. &
                        (y(i) - global_ll(2) > non_absorbing_lw(2) + absorbing_boundary_width(2)) ) then
                        ! High friction around absorbing boundary
                        domain%manning_squared(i,j) = high_friction_manning**2 !1.0_dp * 1.0_dp
                    else
                        ! Regular manning
                        domain%manning_squared(i,j) = low_friction_manning**2 
                    end if
                end if
            end do
        end do

        deallocate(x,y)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        !! Wall boundaries (without boundary conditions)
        !wall = 0.5_dp
        !domain%U(:,1,ELV) = wall
        !domain%U(:,domain%nx(2),ELV) = wall
        !domain%U(domain%nx(1),:,ELV) = wall
        !domain%U(1,:,ELV) = wall

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        !Gauges
        leftmost_x = [5.76_dp, 6.82_dp, 7.56_dp]
        gauges_1_to_4_x_coord = leftmost_x(forcing_case)

        ! Gauges 1-4 share an x coordinate, which varies with forcing-case
        gauge_xy(1, 1:4) = gauges_1_to_4_x_coord
        gauge_xy(2, 1:4) = [16.05_dp, 14.55_dp, 13.05_dp, 11.55_dp] ! y coordinate
        gauge_xy(3, 1:4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp] ! ID
        ! Other gauges have fixed locations
        gauge_xy(1:3, 5) = [ 9.36_dp, 13.8_dp,   6.0_dp]
        gauge_xy(1:3, 6) = [10.36_dp, 13.8_dp,   9.0_dp]
        gauge_xy(1:3, 7) = [12.96_dp, 11.22_dp, 16.0_dp]
        gauge_xy(1:3, 8) = [15.56_dp, 13.80_dp, 22.0_dp]
        call domain%setup_point_gauges(xy_coords = gauge_xy(1:2,:), gauge_ids=gauge_xy(3,:))

    end subroutine

end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program BP06

    use global_mod, only: ip, dp, minimum_allowed_depth
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum
    use timer_mod
    use logging_mod, only: log_output_unit
    use local_routines
    implicit none

    ! Useful misc variables
    integer(ip):: j, i, i0, j0, centoff, nd, lg, forcing_case
    real(dp):: last_write_time, gx(4), gy(4)

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 1.0_dp ! Increase resolution by this amount
    
    real(dp) ::  global_dt = 0.020_dp / mesh_refine

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 0.2_dp
    real(dp) :: final_time = 40._dp

    ! Use this to read command line arguments
    character(len=20) :: tempchar

    !
    ! Key geometric parameters are defined in the 'local routines' module
    !

    ! Grid size (number of x/y cells) in outer domain
    integer(ip), parameter:: global_nx(2) = int(global_lw*10, ip) * mesh_refine
    ! Extent of Inner domain around the island
    integer(ip) :: nest_ratio = 3_ip
    real(dp) :: high_res_ll(2) = island_centre - [1.0_dp, 1.0_dp] * island_radius_base * sqrt(2.0_dp)
    real(dp) :: high_res_ur(2) = island_centre + [1.0_dp, 1.0_dp] * island_radius_base * sqrt(2.0_dp)

    ! Forcing-case = 1, 2, 3
    ! Passed as first command line argument
    call get_command_argument(1, tempchar)  
    read(tempchar, *) forcing_case 

    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 2
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
    md%domains(1)%timestepping_method = 'rk2'

    print*, 1, ' lw: ', md%domains(1)%lw, ' ll: ', md%domains(1)%lower_left, ' dx: ', md%domains(1)%dx, &
        ' nx: ', md%domains(1)%nx

    ! A detailed domain
    call md%domains(2)%match_geometry_to_parent(&
        parent_domain=md%domains(1), &
        lower_left=high_res_ll, &
        upper_right=high_res_ur, &
        dx_refinement_factor=nest_ratio, &
        timestepping_refinement_factor=nest_ratio)
    md%domains(2)%timestepping_method = 'rk2'

    print*, 2, ' lw: ', md%domains(2)%lw, ' ll: ', md%domains(2)%lower_left, ' dx: ', md%domains(2)%dx, &
        ' nx: ', md%domains(2)%nx

     
    ! Set the CFL limit for each model. 
    do j = 1, size(md%domains)
        md%domains(j)%cfl = merge(0.7, 0.99, md%domains(j)%timestepping_method == 'linear')
    end do

    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j), forcing_case)
    end do

    ! Build boundary conditions
    call setup_boundary_information(forcing_case)
    ! This is not really a boundary -- it's an internal wavemaker forcing
    md%domains(1)%boundary_subroutine => forcing_subroutine

    call md%make_initial_conditions_consistent()
    
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        print*, 'domain: ', j, 'ts: ', &
            md%domains(j)%linear_timestep_max()*merge(1.0, 0.5, md%domains(j)%timestepping_method == 'linear')
    end do

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)
        
        ! IO 
        if(md%domains(1)%time - last_write_time >= approximate_writeout_frequency) then
            call program_timer%timer_start('IO')
            call md%print()
            do j = 1, nd
                call md%domains(j)%write_to_output_files()
            end do
            last_write_time = last_write_time + approximate_writeout_frequency
            call program_timer%timer_end('IO')
        end if

        call program_timer%timer_start('IO')
        do j = 1, nd
            call md%domains(j)%write_gauge_time_series()
        end do
        call program_timer%timer_end('IO')

        call md%evolve_one_step(global_dt)

        if (md%domains(1)%time > final_time) exit
    end do

    call program_timer%timer_end('evolve')

    ! Print out timing info for each
    do i = 1, nd
        !lg = md%domains(i)%logfile_unit
        write(log_output_unit, *) ''
        write(log_output_unit, *) 'Timer ', i
        write(log_output_unit, *) ''
        call md%domains(i)%timer%print(output_file_unit=log_output_unit)
        call md%domains(i)%write_max_quantities()
        call md%domains(i)%finalise()
    end do

    print*, ''
    call program_timer%print(output_file_unit=log_output_unit)

end program
