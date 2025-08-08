module local_routines 
    !!
    !! Setup the "wave on a composite beach" problem
    !!

    use global_mod, only: dp, ip, charlen
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: read_gdal_raster
    use which_mod, only: which
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    implicit none

    ! Hold some data for the boundary condition
    type :: boundary_information_type
        character(charlen):: bc_file
        real(dp), allocatable:: boundary_data(:,:)
        type(linear_interpolator_type):: gauge4_ts_function
        real(dp):: boundary_elev
        real(dp):: t0
    end type

    ! This will hold the information -- is seen by other parts of the module
    type(boundary_information_type):: boundary_information

    contains 

    subroutine setup_boundary_information(bc_file, boundary_elev)
        ! Read boundary info from a file and pack into the boundary_information type 
        character(charlen), intent(in):: bc_file
        real(dp), intent(in):: boundary_elev

        integer(ip):: bc_unit, nr, nc, skip, i

        boundary_information%bc_file = bc_file
        boundary_information%boundary_elev = boundary_elev
        open(newunit=bc_unit, file=bc_file)
        nr = count_file_lines(bc_unit)
        nc = 9
        skip = 6
        allocate(boundary_information%boundary_data(nr - skip, nc))
        do i = 1, nr
            if(i > skip) then
                read(bc_unit, *) boundary_information%boundary_data(i - skip,:)
            else
                read(bc_unit, *) 
            end if 
        end do
        close(bc_unit)

        call boundary_information%gauge4_ts_function%initialise(&
                boundary_information%boundary_data(:,1), boundary_information%boundary_data(:,2))
        boundary_information%t0 = boundary_information%boundary_data(1,1)

    end subroutine
    
    function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        ! Function to evaluate the boundary at the domain, passed to model boundary conditions
        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4)

        if(i == 1) then
            call boundary_information%gauge4_ts_function%eval([t + boundary_information%t0], stage_uh_vh_elev(STG:STG))
            stage_uh_vh_elev(UH:VH) = 0.0_dp
            stage_uh_vh_elev(ELV) = domain%U(i,j,ELV) !boundary_information%boundary_elev
            stage_uh_vh_elev(STG) = max(stage_uh_vh_elev(STG), domain%U(i,j,ELV))
        else
            stage_uh_vh_elev(STG) = domain%U(i,j,STG)
            stage_uh_vh_elev(ELV) = domain%U(i,j,ELV)
            stage_uh_vh_elev(UH:VH) = 0.0_dp
        end if
    end function

    subroutine set_initial_conditions_bp2(domain, tank_bases, tank_slopes, tank_width, initial_depth)
        type(domain_type), intent(inout):: domain
        real(dp), intent(in) :: tank_bases(4), tank_slopes(4), tank_width, initial_depth

        real(dp):: tank_x(5)
        real(dp):: x, y, elev
        integer(ip):: j, i, k
        real(dp):: gauge_xy(2,11), wall

        ! Stage
        domain%U(:,:,STG) = 0.0_dp

        domain%msl_linear = 0.0_dp ! Controls radiation boundary at later times

        ! Elevation
        tank_x = 0.0_dp
        do i = 2, 5
            tank_x(i) = sum(tank_bases(1:(i-1)))            
        end do
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                x = domain%lower_left(1) + (i-0.5_dp) * domain%dx(1) 
                elev = -initial_depth
                do k = 2, 4
                     if(x > tank_x(k)) then
                        elev = elev + (min(x, tank_x(k+1)) - tank_x(k))*tank_slopes(k)
                     end if
                end do
            domain%U(i,j,ELV) = elev
            end do
        end do
      
        ! Reflective boundaries on 3 sides
        wall = 1.0_dp
        domain%U(:, 1, ELV) = wall
        domain%U(:, 2, ELV) = wall
        domain%U(:, domain%nx(2), ELV) = wall
        domain%U(:, domain%nx(2)-1, ELV) = wall
        domain%U(domain%nx(1), :, ELV) = wall
        domain%U(domain%nx(1)-1, :, ELV) = wall

        ! Try adding a basin to avoid boundary problems
        !domain%U(1:4, :,ELV) = domain%U(1, 3, ELV)
   
        ! Stage >= bed 
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))
        
        ! Define locations of gauge outputs
        gauge_xy(2,:) = 0.0_dp ! Always y == 0
        gauge_xy(1, 1:4) = 0.0_dp + 1.0e-06_dp ! Nudge x-coordinate inside the domain
        gauge_xy(1,5) = tank_x(2)
        gauge_xy(1,6) = 0.5_dp * (tank_x(3) + tank_x(2))
        gauge_xy(1,7) = tank_x(3)
        gauge_xy(1,8) = 0.5_dp * (tank_x(3) + tank_x(4))
        gauge_xy(1,9) = tank_x(4)
        gauge_xy(1,10) = 0.5_dp * (tank_x(4) + tank_x(5))
        ! Include just before the wall
        gauge_xy(1,11) = domain%x(domain%nx(1) - 2)

        call domain%setup_point_gauges(gauge_xy)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program BP02
    !!
    !! NTHMP benchmark problem 2 (and 5) -- Wave on a composite beach.
    !! This has an analytical solution (for the linear shallow water equations), and
    !! some experimental results. 
    !!
    use global_mod, only: ip, dp, minimum_allowed_depth
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: boundary_stage_transmissive_momentum, flather_boundary, &
        transmissive_boundary, boundary_stage_transmissive_normal_momentum
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    type(multidomain_type):: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.1_dp
    real(dp), parameter :: final_time = 30.0_dp

    ! Domain info
    character(charlen) :: timestepping_method
    
    ! length/width
    real(dp) :: global_lw(2), global_ll(2)
    integer(ip) :: global_nx(2) 

    ! Local variables 
    real(dp) :: timestep, base_l, dx, tank_width, tank_length, initial_depth
    character(charlen):: test_case, bc_file
    real(dp):: tank_bases(4), tank_slopes(4) 

    integer(ip), parameter :: boundary_buffer_cells = 10_ip ! So we can force the boundary with the NLSW
    real(dp) :: boundary_buffer

    ! Get the case. Values should be caseA, caseB, caseC
    call get_command_argument(1, test_case)
    select case(test_case)
        case('caseA')
            base_L = 2.40_dp
            bc_file = '../../nthmp/test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3a_analytical.txt'
        case('caseB')
            base_L = 0.98_dp
            bc_file = '../../nthmp/test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3b_analytical.txt'
        case('caseC')
            base_L = 0.64_dp
            bc_file = '../../nthmp/test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3c_analytical.txt'
        case default
            print*, 'Must specify a test case (one of caseA, caseB, caseC)'
            stop
    end select

    call get_command_argument(2, timestepping_method)  
    
    ! Resolution
    dx = 0.01_dp ! Needs to evenly divide tank_length
    boundary_buffer = dx * boundary_buffer_cells
 
    ! Tank geometry  -- add a little extra at the end so the reflective wall is in the right place
    tank_bases = [base_L, 4.36_dp, 2.93_dp, 0.9_dp + 2.0_dp*dx]
    tank_slopes = [0.0_dp, 1.0_dp/53.0_dp, 1.0_dp/150.0_dp, 1.0_dp/13.0_dp]
    tank_width = dx * 5.0_dp ! This will cause a single cell of flow (with 2dx walls on either side)
    tank_length = sum(tank_bases)
    initial_depth = 0.218_dp

    ! Large scale
    global_lw = [tank_length, tank_width]
    global_ll = [0.0_dp, -tank_width/2.0_dp]
    global_nx = global_lw/dx

    ! Setup model with 2 domains -- a small one near the boundary using NLSW (which the boundary conditions are better suited to),
    ! and a dispersive domain covering most areas.
    allocate(md%domains(2))
    !md%periodic_ys = [global_ll(2), global_ll(2) + global_lw(2)]

    ! Boundary domain with NLSW. Point is to impose the boundary condition consistently.
    md%domains(1)%lw = [boundary_buffer, global_lw(2)]
    md%domains(1)%lower_left =  global_ll 
    md%domains(1)%nx = [boundary_buffer_cells, global_nx(2)]
    md%domains(1)%timestepping_method = timestepping_method
    md%domains(1)%use_dispersion = .false. !
    !md%domains(1)%nc_grid_output%flush_every_n_output_steps = 1_ip !

    ! Main domain with dispersive solver, and halos that cover the NLSW domain
    md%domains(2)%lw = global_lw + [-boundary_buffer, 0.0_dp]
    md%domains(2)%lower_left = global_ll + [boundary_buffer, 0.0_dp]
    md%domains(2)%nx = global_nx - [boundary_buffer_cells, 0]
    md%domains(2)%timestepping_method = timestepping_method
    md%domains(2)%use_dispersion = .true. !
    md%domains(2)%minimum_nesting_layer_thickness = boundary_buffer_cells
    !md%domains(2)%nc_grid_output%flush_every_n_output_steps = 1_ip !

    !! Experiment with more nested iterations
    !md%domains(2)%ds%tridiagonal_inner_iter = 4_ip

    !! Taper off dispersion
    !md%domains(2)%ds%td1 = 0.07_dp
    !md%domains(2)%ds%td2 = 0.05_dp

    !! Experiment with using the dispersive solver without dispersion (does
    !! weird things for 'midpoint' this doesn't turn off the predictor-step
    !! dispersion).
    !md%dispersive_outer_iterations_count = -1

    call md%setup

    ! Call local routine to set initial conditions
    call set_initial_conditions_BP2(md%domains(1), tank_bases, tank_slopes, tank_width, initial_depth)
    call set_initial_conditions_BP2(md%domains(2), tank_bases, tank_slopes, tank_width, initial_depth)

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    call setup_boundary_information(bc_file, -initial_depth)
    md%domains(1)%boundary_function => boundary_function
    ! For nonlinear schemes, this boundary condition causes spurious reflections of the outgoing waves (because depths are too
    ! shallow)
    md%domains(1)%boundary_subroutine => boundary_stage_transmissive_normal_momentum

    call md%make_initial_conditions_consistent() ! Get the initial volume right

    ! Fixed timestep
    timestep = md%stationary_timestep_max() * 0.5_dp
    print*, timestepping_method, ', timestep = ', timestep

    ! Evolve the code
    do while (.TRUE.)

        ! Avoid storing grids often
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = 1_ip) ! 999999_ip)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(timestep)

        if(md%domains(1)%time > 10.0_dp) then
            ! After the initial wave has entered the domain, we change the boundary condition
            ! for nonlinear schemes to reduce reflection of outgoing waves. 
            md%domains(1)%boundary_subroutine => flather_boundary
            md%domains(1)%boundary_function => NULL()
        end if

    END DO

    call md%finalise_and_print_timers

end program
