module local_routines 

    use global_mod, only: dp, ip, charlen, wall_elevation
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
    
    ! Make a function to evaluate the boundary at the domain
    !
    function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4)
        call boundary_information%gauge4_ts_function%eval([t + boundary_information%t0], stage_uh_vh_elev(1:1))
        stage_uh_vh_elev(2:3) = 0.0_dp
        stage_uh_vh_elev(4) = boundary_information%boundary_elev
    end function

    subroutine set_initial_conditions_bp2(domain, tank_bases, tank_slopes, tank_width, initial_depth)
        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: tank_bases(4), tank_slopes(4), tank_width, initial_depth

        real(dp):: tank_x(5)
        real(dp):: x, y, elev
        integer(ip):: j, i, k
        real(dp):: gauge_xy(2,11), wall

        ! Set stage
        domain%U(:,:,STG) = 0.0_dp

        ! Set elevation
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
        wall = 0.5_dp
        domain%U(:, 1, ELV) = wall
        domain%U(:, 2, ELV) = wall
        domain%U(:, domain%nx(2), ELV) = wall
        domain%U(:, domain%nx(2)-1, ELV) = wall
        domain%U(domain%nx(1), :, ELV) = wall
        domain%U(domain%nx(1)-1, :, ELV) = wall
    
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))
        
        ! Get gauge points
        gauge_xy(2,:) = 0.0_dp
        gauge_xy(1, 1:3) = 0.0_dp
        gauge_xy(1, 4) = 0.0_dp
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program bp2
    use global_mod, only: ip, dp, minimum_allowed_depth
    use domain_mod, only: domain_type
    !use boundary_mod, only: boundary_stage_transmissive_normal_momentum
    use boundary_mod, only: boundary_stage_transmissive_momentum
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    integer(ip):: j
    real(dp):: last_write_time, rain_rate
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.1_dp
    real(dp), parameter :: final_time = 30.0_dp

    ! Domain info
    character(charlen) :: timestepping_method != 'linear' !'rk2' !'linear'
    
    !! length/width
    real(dp), dimension(2) :: global_lw, global_ll
    integer(ip), dimension(2) :: global_nx 

    ! Local variables 
    real(dp) :: timestep, base_l, dx, tank_width, tank_length, initial_depth
    character(charlen):: test_case, bc_file
    real(dp):: tank_bases(4), tank_slopes(4) 


    ! Get the case. Values should be caseA, caseB, caseC
    call get_command_argument(1, test_case)
    select case(test_case)
        case('caseA')
            base_L = 2.40_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3a_analytical.txt'

        case('caseB')
            base_L = 0.98_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3b_analytical.txt'

        case('caseC')
            base_L = 0.64_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3c_analytical.txt'
        case default
            print*, 'Must specify a test case (one of caseA, caseB, caseC)'
            stop
    end select

    call get_command_argument(2, timestepping_method)  
    
    ! Resolution
    dx = 0.01_dp
 
    ! Tank geometry  -- add a little extra at the end so the reflective wall is in the right place
    tank_bases = [base_L, 4.36_dp, 2.93_dp, 0.9_dp + 2.0_dp*dx]
    tank_slopes = [0.0_dp, 1.0_dp/53.0_dp, 1.0_dp/150.0_dp, 1.0_dp/13.0_dp]
    tank_width = 1.0_dp
    tank_length = sum(tank_bases)
    initial_depth = 0.218_dp
    

    ! Large scale
    global_lw = [tank_length, tank_width]
    global_ll = [0.0_dp, -tank_width/2.0_dp]
    global_nx = global_lw/dx

    domain%timestepping_method = timestepping_method

    ! Allocate domain -- must have set timestepping method BEFORE this
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    call set_initial_conditions_BP2(domain, tank_bases, tank_slopes, tank_width, initial_depth)

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    call setup_boundary_information(bc_file, -initial_depth)
    domain%boundary_function => boundary_function
    !domain%boundary_subroutine => boundary_stage_transmissive_normal_momentum
    domain%boundary_subroutine => boundary_stage_transmissive_momentum

   
    ! Linear requires a fixed timestep 
    if (.not. domain%adaptive_timestepping) then
        timestep = domain%stationary_timestep_max() * 0.5_dp
    end if


    ! Trick to get the code to write out just after the first timestep
    last_write_time = domain%time - approximate_writeout_frequency

    ! Evolve the code
    do while (.TRUE.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files(time_only=.true.)
            call domain%write_gauge_time_series()
            print*, 'Mass balance: ', domain%mass_balance_interior()

        end if

        if (domain%time > final_time) exit

        ! Suggested to use a transmissive type boundary at this stage
        !IF(domain%time > 10.0_dp) THEN
        !    domain%boundary_type = 'flather_still_water' !'transmissive'
        !    domain%boundary_function => NULL()
        !END IF
        !! Example with fixed timestep
        !CALL domain%evolve_one_step(timestep=2.5_dp)

        ! Variable timestep
        if(.not. domain%adaptive_timestepping) then
            call domain%evolve_one_step(timestep = timestep)
        else
            call domain%evolve_one_step()
        end if

    END DO

    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print()

    call domain%finalise()

end program
