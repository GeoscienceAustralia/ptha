MODULE local_routines 
    USE global_mod, only: dp, ip, charlen, wall_elevation
    USE domain_mod, only: domain_type
    USE read_raster_mod, only: read_gdal_raster
    USE which_mod, only: which
    USE file_io_mod, only: count_file_lines
    USE linear_interpolator_mod, only: linear_interpolator_type
    IMPLICIT NONE

    ! Hold some data for the boundary condition
    TYPE :: boundary_information_type
        CHARACTER(charlen):: bc_file
        REAL(dp), ALLOCATABLE:: boundary_data(:,:)
        TYPE(linear_interpolator_type):: gauge4_ts_function
        REAL(dp):: boundary_elev
        REAL(dp):: t0
    END TYPE

    ! This will hold the information -- is seen by other parts of the module
    TYPE(boundary_information_type):: boundary_information

    CONTAINS 

    SUBROUTINE setup_boundary_information(bc_file, boundary_elev)
        CHARACTER(charlen), INTENT(IN):: bc_file
        REAL(dp), INTENT(IN):: boundary_elev

        INTEGER(ip):: bc_unit, nr, nc, skip, i

        boundary_information%bc_file = bc_file
        boundary_information%boundary_elev = boundary_elev
        OPEN(newunit=bc_unit, file=bc_file)
        nr = count_file_lines(bc_unit)
        nc = 9
        skip = 6
        ALLOCATE(boundary_information%boundary_data(nr - skip, nc))
        DO i = 1, nr
            if(i > skip) then
                READ(bc_unit, *) boundary_information%boundary_data(i - skip,:)
            else
                READ(bc_unit, *) 
            end if 
        END DO
        CLOSE(bc_unit)

        CALL boundary_information%gauge4_ts_function%initialise(&
                boundary_information%boundary_data(:,1), boundary_information%boundary_data(:,2))
        boundary_information%t0 = boundary_information%boundary_data(1,1)

    END SUBROUTINE
    
    ! Make a function to evaluate the boundary at the domain
    !
    FUNCTION boundary_function(domain, t, x, y) result(stage_uh_vh_elev)
        TYPE(domain_type), INTENT(IN):: domain
        REAL(dp), INTENT(IN):: t, x, y
        REAL(dp):: stage_uh_vh_elev(4)
        CALL boundary_information%gauge4_ts_function%eval([t + boundary_information%t0], stage_uh_vh_elev(1:1))
        stage_uh_vh_elev(2:3) = 0.0_dp
        stage_uh_vh_elev(4) = boundary_information%boundary_elev
    END FUNCTION

    SUBROUTINE set_initial_conditions_BP2(domain, tank_bases, tank_slopes, tank_width, initial_depth)
        CLASS(domain_type), TARGET, INTENT(INOUT):: domain
        REAL(dp), INTENT(IN) :: tank_bases(4), tank_slopes(4), tank_width, initial_depth

        REAL(dp):: tank_x(5)
        REAL(dp):: x, y, elev
        INTEGER(ip):: j, i, k
        REAL(dp):: gauge_xy(2,11), wall

        ! Set stage
        domain%U(:,:,1:3) = 0.0_dp

        ! Set elevation
        tank_x = 0.0_dp
        DO i = 2, 5
            tank_x(i) = sum(tank_bases(1:(i-1)))            
        END DO
        DO j = 1, domain%nx(2)
            DO i = 1, domain%nx(1)
                x = domain%lower_left(1) + (i-0.5_dp) * domain%dx(1) 
                elev = -initial_depth
                DO k = 2, 4
                     if(x > tank_x(k)) then
                        elev = elev + (min(x, tank_x(k+1)) - tank_x(k))*tank_slopes(k)
                     end if
                END DO
            domain%U(i,j,4) = elev
            END DO
        END DO
      
        ! Reflective boundaries on 3 sides
        wall = 0.5_dp
        domain%U(:, 1, 4) = wall
        domain%U(:, 2, 4) = wall
        domain%U(:, domain%nx(2), 4) = wall
        domain%U(:, domain%nx(2)-1, 4) = wall
        domain%U(domain%nx(1), :, 4) = wall
        domain%U(domain%nx(1)-1, :, 4) = wall
    
        domain%U(:,:,1) = max(domain%U(:,:,1), domain%U(:,:,4))
        
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

    END SUBROUTINE

END MODULE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM BP2
    USE global_mod, only: ip, dp, minimum_allowed_depth
    USE domain_mod, only: domain_type
    USE boundary_mod, only: boundary_stage_transmissive_normal_momentum
    USE linear_interpolator_mod, only: linear_interpolator_type
    USE local_routines
    IMPLICIT NONE

    INTEGER(ip):: j
    REAL(dp):: last_write_time, rain_rate
    TYPE(domain_type):: domain

    ! Approx timestep between outputs
    REAL(dp), PARAMETER :: approximate_writeout_frequency = 0.1_dp
    REAL(dp), PARAMETER :: final_time = 30.0_dp

    ! Domain info
    CHARACTER(charlen) :: timestepping_method != 'linear' !'rk2' !'linear'
    
    !! length/width
    REAL(dp), DIMENSION(2) :: global_lw, global_ll
    INTEGER(ip), DIMENSION(2) :: global_nx 

    ! Local variables 
    REAL(dp) :: timestep, base_L, dx, tank_width, tank_length, initial_depth
    CHARACTER(charlen):: test_case, bc_file
    REAL(dp):: tank_bases(4), tank_slopes(4) 


    ! Get the case. Values should be caseA, caseB, caseC
    call get_command_argument(1, test_case)
    SELECT CASE(test_case)
        CASE('caseA')
            base_L = 2.40_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3a_analytical.txt'

        CASE('caseB')
            base_L = 0.98_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3b_analytical.txt'

        CASE('caseC')
            base_L = 0.64_dp
            bc_file = '../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3c_analytical.txt'
        CASE DEFAULT
            print*, 'Must specify a test case (one of caseA, caseB, caseC)'
            stop
    END SELECT

    call get_command_argument(2, timestepping_method)  
    
    ! Resolution
    dx = 0.02_dp
 
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
    if(timestepping_method == 'euler') domain%theta = 0.0_dp
    if(timestepping_method == 'rk2') domain%theta = 1.0_dp

    ! Allocate domain -- must have set timestepping method BEFORE this
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_BP2(domain, tank_bases, tank_slopes, tank_width, initial_depth)

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    CALL setup_boundary_information(bc_file, -initial_depth)
    domain%boundary_function => boundary_function
    domain%boundary_subroutine => boundary_stage_transmissive_normal_momentum

   
    ! Linear requires a fixed timestep 
    if (timestepping_method == 'linear') then
        domain%cfl = 0.9_dp
        timestep = domain%linear_timestep_max() 
    end if


    ! Trick to get the code to write out just after the first timestep
    last_write_time = domain%time - approximate_writeout_frequency

    ! Evolve the code
    DO WHILE (.TRUE.)

        IF(domain%time - last_write_time >= approximate_writeout_frequency) THEN

            last_write_time = last_write_time + approximate_writeout_frequency

            ! This avoids any artefacts in the numerical update of the model
            ! which should be overwritten by the boundary condition
            CALL domain%update_boundary()

            CALL domain%print()
            CALL domain%write_to_output_files(time_only=.true.)
            CALL domain%write_gauge_time_series()
            print*, 'Mass balance: ', domain%mass_balance_interior()

        END IF

        IF (domain%time > final_time) THEN
            EXIT 
        END IF

        ! Suggested to use a transmissive type boundary at this stage
        !IF(domain%time > 10.0_dp) THEN
        !    domain%boundary_type = 'flather_still_water' !'transmissive'
        !    domain%boundary_function => NULL()
        !END IF
        !! Example with fixed timestep
        !CALL domain%evolve_one_step(timestep=2.5_dp)

        ! Variable timestep
        if(timestepping_method == 'linear') then
            CALL domain%evolve_one_step(timestep = timestep)
        else
            CALL domain%evolve_one_step()
        end if

    END DO

    CALL domain%write_max_quantities()

    ! Print timing info
    CALL domain%timer%print()

    CALL domain%finalise()

END PROGRAM
