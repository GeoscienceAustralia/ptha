MODULE local_routines 
    USE global_mod, only: dp, ip, wall_elevation
    USE domain_mod, only: domain_type
    USE which_mod, only: which
    IMPLICIT NONE
    
    REAL(dp), PARAMETER:: bed_slope = 0.01_dp, channel_width=20.0_dp, channel_depth=2.0_dp


    CONTAINS 

    SUBROUTINE set_initial_conditions_uniform_channel(domain)            
        CLASS(domain_type), TARGET, INTENT(INOUT):: domain
        INTEGER(ip):: i,j
        REAL(dp):: x, y, cx, cy, slope

        ! Set elevation
        slope = bed_slope 

        cx = (domain%lw(1))*0.5
        cy = (domain%lw(2))*0.5
        DO i = 1,domain%nx(1)
            DO j = 1, domain%nx(2)
                x = (i-0.5_dp)*domain%dx(1) - cx
                y = (j-0.5_dp)*domain%dx(2) - cy 
                domain%U(i,j,4) = y*slope
                IF(x > -channel_width/2.0_dp .AND. x < channel_width/2.0_dp) THEN
                    domain%U(i,j,4) = domain%U(i,j,4) - channel_depth 
                END IF
            END DO
        END DO
        ! Wall boundaries along the sides and top
        domain%U(1,:,4) = wall_elevation 
        domain%U(domain%nx(1),:,4) = wall_elevation 
        domain%U(:,domain%nx(2),4) = wall_elevation 
        domain%U(:,1,4) = wall_elevation 

        ! Stage
        domain%U(:,:,1) = domain%U(:,:,4)

        ! Ensure stage >= elevation
        domain%U(:,:,1 ) = max(domain%U(:,:,1), domain%U(:,:,4))


        domain%manning_squared = 0.03_dp**2

    END SUBROUTINE


END MODULE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM uniform_channel
    USE global_mod, only: ip, dp, charlen
    USE domain_mod, only: domain_type
    USE file_io_mod, only: read_csv_into_array
    USE local_routines
    IMPLICIT NONE

    INTEGER(ip):: i, j
    REAL(dp):: last_write_time
    TYPE(domain_type):: domain

    ! Approx timestep between outputs
    REAL(dp), PARAMETER :: approximate_writeout_frequency = 60.0_dp
    REAL(dp), PARAMETER :: final_time = 600.0_dp

    ! Length/width
    REAL(dp), PARAMETER, DIMENSION(2):: global_lw = [50._dp, 1000._dp] 
    ! Lower-left corner coordinate
    REAL(dp), PARAMETER, DIMENSION(2):: global_ll = [0._dp, 0._dp]
    ! grid size (number of x/y cells)
    INTEGER(ip), PARAMETER, DIMENSION(2):: global_nx = [50, 1000] ! [400, 400] 

    ! Discharge of 1 cubic meter / second per metre width
    REAL(dp) :: discharge_per_unit_width = 1.0_dp
    REAL(dp) :: Qin, model_vd, model_d, model_v, theoretical_d, theoretical_v, theoretical_vol, model_vol

    INTEGER(ip), ALLOCATABLE:: channel_x_indices(:)

    ! This gives the inflow discharge per cell, where it is applied
    Qin = discharge_per_unit_width * (global_lw(1)/global_nx(1))

    domain%theta = 1.0_dp
    domain%timestepping_method = 'rk2'

    ! Prevent dry domain from causing enormous initial step
    domain%maximum_timestep = 5.0_dp

    ! Allocate domain
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_uniform_channel(domain)

    ! Find channel indices 
    ! This corresponds to the lowest points on a horizontal cross-section
    CALL which(domain%U(:,2,4) < (minval(domain%U(:,2,4)) + 1.0e-06), channel_x_indices)
    print*, channel_x_indices


    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    ! Evolve the code
    DO WHILE (.TRUE.)

        IF(domain%time - last_write_time >= approximate_writeout_frequency) THEN

            last_write_time = last_write_time + approximate_writeout_frequency

            CALL domain%print()
            CALL domain%write_to_output_files()

            IF (domain%time > final_time) THEN
                EXIT 
            END IF

        END IF

        CALL domain%evolve_one_step()

        ! Add discharge inflow
        domain%U(channel_x_indices, domain%nx(2)-1, 1) = &
            domain%U(channel_x_indices, domain%nx(2)-1, 1) + Qin*domain%evolve_step_dt/product(domain%dx)

    END DO

    call domain%write_max_quantities()

    call domain%timer%print()


    theoretical_vol = Qin * domain%time * channel_width
    model_vol = sum(domain%U(:,:,1) - domain%U(:,:,4)) * domain%dx(1) * domain%dx(2)

    ! Analytical solution for uniform slope
    ! vd = discharge_per_unit_width
    ! bed_slope = manning_squared * v * abs(v) / d^(4/3)
    ! bed_slope = manning_squared * discharge_per_unit_width**2 / (d^(10/3))
    ! d = (manning_squared * discharge_per_unit_width**2 / bed_slope )**(3/10)

    i = domain%nx(1)/2
    j = domain%nx(2)/2    
    model_vd = -1.0_dp * domain%U(i,j,3)    
    model_d = domain%U(i,j,1) - domain%U(i,j,4)
    model_v = model_vd/model_d

    theoretical_d = (domain%manning_squared(i,j) * discharge_per_unit_width**2 / bed_slope)**(3.0_dp/10.0_dp)

    print*, '## Testing ##'
    print*, ''
    print*, '    Theoretical volume: ', theoretical_vol
    print*, '    Model volume: ', model_vol
    print*, '      error: ', theoretical_vol - model_vol
    IF( abs(theoretical_vol - model_vol) < 1.0e-06_dp * theoretical_vol) THEN
        print*, 'PASS'
    ELSE
        stop 'FAIL: Mass conservation error'
    END IF
    

    print*, ' '
    print*, '   Theoretical vd: ', discharge_per_unit_width
    print*, '   Model vd (near centre): ', model_vd
    IF( abs(discharge_per_unit_width - model_vd) < 1.0e-03_dp*discharge_per_unit_width) THEN
        print*, 'PASS'
    ELSE
        stop 'FAIL: Flux error'
    END IF

    print*, ' '
    print*, '   Theoretical d: ', theoretical_d
    print*, '   Model d: ', model_d
    IF( abs(theoretical_d - model_d) < 1.0e-02_dp*theoretical_d) THEN
        print*, 'PASS'
    ELSE
        stop 'FAIL: Friction or steady flow error'
    END IF

    CALL domain%finalise()    

END PROGRAM
