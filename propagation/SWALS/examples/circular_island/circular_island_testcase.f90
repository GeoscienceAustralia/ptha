!
! Plane wave scattering around a circular conical island.
! This case has an analytical solution, due to: 
! Zhang and Zhu (1994) New solutions for the propagation of long waves over
! variable depth. Journal of fluid mechanics 278: 391-406
!


! Module used to define initial and boundary conditions
MODULE local_routines 

    USE global_mod, only: dp, ip, charlen, wall_elevation, pi, gravity
    USE domain_mod, only: domain_type, STG, UH, VH, ELV
    USE read_raster_mod, only: read_gdal_raster
    USE which_mod, only: which
    USE file_io_mod, only: count_file_lines
    USE linear_interpolator_mod, only: linear_interpolator_type
    IMPLICIT NONE

    ! Hold some data used by the boundary condition. We can set this from
    ! inside the main program.
    TYPE :: boundary_information_type
        REAL(dp) :: offshore_elev
        REAL(dp) :: boundary_wave_period
    END TYPE

    ! The main program will modify this type to set up the boundary condition
    TYPE(boundary_information_type), PUBLIC :: boundary_information

    CONTAINS 

    !
    ! Make a function to evaluate the boundary at the domain. We will use
    ! this in conjunction with a flather type radiation condition
    FUNCTION boundary_function(domain, t, x, y) result(stage_uh_vh_elev)
        TYPE(domain_type), INTENT(IN):: domain
        REAL(dp), INTENT(IN):: t, x, y
        REAL(dp):: stage_uh_vh_elev(4), wavelength, waveperiod, theta


        ! Specify a wave, corresponding to the far-field condition
        waveperiod = boundary_information%boundary_wave_period
        wavelength = (sqrt(-gravity*boundary_information%offshore_elev)*waveperiod)

        ! Make the wave come from the east
        theta = max((t/waveperiod  + (x-domain%x(domain%nx(1)))/wavelength), 0.0_dp)
        stage_uh_vh_elev(STG) = 1.0_dp * sin(2*pi* theta)
        stage_uh_vh_elev(UH) = -wavelength/waveperiod * stage_uh_vh_elev(STG)
        stage_uh_vh_elev(VH) = 0.0_dp
        stage_uh_vh_elev(ELV) = boundary_information%offshore_elev

    END FUNCTION

    ! Initial conditions + locations of gauges
    SUBROUTINE set_initial_conditions_circular_island(domain, offshore_depth, island_radius, slope_radius)
        CLASS(domain_type), TARGET, INTENT(INOUT):: domain
        REAL(dp), INTENT(IN) :: offshore_depth, island_radius, slope_radius

        ! Local variables
        REAL(dp) :: x, y, elev, slope, radius, theta
        REAL(dp), ALLOCATABLE :: gauges(:,:)
        INTEGER(ip) :: i, j

        ! Stage, UH, VH
        domain%U(:,:,[STG, UH, VH]) = 0.0_dp

        ! Set elevation
        slope = offshore_depth/(slope_radius - island_radius)
        DO j = 1, domain%nx(2)
            DO i = 1, domain%nx(1)
                x = domain%x(i)
                y = domain%y(j)
                radius = sqrt(x*x + y*y)
                elev = -offshore_depth + max(slope_radius - radius, 0.0_dp) *  slope
                domain%U(i,j,ELV) = elev
            END DO
        END DO
        ! No negative depth!
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Set the tide gauges 
        ALLOCATE(gauges(2, 150)) 

        ! Gauges around the island (slightly increase the radius to avoid points on dry cells)
        DO i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2,i) = (island_radius + domain%dx(1)) *  [cos(theta), sin(theta)]
        END DO

        ! Gauges a bit further offshore
        DO i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2, i+50) = (island_radius + (slope_radius - island_radius)/2.0_dp) *  [cos(theta), sin(theta)]
        END DO

        ! Gauges yet further offshore 
        DO i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2, i+100) = (slope_radius) *  [cos(theta), sin(theta)]
        END DO

        call domain%setup_point_gauges(gauges)

    END SUBROUTINE

END MODULE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM circular_island
    USE global_mod, only: ip, dp, minimum_allowed_depth
    USE domain_mod, only: domain_type
    USE boundary_mod, only: boundary_stage_transmissive_normal_momentum, flather_boundary
    USE linear_interpolator_mod, only: linear_interpolator_type
    USE local_routines
    IMPLICIT NONE

    INTEGER(ip):: j
    REAL(dp):: last_write_time, rain_rate
    TYPE(domain_type):: domain

    ! Approx timestep between outputs
    REAL(dp), PARAMETER :: approximate_writeout_frequency = 50.0_dp
    REAL(dp), PARAMETER :: final_time = 200000.0_dp

    ! Domain info
    CHARACTER(charlen) :: timestepping_method = 'linear' 
    
    !! length/width
    REAL(dp), DIMENSION(2) :: global_lw, global_ll
    INTEGER(ip), DIMENSION(2) :: global_nx 

    ! Local variables 
    REAL(dp) :: timestep
    REAL(dp) :: offshore_depth, island_radius, slope_radius, dx
    REAL(dp) :: boundary_wave_period
    CHARACTER(charlen):: test_case, bc_file
    REAL(dp):: tank_bases(4), tank_slopes(4) 
    INTEGER(ip) :: full_write_step

    ! Write the stage raster time-series less often than we write at gauges, to avoid
    ! overly large files. 
    INTEGER(ip), PARAMETER :: frequency_full_write_steps = 60
    LOGICAL, PARAMETER :: never_write_grid_time_slices = .true.

    ! Zero the max stage record after this much time has elapsed. The
    ! idea is to allow the transients to pass, then reset the max stage,
    ! so it ultimately only records the 'stationary' part of the run
    REAL(dp), PARAMETER :: reset_max_stage_at_time = 120000.0_dp

    ! Resolution
    dx = 2000.00_dp

 
    ! Large scale domain. Wave comes from east side
    global_lw = [2000.0_dp , 3000.0_dp]*1e+03 
    global_ll = -global_lw/2.0_dp
    global_nx = nint(global_lw/dx)

    ! Geometry
    island_radius = 40000.0_dp
    slope_radius = 160000.0_dp
    offshore_depth = 300.0_dp
    
    ! Incoming wave
    ! wavelength = 2 x slope_radius
    boundary_wave_period = slope_radius * 2.0_dp /sqrt(9.8 * offshore_depth) !12.0_dp * 60.0_dp


    domain%timestepping_method = timestepping_method

    ! Allocate domain -- must have set timestepping method BEFORE this
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_circular_island(domain, offshore_depth, island_radius, slope_radius)

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    domain%boundary_function => boundary_function
    domain%boundary_subroutine => flather_boundary
    boundary_information%offshore_elev = -offshore_depth
    boundary_information%boundary_wave_period = boundary_wave_period

   
    ! Linear requires a fixed timestep 
    if (timestepping_method == 'linear') then
        domain%cfl = 0.7_dp
        timestep = domain%linear_timestep_max() 
    end if


    ! Trick to get the code to write out just after the first timestep
    last_write_time = domain%time - approximate_writeout_frequency

    full_write_step = 0
    ! Evolve the code
    DO WHILE (.TRUE.)

        IF(domain%time - last_write_time >= approximate_writeout_frequency) THEN

            ! Reset the peak stage record once during the simulation, so the final
            ! result reflects the 'stationary' model solution
            if((domain%time >= reset_max_stage_at_time) .AND. &
               (last_write_time <= reset_max_stage_at_time)) then
                domain%max_U(:,:,1) = -HUGE(1.0_dp)
            end if

            last_write_time = last_write_time + approximate_writeout_frequency

            ! This avoids any artefacts in the numerical update of the model
            ! which should be overwritten by the boundary condition
            CALL domain%update_boundary()

            CALL domain%print()

            full_write_step = full_write_step + 1
            if(mod(full_write_step, frequency_full_write_steps) == 0) then 
                CALL domain%write_to_output_files(time_only=never_write_grid_time_slices)
            end if

            CALL domain%write_gauge_time_series()
            print*, 'Mass balance: ', domain%mass_balance_interior()

        END IF

        IF (domain%time > final_time) THEN
            EXIT 
        END IF

        CALL domain%evolve_one_step(timestep = timestep)

    END DO

    CALL domain%write_max_quantities()

    ! Print timing info
    CALL domain%timer%print()

    CALL domain%finalise()

END PROGRAM
