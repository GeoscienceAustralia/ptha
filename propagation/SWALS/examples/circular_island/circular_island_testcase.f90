
module local_routines 
    !! Define initial and boundary conditions

    use global_mod, only: dp, ip, charlen, wall_elevation, pi, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: read_gdal_raster
    use which_mod, only: which
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    implicit none

    type :: boundary_information_type
        !! Hold some data used by the boundary condition. We can set this from
        !! inside the main program.
        real(dp) :: offshore_elev
        real(dp) :: boundary_wave_period
    end type

    type(boundary_information_type), public :: boundary_information
        !! The main program will modify this type to set up the boundary condition

    contains 

    pure function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        !! Make a function to evaluate the boundary at the domain. We will use
        !! this in conjunction with a flather type radiation condition

        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4), wavelength, waveperiod, theta


        ! Specify a wave, corresponding to the far-field condition
        waveperiod = boundary_information%boundary_wave_period
        wavelength = (sqrt(-gravity*boundary_information%offshore_elev)*waveperiod)

        ! Make the wave come from the east
        theta = max((t/waveperiod  + (domain%x(i)-domain%x(domain%nx(1)))/wavelength), 0.0_dp)
        stage_uh_vh_elev(STG) = 1.0_dp * sin(2*pi* theta)
        stage_uh_vh_elev(UH) = -wavelength/waveperiod * stage_uh_vh_elev(STG)
        stage_uh_vh_elev(VH) = 0.0_dp
        stage_uh_vh_elev(ELV) = boundary_information%offshore_elev

    end function

    subroutine set_initial_conditions_circular_island(domain, offshore_depth, island_radius, slope_radius)
        !! Setup initial conditions + locations of gauges

        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in) :: offshore_depth, island_radius, slope_radius

        ! Local variables
        real(dp) :: x, y, elev, slope, radius, theta
        real(dp), allocatable :: gauges(:,:)
        integer(ip) :: i, j

        ! Stage, UH, VH
        domain%U(:,:,[STG, UH, VH]) = 0.0_dp

        ! Set elevation
        slope = offshore_depth/(slope_radius - island_radius)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                x = domain%x(i)
                y = domain%y(j)
                radius = sqrt(x*x + y*y)
                elev = -offshore_depth + max(slope_radius - radius, 0.0_dp) *  slope
                domain%U(i,j,ELV) = elev
            end do
        end do
        ! No negative depth!
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Set the tide gauges 
        allocate(gauges(2, 150)) 

        ! Gauges around the island (slightly increase the radius to avoid points on dry cells)
        do i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2,i) = (island_radius + domain%dx(1)) *  [cos(theta), sin(theta)]
        end do

        ! Gauges a bit further offshore
        do i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2, i+50) = (island_radius + (slope_radius - island_radius)/2.0_dp) *  [cos(theta), sin(theta)]
        end do

        ! Gauges yet further offshore 
        do i = 1, 50
            theta = 2*pi/50.0_dp * (i-1)
            gauges(1:2, i+100) = (slope_radius) *  [cos(theta), sin(theta)]
        end do

        call domain%setup_point_gauges(gauges)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program circular_island
    !! Plane wave scattering around a circular conical island.
    !! This case has an analytical solution, due to: 
    !! Zhang and Zhu (1994) New solutions for the propagation of long waves over
    !! variable depth. Journal of fluid mechanics 278: 391-406

    use global_mod, only: ip, dp, minimum_allowed_depth, default_linear_timestepping_method
    use domain_mod, only: domain_type
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum, flather_boundary
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    integer(ip):: j
    real(dp):: last_write_time, rain_rate
    type(domain_type):: domain

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 50.0_dp
    real(dp), parameter :: final_time = 200000.0_dp

    ! Domain info
    character(charlen) :: timestepping_method = default_linear_timestepping_method
    
    !! length/width
    real(dp), dimension(2) :: global_lw, global_ll
    integer(ip), dimension(2) :: global_nx 

    ! Local variables 
    real(dp) :: timestep
    real(dp) :: offshore_depth, island_radius, slope_radius, dx
    real(dp) :: boundary_wave_period
    character(charlen):: test_case, bc_file
    real(dp):: tank_bases(4), tank_slopes(4) 
    integer(ip) :: full_write_step

    ! Write the stage raster time-series less often than we write at gauges, to avoid
    ! overly large files. 
    integer(ip), parameter :: frequency_full_write_steps = 60
    logical, parameter :: never_write_grid_time_slices = .false.

    ! Zero the max stage record after this much time has elapsed. The
    ! idea is to allow the transients to pass, then reset the max stage,
    ! so it ultimately only records the 'stationary' part of the run
    real(dp), parameter :: reset_max_stage_at_time = 120000.0_dp

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
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    call set_initial_conditions_circular_island(domain, offshore_depth, island_radius, slope_radius)

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    domain%boundary_function => boundary_function
    domain%boundary_subroutine => flather_boundary
    boundary_information%offshore_elev = -offshore_depth
    boundary_information%boundary_wave_period = boundary_wave_period

   
    ! Linear requires a fixed timestep 
    if (.not. domain%adaptive_timestepping) then
        timestep = domain%stationary_timestep_max() 
    end if


    ! Trick to get the code to write out just after the first timestep
    last_write_time = domain%time - approximate_writeout_frequency

    full_write_step = 0
    ! Evolve the code
    do while (.true.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            ! Reset the peak stage record once during the simulation, so the final
            ! result reflects the 'stationary' model solution
            if((domain%time >= reset_max_stage_at_time) .and. &
               (last_write_time <= reset_max_stage_at_time)) then
                domain%max_U(:,:,1) = -huge(1.0_dp)
            end if

            last_write_time = last_write_time + approximate_writeout_frequency

            ! This avoids any artefacts in the numerical update of the model
            ! which should be overwritten by the boundary condition
            call domain%update_boundary()

            full_write_step = full_write_step + 1
            if(mod(full_write_step, frequency_full_write_steps) == 0) then 
                call domain%write_to_output_files(time_only=never_write_grid_time_slices)
                call domain%print()
                print*, 'Mass balance: ', domain%mass_balance_interior()
            end if

            call domain%write_gauge_time_series()

        end if

        if (domain%time > final_time) exit 

        if (.not. domain%adaptive_timestepping) then
            call domain%evolve_one_step(timestep = timestep)
        else
            call domain%evolve_one_step()
        end if

    end do

    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print()

    call domain%finalise()

end program
