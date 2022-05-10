
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

        domain%msl_linear = 0.0_dp

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
            gauges(1:2,i) = (island_radius + 1.01_dp*domain%dx(1)) *  [cos(theta), sin(theta)]
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
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum, flather_boundary
    use linear_interpolator_mod, only: linear_interpolator_type
    use local_routines
    implicit none

    type(multidomain_type) :: md

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 10.0_dp
    real(dp), parameter :: final_time = 200000.0_dp

    ! Timestepping method should be linear, as we compare against a linear analytical solution.
    character(charlen) :: timestepping_method = default_linear_timestepping_method
    
    ! Local variables 
    real(dp) :: timestep, boundary_wave_period
    real(dp) :: offshore_depth, island_radius, slope_radius, dx
    character(charlen):: test_case, bc_file
    integer(ip) :: j

    ! Write the stage raster time-series less often than we write at gauges, to avoid
    ! overly large files. 
    integer(ip), parameter :: frequency_full_write_steps = 300
    logical, parameter :: never_write_grid_time_slices = .false.

    ! Zero the max stage record after this much time has elapsed. The
    ! idea is to allow the transients to pass, then reset the max stage,
    ! so it ultimately only records the 'stationary' part of the run
    real(dp), parameter :: reset_max_stage_at_time = 120000.0_dp
    logical :: have_reset_max_stage = .FALSE.

    ! Scale this to refine the mesh -- e.g. value of 2 will halve the grid cell-size and time-step.
    real(dp), parameter :: mesh_refine = 1.0_dp 

    integer, parameter :: nd = 2
    

    ! Large scale domain. Wave comes from east side
    if(nd == 1) then
        ! Single domain model
        allocate(md%domains(1))
        dx = 2000.00_dp/mesh_refine
        md%domains(1)%lw = [2000.0_dp , 3000.0_dp]*1e+03
        md%domains(1)%lower_left = -md%domains(1)%lw/2.0_dp
        md%domains(1)%nx = nint(md%domains(1)%lw/dx)
        md%domains(1)%timestepping_method = timestepping_method

    else if(nd == 2) then
        ! Nested domain model

        ! Coarser outer domain
        allocate(md%domains(2))
        dx = 4000.00_dp/mesh_refine
        md%domains(1)%lw = [2000.0_dp , 3000.0_dp]*1e+03 
        md%domains(1)%lower_left = -md%domains(1)%lw/2.0_dp
        md%domains(1)%nx = nint(md%domains(1)%lw/dx)
        md%domains(1)%timestepping_method = timestepping_method

        ! Finer nested domain
        call md%domains(2)%match_geometry_to_parent(&
            parent_domain = md%domains(1), &
            lower_left  = [-4.0e+05_dp, -4.0e+05_dp], &
            upper_right = [ 4.0e+05_dp,  4.0e+05_dp], &
            dx_refinement_factor = 3_ip, &
            ! For inner linear domains, timestepping refinement should not be used
            ! (leads to instability).
            timestepping_refinement_factor = 1_ip)
        md%domains(2)%timestepping_method = timestepping_method

    else
        stop 'unsupported nesting setup'
    end if

    call md%setup

    ! Geometry
    island_radius = 40000.0_dp
    slope_radius = 160000.0_dp
    offshore_depth = 300.0_dp
    
    ! Incoming wave
    ! wavelength = 2 x slope_radius
    boundary_wave_period = slope_radius * 2.0_dp /sqrt(9.8_dp * offshore_depth) !12.0_dp * 60.0_dp

    do j = 1, size(md%domains)
        ! Call local routine to set initial conditions
        call set_initial_conditions_circular_island(md%domains(j), offshore_depth, island_radius, slope_radius)
    end do

    ! Get the boundary data and make an interpolation function f(t) for gauge 4
    boundary_information%offshore_elev = -offshore_depth
    boundary_information%boundary_wave_period = boundary_wave_period
    do j = 1, size(md%domains)
        if(any(md%domains(j)%boundary_exterior)) then
            md%domains(j)%boundary_function => boundary_function
            md%domains(j)%boundary_subroutine => flather_boundary
        end if
    end do

    call md%make_initial_conditions_consistent
    call md%set_null_regions_to_dry
   
    ! Use the min timestep on the multidomain (accounting for the fact that wave heights grow somewhat)
    timestep = md%stationary_timestep_max() * offshore_depth/(10.0_dp + offshore_depth)

    ! Evolve the code
    do while (.true.)

        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = frequency_full_write_steps, &
            print_less_often = frequency_full_write_steps)

        ! We are interested in the 'periodic steady-state' solution. To avoid the initial model transients,
        ! reset this after a burn in time
        if(md%domains(1)%time > reset_max_stage_at_time .and. (.not. have_reset_max_stage)) then
            do j = 1, size(md%domains)
                md%domains(j)%max_U = -HUGE(1.0_dp)
            end do
            have_reset_max_stage = .true. 
        end if

        ! End of run
        if(md%domains(1)%time > final_time) exit

        call md%evolve_one_step(timestep)

    end do

    call md%finalise_and_print_timers

end program
