
module local_routines 
    !!
    !! NTHMP tsunami-currents benchmark problem, Seaside OSU Experiment
    !!
    use global_mod, only: dp, ip, charlen, gravity
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use file_io_mod, only: count_file_lines, read_csv_into_array
    use linear_interpolator_mod, only: linear_interpolator_type
    use read_raster_mod, only: multi_raster_type

    implicit none

    !
    ! Parameters defining the domain geometry.
    !
    real(dp), parameter :: initial_sea_level = 0.97_dp

    ! Lower-left corner coordinate based on input bathymetry (slightly inside),
    ! assuming boundary forcing at x=5m. 
    ! This forcing leads to a phase-lag in the model at gauges wg3/wg4, which 
    ! is also reported by other modelling groups.
    real(dp), parameter :: global_ll(2) = [5.0_dp, -13.25_dp] 
    character(len=charlen), parameter :: &
        gauge_boundary_file = './problem_data/ts_5m.txt'

    ! Here we apply a forcing using wg1, located at x=2.068m. This is an 
    ! alternative to using the ideal boundary forcing.
    ! However it still leads to a phase-lag at the model at gauges wg3/wg4
    ! real(dp), parameter :: global_ll(2) = [2.05_dp, -13.25_dp] 
    ! character(len=charlen), parameter :: &
    !      gauge_boundary_file = './problem_data/ts_from_observations_at_wg1.txt'

    ! Length/width = 'upper_right' - 'lower_left' + 'extra-for-boundary-trough'
    real(dp), parameter :: &
        global_lw(2) = [43.63_dp, 8.55_dp] - global_ll + [5.0_dp, 0.0_dp]


    !
    ! Hold some data for the forcing
    !
    type :: boundary_information_type

        character(len=charlen) :: gauge_stage_forcing_file = gauge_boundary_file
        type(linear_interpolator_type) :: gauge_stage_forcing
        real(dp):: t0 = 0.0_dp
        real(dp), allocatable:: stage_data(:,:)
        integer(ip) :: skip_lines = 0_ip

    end type

    ! Hold the data/interpolators to setup the boundary 
    type(boundary_information_type):: boundary_information

    contains 

    !
    ! Read files with wavemaker info
    !
    subroutine setup_boundary_information()

        integer(ip):: bc_unit, nr, nc, skip, i, extra, j

        !
        ! Stage time-series 
        !
        ! Read from file - first column is time, second is stage 
        print*, 'Reading boundary forcing ....'
        open(newunit=bc_unit, file=boundary_information%gauge_stage_forcing_file)
        nr = count_file_lines(bc_unit)
        nc = 2
        skip = boundary_information%skip_lines
        extra = 1 ! One more data point to avoid exceeding time 
        allocate(boundary_information%stage_data(nr - skip + extra, nc))
        do i = 1, nr
            if(i > skip) then
                read(bc_unit, *) boundary_information%stage_data(i - skip,:)
            else
                read(bc_unit, *) 
            end if 
        end do
        close(bc_unit)

        ! Extend the time-series with a constant value, so that time does not
        ! exceed model run time
        boundary_information%stage_data(    nr - skip + extra    , 1:2) = &
            boundary_information%stage_data(nr - skip + extra - 1, 1:2) + &
            [1.0e+06_dp, 0.0_dp]

        ! Convert the data to use the initial_sea_level
        boundary_information%stage_data(:,2) = &
            boundary_information%stage_data(:,2) + initial_sea_level

        ! Ensure the first time is zero.
        boundary_information%stage_data(1,1) = 0.0_dp

        ! Make the linear interpolation function for the correct forcing case
        call boundary_information%gauge_stage_forcing%initialise(&
                boundary_information%stage_data(:,1), &
                boundary_information%stage_data(:,2))

    end subroutine

    !
    ! Make a function to evaluate the boundary at the domain
    !
    function boundary_function(domain, t, i, j) result(stage_uh_vh_elev)
        type(domain_type), intent(in):: domain
        real(dp), intent(in):: t
        integer(ip), intent(in) :: i, j
        real(dp):: stage_uh_vh_elev(4)
        real(dp) :: local_elev

        call boundary_information%gauge_stage_forcing%eval(&
            [t + boundary_information%t0], stage_uh_vh_elev(1:1))
        
        local_elev = domain%U(i,j,4)
        if(local_elev > stage_uh_vh_elev(1)) then
            stage_uh_vh_elev(1) = local_elev
        end if
        ! Semi-transmissive boundary condition will provide uh/vh values
        !stage_uh_vh_elev(2:3) = 0.0_dp 

        ! With flather boundary we need to well approximate the UH value
        ! Estimate from plane wave equation
        !stage_uh_vh_elev(2) = (stage_uh_vh_elev(1) - initial_sea_level) * sqrt(gravity * initial_sea_level)
        !stage_uh_vh_elev(2) = (stage_uh_vh_elev(1) - initial_sea_level) * sqrt(gravity * (stage_uh_vh_elev(1) - local_elev))
        stage_uh_vh_elev(2) = (stage_uh_vh_elev(1) - initial_sea_level) * sqrt(gravity/(initial_sea_level-local_elev)) * &
            (stage_uh_vh_elev(1) - local_elev)

        stage_uh_vh_elev(3) = 0.0_dp

        stage_uh_vh_elev(4) = local_elev

    end function
   
    ! Main geometry setup routine
    subroutine set_initial_conditions(domain)
        class(domain_type), target, intent(inout):: domain

        type(multi_raster_type):: elevation_data
        real(dp), allocatable:: x(:), y(:)
        integer(ip):: j
        
        ! Stage
        domain%U(:,:,STG) = initial_sea_level
        
        ! Elevation (background)
        domain%U(:,:,ELV) = 0.0_dp

        !
        ! Setup elevation
        !
        call elevation_data%initialise(&
            [character(len=charlen) :: './problem_data/bathy_raster.tif'])
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            ! To resolve the buildings, better to NOT use bilinear interpolation.
            call elevation_data%get_xy(x,y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=0_ip)

            ! Add a 'trough' to catch water at the right hand side of the 
            ! domain, and work around the boundary reflection. Park et al. 
            ! (2013) mention their model suffers this boundary reflection, 
            ! which contaminates the comparison with data at later times. It's 
            ! not clear how the experimental setup differs from their model.
            where(domain%U(:,j,ELV) < -100.0_dp )
                domain%U(:,j,ELV) = 0.7_dp
                domain%U(:,j,STG) = 0.7_dp
            end where

            ! Add walls on 3 boundaries -- 20cm thick
            if(domain%y(j) < global_ll(2) + 0.2_dp) &
                domain%U(:,j,ELV) = 10.0_dp
            if(domain%y(j) > global_ll(2) + global_lw(2) - 0.2_dp) &
                domain%U(:,j,ELV) = 10.0_dp
            where(x > global_ll(1) + global_lw(1) - 0.2_dp) &
                    domain%U(:,j,ELV) = 10.0_dp
        end do
        deallocate(x,y)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), &
            maxval(domain%U(:,:,ELV))

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        ! ( Manning coefficient )^2
        domain%manning_squared = 0.005_dp**2

        domain%msl_linear = initial_sea_level

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Seaside_OSU
    !!
    !! NTHMP tsunami-currents benchmark problem, Seaside OSU Experiment
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, &
        default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type, setup_multidomain
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum, boundary_stage_radiation_momentum, flather_boundary
    use timer_mod, only: timer_type
    use logging_mod, only: log_output_unit
    use local_routines
    implicit none

    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 1.0_dp ! Factor to increase resolution 
    
    real(dp), parameter ::  global_dt = 0.0125_dp / mesh_refine 

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 0.05_dp
    real(dp), parameter :: final_time = 40._dp
    integer(ip), parameter :: print_every_nth_step = 4_ip, &
        write_grids_every_nth_step = 4_ip

    ! Most geometric parameters are defined in the 'local routines' module.
    ! Grid size (number of x/y cells) in outer domain.
    integer(ip), parameter:: global_nx(2) = int(global_lw/0.1_dp, ip) * mesh_refine

    ! Useful misc variables
    integer(ip):: j, i, nd, forcing_case
    character(len=20) :: tempchar

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
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method

    ! A detailed domain.
    ! Extend the domain on the right edge to include a 'water catching' 
    ! region, because the reflection off the back-wall occurs at late times and 
    ! contaminates the results. Park et al. (2013) mention their model suffers 
    ! this boundary reflection, which contaminates the comparison with data at 
    ! later times.
    call md%domains(2)%match_geometry_to_parent(&
        parent_domain=md%domains(1), &
        lower_left= [31.0_dp, -7.5_dp], &
        upper_right=[42.25_dp,  7.0_dp], &
        dx_refinement_factor=4_ip, &
        timestepping_refinement_factor=4_ip)
    md%domains(2)%timestepping_method = default_nonlinear_timestepping_method
    md%domains(2)%local_timestepping_scale = 0.9_dp

    ! Allocate domains and prepare comms
    call md%setup()

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
    end do

    call md%set_point_gauges_from_csv(&
        './problem_data/point_gauge_locations.csv', skip_header=1_ip)

    ! Build boundary conditions
    call setup_boundary_information()
    !md%domains(1)%boundary_subroutine => boundary_stage_transmissive_normal_momentum ! OK but slightly attenuates
    !md%domains(1)%boundary_subroutine => boundary_stage_radiation_momentum ! Strongly attenuates
    md%domains(1)%boundary_subroutine => flather_boundary ! Strongly attenuates if UH@boundary is zero
    md%domains(1)%boundary_function => boundary_function

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
