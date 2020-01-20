module local_routines 
    !!
    !! Setup the Tauranga harbour geometry and boundary conditions
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use file_io_mod, only: count_file_lines
    use linear_interpolator_mod, only: linear_interpolator_type
    use points_in_poly_mod, only: points_in_poly

    implicit none

    ! Hold some data for the boundary condition
    type :: boundary_information_type
        character(charlen):: bc_file
        real(dp), allocatable:: boundary_data(:,:)
        type(linear_interpolator_type):: gauge4_ts_function
        real(dp):: boundary_elev
        real(dp):: t0 = 0.0_dp
    end type

    !
    ! Parameter controlling the northern boundary treatment, which substantially
    ! affects the tsunami in this case.
    !
    ! Possible values are:
    !   'boundary_stage_radiation_momentum' -- Works pretty well without tuning for rk2 and related. 
    !   'flather_with_vh_from_continuity' -- about right? But depends on a scale for the boundary 
    !       VH term, which we can estimate from continuity, but needs tuning. OK for leapfrog_nonlinear.
    !   'boundary_stage_transmissive_normal_momentum' -- tsunami waves are over-amplified, probably 
    !       because of reflections from the boundary. Tide is OK
    !   'boundary_stage_transmissive_momentum' -- Similar to previous
    !   'flather_with_uh_equal_zero' -- tsunami waves under-amplified (i.e. too small). Tide is OK
    !   
    character(len=charlen) :: boundary_type = 'boundary_stage_radiation_momentum'
    !character(len=charlen) :: boundary_type = 'flather_with_vh_from_continuity'
    !character(len=charlen) :: boundary_type = 'boundary_stage_transmissive_momentum'

    ! This will hold the information -- is seen by other parts of the module
    type(boundary_information_type):: boundary_information

    contains 

    ! Read files with boundary condition info, and make a BC function
    subroutine setup_boundary_information(bc_file, boundary_elev)
        character(charlen), intent(in):: bc_file
        real(dp), intent(in):: boundary_elev

        integer(ip):: bc_unit, nr, nc, skip, i, extra

        boundary_information%bc_file = bc_file
        boundary_information%boundary_elev = boundary_elev
        open(newunit=bc_unit, file=bc_file)
        nr = count_file_lines(bc_unit)
        nc = 2
        skip = 1
        extra = 1 ! One more data point to avoid exceeding time 
        allocate(boundary_information%boundary_data(nr - skip + extra, nc))
        do i = 1, nr
            if(i > skip) then
                read(bc_unit, *) boundary_information%boundary_data(i - skip,:)
            else
                read(bc_unit, *) 
            end if 
        end do
        close(bc_unit)
        ! Extend the time-series with a constant value, so that time does not
        ! exceed model run time
        boundary_information%boundary_data(nr - skip + extra,1) = 1.0e+06_dp + &
            boundary_information%boundary_data(nr - skip + extra - 1, 1)

        boundary_information%boundary_data(nr - skip + extra,2) = &
            boundary_information%boundary_data(nr - skip + extra - 1, 2)

        !print*, 'ZEROING STAGE '
        !boundary_information$boundary_data(:,2) = 0.0_dp
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
        real(dp) :: local_elev, dhdt(1), dt, next_h(1)

        ! Set the stage
        call boundary_information%gauge4_ts_function%eval([t + boundary_information%t0], stage_uh_vh_elev(1:1))
      
        ! Get the time-derivative of stage. Useful for some approaches
        dt = 1.0e-06
        call boundary_information%gauge4_ts_function%eval([t + boundary_information%t0 + dt], next_h)
        dhdt = (next_h(1)/dt - stage_uh_vh_elev(1)/dt)
       
        ! Set the elevation 
        local_elev = domain%U(i,j,4)
        stage_uh_vh_elev(4) = local_elev

        if(local_elev >= stage_uh_vh_elev(1)) then
            ! Treat dry boundary case
            stage_uh_vh_elev(1) = local_elev
            stage_uh_vh_elev(2:3) = 0.0_dp
        else

            ! Much experimentation was conducted here. This problem is sensitive to
            ! imperfections in our semi-transmissive boundary conditions which allow a 
            ! stage forcing. (Not too surprising, because the boundary is quite close to the
            ! coast, in 'not very deep' water).

            select case(boundary_type)

            case('boundary_stage_transmissive_normal_momentum')
                ! Do nothing, because we do not need uh/vh.
                stage_uh_vh_elev(2:3) = 0.0_dp 

            case('boundary_stage_transmissive_momentum')
                ! Do nothing, because we do not need uh/vh.
                stage_uh_vh_elev(2:3) = 0.0_dp 

            case('boundary_stage_radiation_momentum') 
                ! These will never be used for this boundary
                stage_uh_vh_elev(2:3) = 0.0_dp

            case('flather_with_vh_equal_zero') 
                ! This absorbs, but also distorts the stage at the boundary too much when
                ! wave frequencies are lower
                stage_uh_vh_elev(2:3) = 0.0_dp

            case('flather_with_vh_from_continuity')
                ! Approach 4: Flat-free-surface continuity
                if(j == domain%nx(2)) then
                    stage_uh_vh_elev(2) = 0.0_dp
                    ! Assume flat free surface from offshore to the coast + estuary volume. 
                    ! That gives a rough estimate of -VH. However, in practice the factor needs tuning.
                    stage_uh_vh_elev(3) = 0.0_dp -dhdt(1) * 5000.0_dp
                else
                    stage_uh_vh_elev(2:3) = 0.0_dp
                end if

            case default
                print*, 'boundary_type not recognised'
                stop
            end select

        end if

    end function

    ! Main setup routine
    subroutine set_initial_conditions(domain)
        class(domain_type), target, intent(inout):: domain
        integer(ip):: i, j
        character(len=charlen):: input_elevation, input_stage
        real(dp), allocatable:: x(:), y(:)
        logical, allocatable:: is_inside(:)
        type(multi_raster_type):: elevation_data
        real(dp) :: wall, w
        real(dp) :: gauge_xy(3,6)
        real(dp) :: pol1(4,2), pol2(4,2)
        logical :: flatten_bathymetry_near_boundary

        ! Stage
        domain%U(:,:,STG) = -0.7269_dp !Same as first tide gauge value

        ! Set elevation with the raster
        input_elevation = './bathymetry/TAU_Whole_Harbour_10_m_srf6_rotated.tif' 

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)), is_inside(domain%nx(1)))
        x = domain%x
        call elevation_data%initialise([input_elevation])

        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip)

            flatten_bathymetry_near_boundary = (.not.(&
                boundary_type == 'boundary_stage_transmissive_normal_momentum' .or. &
                boundary_type == 'boundary_stage_transmissive_momentum'))

            if(domain%y(j) > 18500.0_dp .and. flatten_bathymetry_near_boundary) then
                ! Put a 'flat' boundary where we impose a BC, by smoothly merging with
                ! the topography over 500m
                !w = min((domain%y(j) - 18500.0_dp)/500.0_dp, 1.0_dp)
                !domain%U(:,j,ELV) = -25.0_dp * w + (1.0_dp - w) * domain%U(:,j,ELV)

                ! Rapid change
                domain%U(:,j,ELV) = -25.0_dp 
            end if
        end do

        if(domain%timestepping_method == 'cliffs') then
            ! Bathymetry smoothing is essential for CLIFFS, but can also reduce the static-noise in 
            ! the finite volume solvers where the elevation changes rapidly.
            domain%cliffs_minimum_allowed_depth = 0.1_dp
            call domain%smooth_elevation(smooth_method='cliffs')
        end if

        !call domain%smooth_elevation(smooth_method='9pt_average')

        ! The DEM needs to be 'fixed' in a few places where bridges remain. Google earth    
        ! suggests the bridges should not strongly impede the flow. So based on checks of the DEM, ....
        ! This doesn't really affect the model results near the gauges, but....
        !
        ! pol1 should have a value about -6.0
        !
        pol1(1,:) = [34259.0_dp, 10988.0_dp]
        pol1(2,:) = [34545.0_dp, 11233.0_dp]
        pol1(3,:) = [34366.0_dp, 11386.0_dp]
        pol1(4,:) = [34070.0_dp, 11161.0_dp]
        !
        ! pol2 should have a value about -0.60 
        !
        pol2(1,:) =  [34545.0_dp, 11233.0_dp] 
        pol2(2,:) =  [34892.0_dp, 11596.0_dp]
        pol2(3,:) =  [34716.0_dp, 11736.0_dp]
        pol2(4,:) =  [34366.0_dp, 11386.0_dp]
        ! apply the patch to the dem
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call points_in_poly(pol1(:,1), pol1(:,2), x, y, is_inside)
            where(is_inside) domain%U(:,j,ELV) = -6.0_dp
            call points_in_poly(pol2(:,1), pol2(:,2), x, y, is_inside)
            where(is_inside) domain%U(:,j,ELV) = -0.6_dp
        end do

        deallocate(x,y, is_inside)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        ! Wall boundaries (without boundary conditions)
        ! In interior domains these will be overwritten
        wall = 20._dp
        domain%U(:,1,ELV) = wall
        domain%U(domain%nx(1),:,ELV) = wall
        domain%U(1,:,ELV) = wall

        ! Friction 
        if(domain%timestepping_method /= 'linear') then
            domain%manning_squared = 0.025_dp * 0.025_dp
        end if

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Tauranga
    !!
    !! NTHMP currents test problem -- Tohoku tsunami observed at Tauganga Harbour, New Zealand.
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: boundary_stage_transmissive_normal_momentum, flather_boundary, &
        boundary_stage_transmissive_momentum, boundary_stage_radiation_momentum
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit
    implicit none

    ! Useful misc variables
    integer(ip):: j, nd

    ! Type holding all domains 
    type(multidomain_type) :: md

    type(timer_type) :: program_timer

    real(dp), parameter :: mesh_refine = 0.06_dp ! Increase resolution by this amount.  {e.g. 1.0 = 10m; 2 = 5m; 0.1 = 100m; etc}
    !real(dp), parameter :: mesh_refine = 0.12_dp ! Increase resolution by this amount.  {e.g. 1.0 = 10m; 2 = 5m; 0.1 = 100m; etc}
    
    real(dp) ::  global_dt = 0.22_dp / mesh_refine
    real(dp), parameter :: final_time = 3600.0_dp * 24.0_dp * 2.0_dp !6.8_dp

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 30.0_dp
    integer(ip), parameter :: only_write_grids_every_nth_output_step = 20_ip ! Write grids less often than gauges

    character(len=charlen) ::  bc_file = './boundary/ABeacon_stage_timeseries.csv'
    real(dp) :: bc_elev

    ! Length/width
    real(dp), dimension(2):: global_lw = [41000.0_dp, 18600.0_dp] ![41000.0_dp, 22400.0_dp]
    ! Lower-left corner coordinate
    real(dp), dimension(2):: global_ll = [0.0_dp, 0.0_dp]
    ! grid size (number of x/y cells)
    integer(ip), dimension(2):: global_nx = nint([4100_ip, 1860_ip] * mesh_refine) !nint([4100_ip, 2240_ip] * mesh_refine) !
    integer(ip), parameter :: nest_ratio = 3_ip
    integer(ip), parameter :: boundary_domain_thickness = 0_ip

    call program_timer%timer_start('setup')

    ! nd domains in this model
    nd = 2
    allocate(md%domains(nd))
    !md%load_balance_file = 'load_balance_partition.txt'

    !
    ! Setup basic metadata
    !

    ! Main domain, with the northern-end optionally shorn off and replaced with a
    ! boundary_domain
    md%domains(1)%lower_left = global_ll
    md%domains(1)%nx = global_nx - [0_ip, boundary_domain_thickness]
    md%domains(1)%lw = global_lw * ( ( 1.0_dp * md%domains(1)%nx ) / (1.0_dp * global_nx) )
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method
    !md%domains(1)%compute_fluxes_inner_method = 'DE1'
    !md%domains(1)%theta = 4.0_dp !! No big impact for this problem.
    !md%domains(1)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
    !md%domains(1)%linear_solver_is_truely_linear = .false.


    print*, 1, ' lw: ', md%domains(1)%lw, ' ll: ', md%domains(1)%lower_left, ' dx: ', md%domains(1)%dx, &
        ' nx: ', md%domains(1)%nx

    ! A detailed domain [Cannot partially share a physical boundary with the outer domain]
    call md%domains(2)%match_geometry_to_parent(&
        parent_domain=md%domains(1), &
        lower_left=[28500.0_dp, 11500.0_dp], &
        upper_right=[32200.0_dp, 16900.0_dp], &
        dx_refinement_factor=nest_ratio, &
        timestepping_refinement_factor=nest_ratio,&
        rounding_method='nearest')
    md%domains(2)%timestepping_method = default_nonlinear_timestepping_method !'cliffs' !'midpoint' !'rk2'
    !md%domains(2)%theta = 4.0_dp ! No big impact for this problem.

    print*, 2, ' lw: ', md%domains(2)%lw, ' ll: ', md%domains(2)%lower_left, ' dx: ', md%domains(2)%dx, &
        ' nx: ', md%domains(2)%nx

    ! Allocate domains and prepare comms
    call md%setup()

    if(md%domains(1)%is_staggered_grid) then
        ! The 'boundary_stage_radiation_momentum' boundary is not sufficiently radiative with the
        ! leapfrog nonlinear solver.
        boundary_type = 'flather_with_vh_from_continuity'
    else
        ! Good for rk2 and related solvers, and OK for cliffs.
        boundary_type = 'boundary_stage_radiation_momentum'
    end if
        

    ! Initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j))
        ! Even for linear solver, allow the (g depth dstage/dx) term to have depth varying
        md%domains(j)%linear_solver_is_truely_linear = .false.
    end do
    call md%make_initial_conditions_consistent 
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()

    ! Setup hazard points
    call md%set_point_gauges_from_csv("point_gauge_locations.csv", skip_header=1_ip)

    ! Build boundary conditions
    bc_elev = minval(md%domains(1)%U(:,:,4)) 
    call setup_boundary_information(bc_file, bc_elev)
    
    ! Boundary. Care is required in this problem, because the boundary is so close to the
    ! coast -- reflections can be a problem. A number of approaches can be tested by 
    ! changing 'boundary_type' in the local_routines module
    select case(boundary_type)
    case('boundary_stage_transmissive_normal_momentum')
        md%domains(1)%boundary_subroutine => boundary_stage_transmissive_normal_momentum
    case('boundary_stage_transmissive_momentum')
        md%domains(1)%boundary_subroutine => boundary_stage_transmissive_momentum
    case('flather_with_vh_from_continuity')
        md%domains(1)%boundary_subroutine => flather_boundary
    case('flather_with_vh_equal_zero')
        md%domains(1)%boundary_subroutine => flather_boundary
    case('boundary_stage_radiation_momentum')
        md%domains(1)%boundary_subroutine => boundary_stage_radiation_momentum
    case default
        stop "Invalid boundary_type value"
    end select
    md%domains(1)%boundary_function => boundary_function
   
    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        write(log_output_unit, *) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    print*, 'End setup'
    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    ! Evolve the code
    do while (.true.)
       
        ! Write gauges and print after 'approximate_writeout_frequency' time has passed 
        ! Don't write gauges every time
        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency, &
            write_grids_less_often = only_write_grids_every_nth_output_step, &
            write_gauges_less_often = 1_ip, &
            print_less_often = 1_ip, &
            timing_tol = 1.0e-06_dp)
        call program_timer%timer_end('IO')

        ! Finish at some point
        if (md%domains(1)%time > final_time) exit

        ! Main evolve
        call md%evolve_one_step(global_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)

end program
