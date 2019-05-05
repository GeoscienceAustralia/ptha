#ifdef TIMER
#   define TIMER_START(tname) call domain%timer%timer_start(tname)
#   define TIMER_STOP(tname)  call domain%timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

module local_routines 

    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: gdal_raster_dataset_type
    use which_mod, only: which
    use file_io_mod, only: read_csv_into_array
    use stop_mod, only: generic_stop
    implicit none

    contains 

    subroutine set_initial_conditions_generic(domain, input_elevation_raster, &
        input_stage_raster, hazard_points_file, skip_header,&
        adaptive_computational_extents, negative_elevation_raster)

        class(domain_type), target, intent(inout):: domain
        character(len=charlen), intent(in):: input_elevation_raster, &
            input_stage_raster, hazard_points_file
        integer(ip), intent(in):: skip_header
        logical, intent(in):: adaptive_computational_extents, negative_elevation_raster

        integer(ip):: i, j
        real(dp), allocatable:: x(:), y(:), xy_coords(:,:)
        integer(ip):: xl, xU, yl, yU
        type(gdal_raster_dataset_type):: elevation_data, stage_data

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))

        x = domain%x

        print*, "Setting elevation ..."

        call elevation_data%initialise(input_elevation_raster)

        print*, '    bounding box of input elevation: ' 
        print*, '    ', elevation_data%lowerleft
        print*, '    ', elevation_data%upperright

        ! Allow periodic elevation data
        where( x < elevation_data%lowerleft(1) )
            x = x + 360.0_dp
        end where
        where( x > elevation_data%upperright(1) )
            x = x - 360.0_dp
        end where

        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip)
            if(negative_elevation_raster) domain%U(:,j,ELV) = -domain%U(:,j,ELV)
        end do
        call elevation_data%finalise()
    
        print*, "Setting stage ..."
        ! Set stage -- zero outside of initial condition file range
        domain%U(:,:,[STG,UH,VH]) = 0.0_dp
        CALL stage_data%initialise(input_stage_raster)

        ! Do not allow periodic stage data (unlike elevation, since the latter has to be set
        ! everywhere, whereas the former will just be set inside the provided data and use
        ! the default 0 everywhere else)
        x = domain%x

        ! Get the x indices which are inside the stage raster 
        ! We permit this to only cover a small part of the domain
        xl = count(domain%x < stage_data%lowerleft(1)) + 1
        xU = count(domain%x < stage_data%upperright(1))
        yl = count(domain%y < stage_data%lowerleft(2)) + 1
        yU = count(domain%y < stage_data%upperright(2))

        print*, '    bounding box of input stage: ' 
        print*, '    ', stage_data%lowerleft
        print*, '    ', stage_data%upperright
        print*, '    xl: ', xl, ' xU: ', xU
        print*, '    yl: ', yl, ' yU: ', yU

        if(xl > 0 .AND. xU > 0) then
            do j = 1, domain%nx(2)
                if((domain%y(j) > stage_data%lowerleft(2)) .AND. &
                    (domain%y(j) < stage_data%upperright(2))) then
                    y(xl:xU) = domain%y(j)
                    CALL stage_data%get_xy(x(xl:xU), y(xl:xU), &
                        domain%U(xl:xU, j, STG), (xU-xl+1))
                end if
            end do
        end if
        call stage_data%finalise() 

        print*, '    max stage: ', maxval(domain%U(:,:,STG))
        print*, '    min stage: ', minval(domain%U(:,:,STG))

        ! By updating the boundary, we prevent changes to the elevation that can occur
        ! e.g. if the elevation could be updated at edges by the boundary, 
        call domain%update_boundary()

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV)+1.0e-06)


        ! Setup gauges
        if(hazard_points_file /= '') then
            print*, ''
            print*, 'Reading hazard points'
            call read_csv_into_array(xy_coords, hazard_points_file, skip_header)
            print*, '    Setting up gauges'
            print*, domain%lower_left, domain%lw
            call domain%setup_point_gauges(xy_coords(1:2,:), time_series_var=[STG], static_var=[ELV], gauge_ids=xy_coords(3,:))
            if(allocated(domain%point_gauges%xy)) then
                print*, '    The number of points is', domain%point_gauges%n_gauges
                print*, '    The first point is (lon,lat): ', domain%point_gauges%xy(1:2,1)
                print*, '    The last  point is (lon,lat): ', domain%point_gauges%xy(1:2, domain%point_gauges%n_gauges)
                print*, ''
            end if
        else
            print*, ''
            print*, 'No hazard points'
            print*, ''
        end if
       
        if(adaptive_computational_extents) then
#ifdef COARRAY
            print*, 'Cannot have adaptive computational extents with coarrays'
            call generic_stop()
#endif
            ! Make up some xL/xU limits to do preliminary checks on allowing that
            ! to evolve
            domain%xL = xl
            domain%xU = xU
            domain%yL = yL
            domain%yU = yU
        end if

        print*, 'Initial conditions set'
        print*, ''

        
    end subroutine

end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program coarray_example

    use global_mod, only: ip, dp, minimum_allowed_depth
    use domain_mod, only: domain_type
    use boundary_mod, only: flather_boundary, periodic_EW_reflective_NS
#ifdef COARRAY
    use coarray_point2point_comms_mod, only: allocate_p2p_comms
#endif
    use local_routines
    implicit none

    integer(ip) :: i, j, ti, ni, nx_ca, ny_ca
    real(dp) :: last_write_time
    type(domain_type) :: domain
    character(charlen) :: timestepping_method, input_parameter_file
    real(dp) :: approximate_writeout_frequency, final_time 

    real(dp) :: global_ur(2), global_ll(2), global_lw(2), cfl, dx(2)
    integer(ip) :: global_nx(2)
    integer(ip) :: skip_header_hazard_points_file
    character(len=charlen) :: input_elevation_raster, input_stage_raster, &
        hazard_points_file, output_basedir
    logical :: record_max_u, output_grid_timeseries, adaptive_computational_extents, &
        negative_elevation_raster, ew_periodic, ns_periodic

    real(dp) :: timestep
    
    ! Input data
    namelist /MODELCONFIG/ &
        input_elevation_raster, input_stage_raster, global_ll, &
        global_ur, global_nx, approximate_writeout_frequency, output_basedir, &
        final_time, timestepping_method,&
        cfl, hazard_points_file, skip_header_hazard_points_file, record_max_U,&
        output_grid_timeseries, adaptive_computational_extents, negative_elevation_raster, &
        nx_ca, ny_ca

#ifdef COARRAY
    sync all
#endif
    TIMER_START('SETUP')

#ifdef COARRAY
    ni = num_images()
    ti = this_image()
#else
    ti = 1
    ni = 1
#endif

    ! Read input namelist 
    call get_command_argument(1, input_parameter_file)
    open(newunit=j, file=input_parameter_file, action='read')
    read(j, nml=MODELCONFIG)
    close(j)

    domain%myid = 1000 + ti

    print*, ''
    print*, '#### INPUT FROM MODELCONFIG ####'
    print*, ''
    print*, 'input_elevation_raster: ', trim(input_elevation_raster)
    print*, 'input_stage_raster: ', trim(input_stage_raster)
    print*, 'global_ll: ', global_ll
    print*, 'global_ur: ', global_ur
    print*, 'global_nx: ', global_nx
    print*, 'approximate_writeout_frequency: ', approximate_writeout_frequency
    print*, 'output_basedir: ', trim(output_basedir)
    print*, 'final_time: ', final_time
    print*, 'timestepping_method: ', trim(timestepping_method)
    print*, 'cfl: ', cfl
    print*, 'hazard_points_file: ', trim(hazard_points_file)
    print*, 'skip_header_hazard_points_file: ', skip_header_hazard_points_file
    print*, 'record_max_U: ', record_max_U
    print*, 'output_grid_timeseries: ', output_grid_timeseries
    print*, 'adaptive_computational_extents: ', adaptive_computational_extents
    print*, 'negative_elevation_raster: ', negative_elevation_raster
    print*, ''
    print*, '#### END INPUT ####'

    print*, ''
    print*, 'Other key variables: '
    print*, 'dp: ', dp
    print*, 'ip: ', ip
    print*, 'charlen: ', charlen

    global_lw = global_ur - global_ll
    dx = global_lw/(1.0_dp*global_nx)

    ! Use either ew-periodic or flather boundaries
    if(abs(global_lw(1) - 360.0_dp) < (0.1_dp*dx(1))) then
        ! Global model in EW direction.
        ! Case with periodic EW boundaries
        print*, ''
        print*, 'Assuming global model with periodic boundaries'
        print*, ''
        domain%boundary_subroutine => periodic_EW_reflective_NS 

#ifndef COARRAY
        global_nx = global_nx + [4,0]
        global_ur = global_ur + [2*dx(1), 0.0_dp]
        global_ll = global_ll - [2*dx(1), 0.0_dp]
        global_lw = global_ur - global_ll
#endif

        ew_periodic = .TRUE.
        ns_periodic = .FALSE.

    else
        domain%boundary_subroutine => flather_boundary

        ew_periodic = .FALSE.
        ns_periodic = .FALSE.

    end if

    domain%timestepping_method = timestepping_method
    domain%record_max_U = record_max_U
    domain%output_basedir = output_basedir

#ifndef COARRAY    
    call domain%allocate_quantities(global_lw, global_nx, global_ll)
#else
    print*, 'Allocating image ', ti
    ! Allocate domain with coarrays
    call domain%allocate_quantities(global_lw, global_nx, global_ll, co_size_xy = [nx_ca, ny_ca], &
        ew_periodic = ew_periodic, ns_periodic = ns_periodic)
    sync all
    call allocate_p2p_comms
#endif
    domain%cfl = cfl

    call domain%log_outputs()
    write(domain%logfile_unit, MODELCONFIG)

    ! Call local routine to set initial conditions
    call set_initial_conditions_generic(domain, input_elevation_raster,&
        input_stage_raster, hazard_points_file, skip_header_hazard_points_file,&
        adaptive_computational_extents, negative_elevation_raster)


#ifdef COARRAY
    ! Also print information about whether boundary conditions are applied
    write(domain%logfile_unit, *) 'ti: ', ti, 'boundary_exterior: ', domain%boundary_exterior
#endif

    timestep = domain%linear_timestep_max()
    if(domain%timestepping_method /= 'linear') then
        ! The Finite volume methods have 1/2 as long a timestep for the same cfl
        timestep = 0.5_dp * timestep
    end if

#ifdef COARRAY
    call co_min(timestep)
    write(domain%logfile_unit, *) 'reduced ts: ', timestep
#else
    write(domain%logfile_unit, *) 'ts: ', timestep
#endif

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

#ifdef COARRAY
    sync all
#endif

    TIMER_STOP('SETUP')

    ! Evolve the code
    do while (.true.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files(time_only = (output_grid_timeseries .eqv. .false.))
            call domain%write_gauge_time_series()

        end if

        if (domain%time > final_time) exit

        !! Example with fixed timestep
        call domain%evolve_one_step(timestep=timestep)

        ! Evolve the active domain? NOT WITH COARRAYS
        domain%xL = max(domain%xL - 1, 1)
        domain%xU = min(domain%xU + 1, domain%nx(1))
        domain%yL = max(domain%yL - 1, 1)
        domain%yU = min(domain%yU + 1, domain%nx(2))

        ! Treatment of spherical models with periodic EW conditions
        ! The BC region is within 4 cells of the boundaries (considering
        ! cells that copy-out, as well as cells that copy-in)
        if((domain%xl <= 4).or.(domain%xU >= domain%nx(1) - 3)) then
            domain%xl = 1_ip
            domain%xU = domain%nx(1)
        end if

    end do

    print*, ''
    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print(output_file_unit = domain%logfile_unit)

    call domain%finalise()

end program
