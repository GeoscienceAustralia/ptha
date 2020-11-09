module local_routines 
    !! 
    !! Subroutines and data used to setup the model scenario 
    !!

    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type 
    use which_mod, only: which
    use file_io_mod, only: read_csv_into_array
    use forcing_mod, only: forcing_patch_type, apply_forcing_patch
    use iso_c_binding, only: c_loc
    implicit none

    contains 

    subroutine set_initial_conditions_generic_model(domain, input_elevation_raster, &
        input_stage_raster, hazard_points_file, skip_header,&
        adaptive_computational_extents, negative_elevation_raster, manning_n, rise_time)

        class(domain_type), target, intent(inout):: domain
        character(len=charlen), intent(in):: input_elevation_raster, &
            input_stage_raster, hazard_points_file
        integer(ip), intent(in):: skip_header
        logical, intent(in):: adaptive_computational_extents, negative_elevation_raster
        real(dp), intent(in) :: manning_n, rise_time

        integer(ip):: i, j, extra_buffer
        real(dp), allocatable:: x(:), y(:), xy_coords(:,:)
        integer(ip):: stage_raster_dim(2), xl, xu, yl, yu
        real(dp) :: stage_raster_ll(2), stage_raster_ur(2)
        type(multi_raster_type):: elevation_data, stage_data
        real(dp) :: delay_forcing_time

        character(charlen):: attribute_names(8), attribute_values(8)

        type(forcing_patch_type), pointer :: forcing_patch
        
        ! Attributes to be stored in the netcdf files as attributes
        attribute_names(1) = 'input_elevation_raster'
        attribute_values(1) = input_elevation_raster
        attribute_names(2) = 'input_stage_raster'
        attribute_values(2) = input_stage_raster
        attribute_names(3) = 'hazard_points_file'
        attribute_values(3) = hazard_points_file
        attribute_names(4) = 'timestepping_method'
        attribute_values(4) = domain%timestepping_method
        attribute_names(5) = 'metadata_ascii_filename'
        attribute_values(5) = domain%metadata_ascii_filename

        ! Include the git version
        attribute_names(6) = 'sourcecode_git_version_number'
#ifdef GITVERSION        
        attribute_values(6) = & ! Continuation to reduce the chance of > 132 characters in code
GITVERSION 
#else
        attribute_values(6) = 'not provided'
#endif

        ! Include the original source code directory
        attribute_names(7) = 'sourcecode_directory'
#ifdef SOURCEDIR       
        attribute_values(7) = & ! Continuation to reduce the chance of > 132 characters in code
SOURCEDIR 
#else
        attribute_values(7) = 'not provided'
#endif

        attribute_names(8) = 'model_run_directory'
        call get_environment_variable("PWD", attribute_values(8))


        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))

        x = domain%x

        print*, "Setting elevation ..."

        call elevation_data%initialise([input_elevation_raster])

        print*, '    bounding box of input elevation: ' 
        print*, '    ', elevation_data%lowerleft
        print*, '    ', elevation_data%upperright

        ! Treat periodic case, where there will be 2 points on either side
        ! which exceed the input raster extent
        if(x(domain%nx(1)) - x(1) > 360.0_dp) then
            x(1:2) = x(1:2) + 360.0_dp
            x((domain%nx(1) - 1):domain%nx(1)) = x((domain%nx(1) - 1):domain%nx(1)) - 360.0_dp
        end if

        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip)
            if(negative_elevation_raster) domain%U(:,j,ELV) = -domain%U(:,j,ELV)

        end do
        call elevation_data%finalise()

        ! This makes cliffs stable
        if(domain%timestepping_method == 'cliffs') then
            call domain%smooth_elevation(smooth_method='cliffs')
        end if

        print*, "Setting stage ..."

        ! Set stage -- zero outside of initial condition file range
        domain%U(:,:,[STG,UH,VH]) = 0.0_dp
       
        call stage_data%initialise([input_stage_raster])

        ! Get the x indices which are inside the stage raster 
        ! We permit this to only cover a small part of the domain
        extra_buffer = 0
        xl = count(domain%x < stage_data%lowerleft(1)) + 1  - extra_buffer
        xU = count(domain%x < stage_data%upperright(1)) + extra_buffer
        yl = count(domain%y < stage_data%lowerleft(2)) + 1 - extra_buffer
        yU = count(domain%y < stage_data%upperright(2)) + extra_buffer

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
                    call stage_data%get_xy(x(xl:xU), y(xl:xU), &
                        domain%U(xl:xU, j, STG), (xU-xl+1))
                end if
            end do
        end if
        CALL stage_data%finalise() 


        if(rise_time > 0.0_dp .and. xl > 0 .and. xU > 0 .and. yl > 0 .and. yU > 0) then
            ! Apply the source over some finite rise time, rather than as the initial stage.

            ! For rk2, we should not start the forcing right away to ensure the full forcing
            ! is applied. This is because of rk2's approach:
            !     [start, advance 2-time-steps to end, then-average(start,end)].
            ! If we don't have timestepping before the forcing start_time, we miss out on part of
            ! of the forcing that should have been obtained (hypothetically if we evolved from time=-dt to 
            ! time=0, we would have included some of the forcing -- this is the part we miss).
            if(any(domain%timestepping_method == [character(len=charlen) :: 'rk2', 'rk2n'])) then
                delay_forcing_time = 30.0_dp
            else
                delay_forcing_time = 0.0_dp
            end if

            
            ! Make a forcing_patch that applies the source 
            allocate(forcing_patch)
            call forcing_patch%setup(start_time=delay_forcing_time, end_time=rise_time+delay_forcing_time, &
                i0=xl, i1=xU, j0=yl, j1=yU, k0=STG, k1=STG)

            ! Set the forcing to be equal to the stage that we just created
            forcing_patch%forcing_work(xl:xU, yl:yU, STG) = domain%U(xl:xU, yl:yU, STG)

            ! Move the forcing patch into the domain, and ensure it is used in a forcing term
            domain%forcing_context_cptr = c_loc(forcing_patch)
            forcing_patch => NULL()
            domain%forcing_subroutine => apply_forcing_patch

            ! Re-set the stage deformation to zero, as it will be applied in a forcing instead.
            domain%U(:,:,STG) = 0.0_dp
        end if


        print*, '    max stage: ', maxval(domain%U(:,:,STG))
        print*, '    min stage: ', minval(domain%U(:,:,STG))
        print*, '    max elev: ', maxval(domain%U(:,:,ELV))
        print*, '    min elev: ', minval(domain%U(:,:,ELV))

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV))

        ! Setup gauges
        if(hazard_points_file /= '') then
            print*, ''
            print*, "Reading hazard points and ID's"
            call read_csv_into_array(xy_coords, hazard_points_file, skip_header)
            if( (any(xy_coords(1,:) < domain%lower_left(1))).OR. &
                (any(xy_coords(1,:) > domain%lower_left(1) + domain%lw(1))).OR. &
                (any(xy_coords(2,:) > domain%lower_left(2) + domain%lw(2))).OR. &
                (any(xy_coords(2,:) < domain%lower_left(2))) ) then
                print*, '    # WARNING: Some hazard points outside domain extents are being discarded'
                print*, '    #          (Could reflect running in parallel, or round-off error)'
            end if

            print*, '    Setting up gauges'
            call domain%setup_point_gauges(xy_coords(1:2,:), time_series_var=[STG, UH, VH], static_var=[ELV], &
                gauge_ids = xy_coords(3,:), attribute_names=attribute_names, attribute_values=attribute_values)
            print*, '    The number of points on this domain is', domain%point_gauges%n_gauges
            print*, '    The first point is (lon,lat): ', domain%point_gauges%xy(1:2, 1)
            print*, '    The last  point is (lon,lat): ', domain%point_gauges%xy(1:2, domain%point_gauges%n_gauges)
            print*, ''
        else
            print*, ''
            print*, 'No hazard points'
            print*, ''
        end if
       
        if(adaptive_computational_extents) then
            ! Make up some xL/xU limits to do preliminary checks on allowing that
            ! to evolve
            domain%xL = xl
            domain%xU = xU
            domain%yL = yL
            domain%yU = yU
        end if

        deallocate(x,y)
        if(allocated(xy_coords)) then
            deallocate(xy_coords)
        end if

        if(allocated(domain%manning_squared)) then
            print*, 'Setting manning friction'
            domain%manning_squared = manning_n * manning_n
        end if
        
        print*, 'Initial conditions set'
        print*, ''

    end subroutine
end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program generic_model
    !! Generic single-grid spherical coordinate model. Useful for quickly running ocean-scale tsunami propagation problems.

    use global_mod, only: ip, dp, minimum_allowed_depth
    use domain_mod, only: domain_type
    use boundary_mod, only: flather_boundary, periodic_EW_reflective_NS
    use local_routines
    implicit none

    integer(ip):: j
    real(dp):: last_write_time, rain_rate
    type(domain_type):: domain

    character(charlen) :: timestepping_method, input_parameter_file, namelist_output_filename

    real(dp) :: approximate_writeout_frequency, final_time, manning_n, cliffs_minimum_allowed_depth
    real(dp):: timestep, global_ur(2), global_ll(2), global_lw(2), cfl, dx(2), linear_friction_coeff, rise_time
    integer(ip):: global_nx(2), skip_header_hazard_points_file, file_unit_temp, grid_output_spatial_stride
    character(len=charlen):: input_elevation_raster, input_stage_raster, &
        hazard_points_file, output_basedir
    logical:: record_max_U, output_grid_timeseries, adaptive_computational_extents, &
        negative_elevation_raster, linear_solver_is_truely_linear

    ! Key input data
    namelist /MODELCONFIG/ &
        input_elevation_raster, input_stage_raster, global_ll, &
        global_ur, global_nx, approximate_writeout_frequency, output_basedir, &
        final_time, timestepping_method, manning_n, linear_friction_coeff, rise_time, &
        cfl, hazard_points_file, skip_header_hazard_points_file, record_max_U,&
        output_grid_timeseries, adaptive_computational_extents, negative_elevation_raster, &
        linear_solver_is_truely_linear, grid_output_spatial_stride, cliffs_minimum_allowed_depth

    ! Predefine some variables that might not be in the input file
    manning_n = 0.0_dp
    linear_friction_coeff = 0.0_dp
    linear_solver_is_truely_linear = .true.
    grid_output_spatial_stride = 1
    cliffs_minimum_allowed_depth = 1.0e-03_dp
    rise_time = 0.0_dp

    ! Read the input file -- the name of this file should be the first
    ! commandline argument
    call get_command_argument(1, input_parameter_file)
    open(newunit=file_unit_temp, file=input_parameter_file)
    read(file_unit_temp, nml=MODELCONFIG)
    close(file_unit_temp)


    ! Logical check
#ifndef SPHERICAL
    print*, 'This code is setup assuming spherical coordinates'
    print*, 'but -DSPHERICAL was not passed to the compiler.'
    stop
#endif
  
    ! Report key output information 
    print*, ''
    print*, '#### INPUT FROM MODELCONFIG ####'
    print*, ''
    print*, 'input_elevation_raster: ', TRIM(input_elevation_raster)
    print*, 'input_stage_raster: ', TRIM(input_stage_raster)
    print*, 'global_ll: ', global_ll
    print*, 'global_ur: ', global_ur
    print*, 'global_nx: ', global_nx
    print*, 'approximate_writeout_frequency: ', approximate_writeout_frequency
    print*, 'output_basedir: ', TRIM(output_basedir)
    print*, 'final_time: ', final_time
    print*, 'timestepping_method: ', TRIM(timestepping_method)
    print*, 'linear_solver_is_truely_linear: ', linear_solver_is_truely_linear
    print*, 'manning n: ', manning_n
    print*, 'linear friction coefficient: ', linear_friction_coeff
    print*, 'cfl: ', cfl
    print*, 'hazard_points_file: ', TRIM(hazard_points_file)
    print*, 'skip_header_hazard_points_file: ', skip_header_hazard_points_file
    print*, 'record_max_U: ', record_max_U
    print*, 'output_grid_timeseries: ', output_grid_timeseries
    print*, 'adaptive_computational_extents: ', adaptive_computational_extents
    print*, 'negative_elevation_raster: ', negative_elevation_raster
    print*, 'grid_output_spatial_stride: ', grid_output_spatial_stride
    print*, ''
    print*, '#### END INPUT ####'

    print*, ''
    print*, 'Other key variables: '
    print*, 'dp: ', dp
    print*, 'ip: ', ip
    print*, 'charlen: ', charlen

    global_lw = global_ur - global_ll
    dx = global_lw/(1.0_dp*global_nx)

    if(abs(global_lw(1) - 360.0_dp) < (0.1_dp*dx(1))) then
        ! Global model in EW direction.
        ! Case with periodic EW boundaries
        print*, ''
        print*, 'Assuming global model with periodic boundaries: Appending cells to domain'
        print*, ''
        domain%boundary_subroutine => periodic_EW_reflective_NS 
        ! The periodic boundaries involve 2 cells on the model exterior
        domain%exterior_cells_width = 2

        global_nx = global_nx + [4_ip, 0_ip]
        global_ur = global_ur + [2*dx(1), 0.0_dp]
        global_ll = global_ll - [2*dx(1), 0.0_dp]
        global_lw = global_ur - global_ll
    else
        domain%boundary_subroutine => flather_boundary
    end if

    domain%timestepping_method = timestepping_method
    domain%record_max_U = record_max_U
    domain%output_basedir = output_basedir
    domain%linear_solver_is_truely_linear = linear_solver_is_truely_linear
    domain%nc_grid_output%spatial_stride = grid_output_spatial_stride

    !domain%friction_type = 'chezy'
    domain%linear_friction_coeff = linear_friction_coeff

    ! !Optionally suppress limiting with rk2
    !domain%theta = 4.0_dp
    
    ! Allocate domain -- must have set timestepping method BEFORE this
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    domain%cfl = cfl

    ! Append the input namelist to the metadata
    namelist_output_filename = trim(domain%output_folder_name) // '/modelconfig.in'
    open(newunit=file_unit_temp, file=namelist_output_filename) 
    write(file_unit_temp, MODELCONFIG)
    close(file_unit_temp)
    call domain%log_outputs()
    write(domain%logfile_unit, MODELCONFIG)

    if(domain%timestepping_method == 'cliffs') then
        domain%cliffs_minimum_allowed_depth = cliffs_minimum_allowed_depth
    end if
    ! Call local routine to set initial conditions
    call set_initial_conditions_generic_model(domain, input_elevation_raster,&
        input_stage_raster, hazard_points_file, skip_header_hazard_points_file,&
        adaptive_computational_extents, negative_elevation_raster, manning_n, rise_time)

    timestep = domain%stationary_timestep_max()
    write(domain%logfile_unit, *) 'ts: ', timestep

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    ! Evolve the code
    do while (.true.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files(time_only = (output_grid_timeseries .EQV. .FALSE.))
            call domain%write_gauge_time_series()
            write(domain%logfile_unit, *) 'Mass balance: ', domain%mass_balance_interior()

        end if

        if (domain%time > final_time) exit 

        if(.not. domain%adaptive_timestepping) then
            ! Example with fixed timestep
            call domain%evolve_one_step(timestep=timestep)
        else
            call domain%evolve_one_step()
        end if

        ! Evolve the active domain?
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

    END DO

    write(domain%logfile_unit, *) ''
    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print(output_file_unit=domain%logfile_unit)

    call domain%finalise()

end program
