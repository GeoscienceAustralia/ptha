MODULE local_routines 
    ! 
    ! Subroutines & data used to setup the model scenario (comparable to
    ! project.py in ANUGA)
    !
    USE global_mod, only: dp, ip, charlen, wall_elevation
    USE domain_mod, only: domain_type, STG, UH, VH, ELV
    USE read_raster_mod, only: gdal_raster_dataset_type
    USE which_mod, only: which
    USE file_io_mod, only: read_csv_into_array
    IMPLICIT NONE

    CONTAINS 

    SUBROUTINE set_initial_conditions_generic_model(domain, input_elevation_raster, &
        input_stage_raster, hazard_points_file, skip_header,&
        adaptive_computational_extents, negative_elevation_raster)

        CLASS(domain_type), TARGET, INTENT(INOUT):: domain
        CHARACTER(len=charlen), INTENT(IN):: input_elevation_raster, &
            input_stage_raster, hazard_points_file
        INTEGER(ip), INTENT(IN):: skip_header
        LOGICAL, INTENT(IN):: adaptive_computational_extents, negative_elevation_raster

        INTEGER(ip):: i, j
        REAL(dp), ALLOCATABLE:: x(:), y(:), xy_coords(:,:)
        INTEGER(ip):: stage_raster_dim(2), xl, xU, yl, yU
        REAL(dp) :: stage_raster_ll(2), stage_raster_ur(2)
        TYPE(gdal_raster_dataset_type):: elevation_data, stage_data

        CHARACTER(charlen):: attribute_names(8), attribute_values(8)
        
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
        ALLOCATE(x(domain%nx(1)), y(domain%nx(1)))

        x = domain%x

        print*, "Setting elevation ..."

        CALL elevation_data%initialise(input_elevation_raster)

        print*, '    bounding box of input elevation: ' 
        print*, '    ', elevation_data%lowerleft
        print*, '    ', elevation_data%upperright

        ! Treat periodic case, where there will be 2 points on either side
        ! which exceed the input raster extent
        if(x(domain%nx(1)) - x(1) > 360.0_dp) then
            x(1:2) = x(1:2) + 360.0_dp
            x((domain%nx(1) - 1):domain%nx(1)) = x((domain%nx(1) - 1):domain%nx(1)) - 360.0_dp
        end if

        DO j = 1, domain%nx(2)
            y = domain%y(j)
            CALL elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip)
            if(negative_elevation_raster) domain%U(:,j,ELV) = -domain%U(:,j,ELV)
        END DO
        CALL elevation_data%finalise()
    

        print*, "Setting stage ..."

        ! Set stage -- zero outside of initial condition file range
        domain%U(:,:,[STG,UH,VH]) = 0.0_dp
       
        CALL stage_data%initialise(input_stage_raster)

        ! Get the x indices which are inside the stage raster 
        ! We permit this to only cover a small part of the domain
        xl = COUNT(domain%x < stage_data%lowerleft(1)) + 1
        xU = COUNT(domain%x < stage_data%upperright(1))
        yl = COUNT(domain%y < stage_data%lowerleft(2)) + 1
        yU = COUNT(domain%y < stage_data%upperright(2))

        print*, '    bounding box of input stage: ' 
        print*, '    ', stage_data%lowerleft
        print*, '    ', stage_data%upperright
        print*, '    xl: ', xl, ' xU: ', xU
        print*, '    yl: ', yl, ' yU: ', yU

        if(xl > 0 .AND. xU > 0) then
            DO j = 1, domain%nx(2)
                if((domain%y(j) > stage_data%lowerleft(2)) .AND. &
                    (domain%y(j) < stage_data%upperright(2))) then
                    y(xl:xU) = domain%y(j)
                    CALL stage_data%get_xy(x(xl:xU), y(xl:xU), &
                        domain%U(xl:xU, j, STG), (xU-xl+1))
                end if
            END DO
        end if
        CALL stage_data%finalise() 

        print*, '    max stage: ', maxval(domain%U(:,:,STG))
        print*, '    min stage: ', minval(domain%U(:,:,STG))

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
                print*, '    # WARNING: Some hazard points outside domain extents are being clipped'
                print*, '    #          (Could include points on the domain edge, due to round-off error)'
                print*, '    # The original hazard point range:'
                print*, '    #     lon: ', minval(xy_coords(1,:)), maxval(xy_coords(1,:))
                print*, '    #     lat: ', minval(xy_coords(2,:)), maxval(xy_coords(2,:))
                print*, '    # will be truncated to the range of the cell midpoints'
                xy_coords(1,:) = min(xy_coords(1,:), maxval(domain%x))
                xy_coords(1,:) = max(xy_coords(1,:), minval(domain%x))
                xy_coords(2,:) = min(xy_coords(2,:), maxval(domain%y))
                xy_coords(2,:) = max(xy_coords(2,:), minval(domain%y))
            end if

            print*, '    Setting up gauges'
            call domain%setup_point_gauges(xy_coords(1:2,:), time_series_var=[STG, UH, VH], static_var=[ELV], &
                gauge_ids = xy_coords(3,:), attribute_names=attribute_names, attribute_values=attribute_values)
            print*, '    The number of points is', domain%point_gauges%n_gauges
            print*, '    The first point is (lon,lat): ', xy_coords(1,1), xy_coords(2,1)
            print*, '    The last  point is (lon,lat): ', xy_coords(1:2, domain%point_gauges%n_gauges)
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

        DEALLOCATE(x,y)
        if(allocated(xy_coords)) then
            DEALLOCATE(xy_coords)
        end if
        
        print*, 'Initial conditions set'
        print*, ''
        
    END SUBROUTINE
END MODULE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM generic_model
    USE global_mod, only: ip, dp, minimum_allowed_depth
    USE domain_mod, only: domain_type
    USE boundary_mod, only: flather_boundary, periodic_EW_reflective_NS
    USE local_routines
    IMPLICIT NONE

    INTEGER(ip):: j
    REAL(dp):: last_write_time, rain_rate
    TYPE(domain_type):: domain

    CHARACTER(charlen) :: timestepping_method, input_parameter_file, namelist_output_filename

    REAL(dp) :: approximate_writeout_frequency, final_time 
    REAL(dp):: timestep, global_ur(2), global_ll(2), global_lw(2), cfl, dx(2)
    INTEGER(ip):: global_nx(2), skip_header_hazard_points_file, file_unit_temp
    CHARACTER(len=charlen):: input_elevation_raster, input_stage_raster, &
        hazard_points_file, output_basedir
    LOGICAL:: record_max_U, output_grid_timeseries, adaptive_computational_extents, &
        negative_elevation_raster

    ! Key input data
    namelist /MODELCONFIG/ &
        input_elevation_raster, input_stage_raster, global_ll, &
        global_ur, global_nx, approximate_writeout_frequency, output_basedir, &
        final_time, timestepping_method,&
        cfl, hazard_points_file, skip_header_hazard_points_file, record_max_U,&
        output_grid_timeseries, adaptive_computational_extents, negative_elevation_raster

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
    print*, 'cfl: ', cfl
    print*, 'hazard_points_file: ', TRIM(hazard_points_file)
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

    if(abs(global_lw(1) - 360.0_dp) < (0.1_dp*dx(1))) then
        ! Global model in EW direction.
        ! Case with periodic EW boundaries
        print*, ''
        print*, 'Assuming global model with periodic boundaries: Appending cells to domain'
        print*, ''
        domain%boundary_subroutine => periodic_EW_reflective_NS 
        ! The periodic boundaries involve 2 cells on the model exterior
        domain%exterior_cells_width = 2

        global_nx = global_nx + [4,0]
        global_ur = global_ur + [2*dx(1), 0.0]
        global_ll = global_ll - [2*dx(1), 0.0]
        global_lw = global_ur - global_ll
    else
        domain%boundary_subroutine => flather_boundary
    end if

    domain%timestepping_method = timestepping_method
    domain%cfl = cfl
    domain%record_max_U = record_max_U
    domain%output_basedir = output_basedir
    
    ! Allocate domain -- must have set timestepping method BEFORE this
    CALL domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Append the input namelist to the metadata
    namelist_output_filename = TRIM(domain%output_folder_name) // '/modelconfig.in'
    open(newunit=file_unit_temp, file=namelist_output_filename) 
    write(file_unit_temp, MODELCONFIG)
    close(file_unit_temp)
    CALL domain%log_outputs()
    write(domain%logfile_unit, MODELCONFIG)

    ! Call local routine to set initial conditions
    CALL set_initial_conditions_generic_model(domain, input_elevation_raster,&
        input_stage_raster, hazard_points_file, skip_header_hazard_points_file,&
        adaptive_computational_extents, negative_elevation_raster)

    timestep = domain%linear_timestep_max()
    write(domain%logfile_unit, *) 'ts: ', timestep

    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    ! Evolve the code
    DO WHILE (.TRUE.)

        IF(domain%time - last_write_time >= approximate_writeout_frequency) THEN

            last_write_time = last_write_time + approximate_writeout_frequency

            CALL domain%print()
            CALL domain%write_to_output_files(time_only = (output_grid_timeseries .EQV. .FALSE.))
            CALL domain%write_gauge_time_series()
            write(domain%logfile_unit, *) 'Mass balance: ', domain%mass_balance_interior()

        END IF

        IF (domain%time > final_time) THEN
            EXIT 
        END IF

        !! Example with fixed timestep
        if(domain%timestepping_method == 'linear') then
            CALL domain%evolve_one_step(timestep=timestep)
        else
            CALL domain%evolve_one_step()
        end if

        !! Evolve the active domain?
        domain%xL = max(domain%xL - 1, 1)
        domain%xU = min(domain%xU + 1, domain%nx(1))
        domain%yL = max(domain%yL - 1, 1)
        domain%yU = min(domain%yU + 1, domain%nx(2))

        ! Treatment of spherical models with periodic EW conditions
        ! The BC region is within 4 cells of the boundaries (considering
        ! cells that copy-out, as well as cells that copy-in)
        if((domain%xl <= 4).OR.(domain%xU >= domain%nx(1) - 3)) then
            domain%xl = 1_ip
            domain%xU = domain%nx(1)
        end if

    END DO

    write(domain%logfile_unit, *) ''
    CALL domain%write_max_quantities()

    ! Print timing info
    CALL domain%timer%print(output_file_unit=domain%logfile_unit)

    CALL domain%finalise()

END PROGRAM
