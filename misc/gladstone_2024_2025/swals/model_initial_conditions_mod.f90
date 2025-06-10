module initial_conditions_mod
    !!
    !! Routines for initial conditions and other conveniences, for a
    !! global-scale model with regional high-resolution domains
    !!
    
    use iso_c_binding, only: c_loc
    
    ! Default integer / real precision and character length
    use global_mod, only: dp, ip, charlen, minimum_allowed_depth
    ! domain type and indices of key variables in the array domain%U:
    !     stage (domain%U(:,:,STG))
    !     depth-integrated-velocities (domain%U(:,:,UH:VH))
    !     elevation (domain%U(:,:,ELV) )
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    ! Interpolation from a collection of rasters (with an order of preference)
    use read_raster_mod, only: multi_raster_type
    ! File unit for logging
    use logging_mod, only: log_output_unit
    ! burn lines  into the elevation grid
    use burn_into_grid_mod, only: xyz_lines_type, polygons_values_type
    ! Throw errors gracefully in parallel or serial
    use stop_mod, only: generic_stop
    ! Used to apply an Okada source over a finite time.
    use forcing_mod, only: forcing_patch_type, apply_forcing_patch
    ! Used to read the elevation data filenames in preference order
    use file_io_mod, only: read_character_file, count_file_lines, read_csv_into_array
    
    ! Multidomain variables defining this specific case
    use multidomain_design_mod, only: &
        global_lw, global_ll, use_periodic_EW_multidomain, &
        swals_elevation_files_in_preference_order, &
        open_estuary_entrances_polygons_values_file, &
        tidal_adjustment_file, &
        override_initial_stage_polygons_values_file, &
        raster_na_below_limit, smooth_elevation_along_nesting_boundaries, &
        mannning_n_file_list_in_preference_order, &
        default_manning_staggered_grid, default_manning_non_staggered_grid, &
        breakwalls_file_list, inverts_file_list, &
        final_time_full_runs, final_time_test_runs, &
        NNL, read_multidomain_design_variables, &
        real_table_type ! = type used for metadata about domains in each nesting level
    
    implicit none
    
    private
    public :: set_initial_conditions, rise_time, stage_file_is_raster, ambient_sea_level, parse_inputs
    
    logical :: stage_file_is_raster
        !! Is the input stage file is a raster (.true.) or a .csv (.false.).
        !! In the latter case, the csv file will contain the time-varying source info.
    
    real(dp), parameter :: rise_time = 0.0_dp ! 180.0_dp
        !! Okada deformation is applied over rise_time seconds.
        !! If using non-zero rise time, beware interpretation of max-stage -
        !! in areas of subsidence it could reflect the pre-subsidence topography.
        !! This could be worked-around by resetting domain%max_U at some time after the
        !! forcing has finished, but I haven't done that.
        !! If stage_forcing is a CSV file defining a time-varying source inversion,
        !! then the rise_time here is not used (because it is provided in the CSV file)
    
    real(dp) :: ambient_sea_level
        !! Background sea level - taken as a commandline argument.
    
    character(len=charlen), allocatable :: input_elevation_files(:)
        !! Elevation raster files

    character(len=charlen), allocatable :: input_manning_n_files(:)
        !! Manning's n raster files
    
    type(xyz_lines_type) :: breakwalls_forced, inverts_forced
        !! Type with lon/lat/elev representing the tops of breakwalls (or channel inverts, for channel_inverts_forced).
        !! These features are "burned" into the elevation so that on coarser
        !! grids, we don't accidently miss key flow features. 
    
    type(polygons_values_type) :: override_initial_stage_polygons, open_estuary_entrances_polygons
        !! Used to set the initial stage in polygons. This will override the 'ambient_sea_level' and any
        !! static initial condition. 
        !! Also limit the maximum elevation inside estuary entrances, to open them
    !
    !logical, parameter :: close_bunbury_floodgate = .FALSE.
    !    !! If .TRUE., this will insert "closed-gate" elevations along the Bunbury "plug" floodgate.
    !    !! 12/2022, Adrian Brannigan indicated we should run with the floodgate open, as it is not always
    !    !! monitored and reaching the end of its design life.
    !logical, parameter :: burn_gap_in_bunbury_floodgate = .TRUE.
    !    !! If .TRUE., then burn the channel bed (e.g. gap in the floodgate) into the elevation. This is to 
    !    !! ensure the floodgate is 'cleanly' represented in the model, despite the regular grid. It will be applied
    !    !! BEFORE we enforce elevation maxima (e.g. breakwalls)
    
    !character(len=charlen), parameter :: vasse_estuary_polygon = &
    !    '../elevation/initial_stage_40cmAHD/initial_stage_40cmAHD.csv'
    !real(dp), parameter :: initial_stage_40cm_AHD = 0.4_dp
    !    !! Due to floodgates, the Vasse estuary does not have the same sea-level as offshore areas.
    !    !! It is kept at 0.4m AHD (see Shane Martin's 2014 report on storm surge in Busselton)
    
    contains
    
    subroutine parse_inputs(stage_file, multidomain_design_namelists, &
        model_name, run_type, &
        ambient_sea_level, final_time, output_basedir, &
        num_domains_per_nesting_level, nesting_domain_extents_file, &
        nesting_domain_extents)
    
        character(len=charlen), intent(inout) :: stage_file, model_name, run_type, &
            output_basedir, multidomain_design_namelists
        character(len=charlen), allocatable, intent(inout) :: nesting_domain_extents_file(:)
        real(dp), intent(out) :: ambient_sea_level, final_time
        integer(ip), allocatable, intent(inout) :: num_domains_per_nesting_level(:)
        ! Metadata for each nesting level
        type(real_table_type), allocatable, intent(inout) :: nesting_domain_extents(:)
    
        integer(ip), parameter :: narg = 5 ! Expected number of input arguments
        integer(ip) :: n, j
        character(len=charlen) :: ambient_sea_level_char
    
#ifndef SPHERICAL
        ! Ensure spherical coordinates are in use
        write(log_output_unit,*) &
            'Code assumes spherical coordinates, but SPHERICAL is not defined. Compile with -DSPHERICAL'
        call generic_stop
#endif
    
        if(command_argument_count() /= narg) then
            ! Error message
            write(log_output_unit, *) "Incorrect number of commandline arguments - exactly ", narg, " are required:"
            write(log_output_unit, *) "  stage_file (raster with stage perturbation, or csv with time-varying source)"
            write(log_output_unit, *) "  multidomain_design_namelists (file with namelists controlling high-level setup)"
            write(log_output_unit, *) "  model_name (used within the output_folder name) "
            write(log_output_unit, *) "  run_type (either 'test' or 'test_load_balance' or 'full') "
            write(log_output_unit, *) "  ambient_sea_level (background sea-level in m, e.g. 0.0 or 0.95)"
            write(log_output_unit, *) ""
            call generic_stop
        end if
    
        call get_command_argument(1, stage_file)
        call get_command_argument(2, multidomain_design_namelists)
        call get_command_argument(3, model_name)
        call get_command_argument(4, run_type)
        call get_command_argument(5, ambient_sea_level_char)
        read(ambient_sea_level_char, *) ambient_sea_level
    
        write(log_output_unit, *) 'stage_file: ', trim(stage_file)
        write(log_output_unit, *) 'multidomain_design_namelists: ', trim(multidomain_design_namelists)
        write(log_output_unit, *) 'run_type: ', trim(run_type)
        write(log_output_unit, *) 'model_name: ', trim(model_name)
        write(log_output_unit, *) 'ambient_sea_level: ', ambient_sea_level
    
        ! Set non-default values from the namelist if needed.
        call read_multidomain_design_variables(multidomain_design_namelists)
    
        ! If stage_file ends in csv, assume it contains data for a 
        ! multi-unit-source inversion. Otherwise assume it is a raster
        n = len_trim(stage_file)
        stage_file_is_raster = (stage_file( (n-3):n) /= '.csv')
    
        ! Set final_time based on the run_type
        if(run_type == 'test' .or. run_type == 'test_load_balance') then
            final_time = final_time_test_runs
        else if(run_type == 'full') then
            final_time = final_time_full_runs
        else
            write(log_output_unit,*) &
                'run_type must be either "test" or "test_load_balance" or "full"'
            call generic_stop
        end if
    
        ! Informative name for output folder
        output_basedir = './OUTPUTS/' // &
            trim(model_name) // '-' // &
            trim(run_type) // '-' // &
            'ambient_sea_level_' // trim(ambient_sea_level_char)
        write(log_output_unit, *) 'output_basedir: ', trim(output_basedir)
    
        ! Read metadata on each nesting level, and count the domains
        num_domains_per_nesting_level(1) = 1_ip ! Global domain
        do j = 2, NNL ! Avoid j=1 (global domain, treated separately)
            call read_csv_into_array(nesting_domain_extents(j)%metadata, &
                nesting_domain_extents_file(j), skip_header=1_ip)
            num_domains_per_nesting_level(j) = size(nesting_domain_extents(j)%metadata, 2)
        end do
        ! Sanity check. We should have >=1 domain on all nesting levels,
        ! and exactly 1 global domain
        if(any(num_domains_per_nesting_level <= 0) .or. &
           (num_domains_per_nesting_level(1) /= 1 ) ) then
            write(log_output_unit, *) 'ERROR: Invalid value of num_domains_per_nesting_level'
            write(log_output_unit, *) num_domains_per_nesting_level
            call generic_stop
        end if
    
        ! Sanity check periodic EW multidomain on earth
        if(use_periodic_EW_multidomain) then
            if( abs(global_lw(1) - 360.0_dp) > (spacing(360.0_dp)*100.0_dp) ) then
                write(log_output_unit, *) 'ERROR: Requested periodic EW boundary, but global_lw(1) is not close to 360 degrees'
                call generic_stop
            end if 
        end if
    
    end subroutine
    
    subroutine set_initial_conditions(domain, stage_file, global_dt, all_dx_md)
        !! Setup initial conditions and earthquake forcing
    
        class(domain_type), intent(inout):: domain
            ! We initialise domain%U and domain%manning_squared, and
            ! setup the forcing term if rise_time > 0
        character(len=charlen), intent(in) :: stage_file
            ! Raster filename with the Okada stage perturbation, or CSV file with time-varying forcing.
        real(dp), intent(in) :: global_dt
            ! The global time-step. This is used in-case we need to use a forcing rise-time > 0.
            ! In that case, the rk2 timestepping will not include the full forcing if
            ! we begin the forcing at time=0. However if we begin the forcing after global_dt,
            ! it will include the full forcing.
        real(dp), intent(in) :: all_dx_md(:,:,:)
            ! Value of md%all_dx_md -- gives [dx,dy] for all domains in the multidomain on all images.
            ! Needed if we smooth elevation near the domain nesting boundaries
    
        ! Basic setting of elevation -- modified again later by the earthquake forcing.
        call setup_elevation(domain, all_dx_md)
    
        ! Set the earthquake initial condition [stage + elevation perturbation]
        call setup_stage_and_forcing(domain, stage_file, global_dt)

        ! Set the Manning coefficient
        call setup_manning(domain)
    
    end subroutine

    subroutine setup_manning(domain)
        !! Set the Manning coefficient for the domain
        !! Read from a raster file, otherwise set to default value.
        type(domain_type), intent(inout) :: domain

        ! Objects to interpolate from manning files
        type(multi_raster_type):: manning_data
        ! Coordinates for stage lookup
        real(dp), allocatable:: x(:), y(:)
        integer(ip) :: i, j
        character(len=charlen), parameter :: char_format = '(A)'
        
        allocate(x(domain%nx(1)), y(domain%nx(2)))
        x = domain%x

        ! Set the manning value row-by-row, to save memory compared to doing it all at once.
        if(use_periodic_EW_multidomain) then
            ! To avoid x values going outside the rasters, enforce the periodicity of the multidomain.
            where(x < global_ll(1)) x = x + 360.0_dp
            where(x > (global_ll(1) + global_lw(1))) x = x - 360.0_dp
        end if

        ! Set the Manning coefficient from file if available
        if(len_trim(mannning_n_file_list_in_preference_order) > 0) then
            ! Read the manning data files in preference order
            if(.not. allocated(input_manning_n_files)) then
                call read_character_file(mannning_n_file_list_in_preference_order, input_manning_n_files, char_format)
            end if
            call manning_data%initialise(input_manning_n_files)
            write(log_output_unit, *) 'Initialised input manning files:', input_manning_n_files 
           
            ! Interpolate the Manning coefficient
            do j = 1, domain%nx(2)
                y = domain%y(j)

                call manning_data%get_xy(x, &
                    y, &
                    domain%manning_squared(:,j), &
                    domain%nx(1), &
                    bilinear=1_ip, &
                    na_below_limit=raster_na_below_limit)

                ! Square the Manning coefficient
                where (domain%manning_squared(:,j) > 0) &
                    domain%manning_squared(:,j) = domain%manning_squared(:,j)**2
              
                ! Set the Manning coefficient to the defaults where the raster is NA
                if(domain%is_staggered_grid) then
                    where(domain%manning_squared(:,j) <= raster_na_below_limit ) &
                        domain%manning_squared(:,j) = default_manning_staggered_grid**2
                else
                    where(domain%manning_squared(:,j) <= raster_na_below_limit ) &
                        domain%manning_squared(:,j) = default_manning_non_staggered_grid**2
                end if
     
            end do
            call manning_data%finalise()
        end if

        ! Sanity check that we don't have NA values in the Manning coefficient
        if(any(domain%manning_squared /= domain%manning_squared)) then
            write(log_output_unit, *) 'NA values in Manning coefficient at:' 
            do i = 1, domain%nx(1)
                if(domain%manning_squared(i,j) /= domain%manning_squared(i,j)) write(log_output_unit,*) '  ', x(i), y(i)
            end do
            call generic_stop
        end if
    
        write(log_output_unit, *) 'Manning^2 range:'
        write(log_output_unit, *) minval(domain%manning_squared), maxval(domain%manning_squared)
        flush(log_output_unit)

    end subroutine
    
    subroutine setup_elevation(domain, all_dx_md)
        !! Set the elevation of the domain from a set of rasters.
        !! Burn breakwalls into the elevation.
        !! Enforce NS walls, so that the domain is completely closed and we can
        !! reason clearly about energy and mass conservation.
        type(domain_type), intent(inout) :: domain
        !! The domain for which domain%U(:,:,ELV) is set
        real(dp), intent(in) :: all_dx_md(:,:,:)
        !! Copy of md%all_dx_md(:,:,:), giving dx/dy on all domains in multidomain.
        !! If we smooth the elevation near the domain nesting boundaries, then this 
        !! is used to define those boundaries
        type(multi_raster_type):: elevation_data, tidal_adjustment
        !! Objects to interpolate from elevation files
        real(dp), allocatable:: x(:), y(:), adjustment(:)
        !! Coordinates for elevation lookup
        integer(ip) :: j, num_smooth, i
        character(len=charlen), parameter :: char_format = '(A)'
    
        if(.not. allocated(input_elevation_files)) then
            ! Read the elevation data files in preference order
            call read_character_file(swals_elevation_files_in_preference_order, &
                input_elevation_files, char_format)
        end if
        call elevation_data%initialise(input_elevation_files)

        !
        ! Set elevation row-by-row, to save memory compared to doing it all at once.
        !
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        if(use_periodic_EW_multidomain) then
            ! To avoid x values going outside the rasters, enforce the periodicity
            ! of the multidomain.
            where(x <  global_ll(1)                ) x = x + 360.0_dp
            where(x > (global_ll(1) + global_lw(1))) x = x - 360.0_dp
        end if
    
        do j = 1, domain%nx(2)
            y = domain%y(j)
    
            if(allocated(domain%elevation_source_file_index)) then
                ! Set the elevation on a row
                ! Also store the elevation source file index
                call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip, &
                    ! Cells below a threshold are treated as NA
                    na_below_limit=raster_na_below_limit, &
                    raster_index=domain%elevation_source_file_index(:,j))
            else
                ! Set the elevation on a row
                call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip, &
                    ! Cells below a threshold are treated as NA
                    na_below_limit=raster_na_below_limit)
            end if
    
            ! Sanity check
            if(any(domain%U(:,j,ELV) /= domain%U(:,j,ELV))) then
                write(log_output_unit, *) 'NA Elevation at points:' 
                do i = 1, domain%nx(1)
                    if(domain%U(i,j,ELV) /= domain%U(i,j,ELV)) write(log_output_unit,*) '  ', x(i), y(i)
                end do
                call generic_stop
            end if
        end do
        call elevation_data%finalise()
        deallocate(input_elevation_files)

        write(log_output_unit, *) 'Elevation ranges read from rasters:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
    
        if(smooth_elevation_along_nesting_boundaries) then
            call domain%smooth_elevation_near_nesting_fine2coarse_boundaries(all_dx_md)
    
            !! Try smoothing the elevation near a specific site (where a nesting artefact had occurred)
            !! Because of the diffusive nature of local smoothing, we might smooth more times on fine
            !! grids than on coarse grids. This would also prevent smoothing on the coarser domains
            !! num_smooth = nint(1.0_dp/(60.0_dp * 9.0_dp * 6.0_dp) * 1.0_dp/domain%dx(1))**2
            !! Site location and smoothing radius chosen from inspection of nesting artefact.
            ! Kolan River: Checkerboard speed on boundary of two 2nd level nesting cells 
            call domain%smooth_elevation_near_point(&
                number_of_9pt_smooths = 2_ip, &
                pt = [152.200815_dp, -24.666313_dp], &
                smooth_region_radius_meters = 800.0_dp, &
                transition_region_radius_meters = 1000.0_dp)

            write(log_output_unit, *) 'Elevation ranges after smoothing specified points and nesting boundaries:'
            write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        end if

        call setup_breakwalls_and_inverts(domain)

        write(log_output_unit, *) 'Elevation ranges after setting breakwalls and inverts:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        
        ! Add in the tidal adjustment. It must match the model domain's extent with a buffer for get_xy.
        if(tidal_adjustment_file /= "") then
            allocate(adjustment(domain%nx(1)))
            call tidal_adjustment%initialise([tidal_adjustment_file])
            ! Set the adjustment row-by-row, to save memory compared to doing it all at once.
            do j = 1, domain%nx(2)
                y = domain%y(j)
                call tidal_adjustment%get_xy(&
                    x, y, adjustment, domain%nx(1), bilinear=1_ip, na_below_limit=raster_na_below_limit &
                )
                domain%U(:,j,ELV) = domain%U(:,j,ELV) + adjustment
                
                ! Sanity check that we don't have values below the NA threshold in the domain
                if(any(domain%U(:,j,ELV) < raster_na_below_limit)) then
                    write(log_output_unit, *) 'Elevation below NA threshold from tidal adjustment at points:' 
                    do i = 1, domain%nx(1)
                        if(domain%U(i,j,ELV) < raster_na_below_limit) write(log_output_unit,*) '  ', x(i), y(i)
                    end do
                    write(log_output_unit, *) 'Possibly inconsistent extent in tidal adjustment and model domain.', &
                        'The tidal extent may need to be larger than the model domain to get_xy for the model.'
                    call generic_stop
                end if

                ! Sanity check that we don't have NA values in the domain elevation
                if(any(domain%U(:,j,ELV) /= domain%U(:,j,ELV))) then
                    write(log_output_unit, *) 'NA Elevation from tidal adjustment at points:' 
                    do i = 1, domain%nx(1)
                        if(domain%U(i,j,ELV) /= domain%U(i,j,ELV)) write(log_output_unit,*) '  ', x(i), y(i)
                    end do
                    call generic_stop
                end if

            end do
            call tidal_adjustment%finalise()
            deallocate(adjustment)
        end if
        deallocate(x, y)

        write(log_output_unit, *) 'Elevation ranges after setting tidal adjustment:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
          
        if( any(domain%timestepping_method == [character(len=charlen) :: &
                'linear', 'leapfrog_linear_plus_nonlinear_friction']) ) then
            ! Avoid 'very shallow' cells in linear domain, because when linear
            ! sends to the nonlinear domain one can have wet-dry issues
            where(domain%U(:,:,ELV) < ambient_sea_level .and. domain%U(:,:,ELV) > ambient_sea_level-1.0_dp) &
                    domain%U(:,:,ELV) = ambient_sea_level - 1.0_dp
        end if
            
        ! Optionally limit the maximum elevation in user-provided polygons
        if(open_estuary_entrances_polygons_values_file /= "") then
    
            ! Setup the polygons_values data structure the first time this is called
            if(.not. allocated(open_estuary_entrances_polygons%polyvalue)) then
    
                call open_estuary_entrances_polygons%setup_from_csv_with_file_value_columns(&
                    open_estuary_entrances_polygons_values_file, &
                    skip_header_file_value_csv=0_ip, &  ! Assume no header in the main csv
                    skip_header_polygons=1_ip, & ! Assume a single header in the polygon geometries
                    verbose=.TRUE.)
                flush(log_output_unit)
    
            end if
    
            ! Ensure the elevation in the estuary-entrance polygons doesn't exceed the specified value. 
            call open_estuary_entrances_polygons%burn_into_grid(domain%U(:,:,ELV), &
                domain%lower_left, domain%lower_left + domain%lw, burn_type='min')
    
        end if
    
        write(log_output_unit, *) 'Elevation ranges after clipping artefacts:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        flush(log_output_unit)
    
    end subroutine
    
    subroutine setup_breakwalls_and_inverts(domain)
        !! Burn breakwalls (defining elevation maxima) and inverts 
        !! (defining elevation minimum) into the domain's elevation grid, i.e.
        !! domain%U(:,:,ELV)
        type(domain_type), intent(inout) :: domain
            !! domain%U(:,:,ELV) is modified
    
        !! Filenames containing breakwall and invert x/y/z data
        character(len=charlen), allocatable :: breakwall_csv_lon_lat_z(:), &
            invert_csv_lon_lat_z(:)
        character(len=charlen), parameter :: char_format = '(A)'
    
        if(breakwalls_file_list /= "" .and. &
           (.not. allocated(breakwalls_forced%lines))) then
            ! On the first call, setup the breakwalls_forced object.
    
            ! Read the breakwall files
            call read_character_file(breakwalls_file_list, &
                breakwall_csv_lon_lat_z, char_format)
    
            ! Setup breakwalls interpolation object
            call breakwalls_forced%read_from_csv(breakwall_csv_lon_lat_z, &
                skip_header=1_ip)
    
        end if
        
        ! Burn breakwalls into the grid
        if(allocated(breakwalls_forced%lines)) &
            call breakwalls_forced%burn_into_grid(domain%U(:,:,ELV), &
                domain%lower_left, (domain%lower_left + domain%lw), &
                burn_type='max')
    
        if(inverts_file_list /= "" .and. &
           (.not. allocated(inverts_forced%lines))) then
            ! On the first call, setup the inverts_forced object.
    
            ! Read the invert files
            call read_character_file(inverts_file_list, &
                invert_csv_lon_lat_z, char_format)
    
            ! Setup inverts interpolation object
            call inverts_forced%read_from_csv(invert_csv_lon_lat_z, &
                skip_header=1_ip)
    
        end if
    
        ! Burn inverts into the grid
        if(allocated(inverts_forced%lines)) &
            call inverts_forced%burn_into_grid(domain%U(:,:,ELV), &
                domain%lower_left, (domain%lower_left + domain%lw), &
                burn_type='min')
    
    end subroutine
    
    subroutine setup_stage_and_forcing(domain, stage_file, global_dt)
        !! Set the domain's stage in domain%U(:,:,STG). If the rise_time is
        !! greater than zero, then also setup the rise_time forcing.
    
        class(domain_type), intent(inout):: domain
            !! This routine sets domain%U(:,:,STG), and may set
            !! domain%forcing_context_cptr and domain%forcing_subroutine
        character(len=charlen), intent(in) :: stage_file
            !! Raster filename with the Okada stage perturbation, OR csv file with
            !! information on time-varying source inversion.
        real(dp), intent(in) :: global_dt
            !! Model time-step. Useful in case we want to apply a rise-time forcing,
            !! [we should not start the forcing right away, because it can cause
            !!  conservation issues for rk2 -- rather we apply the forcing after a
            !! time-step].
    
        type(multi_raster_type):: stage_data
            ! Object to interpolate from stage file
        real(dp), allocatable:: x(:), y(:)
            ! Coordinates for stage lookup
        integer(ip):: j, n, fid
            ! Convenience integers
        character(len=charlen), allocatable :: stage_files(:)
        real(dp), allocatable :: slips(:), start_times(:), end_times(:)
        integer, parameter :: csv_header_size = 1
    
        if(stage_file /= "") then
            if(stage_file_is_raster .and. rise_time == 0.0_dp) then
                ! Set stage and elevation directly from the raster
                call stage_data%initialise([stage_file])
    
                ! Stage will be corrected later for the tsunami initial condition
                ! and sea-level
                domain%U(:,:,STG) = 0.0_dp
    
                ! Set stage PERTURBATION row-by-row.
                ! This saves memory compared to doing it all at once.
                allocate(x(domain%nx(1)), y(domain%nx(1)))
                x = domain%x
                do j = 1, domain%nx(2)
                    y = domain%y(j)
    
                    ! Set stage PERTURBATION
                    call stage_data%get_xy(x,y, domain%U(:,j,STG), domain%nx(1), &
                        bilinear=1_ip, na_below_limit=raster_na_below_limit)
    
                    ! Clip 'NA' regions (since the stage raster does not cover the
                    ! entire domain). In this case revert to the "0" value, later
                    ! we will add an ambient_sea_level offset
                    where(domain%U(:,j,STG) < (raster_na_below_limit) ) &
                        domain%U(:,j,STG) = 0.0_dp
                end do
                deallocate(x,y)
                call stage_data%finalise()
    
                !write(log_output_unit, *) 'INFLATING STAGE PERTURBATION'
                !domain%U(:,:,STG) = 5.0_dp * domain%U(:,:,STG)
    
                ! Add the stage perturbation to the elevation as well
                domain%U(:,:,ELV) = domain%U(:,:,ELV) + domain%U(:,:,STG)
    
                write(log_output_unit, *) 'Stage perturbation range:'
                write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
    
            else
                ! Here we have a rise-time treatment, implemented with forcing terms.
                ! There are 2 cases
                ! - The stage_file is a raster
                ! - The stage_file is a csv with metadata describing a unit-source inversion,
                !   where the rows are:
                !     "water-surface-unit-source rasters", "slip", 
                !     "forcing_start_time", "forcing_end_time"
                !
                ! Do not start the forcing at start_time=0. Why not?
                ! For 'rk2', it will mean the full forcing is not applied. This
                ! is because of rk2's approach:
                !     [start, advance 2-time-steps to end, then-average(start,end)].
                ! If we don't have timestepping before the forcing start_time, we miss
                ! out on part of the forcing that should have been obtained
                ! (hypothetically if we had evolved from time=-dt to
                ! time=0, we would have included the missing part of the forcing).
                ! A safe workaround is to only start forcing when time >= global_dt,
                ! because no domain will have a larger time-step.
                !
                if(stage_file_is_raster) then
                    ! Mimic the data we would need in the unit source case
                    n = 1
                    stage_files = [stage_file]
                    start_times = [global_dt] ! Avoid forcing in first time-step
                    end_times   = [global_dt + rise_time]! Avoid forcing in first time-step
                    slips       = [1.0_dp] ! Equivalent to using the raw raster values
                else
                    ! Read the data to perform the unit-source summation
                    open(newunit=fid, file=stage_file)
                    n = count_file_lines(fid) - csv_header_size
                    allocate(stage_files(n), slips(n), start_times(n), end_times(n))
                    ! Skip header
                    do j = 1, csv_header_size
                        read(fid, *)
                    end do
                    ! Read data
                    do j = 1, n
                        read(fid, *) stage_files(j), slips(j), start_times(j), end_times(j)
                    end do
                    close(fid)
                    start_times = start_times + global_dt ! Avoid forcing in first time-step
                    end_times = end_times + global_dt ! Avoid forcing in first time-step
                end if
    
                do j = 1, n
                    !write(log_output_unit, *) 'Forcing metadata: ', trim(stage_files(j)), &
                    !    slips(j), start_times(j), end_times(j)
                    call setup_forcing_with_rise_time(domain, stage_files(j), slips(j), &
                        start_times(j), end_times(j))
                end do
    
                deallocate(stage_files, slips, start_times, end_times)
            end if
        end if
    
        ! Adjust for ambient sea-level
        domain%U(:,:,STG) = domain%U(:,:,STG) + ambient_sea_level
        domain%msl_linear = ambient_sea_level ! Influences linear solvers & potential-energy calculation.
    
        ! Optionally override the initial stage in user-provided polygons
        if(override_initial_stage_polygons_values_file /= "") then
    
            ! Setup the polygons_values data structure the first time this is called
            if(.not. allocated(override_initial_stage_polygons%polyvalue)) then
    
                call override_initial_stage_polygons%setup_from_csv_with_file_value_columns(&
                    override_initial_stage_polygons_values_file, &
                    skip_header_file_value_csv=0_ip, &  ! Assume no header in the main csv
                    skip_header_polygons=1_ip, & ! Assume a single header in the polygon geometries
                    verbose=.TRUE.)
                flush(log_output_unit)
    
            end if
    
            ! Set stage
            call override_initial_stage_polygons%burn_into_grid(domain%U(:,:,STG), &
                domain%lower_left, domain%lower_left + domain%lw)
    
        end if
    
        ! Alway need stage >= elevation.
        domain%U(:,:,STG) = max(domain%U(:,:,STG), &
                                domain%U(:,:,ELV) + (minimum_allowed_depth/100.0_dp) )
    
        write(log_output_unit, *) 'Stage range:'
        write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        flush(log_output_unit)
    
    end subroutine
    
    
    subroutine setup_forcing_with_rise_time(domain, stage_file, slip, start_time, end_time)
        !! If the stage-perturbation is applied over a non-zero rise-time
        !! (i.e. not instantaneous), then this routine will set up the forcing
        !! terms.
        class(domain_type), intent(inout):: domain
            !! If the stage-perturbation affects the doman, then this routine
            !! will append a forcing term to that domain.
        character(len=charlen), intent(in) :: stage_file
            !! The unit-source ocean-surface perturbation (as a raster file)
        real(dp), intent(in) :: slip
            !! The slip on the unit source (m). A value of 1.0 will use the raster unchanged.
        real(dp), intent(in) :: start_time
            !! The start time over which slip occurs (seconds)
        real(dp), intent(in) :: end_time
            !! The end time of slip (seconds), with end_time >= start_time
    
        type(forcing_patch_type), pointer :: forcing_context
            ! Use this to apply the Okada forcing with a prescribed rise-time.
            ! It will be attached to the domain via a c_ptr.
        integer(ip):: i0, i1, j0, j1, k0, k1, j
        real(dp), allocatable :: y(:)
            ! Convenience integers
        type(multi_raster_type):: stage_data
    
        call stage_data%initialise([stage_file])
    
        ! Only need to apply the forcing over the stage file spatial domain
        i0 = count(domain%x < stage_data%lowerleft(1)) + 1
        i1 = count(domain%x <= stage_data%upperright(1))
        j0 = count(domain%y < stage_data%lowerleft(2)) + 1
        j1 = count(domain%y <= stage_data%upperright(2))
        k0 = STG
        k1 = ELV ! Forcing the elevation will fail for some FV timestepping methods.
    
        if(i0 <= i1 .and. j0 <= j1) then
            ! The stage data overlaps with this domain, so setup the forcing
    
            if(start_time > end_time) then
                write(log_output_unit, *) "Cannot have start_time > end_time"
                call generic_stop
            end if
    
            allocate(forcing_context)
            call forcing_context%setup(&
                start_time = start_time, end_time = end_time, &
                i0=i0, i1=i1, j0=j0, j1=j1, k0=k0, k1=k1)
    
            ! Read the 'unit-source' forcing row-by-row
            allocate(y(i0:i1))
            do j = j0, j1
                y = domain%y(j) ! Single y coordinate for each row
                call stage_data%get_xy(domain%x(i0:i1), y, &
                    forcing_context%forcing_work(i0:i1, j, STG), (i1-i0+1_ip),&
                    bilinear = 1_ip, na_below_limit=raster_na_below_limit)
                ! Clip NA regions
                where(forcing_context%forcing_work(:,j,STG) < raster_na_below_limit) &
                        forcing_context%forcing_work(:,j,STG) = 0.0_dp
            end do
            deallocate(y)
    
            ! Multiply the 'unit-source' stage forcing by the slip
            forcing_context%forcing_work(:,:,STG) = slip * forcing_context%forcing_work(:,:,STG)
    
            ! No effect on UH/VH terms
            forcing_context%forcing_work(:, :, UH:VH) = 0.0_dp
    
            ! Force the elevation with the stage perturbation. This will
            ! 'fail deliberately' for some finite-volume methods which are not
            ! designed to evolve stage.
            forcing_context%forcing_work(:, :, ELV) = &
                forcing_context%forcing_work(:, :, STG)
    
            ! Ensure any previously defined forcing terms have been stored
            call domain%store_forcing()
            ! Move the forcing_context into the domain
            domain%forcing_context_cptr = c_loc(forcing_context)
            forcing_context => NULL() ! Not deallocated -- the memory is now in the domain
            ! The domain will call apply_forcing_patch during each step
            domain%forcing_subroutine => apply_forcing_patch
            ! Store the forcing (so we can make new forcing terms by redefining domain%forcing_context_cptr and
            ! domain%forcing_subroutine)
            call domain%store_forcing()
    
        end if
    
        call stage_data%finalise()
    
    end subroutine
    
    end module
    