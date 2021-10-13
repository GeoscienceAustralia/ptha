module local_routines 
    !!
    !! Routines for initial conditions and other conveniences, for a 
    !! global-scale model with regional high-resolution domains at various 
    !! sites around Tonga.
    !!

    use iso_c_binding, only: c_loc, c_double

    !
    ! SWALS objects below
    !
    use global_mod, only: dp, ip, charlen
        ! Default integer / real precision and character length
    use domain_mod, only: domain_type, STG, UH, VH, ELV
        ! domain type and indices of key variables in the array domain%U:
        !     stage (domain%U(:,:,STG))
        !     depth-integrated-velocities (domain%U(:,:,UH:VH))
        !     elevation (domain%U(:,:,ELV) )
    use read_raster_mod, only: multi_raster_type
        ! Interpolation from a collection of rasters (with an order of preference)
    use logging_mod, only: log_output_unit
        ! File unit for logging
    use burn_into_grid_mod, only: xyz_lines_type
        ! burn linear breakwalls into the elevation grid
    use stop_mod, only: generic_stop
        ! Throw errors gracefully in parallel or serial
    use forcing_mod, only: forcing_patch_type, apply_forcing_patch
        ! Used to apply an Okada source over a finite time.

    implicit none

    private

    public :: set_initial_conditions, parse_commandline_args, rise_time
        !! These are called by the main program

    type(xyz_lines_type) :: breakwalls_forced
        !! Type with lon/lat/elev representing the tops of breakwalls. 
        !! These features are "burned" into the elevation so that on coarser 
        !! meshes, we don't accidently miss barriers to flow.

    real(dp) :: rise_time
        !! Okada deformation is applied at a constant rate over the rise_time 
        !! (seconds), starting at time=0.0_dp. Defined with a commandline argument. 
    real(dp) :: offshore_manning
        !! Manning friction for offshore domains (if they use manning 
        !! friction - not all solvers do). Defined with a commandline argument
    real(dp), parameter :: nearshore_manning = 0.03_dp
        !! Manning friction for nonlinear domains. Not changed in this study
    real(dp) :: ambient_sea_level
        !! Background sea level. Allowed to vary for this study.

    contains 

    subroutine parse_commandline_args(stage_file, run_type, final_time, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir, highres_regions, outer_domain_extent, output_style)
        !!
        !! Convenience routine to read the commandline arguments. As well as
        !! setting the input arguments, this routine sets "rise_time" and "offshore_manning".
        !!

        character(len=charlen), intent(inout) :: stage_file, run_type, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir, highres_regions, outer_domain_extent, output_style
        real(dp), intent(inout) :: final_time
        
        character(len=charlen) :: rise_time_char, ambient_sea_level_char, offshore_manning_char
        integer(ip), parameter :: narg = 11

        if(command_argument_count() /= narg) then
            !! Error message
            write(log_output_unit, *) "Incorrect number of commandline arguments - exactly ", narg, " are required:"
            write(log_output_unit, *) "  stage_file (raster filename with stage perturbation)"
            write(log_output_unit, *) "  rise_time (time in seconds over which stage perturbation is applied)"
            write(log_output_unit, *) "  ambient_sea_level (static sea level, e.g. 0.0 for MSL)"
            write(log_output_unit, *) "  model_name (used within the output_folder name) "
            write(log_output_unit, *) "  run_type (either 'test' or 'test_load_balance' or 'full') "
            write(log_output_unit, *) "  load_balance_file (file with load balancing metadata, or '' to use crude defaults)"
            write(log_output_unit, *) "  offshore_solver_type ( 'linear_with_manning' or 'leapfrog_nonlinear') )"
            write(log_output_unit, *) &
                "  offshore_manning (manning coefficient for offshore solver if using 'linear_with_manning'/'leapfrog_nonlinear')"
            write(log_output_unit, *) "  highres_regions (either 'none' [no highres domains] or 'tonga' [all highres domains] "
            write(log_output_unit, *) "  outer_domain_extent (either 'global' or 'regional' or 'pacific')"
            write(log_output_unit, *) &
                "  output_style (either 'animation' or 'few_grids'). If 'animation' we write grids quite often, otherwise rarely."
            write(log_output_unit, *) ""
            call generic_stop
        end if

        call get_command_argument(1, stage_file)
        write(log_output_unit, *) 'stage_file: ', trim(stage_file)

        call get_command_argument(2, rise_time_char)
        read(rise_time_char, *) rise_time
        write(log_output_unit, *) 'rise_time: ', rise_time

        call get_command_argument(3, ambient_sea_level_char)
        read(ambient_sea_level_char, *) ambient_sea_level
        write(log_output_unit, *) 'ambient_sea_level: ', ambient_sea_level

        call get_command_argument(4, model_name)
            ! Outputs will go in a directory like './OUTPUTS/model_name/....'
        write(log_output_unit, *) 'model_name: ', trim(model_name)

        call get_command_argument(5, run_type)
        write(log_output_unit, *) 'run_type: ', trim(run_type)
            ! run_type will control whether we do a long run, or a short run 
            ! without load balancing, or a short run with load balancing
        if(.not. any(run_type == [character(len=charlen) :: &
                     'test', 'test_load_balance', 'full']) ) then
            ! An invalid value of run_type was passed
            write(log_output_unit,*) &
                'run_type must be either "test" or "test_load_balance" or "full"'
            call generic_stop
        end if

        if(run_type == 'test') then
            ! Use default load balancing when run_type=='test'. Main purpose is 
            ! for short runs that are then used to create a load_balance_file
            load_balance_file = ''
        else
            call get_command_argument(6, load_balance_file)
        end if
        write(log_output_unit, *) 'load_balance_file: ', trim(load_balance_file)

        call get_command_argument(7, offshore_solver_type)
        write(log_output_unit, *) 'offshore_solver_type: ', trim(offshore_solver_type)

        call get_command_argument(8, offshore_manning_char)
            ! The manning coefficient used in the offshore solver. Only used
            ! when offshore_solver_type == 'linear_with_manning'
        read(offshore_manning_char, *) offshore_manning
        write(log_output_unit, *) 'offshore_manning: ', offshore_manning

        call get_command_argument(9, highres_regions)
        write(log_output_unit, *) 'highres_regions: ', highres_regions
        if(.not. any(highres_regions == [character(len=charlen) :: &
            'none', 'tonga'])) then
            write(log_output_unit, *) 'highres_regions should be either "none" or "tonga"'
            call generic_stop
        end if

        call get_command_argument(10, outer_domain_extent)
        write(log_output_unit, *) 'outer_domain_extent: ', outer_domain_extent
        if(.not. any(outer_domain_extent == [character(len=charlen) :: &
            'global', 'regional', 'pacific'])) then
            write(log_output_unit, *) 'outer_domain_extent should be one of "global" or "regional" or "pacific"'
            call generic_stop
        end if

        call get_command_argument(11, output_style)
        write(log_output_unit, *) 'output_style: ', output_style
        if(.not. any(output_style == [character(len=charlen) :: 'animation', 'few_grids'])) then
            write(log_output_unit, *) 'output_style should be one of "animation" or "few_grids"'
            call generic_stop
        end if


        if(run_type == 'test' .or. run_type == 'test_load_balance') then
            ! evolve for a short time for testing purposes, or to create a run 
            ! that can be used to create a load_balance_file.
            final_time = 3600.0_dp  * 1.0_dp/12.0_dp

        else if(run_type == 'full') then
            if(outer_domain_extent == 'regional') then
                ! The regional domain does not need to evolve for so long
                final_time = 3600.0_dp * 5.0_dp 
                !final_time = 3600.0_dp * 1.0_dp 
            else
                ! For Pacific-wide or global domains, give it some time.
                final_time = 3600.0_dp * 24.0_dp 
            end if

        end if

        write(log_output_unit, *) 'final_time: ', final_time

        output_basedir = './OUTPUTS/' // &
            trim(model_name) // '-' // &
            'risetime_' // trim(rise_time_char) // '-'  // &
            'ambientsealevel_' // trim(ambient_sea_level_char) // '-'  // &
            trim(run_type)   // '-' // &
            trim(offshore_solver_type) // '-' // &
            trim(offshore_manning_char) // '-' // &
            'highres_' // trim(highres_regions)
            ! Informative name for output folder
        write(log_output_unit, *) 'output_basedir: ', trim(output_basedir)

        !    
        ! Sanity checks
        !

        if(abs(ambient_sea_level) > 1.0_dp) then
            write(log_output_unit, *) &
                'DELIBERATE STOP (abs(ambient_sea_level) exceeds 1m -- not expected in this application) ', &
                ambient_sea_level, __LINE__
            call generic_stop
        end if

        if(offshore_manning < 0.0_dp .or. offshore_manning > 0.04_dp) then
            write(log_output_unit, *) &
                'DELIBERATE STOP (MANNING COEFF IS NEGATIVE OR VERY HIGH - typo?) ', &
                offshore_manning, __LINE__
            call generic_stop
        end if

        if(offshore_solver_type == 'linear_with_manning' .or. offshore_solver_type == 'leapfrog_nonlinear') then
            if(offshore_manning == 0.0_dp) then
                write(log_output_unit, *) 'DELIBERATE_STOP (likely input error)'
                write(log_output_unit, *) &
                    'In this study the linear_with_manning/leapfrog_nonlinear cases are intended to have manning coefficient > 0'
                write(log_output_unit, *) &
                    'Although 0.0 manning is not logically wrong, for this study it suggests a mistake'
                call generic_stop
            end if
        end if

    end subroutine

    subroutine setup_breakwalls(domain)
        !! Burn breakwalls into the domain's elevation grid (i.e. 
        !! domain%U(:,:,ELV))
        type(domain_type), intent(inout) :: domain
            !! domain%U(:,:,ELV) is modified

        character(len=charlen), allocatable :: breakwall_csv_lon_lat_z(:)
            !! Filenames containing breakwall x/y/z data, used to setup
            !! 'breakwalls_forced' if it is not already setup

        if(.not. allocated(breakwalls_forced%lines)) then
            ! On the first call, setup the breakwalls_forced object.

            breakwall_csv_lon_lat_z = [character(len=charlen) :: & 
                '../elevation/walls/coast01/coast01.csv', &
                '../elevation/walls/coast02/coast02.csv', &
                '../elevation/walls/coast03/coast03.csv', &
                '../elevation/walls/coast04/coast04.csv', &
                '../elevation/walls/coast05/coast05.csv', &
                '../elevation/walls/coast06/coast06.csv', &
                '../elevation/walls/coast07/coast07.csv', &
                '../elevation/walls/coast08/coast08.csv', &
                '../elevation/walls/coast09/coast09.csv', &
                '../elevation/walls/coast10/coast10.csv', &
                '../elevation/walls/coast11/coast11.csv', &
                '../elevation/walls/coast12/coast12.csv', &
                '../elevation/walls/coast13/coast13.csv', &
                '../elevation/walls/coast14/coast14.csv', &
                '../elevation/walls/ridge01/ridge01.csv', &
                '../elevation/walls/road01/road01.csv', &
                '../elevation/walls/road02/road02.csv', &
                '../elevation/walls/road03/road03.csv', &
                '../elevation/walls/road04/road04.csv', &
                '../elevation/walls/road05/road05.csv', &
                '../elevation/walls/road06/road06.csv', &
                '../elevation/walls/road07/road07.csv', &
                '../elevation/walls/road08/road08.csv', &
                '../elevation/walls/road09/road09.csv', &
                '../elevation/walls/road10/road10.csv', &
                '../elevation/walls/road11/road11.csv', &
                '../elevation/walls/road12/road12.csv', &
                '../elevation/walls/road13/road13.csv', &
                '../elevation/walls/road14/road14.csv', &
                '../elevation/walls/road15/road15.csv']

            call breakwalls_forced%read_from_csv(breakwall_csv_lon_lat_z, &
                skip_header=1_ip)

        end if

        call breakwalls_forced%burn_into_grid(domain%U(:,:,ELV), &
            domain%lower_left, (domain%lower_left + domain%lw), &
            burn_type='max')

    end subroutine

    subroutine setup_elevation(domain)
        !! Set the elevation of the domain from a set of rasters. 
        !! Burn breakwalls into the elevation. 
        !! Enforce NS walls, so that the domain is completely closed and we can 
        !! reason clearly about energy and mass conservation.
        type(domain_type), intent(inout) :: domain
            !! The domain for which domain%U(:,:,ELV) is set

        character(len=charlen), allocatable :: input_elevation_files(:)
            !! Names of elevation files
        type(multi_raster_type):: elevation_data
            !! Objects to interpolate from elevation files
        real(dp), allocatable:: x(:), y(:)
            !! Coordinates for elevation lookup
        integer(ip) :: j

        input_elevation_files = [character(len=charlen) :: &
            ! Initialised in order 'low-preferece to high-preference'
            ! We have to reverse it later as required for the multi_raster_type
             
            !
            ! Large scale background raster data
            !

            ! The PTHA18 DEM
            "../elevation/for_model/merged_gebco_ga250_dem_patched.tif", &
            ! GEBCO2014, with full 30s resolution
            "../elevation/for_model/gebco2014_cropped.tif", &

            ! The regional tongatapu DEM
            "../elevation/for_model/tongatapu_and_offshore_patched.tif", &

            ! Various transition DEMs
            "../elevation/for_model/0.0008mgrid_patched.tif", &
            "../elevation/for_model/0.0002mgrid_patched.tif", &

            ! Highest-res bathy and onshore lidar around Tonga
            "../elevation/for_model/tonga_bathy_moasaic_patched.tif", &
            "../elevation/for_model/tonga_lidar_mosaic_dem_patched.tif", &

            ! Highest-res bathy and onshore lidar around Lifuka. Beware, we do not have good
            ! "transition to offshore" data in this region, so I would not assume the results are
            ! good quality
            "../elevation/for_model/LifukaBathy2011_wgs84_shifted_by_360lon.tif", &
            "../elevation/for_model/LifukaLidar2011_wgs84_shifted_by_360lon.tif" &
            ]
        
        input_elevation_files = input_elevation_files(size(input_elevation_files):1:-1)
            ! Reverse file order (so high-preference files are first) as required by the
            ! multi_raster type
        call elevation_data%initialise(input_elevation_files)

        !
        ! Set elevation row-by-row.
        ! This saves memory compared to doing it all at once.
        !
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), &
                bilinear=1_ip, na_below_limit=-1.0e+20_dp)
                ! On occasion the raster NA cell values are very slightly 
                ! different to the stated NA value in the raster metadata. I 
                ! presume this is because of precision adjustments to the file 
                ! during processing.
                !
                ! To catch those cases we interpret (numbers < -1.0e+20_dp) as 
                ! NA. It works for this data because the offending NA values 
                ! are denoted by large negative numbers.
                ! 
                ! If the NA value is smaller than -1.0e+20, then there is no 
                ! problem so long as the NA value is not corrupted.
        end do
        call elevation_data%finalise()
        deallocate(input_elevation_files)
        deallocate(x, y)

        if( any(domain%timestepping_method == [character(len=charlen) :: &
                'linear', 'leapfrog_linear_plus_nonlinear_friction']) ) then
            ! Avoid 'very shallow' cells in linear domain, because when linear 
            ! sends to the nonlinear domain one can have wet-dry issues
            where(domain%U(:,:,ELV) < 0.0_dp .and. domain%U(:,:,ELV) > -1.0_dp) &
                    domain%U(:,:,ELV) = -1.0_dp
        end if

        write(log_output_unit, *) 'Elevation ranges before setting breakwals:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        !domain%U(:,1:2,ELV) = 100.0_dp
        !domain%U(:,(domain%nx(2)-1):domain%nx(2), ELV) = 100.0_dp
        !    ! Wall boundaries south/north -- for this model the expected radiation in 
        !    ! North-Pacific/Atlantic will be small, given we are running
        !    ! for sources in Pacific/Indian Ocean. We close the domain 
        !    ! so it is easier to reason about the results (e.g. we expect 
        !    ! mass/energy conservation, noting any friction or rise_time source). 
        !    ! Note this would be no good if we relied on the North Atlantic (then 
        !    ! we should use a radiation condition instead).

        call setup_breakwalls(domain)

        write(log_output_unit, *) 'Elevation ranges AFTER setting breakwals:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        flush(log_output_unit)

    end subroutine

    subroutine setup_forcing_with_rise_time(domain, stage_data_lowerleft, &
            stage_data_upperright, global_dt)
        !! If the stage-perturbation is applied over a non-zero rise-time 
        !! (i.e. not instantaneous), then this routine will set up the forcing
        !! terms.
        class(domain_type), intent(inout):: domain
            !! If the stage-perturbation affects the doman, then this routine 
            !! will affect domain%U(:,:,STG), domain%forcing_context_cptr and 
            !! domain%forcing_subroutine
        real(c_double), intent(in) :: stage_data_lowerleft(2), &
            stage_data_upperright(2)
            !! xlimit and ylimit of the stage data: the forcing will only be
            !! applied within this region. 
        real(dp), intent(in) :: global_dt
            !! The outer grid timestep -- don't start forcing until we have
            !! taken one time-step, to avoid any conservation issues with rk2

        type(forcing_patch_type), pointer :: forcing_context
            ! Use this to apply the Okada forcing with a prescribed rise-time.
            ! It will be attached to the domain via a c_ptr.
        integer(ip):: i0, i1, j0, j1, k0, k1
            ! Convenience integers

        ! Only need to apply the forcing over the stage file spatial domain 
        i0 = count(domain%x < stage_data_lowerleft(1)) + 1
        i1 = count(domain%x <= stage_data_upperright(1))
        j0 = count(domain%y < stage_data_lowerleft(2)) + 1
        j1 = count(domain%y <= stage_data_upperright(2))
        k0 = STG
        k1 = ELV ! Forcing the elevation will fail for some FV timestepping methods.

        if(i0 <= i1 .and. j0 <= j1) then

            allocate(forcing_context)
            ! Do not start the forcing at start_time=0. Why not? 
            ! For 'rk2', it will mean the full forcing is not applied. This
            ! is because of rk2's approach [start, advance 2-time-steps to end, then-average(start,end)].
            ! If we don't have timestepping before the forcing start_time, we miss out on part of
            ! of the forcing that should have been obtained (hypothetically if we evolved from time=-dt to 
            ! time=0, we would have included some of the forcing -- this is the part we miss).
            call forcing_context%setup(&
                start_time = global_dt, end_time = rise_time+global_dt, &
                i0=i0, i1=i1, j0=j0, j1=j1, k0=k0, k1=k1)

            forcing_context%forcing_work(i0:i1, j0:j1, STG) = &
                domain%U(i0:i1, j0:j1, STG)
                ! Forcing the stage with the stage perturbation.
            forcing_context%forcing_work(i0:i1, j0:j1, UH:VH) = 0.0_dp
                ! Do not force the UH/VH terms

            forcing_context%forcing_work(i0:i1, j0:j1, ELV) = &
                domain%U(i0:i1, j0:j1, STG)
                ! Force the elevation with the stage perturbation. This will
                ! 'fail deliberately' for some finite-volume methods which are not
                ! designed to evolve stage.
            domain%U(i0:i1, j0:j1, STG) = 0.0_dp
                ! Zero the initial stage, since now we apply it over a 
                ! rise-time

            domain%forcing_context_cptr = c_loc(forcing_context)
            forcing_context => NULL()
                ! Move the forcing_context into the domain
            domain%forcing_subroutine => apply_forcing_patch
                ! The domain will call apply_forcing_patch during each step
        end if


    end subroutine

    subroutine setup_stage(domain, stage_file, global_dt)
        !! Set the domain's stage in domain%U(:,:,STG). If the rise_time is 
        !! greater than zero, then also setup the rise_time forcing.

        class(domain_type), intent(inout):: domain
            !! This routine sets domain%U(:,:,STG), and may set 
            !! domain%forcing_context_cptr and domain%forcing_subroutine
        character(len=charlen), intent(in) :: stage_file
            !! Raster filename with the Okada stage perturbation
        real(dp), intent(in) :: global_dt
            !! Model time-step. Useful in case we want to apply a rise-time forcing,
            !! [we should not start the forcing right away, because it can cause
            !!  conservation issues for rk2 -- rather we apply the forcing after a
            !! time-step].

        type(multi_raster_type):: stage_data
            ! Object to interpolate from stage file
        real(dp), allocatable:: x(:), y(:)
            ! Coordinates for stage lookup
        type(forcing_patch_type), pointer :: forcing_context
            ! Use this to apply the Okada forcing with a prescribed rise-time
        integer(ip):: j
            ! Convenience integers

        domain%U(:,:,STG) = 0.0_dp
            ! This will be corrected later for the tsunami initial condition 
            ! and sea-level

        call stage_data%initialise([stage_file])

        ! Set stage PERTURBATION row-by-row.
        ! This saves memory compared to doing it all at once.
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call stage_data%get_xy(x,y, domain%U(:,j,STG), domain%nx(1), &
                bilinear=1_ip, na_below_limit=-1.0e+20_dp)
                ! Set stage PERTURBATION
            where(domain%U(:,j,STG) < (-1.0e+20_dp) ) domain%U(:,j,STG) = 0.0_dp
                ! Clip 'NA' regions (since the stage raster does not cover the 
                ! entire domain). In this case revert to the "0" value, later 
                ! we will add an ambient_sea_level offset
        end do
        deallocate(x,y)
        call stage_data%finalise()

        write(log_output_unit, *) 'Stage perturbation range:'
        write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))

        ! rise-time treatment here -- only applied if we actually need to.
        if( (rise_time > 0.0_dp) .and. &
            ( (maxval(domain%U(:,:,STG)) > 0.0_dp) .or. &
              (minval(domain%U(:,:,STG)) < 0.0_dp) )) then
            ! This will setup a time-forcing of stage and elevation, and
            ! then set the stage perturbation to zero
            call setup_forcing_with_rise_time(domain, stage_data%lowerleft, &
                stage_data%upperright, global_dt)
        else
            ! Add the stage perturbation to the elevation as well
            domain%U(:,:,ELV) = domain%U(:,:,ELV) + domain%U(:,:,STG)
        end if

        domain%U(:,:,STG) = domain%U(:,:,STG) + ambient_sea_level
            ! Adjust for ambient sea-level 
        domain%msl_linear = ambient_sea_level
            ! This is required for some linear solvers, and provides a baseline for potential energy
            ! calculation.

        domain%U(:,:,STG) = max(domain%U(:,:,STG), &
                                domain%U(:,:,ELV) + 1.0e-07_dp)
            ! Alway need stage >= elevation. 'Dry' if (depth < 1.0e-05_dp)

        write(log_output_unit, *) 'Stage range:'
        write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        flush(log_output_unit)

    end subroutine

    subroutine set_initial_conditions(domain, stage_file, global_dt)
        !! Setup initial conditions. Also setup the stage-forcing if 
        !! rise_time > 0

        class(domain_type), intent(inout):: domain
            ! We initialise domain%U and domain%manning_squared, and 
            ! setup the forcing term if rise_time > 0
        character(len=charlen), intent(in) :: stage_file
            ! Raster filename with the Okada stage perturbation
        real(dp), intent(in) :: global_dt
            ! The global time-step. This is used in-case we choose to use a rise-time >0.
            ! In that case, the rk2 timestepping will not include the full forcing if
            ! we begin the forcing at time=0. However if we begin the forcing after some time,
            ! it will include the full forcing.

        call setup_elevation(domain)

        call setup_stage(domain, stage_file, global_dt)

        if(allocated(domain%manning_squared)) then
            ! Friction 
            if(any(domain%timestepping_method == [character(len=charlen):: &
                'leapfrog_nonlinear', &
                'leapfrog_linear_plus_nonlinear_friction'])) then
                ! Alternative friction for offshore domains.
                ! Note for other types of offshore domains, the friction coefficient
                ! is set in the main program
                domain%manning_squared = offshore_manning**2
            else 
                domain%manning_squared = nearshore_manning**2 
            end if
        end if

    end subroutine

end module 


