module local_routines 
    !!
    !! Routines for initial conditions and other conveniences, for a 
    !! global-scale model with regional high-resolution domains at various 
    !! sites around Australia.
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

    public :: set_initial_conditions, parse_commandline_args
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
    real(dp), parameter :: ambient_sea_level = 0.0_dp
        !! Background sea level. Not changed in this study

    contains 

    subroutine parse_commandline_args(stage_file, run_type, final_time, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir, highres_regions)
        !!
        !! Convenience routine to read the commandline arguments. As well as
        !! setting the input arguments, this routine sets "rise_time" and "offshore_manning".
        !!

        character(len=charlen), intent(inout) :: stage_file, run_type, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir, highres_regions
        real(dp), intent(inout) :: final_time
        
        character(len=charlen) :: rise_time_char, offshore_manning_char
        integer(ip), parameter :: narg = 8

        if(command_argument_count() /= narg) then
            !! Error message
            write(log_output_unit, *) "Incorrect number of commandline arguments - exactly ", narg, " are required:"
            write(log_output_unit, *) "  stage_file (raster filename with stage perturbation)"
            write(log_output_unit, *) "  rise_time (time in seconds over which stage perturbation is applied)"
            write(log_output_unit, *) "  model_name (used within the output_folder name) "
            write(log_output_unit, *) "  run_type (either 'test' or 'test_load_balance' or 'full') "
            write(log_output_unit, *) "  load_balance_file (file with load balancing metadata, or '' to use crude defaults)"
            write(log_output_unit, *) "  offshore_solver_type ( 'linear_with_manning' or 'linear_with_linear_friction' or "
            write(log_output_unit, *) "    'linear_with_reduced_linear_friction' or 'linear_with_delayed_linear_friction' or"
            write(log_output_unit, *) "    'linear_with_no_friction' or 'leapfrog_nonlinear') )"
            write(log_output_unit, *) &
                "  offshore_manning (manning coefficient for offshore solver if using 'linear_with_manning'/'leapfrog_nonlinear')"
            write(log_output_unit, *) "  highres_regions (either 'none' [no highres domains] or 'australia' [all highres domains] "
            write(log_output_unit, *) "    or 'NSW' [NSW high-res domains] )"
            write(log_output_unit, *) ""
            call generic_stop
        end if

        call get_command_argument(1, stage_file)
        write(log_output_unit, *) 'stage_file: ', trim(stage_file)

        call get_command_argument(2, rise_time_char)
        read(rise_time_char, *) rise_time
        write(log_output_unit, *) 'rise_time: ', rise_time

        call get_command_argument(3, model_name)
            ! Outputs will go in a directory like './OUTPUTS/model_name/....'
        write(log_output_unit, *) 'model_name: ', trim(model_name)

        call get_command_argument(4, run_type)
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
            call get_command_argument(5, load_balance_file)
        end if
        write(log_output_unit, *) 'load_balance_file: ', trim(load_balance_file)

        call get_command_argument(6, offshore_solver_type)
        write(log_output_unit, *) 'offshore_solver_type: ', trim(offshore_solver_type)

        call get_command_argument(7, offshore_manning_char)
            ! The manning coefficient used in the offshore solver. Only used
            ! when offshore_solver_type == 'linear_with_manning'
        read(offshore_manning_char, *) offshore_manning
        write(log_output_unit, *) 'offshore_manning: ', offshore_manning

        call get_command_argument(8, highres_regions)
        write(log_output_unit, *) 'highres_regions: ', highres_regions
        if(.not. any(highres_regions == [character(len=charlen) :: &
            'none', 'australia', 'NSW'])) then
            write(log_output_unit, *) 'highres_regions should be either "none" or "australia" or "NSW"'
            call generic_stop
        end if

        if(run_type == 'test' .or. run_type == 'test_load_balance') then
            ! evolve for a short time for testing purposes, or to create a run 
            ! that can be used to create a load_balance_file.
            final_time = 3600.0_dp  * 0.5_dp
        else if(run_type == 'full') then
            ! evolve for 2.5 days
            final_time = 3600.0_dp * 24.0_dp * 2.5_dp
        end if
        write(log_output_unit, *) 'final_time: ', final_time

        output_basedir = './OUTPUTS/' // &
            trim(model_name) // '-' // &
            'risetime_' // trim(rise_time_char) // '-'  // &
            trim(run_type)   // '-' // &
            trim(offshore_solver_type) // '-' // &
            trim(offshore_manning_char) // '-' // &
            'highres_' // trim(highres_regions)
            ! Informative name for output folder
        write(log_output_unit, *) 'output_basedir: ', trim(output_basedir)

        !    
        ! Sanity checks
        !

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
        else
            if(offshore_manning /= 0.0_dp) then
                write(log_output_unit, *) 'DELIBERATE_STOP (likely input error)'
                write(log_output_unit, *) &
                    'In this study only the linear_with_manning/leapfrog_nonlinear cases have manning coefficient > 0'
                write(log_output_unit, *) &
                    'Although manning > 0 should have no effect, the fact it was specified suggests a mistake'
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
                ! All breakwall files
                '../breakwalls/perth/hillaries2.csv', &
                '../breakwalls/perth/hillarys1.csv', &
                '../breakwalls/perth/southperth2.csv', &
                '../breakwalls/perth/southperth3.csv', &
                '../breakwalls/perth/southperth4.csv', &
                '../breakwalls/perth/southperth5.csv', &
                '../breakwalls/perth/southperth6.csv', &
                '../breakwalls/perth/southperth7.csv', &
                '../breakwalls/perth/southperth.csv', &
                '../breakwalls/perth/swanCanning1.csv', &
                '../breakwalls/perth/swanCanning2.csv', &
                '../breakwalls/perth/swanCanning3.csv', &
                '../breakwalls/perth/swanCanning4.csv', &
                '../breakwalls/perth/swanCanning5.csv', &
                '../breakwalls/portKembla/portKembla2.csv', &
                '../breakwalls/portKembla/portKembla.csv', &
                '../breakwalls/portland/portland1.csv', &
                '../breakwalls/portland/portland2.csv', &
                '../breakwalls/ulladulla/ulladulla2.csv', &
                '../breakwalls/ulladulla/ulladulla.csv' ]

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
            ! !"/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif", &
            "../elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif", &
            ! The GA250m data (It is clipped to nearer Australia, given various 
            ! artefacts in the Pacific discussed in PTHA18 report)
            !"/g/data/w85/tsunami/DATA/ELEVATION/GA250m/ER_Mapper_ers/ausbath_09_v4_ex_ex_106_157_-47_-8.tif", &
            "../elevation/orig/GA250/ausbath_09_v4_ex_ex_106_157_-47_-8.tif", &

            ! 
            ! Offshore DEMs derived for this study, to transition smoothly between the global-scale 
            ! data and the nearshore data
            ! 

            ! Victoria
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/Victoria_smooth_between_GA250m_and_HighresCoastalDem/" // &
            !    "Victoria_smooth_between_GA250_and_Vic10m.tif", &
            "../elevation/derived_for_model/vic/smooth_between_GA250m_and_coast/Victoria_smooth_between_GA250_and_Vic10m.tif", &
            ! NSW
            !'/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/ocean_smooth_transition/' // &
            !    'transition_DEM_149_153_-38_-33.25.tif', &
            '../elevation/derived_for_model/nsw/smooth_between_GA250m_and_coast/transition_DEM_149_153_-38_-33.25.tif', &
            ! SW WA
            !'/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/GA250m_contour_adjusted/' // &
            !    'WA_smooth_between_GA250_and_Multibeam_LADS.tif', & 
            '../elevation/derived_for_model/wa/smooth_between_GA250m_and_coast/WA_smooth_between_GA250_and_Multibeam_LADS.tif', & 

            !
            ! Good quality / high-res data below here
            !

            !
            ! Southern half of NSW Onshore LIDAR.
            ! This is good on land but not good in water areas - so is often 
            ! our lowest-preference good-quality dataset.
            !
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/coastal_lidar/NSW_5m_near_coast/tifs_masked_by_WOFS/CLIP_697266.tif', &
            '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_697266.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/coastal_lidar/NSW_5m_near_coast/tifs_masked_by_WOFS/CLIP_588162.tif', &
            '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588162.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/coastal_lidar/NSW_5m_near_coast/tifs_masked_by_WOFS/CLIP_588152.tif', &
            '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588152.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/coastal_lidar/NSW_5m_near_coast/tifs_masked_by_WOFS/CLIP_588113.tif', &
            '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588113.tif', &

            ! SYDNEY REGION.
            ! Bathymetry from Wilson and Power (Scientific Data paper).
            ! Beware in some regions this is lower resolution than the onshore 
            ! lidar, and we could improve the use of the latter. But for our purposes, 
            ! mostly the nearshore bathy-topo provides good quality data in such areas. 
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/botany_0.0001_gcs.txt.tif', &
            '../elevation/orig/nsw/Sydney_Wilson_and_Power/botany_0.0001_gcs.txt.tif', &
            !"/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/syd_0.0001_gcs.txt.tif", &
            '../elevation/orig/nsw/Sydney_Wilson_and_Power/syd_0.0001_gcs.txt.tif', &
            !"/g/data/w85/tsunami/DATA/ELEVATION/NSW/SydneyRegion/hawkesbury_0.0005_gcs.txt.tif", & 
            '../elevation/orig/nsw/Sydney_Wilson_and_Power/hawkesbury_0.0005_gcs.txt.tif', &

            ! Port Kembla -- will add in breakwalls
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/PortKembla_high_res.tif", & 
            "../elevation/derived_for_model/nsw/highres_coastal_merge/PortKembla_high_res.tif", & 

            ! Jervis Bay
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/JervisBay_high_res.tif", & 
            "../elevation/derived_for_model/nsw/highres_coastal_merge/JervisBay_high_res.tif", & 

            ! Ulladullah -- will add in breakwalls
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/Ulladullah_high_res.tif", & 
            "../elevation/derived_for_model/nsw/highres_coastal_merge/Ulladullah_high_res.tif", & 

            ! Batemans Bay
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/Batemans_high_res.tif", &
            "../elevation/derived_for_model/nsw/highres_coastal_merge/Batemans_high_res.tif", & 

            ! Eden
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/NSW/highres_coastal_merge/Eden_high_res.tif", &
            "../elevation/derived_for_model/nsw/highres_coastal_merge/Eden_high_res.tif", & 

            ! 2018 Bathy-topo -- high-quality offshore/onshore LADS.
            ! In general this is our first-preference NSW dataset
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338878/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338878/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_586366/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_586366/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338874/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338874/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338859/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338859/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_340181/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_340181/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338844/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338844/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338896/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338896/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_338901/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_338901/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_466210/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_466210/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            !'/g/data/w85/tsunami/DATA/ELEVATION/NSW/bathytopo/DATA_466239/NSW Government - DPIE/DEMs/' // &
            !    '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
            '../elevation/orig/nsw/bathytopo/DATA_466239/NSW Government - DPIE/DEMs/' // &
                '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &

            ! VICTORIA -- will add in breakwalls for Portland harbour
            !'/g/data/w85/tsunami/DATA/ELEVATION/Victoria/Vic/Victorian-coast_Bathy_10m_EPSG4326.tif', &
            '../elevation/orig/vic/vcdem2017/Victorian-coast_Bathy_10m_EPSG4326.tif', &

            ! WA -- will add in breakwalls for Hillarys harbour
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_1.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_1.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_2.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_2.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_3.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_3.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_4.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_4.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_5.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_5.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_6.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_6.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_7.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_7.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_8.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_8.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_9.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_9.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_10.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_10.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_11.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_11.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_12.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_12.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_13.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_13.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_14.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_14.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_15.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_15.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_16.tif", &
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_16.tif", &
            !"/g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/WA/Perth_nearshore_tifs/merged_rasters/tile_17.tif" ]
            "../elevation/derived_for_model/wa/nearshore_tifs/tile_17.tif" ]
        
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

        domain%U(:,1:2,ELV) = 100.0_dp
        domain%U(:,(domain%nx(2)-1):domain%nx(2), ELV) = 100.0_dp
            ! Wall boundaries south/north -- for this model the expected radiation in 
            ! North-Pacific/Atlantic will be small, given we are running
            ! for sources in Pacific/Indian Ocean. We close the domain 
            ! so it is easier to reason about the results (e.g. we expect 
            ! mass/energy conservation, noting any friction or rise_time source). 
            ! Note this would be no good if we relied on the North Atlantic (then 
            ! we should use a radiation condition instead).

        call setup_breakwalls(domain)

    end subroutine

    subroutine setup_forcing_with_rise_time(domain, stage_data_lowerleft, &
            stage_data_upperright)
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
        k0 = STG ! Stage only
        k1 = STG

        if(i0 <= i1 .and. j0 <= j1) then
            allocate(forcing_context)
            call forcing_context%setup(&
                start_time = 0.0_dp, end_time = rise_time, &
                i0=i0, i1=i1, j0=j0, j1=j1, k0=k0, k1=k1)

            forcing_context%forcing_work(i0:i1, j0:j1, STG) = &
                domain%U(i0:i1, j0:j1, STG)
                ! Forcing work = stage perturbation.
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

    subroutine setup_stage(domain, stage_file)
        !! Set the domain's stage in domain%U(:,:,STG). If the rise_time is 
        !! greater than zero, then also setup the rise_time forcing.

        class(domain_type), intent(inout):: domain
            !! This routine sets domain%U(:,:,STG), and may set 
            !! domain%forcing_context_cptr and domain%forcing_subroutine
        character(len=charlen), intent(in) :: stage_file
            !! Raster filename with the Okada stage perturbation

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

        ! rise-time treatment here -- only applied if we actually need to.
        if( (rise_time > 0.0_dp) .and. &
            ( (maxval(domain%U(:,:,STG)) > 0.0_dp) .or. &
              (minval(domain%U(:,:,STG)) < 0.0_dp) )) then
            call setup_forcing_with_rise_time(domain, stage_data%lowerleft, &
                stage_data%upperright)
        end if

        call stage_data%finalise()

        domain%U(:,:,STG) = domain%U(:,:,STG) + ambient_sea_level
            ! Adjust for ambient sea-level 
        domain%msl_linear = ambient_sea_level
            ! This is required for some linear solvers, and provides a baseline for potential energy
            ! calculation.

        domain%U(:,:,STG) = max(domain%U(:,:,STG), &
                                domain%U(:,:,ELV) + 1.0e-07_dp)
            ! Alway need stage >= elevation. 'Dry' if (depth < 1.0e-05_dp)

    end subroutine

    subroutine set_initial_conditions(domain, stage_file)
        !! Setup initial conditions. Also setup the stage-forcing if 
        !! rise_time > 0

        class(domain_type), intent(inout):: domain
            ! We initialise domain%U and domain%manning_squared, and 
            ! setup the forcing term if rise_time > 0
        character(len=charlen), intent(in) :: stage_file
            ! Raster filename with the Okada stage perturbation

        call setup_elevation(domain)

        call setup_stage(domain, stage_file)

        if(allocated(domain%manning_squared)) then
            ! Friction 
            domain%manning_squared = nearshore_manning**2 
        end if

        if(any(domain%timestepping_method == [character(len=charlen):: &
            'leapfrog_nonlinear', &
            'leapfrog_linear_plus_nonlinear_friction'])) then
            ! Alternative friction for offshore domains.
            ! Note for other types of offshore domains, the friction coefficient
            ! is set in the main program
            domain%manning_squared = offshore_manning**2
        end if

    end subroutine

end module 


