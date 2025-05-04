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
    use global_mod, only: dp, ip, charlen, minimum_allowed_depth
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
    ! Used to read the elevation data filenames in preference order
    use file_io_mod, only: count_file_lines


    implicit none

    private

    public :: set_initial_conditions, parse_commandline_args, highres_regions, included_regions
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
    character(len=charlen) :: highres_regions
        !! Denotes the configuration of high-res regions in the model. Either 'none' or 'australia' or 'perth' or 'australiaSWWA' or 'SWWA'
    real(dp), parameter :: raster_na_below_limit = -1.0e+20_dp
        !! Interpret raster values below this as NA

    character(len=charlen), allocatable :: included_regions(:)
        ! Define which high-res regions should be included in a more convenient way


    contains 

    subroutine parse_commandline_args(stage_file, run_type, final_time, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir)
        !!
        !! Convenience routine to read the commandline arguments. As well as
        !! setting the input arguments, this routine sets "rise_time" and "offshore_manning".
        !!

        character(len=charlen), intent(inout) :: stage_file, run_type, &
            model_name, load_balance_file, offshore_solver_type, &
            output_basedir
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
            write(log_output_unit, *) "  highres_regions (either 'none' [no highres domains] or 'australia' [like nsw+perth] "
            write(log_output_unit, *) "    or 'NSW' [NSW high-res domains] or 'perth' [Perth high-res domains] or "
            write(log_output_unit, *) "    'SWWA' [South West WA highres domains] or 'australiaSWWA' [like nsw+SWWA] or"
            write(log_output_unit, *) "    'NWWA' [Geraldton + Exmouth to just north of Dampier] or 'WA' [all WA] or 'australiaWA'[like nsw+WA]"
            write(log_output_unit, *) ""
            call generic_stop
        end if

        call get_command_argument(1, stage_file)
        write(log_output_unit, *) 'stage_file: ', trim(stage_file)

        call get_command_argument(2, rise_time_char)
        read(rise_time_char, *) rise_time
        write(log_output_unit, *) 'rise_time: ', rise_time

        call get_command_argument(3, model_name)
            ! Outputs will go in a directory like './OUTPUTS_SCRATCH/model_name/....'
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
            'none', 'australia', 'NSW', 'perth', 'SWWA', 'australiaSWWA', 'NWWA', 'WA', 'australiaWA'])) then
            write(log_output_unit, *) 'highres_regions should be either "none" or "australia" or '
            write(log_output_unit, *) '  "NSW" or "perth" or "SWWA" or "australiaSWWA" or "NWWA" or "WA" or "australiaWA"'
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

        !! Note this was changed over time when doing the runs
        !output_basedir = './OUTPUTS/' // &
        !output_basedir = './OUTPUTS_SCRATCH/' // &
        !output_basedir = './OUTPUTS_2021_march_sources/' // &
        !output_basedir = './OUTPUTS_2022_new_events/' // &
        !output_basedir = './OUTPUTS_2023_new_events/' // &
        !output_basedir = './OUTPUTS_2024_extend_SWWA/' // &
        !output_basedir = './OUTPUTS_new_validation_events/' // &
        !output_basedir = './OUTPUTS_2025_testing/' // &
        !output_basedir = './OUTPUTS_2025_extend_WA/' // &
        output_basedir = './OUTPUTS_2025_NWWA/' // &
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

            !if(highres_regions == 'SWWA' .or. highres_regions == 'australiaSWWA' .or. &
            !    highres_regions == 'WA' .or. highres_regions == 'australiaWA') then
            if(any(included_regions == 'swwa')) then
                ! New breakwall near Port Geographe, only in the new models for backward compatability
                breakwall_csv_lon_lat_z = [ character(len=charlen) :: breakwall_csv_lon_lat_z, &
                    '../breakwalls/portgeographe_entrance/portgeographe_entrance.csv']
            end if

            !if(highres_regions == 'NWWA' .or. highres_regions == 'WA' .or. highres_regions == 'australiaWA') then
            if(any(included_regions == 'nwwa')) then
                ! Here the elevation data isn't as good, and we help give definition to key breakwalls
                breakwall_csv_lon_lat_z = [character(len=charlen):: breakwall_csv_lon_lat_z, &
                    '../breakwalls/exmouth_entrance_estimate/exmouth_entrance_breakwall_estimate.csv', &
                    '../breakwalls/onslow_creek_entrance_estimate/onslow_creek_entrance_breakwall_estimate.csv', &
                    '../breakwalls/dampier_causeway/dampier_causeway_breakwall_estimate.csv']
            end if

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

        ! For backwards compatability we use different datasets when
        ! highres_regions includes the more extensive SWWA domains.
        if(any(included_regions == 'swwa' .or. included_regions == 'nwwa')) then
            !
            ! These datasets are similar to data in Davies et al 2020, with updates in south-west WA, and in some cases NWWA
            !
            input_elevation_files = [character(len=charlen) :: &
              ! Initialised in order 'low-preferece to high-preference'
              ! We have to reverse it later as required for the multi_raster_type
               
              !
              ! Large scale background raster data
              !

              ! The PTHA18 DEM
              "../elevation/derived_for_model/global/ptha18/merged_gebco_ga250_dem_patched.tif", &
              ! The GA250m data (It is clipped to nearer Australia, given various 
              ! artefacts in the Pacific discussed in PTHA18 report)
              "../elevation/orig/GA250/ausbath_09_v4_ex_ex_106_157_-47_-8.tif", &

              ! 
              ! Offshore DEMs derived for this study, to transition smoothly between the global-scale 
              ! data and the nearshore data
              ! 

              ! Victoria
              "../elevation/derived_for_model/vic/smooth_between_GA250m_and_coast/Victoria_smooth_between_GA250_and_Vic10m.tif", &
              ! NSW
              '../elevation/derived_for_model/nsw/smooth_between_GA250m_and_coast/transition_DEM_149_153_-38_-33.25.tif', &

              ! Updated transition DEM for SW WA
              '../elevation/derived_for_model_SWWA_update/' // &
                  'Transition_DEM_near_Perth_Oct2021/WA_smooth_between_GA250_and_MultibeamNearshore.tif', &

              !
              ! Good quality / high-res data below here
              !

              !
              ! Southern half of NSW Onshore LIDAR.
              ! This is good on land but not good in water areas - so is often 
              ! our lowest-preference good-quality dataset.
              !
              '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_697266.tif', &
              '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588162.tif', &
              '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588152.tif', &
              '../elevation/derived_for_model/nsw/onshore_lidar_masked_by_WOFS/CLIP_588113.tif', &

              ! SYDNEY REGION.
              ! Bathymetry from Wilson and Power (Scientific Data paper).
              ! Beware in some regions this is lower resolution than the onshore 
              ! lidar, and we could improve the use of the latter. But for our purposes, 
              ! mostly the nearshore bathy-topo provides good quality data in such areas. 
              '../elevation/orig/nsw/Sydney_Wilson_and_Power/botany_0.0001_gcs.txt.tif', &
              '../elevation/orig/nsw/Sydney_Wilson_and_Power/syd_0.0001_gcs.txt.tif', &
              '../elevation/orig/nsw/Sydney_Wilson_and_Power/hawkesbury_0.0005_gcs.txt.tif', &

              ! Port Kembla -- will add in breakwalls
              "../elevation/derived_for_model/nsw/highres_coastal_merge/PortKembla_high_res.tif", & 

              ! Jervis Bay
              "../elevation/derived_for_model/nsw/highres_coastal_merge/JervisBay_high_res.tif", & 

              ! Ulladullah -- will add in breakwalls
              "../elevation/derived_for_model/nsw/highres_coastal_merge/Ulladullah_high_res.tif", & 

              ! Batemans Bay
              "../elevation/derived_for_model/nsw/highres_coastal_merge/Batemans_high_res.tif", & 

              ! Eden
              "../elevation/derived_for_model/nsw/highres_coastal_merge/Eden_high_res.tif", & 

              ! 2018 Bathy-topo -- high-quality offshore/onshore LADS.
              ! In general this is our first-preference NSW dataset
              '../elevation/orig/nsw/bathytopo/DATA_338878/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_586366/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_338874/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_338859/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_340181/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_338844/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_338896/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_338901/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_466210/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &
              '../elevation/orig/nsw/bathytopo/DATA_466239/NSW Government - DPIE/DEMs/' // &
                  '5 Metre/NSW_Marine_5m_TopoBathy_DEM.tif', &

              ! VICTORIA -- will add in breakwalls for Portland harbour
              '../elevation/orig/vic/vcdem2017/Victorian-coast_Bathy_10m_EPSG4326.tif', &
              ! Updated and more extensive SWWA merged raster tiles
              '../elevation/derived_for_model_SWWA_update/SWWA_nearshore_tifs_2021/merged_rasters/all_tiles.vrt', &
              ! Patch at Bunbury
              '../elevation/derived_for_model_SWWA_update/WA_Bunbury_revised_November2022/Bunbury_patch_tile.tif', &
              ! Patch at Port Geographe
              '../elevation/derived_for_model_SWWA_update/WA_Busselton_PortGeographe_merged_tile15/patch_near_portgeographe.tif']

            if(any(included_regions == 'nwwa')) then
                ! Add in NWWA datasets in a backward-compatible way
                input_elevation_files = [character(len=charlen) :: input_elevation_files, &
                    '../elevation/derived_for_model_NWWA_update/North_West_Shelf_DEM_v2_Bathymetry_2020_30m_MSL_cog_WGS84.tif', &
                    '../elevation/derived_for_model_NWWA_update/Exmouth_5m_wgs84.tif', &
                    '../elevation/derived_for_model_NWWA_update/Dampier_5m_WGS84.tif', &
                    '../elevation/derived_for_model_NWWA_update/Onslow_5m_wgs84.tif', &
                    '../elevation/derived_for_model_NWWA_update/PointSamson_5m_WGS84.tif']

            end if              
        else 
            !
            ! Original datasets from Davies et al. (2020)
            !
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
        end if
    
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

    subroutine setup_stage_and_forcing(domain, stage_file, global_dt)
        !! Set the domain's stage in domain%U(:,:,STG). If the rise_time is
        !! greater than zero, then also setup the rise_time forcing.

        type(domain_type), intent(inout):: domain
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
        logical :: stage_file_is_raster

        ! If stage_file ends in csv, assume it contains data for a 
        ! multi-unit-source inversion. Otherwise assume it is a raster
        n = len_trim(stage_file)
        stage_file_is_raster = (stage_file( (n-3):n) /= '.csv')

        if(stage_file /= "") then
            if(stage_file_is_raster .and. rise_time == 0.0_dp) then
                ! Set stage directly from the raster. We don't adjust elevation
                ! so as to maintain backward compatibility with already run simulations.
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

                !!
                !! Can add the stage perturbation to the elevation as well, like
                !domain%U(:,:,ELV) = domain%U(:,:,ELV) + domain%U(:,:,STG)
                !! That is good practice but not done here for backward compatibility with earlier runs.

                !write(log_output_unit, *) 'Stage perturbation range:'
                !write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))

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

        ! Alway need stage >= elevation.
        domain%U(:,:,STG) = max(domain%U(:,:,STG), &
                                domain%U(:,:,ELV) + (minimum_allowed_depth/100.0_dp) )

        !write(log_output_unit, *) 'Stage range:'
        !write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        !flush(log_output_unit)

    end subroutine


    subroutine setup_forcing_with_rise_time(domain, stage_file, slip, start_time, end_time)
        !! If the stage-perturbation is applied over a non-zero rise-time
        !! (i.e. not instantaneous), then this routine will set up the forcing
        !! terms.
        !! In this case we also force the elevation. Note the elevation is not perturbed for instantaneous
        !! initial conditions, purely for backward compatibility with earlier runs. So there will be some
        !! difference between using a rise time of zero, versus a vanishingly small rise time.
        type(domain_type), intent(inout):: domain
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
            ! designed to do it.
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

    subroutine set_initial_conditions(domain, stage_file, global_dt)
        !! Setup initial conditions. Also setup the stage-forcing if 
        !! rise_time > 0

        type(domain_type), intent(inout):: domain
            ! We initialise domain%U and domain%manning_squared, and 
            ! setup the forcing term if rise_time > 0
        character(len=charlen), intent(in) :: stage_file
            ! Raster filename with the Okada stage perturbation
        real(dp), intent(in) :: global_dt

        call setup_elevation(domain)

        call setup_stage_and_forcing(domain, stage_file, global_dt)

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


