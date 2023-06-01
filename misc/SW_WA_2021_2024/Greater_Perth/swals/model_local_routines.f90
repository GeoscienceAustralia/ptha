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
    use file_io_mod, only: read_character_file, count_file_lines
        ! Used to read the elevation data filenames in preference order

    implicit none

    private

    public :: set_initial_conditions, parse_commandline_args, rise_time
        !! These are called by the main program

    character(len=charlen), allocatable :: input_elevation_files(:)
        !! Names of elevation files

    type(xyz_lines_type) :: breakwalls_forced, linear_patches
        !! Type with lon/lat/elev representing the tops of breakwalls. 
        !! These features are "burned" into the elevation so that on coarser 
        !! meshes, we don't accidently miss barriers to flow.

    real(dp), parameter :: rise_time = 0.0_dp !180.0_dp ! Fujii and Satake 2007 tried rise-time of 1,2,3 min.
        !! Okada deformation is applied instantaneously in this model, even though
        !! the code supports use of a non-zero rise-time. If using non-zero rise time,
        !! beware interpretation of max-stage -- in areas of subsidence, it could simply
        !! reflect the pre-subsidence topography. This could be worked-around by resetting
        !! domain%max_U at some time after the forcing has finished, but I haven't done that.
    real(dp), parameter :: offshore_manning = 0.03_dp 
        !! Manning friction for offshore domains (if they use manning 
        !! friction - not all solvers do). 
        !! Slightly smaller than in Davies et al. (2020) for consistency with nearshore sites.
    real(dp), parameter :: nearshore_manning = 0.03_dp
        !! Manning friction for nonlinear domains. Not changed in this study
    real(dp) :: ambient_sea_level
        !! Background sea level - taken as a commandline argument.

    contains 

    subroutine parse_commandline_args(stage_file, run_type, final_time, &
            model_name, coarse_nonlinear_domain_extents_file, &
            fine_nonlinear_domain_extents_file, &
            highres_nonlinear_domain_extents_file, load_balance_file, &
            output_basedir)
        !!
        !! Convenience routine to read the commandline arguments. 
        !!

        character(len=charlen), intent(inout) :: stage_file, run_type, &
            model_name, coarse_nonlinear_domain_extents_file, &
            fine_nonlinear_domain_extents_file, &
            highres_nonlinear_domain_extents_file, &
            load_balance_file, output_basedir
        real(dp), intent(inout) :: final_time
        
        character(len=charlen) :: ambient_sea_level_char
        integer(ip), parameter :: narg = 8

        if(command_argument_count() /= narg) then
            !! Error message
            write(log_output_unit, *) "Incorrect number of commandline arguments - exactly ", narg, " are required:"
            write(log_output_unit, *) "  stage_file (raster filename with stage perturbation)"
            write(log_output_unit, *) "  model_name (used within the output_folder name) "
            write(log_output_unit, *) "  run_type (either 'test' or 'test_load_balance' or 'full') "
            write(log_output_unit, *) "  coarse_nonlinear_domain_extents_file (file with extents of coarse nonlinear domains)" 
            write(log_output_unit, *) "  fine_nonlinear_domain_extents_file (file with extents of fine nonlinear domains)" 
            write(log_output_unit, *) "  highres_nonlinear_domain_extents_file (file with extents of highres nonlinear domains)" 
            write(log_output_unit, *) "  load_balance_file (file with load balancing metadata, or '' to use crude defaults)"
            write(log_output_unit, *) "  ambient_sea_level (background sea-level in m, e.g. 0.0 or 0.95)"
            write(log_output_unit, *) ""
            call generic_stop
        end if

        call get_command_argument(1, stage_file)
        write(log_output_unit, *) 'stage_file: ', trim(stage_file)

        call get_command_argument(2, model_name)
            ! Outputs will go in a directory like './OUTPUTS/model_name_EXTRA_INFO_HERE/....'
        write(log_output_unit, *) 'model_name: ', trim(model_name)

        call get_command_argument(3, run_type)
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
        
        call get_command_argument(4, coarse_nonlinear_domain_extents_file)
        write(log_output_unit, *) 'coarse_nonlinear_domain_extents_file: ', trim(coarse_nonlinear_domain_extents_file)
        call get_command_argument(5, fine_nonlinear_domain_extents_file)
        write(log_output_unit, *) 'fine_nonlinear_domain_extents_file: ', trim(fine_nonlinear_domain_extents_file)
        call get_command_argument(6, highres_nonlinear_domain_extents_file)
        write(log_output_unit, *) 'highres_nonlinear_domain_extents_file: ', trim(highres_nonlinear_domain_extents_file)

        call get_command_argument(7, load_balance_file)
        write(log_output_unit, *) 'load_balance_file: ', trim(load_balance_file)

        if(run_type == 'test' .or. run_type == 'test_load_balance') then
            ! evolve for a short time for testing purposes, or to create a run 
            ! that can be used to create a load_balance_file.
            final_time = 3600.0_dp  * 0.1_dp
        else if(run_type == 'full') then
            ! evolve for 24 hours 
            final_time = 3600.0_dp * 24.0_dp
            ! evolve for 60 hours 
            !final_time = 3600.0_dp * 60.0_dp
        end if
        write(log_output_unit, *) 'final_time: ', final_time

        call get_command_argument(8, ambient_sea_level_char)
        read(ambient_sea_level_char, *) ambient_sea_level
        write(log_output_unit, *) 'ambient_sea_level: ', ambient_sea_level

        output_basedir = './OUTPUTS/' // &
            trim(model_name) // '-' // &
            trim(run_type) // '-' // &
            'ambient_sea_level_' // trim(ambient_sea_level_char)
            ! Informative name for output folder
        write(log_output_unit, *) 'output_basedir: ', trim(output_basedir)


    end subroutine

    subroutine setup_breakwalls(domain)
        !! Burn breakwalls into the domain's elevation grid (i.e. 
        !! domain%U(:,:,ELV))
        type(domain_type), intent(inout) :: domain
            !! domain%U(:,:,ELV) is modified

        !! Filenames containing breakwall x/y/z data, used to setup
        !! 'breakwalls_forced' if it is not already setup
        character(len=charlen), allocatable :: breakwall_csv_lon_lat_z(:), &
            linear_patches_csv_lon_lat_z(:)
        integer(ip) :: file_unit
        character(len=charlen), parameter :: char_format = '(A)'

        if(.not. allocated(breakwalls_forced%lines)) then
            ! On the first call, setup the breakwalls_forced object.

            ! Read the breakwall files
            open(newunit=file_unit, file="../breakwalls/swals_breakwall_files.txt", &
                status="old", action="read")
            call read_character_file(file_unit, breakwall_csv_lon_lat_z, char_format)
            close(file_unit)
            ! Add in the relative directory path
            breakwall_csv_lon_lat_z = "../breakwalls/" // breakwall_csv_lon_lat_z

            ! Setup breakwalls interpolation object
            call breakwalls_forced%read_from_csv(breakwall_csv_lon_lat_z, &
                skip_header=1_ip)

        end if

        call breakwalls_forced%burn_into_grid(domain%U(:,:,ELV), &
            domain%lower_left, (domain%lower_left + domain%lw), &
            burn_type='max')

        !! Now repeat the above for 'linear_patches'. These are lines where
        !! we do not want burn_type='max' -- we just want to burn the value.

        !if(.not. allocated(linear_patches%lines)) then
        !    ! Lines with 'patches' 
        !    linear_patches_csv_lon_lat_z = [character(len=charlen) :: &
        !        '../breakwalls/gappatch/gappatch1.csv' &
        !        ]

        !    call linear_patches%read_from_csv(linear_patches_csv_lon_lat_z, &
        !        skip_header=1_ip)
        !end if

        !call linear_patches%burn_into_grid(domain%U(:,:,ELV), &
        !    domain%lower_left, (domain%lower_left + domain%lw), &
        !    burn_type='point_value')

    end subroutine

    subroutine setup_elevation(domain)
        !! Set the elevation of the domain from a set of rasters. 
        !! Burn breakwalls into the elevation. 
        !! Enforce NS walls, so that the domain is completely closed and we can 
        !! reason clearly about energy and mass conservation.
        type(domain_type), intent(inout) :: domain
            !! The domain for which domain%U(:,:,ELV) is set

        type(multi_raster_type):: elevation_data
            !! Objects to interpolate from elevation files
        real(dp), allocatable:: x(:), y(:)
            !! Coordinates for elevation lookup
        integer(ip) :: j, file_unit, num_smooth
        character(len=charlen), parameter :: char_format = '(A)'

        if(.not. allocated(input_elevation_files)) then
            ! Read the elevation data files in preference order
            open(newunit=file_unit, file="../elevation/swals_elevation_files_in_preference_order.txt", &
                status="old", action="read")
            call read_character_file(file_unit, input_elevation_files, char_format)
            close(file_unit)
        end if

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
                ! are denoted by large negative numbers (e.g. approx -3.0e+38)
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

        ! Try smoothing the elevation near a specific site (where a nesting artefact had occurred)
        !
        ! Because of the diffusive nature of local smoothing, we might smooth more times on fine 
        ! grids than on coarse grids. This would also prevent smoothing on the coarser domains
        ! num_smooth = nint(1.0_dp/(60.0_dp * 9.0_dp * 6.0_dp) * 1.0_dp/domain%dx(1))**2
        ! Site location and smoothing radius chosen from inspection of nesting artefact.
        call domain%smooth_elevation_near_point(number_of_9pt_smooths = 1_ip, & ! num_smooth,
            pt = [115.633333_dp, -32.097179_dp], smooth_region_radius_meters = 400.0_dp, &
            transition_region_radius_meters = 800.0_dp)


        write(log_output_unit, *) 'Elevation ranges before setting breakwals:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        call setup_breakwalls(domain)

        write(log_output_unit, *) 'Elevation ranges AFTER setting breakwals:'
        write(log_output_unit, *) minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
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
                    bilinear = 1_ip, na_below_limit=-1.0e+20_dp) 
                ! Clip NA regions
                where(forcing_context%forcing_work(:,j,STG) < (-1.0e+20_dp)) &
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

    subroutine setup_stage_and_forcing(domain, stage_file, global_dt)
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
        integer(ip):: j, n, fid
            ! Convenience integers
        logical :: stage_file_is_raster
        character(len=charlen), allocatable :: stage_files(:)
        real(dp), allocatable :: slips(:), start_times(:), end_times(:)
        integer, parameter :: csv_header_size = 1


        ! If stage_file ends in csv, assume it contains data for a multi-unit-source inversion
        ! Otherwise we assume it is a raster
        n = len_trim(stage_file)
        stage_file_is_raster = (stage_file( (n-3):n) /= '.csv')

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
                    bilinear=1_ip, na_below_limit=-1.0e+20_dp)

                ! Clip 'NA' regions (since the stage raster does not cover the 
                ! entire domain). In this case revert to the "0" value, later 
                ! we will add an ambient_sea_level offset
                where(domain%U(:,j,STG) < (-1.0e+20_dp) ) domain%U(:,j,STG) = 0.0_dp
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
            !     "water-surface-unit-source rasters", "slip", "forcing_start_time", "forcing_end_time"
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

        ! Adjust for ambient sea-level 
        domain%U(:,:,STG) = domain%U(:,:,STG) + ambient_sea_level
        domain%msl_linear = ambient_sea_level ! Influences linear solvers & potential-energy calculation.

        ! Alway need stage >= elevation.
        domain%U(:,:,STG) = max(domain%U(:,:,STG), &
                                domain%U(:,:,ELV) + (minimum_allowed_depth/100.0_dp) )

        write(log_output_unit, *) 'Stage range:'
        write(log_output_unit, *) minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        flush(log_output_unit)

    end subroutine

    subroutine set_initial_conditions(domain, stage_file, global_dt)
        !! Setup initial conditions and earthquake forcing

        class(domain_type), intent(inout):: domain
            ! We initialise domain%U and domain%manning_squared, and 
            ! setup the forcing term if rise_time > 0
        character(len=charlen), intent(in) :: stage_file
            ! Raster filename with the Okada stage perturbation
        real(dp), intent(in) :: global_dt
            ! The global time-step. This is used in-case we choose to use a rise-time >0.
            ! In that case, the rk2 timestepping will not include the full forcing if
            ! we begin the forcing at time=0. However if we begin the forcing after global_dt,
            ! it will include the full forcing.

        ! Basic setting of elevation -- modified again later by the earthquake forcing.
        call setup_elevation(domain)

        ! Set the earthquake initial condition [stage + elevation perturbation]
        call setup_stage_and_forcing(domain, stage_file, global_dt)

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


