module local_routines 
    !!
    !! Set initial conditions for periodic_multidomain
    !!
    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use logging_mod, only: log_output_unit
    use forcing_mod, only: forcing_patch_type, apply_forcing_patch
    use iso_c_binding, only: c_loc, c_f_pointer
    implicit none

    contains 

    subroutine set_initial_conditions(domain, stage_file, rise_time, global_dt)
        class(domain_type), intent(inout):: domain
        character(len=charlen), intent(in) :: stage_file
        real(dp), intent(in) :: rise_time
        real(dp), intent(in) :: global_dt

        integer(ip):: i, j
        character(len=charlen):: input_elevation(1), input_stage(1)
        type(multi_raster_type):: elevation_data, stage_data
        real(dp), allocatable:: x(:), y(:)
        type(forcing_patch_type), pointer :: forcing_context
            ! Use the forcing_context to optionally apply the stage deformation over time
        integer :: i0, i1, j0, j1

        ! Stage
        domain%U(:,:,STG) = 0.0e-0_dp

        ! Input rasters
        !input_stage(1) = "../generic_example/model_scenario_114_similar_to_tohoku_tsunami.tif"
        input_stage(1) = stage_file
        call stage_data%initialise(input_stage)

        input_elevation(1) = '../generic_example/merged_gebco_ga250_dem.tif'
        call elevation_data%initialise(input_elevation)

        ! Make space for x/y coordinates, at which we will look-up the rasters
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x
        
        ! Set stage and elevation row-by-row.
        ! This saves memory compared to doing it all at once.
        do j = 1, domain%nx(2)
            y = domain%y(j)

            ! Set elevation
            call elevation_data%get_xy(x,y, domain%U(:,j,ELV), domain%nx(1), bilinear=1_ip)

            ! Set stage
            call stage_data%get_xy(x,y, domain%U(:,j,STG), domain%nx(1), bilinear=1_ip)
            ! Clip 'NA' regions (since the stage raster does no cover the entire domain)
            do i = 1, domain%nx(1)
                if(domain%U(i,j,STG) < (-HUGE(1.0_dp)*0.99_dp) ) domain%U(i,j,STG) = 0.0_dp
            end do

            !do i = 1, domain%nx(1)
            !    if(domain%U(i,j,ELV) < (-HUGE(1.0_dp)*0.99_dp) ) print*, i,j,x(i), y(j)
            !end do

        end do
        call elevation_data%finalise()
        call stage_data%finalise()

        if(rise_time > 0.0_dp) then
            ! If a non-zero rise time was used, then create the forcing context to apply the stage perturbation over time.
            ! Find the region where the stage deformation applies
            i0 = count(domain%x < stage_data%lowerleft(1)) + 1
            i1 = count(domain%x <= stage_data%upperright(1))
            j0 = count(domain%y < stage_data%lowerleft(2)) + 1
            j1 = count(domain%y <= stage_data%upperright(2))


            if(i0 < i1 .and. j0 < j1) then
                ! The domain overlaps with the stage raster, so we need to apply the forcing.

                ! This is an object which contains the stage perturbations (and optionally ELV/UH/VH perturbations),
                ! and the start-time and end-time over which it is applied.
                allocate(forcing_context)

                ! Start forcing only after the model has evolved for one timestep. Otherwise the rk2 solver may
                ! not apply the full forcing [since in theory, some of the forcing should have been applied in evolving
                ! juts before time 0.0, because of how rk2 works].

                ! ! We can either force all of domain%U -- or if we provide k0,k1 we can only force some slice - which
                ! ! requires less memory in this case. 
                call forcing_context%setup(start_time=global_dt, end_time=rise_time+global_dt, &
                    i0=i0, i1=i1, j0=j0, j1=j1, k0=STG, k1=STG)
                forcing_context%forcing_work(i0:i1,j0:j1,STG) = domain%U(i0:i1, j0:j1, STG)
                ! ! Alternative where we can force everything
                !call forcing_context%setup(start_time=global_dt, end_time=rise_time+global_dt, &
                !    i0=i0, i1=i1, j0=j0, j1=j1)
                ! Set the forcing work to be equal to the stage deformation. 
                ! forcing_context%forcing_work(i0:i1,j0:j1,STG) = domain%U(i0:i1, j0:j1, STG)
                !forcing_context%forcing_work(i0:i1,j0:j1,ELV) = 0.0_dp !domain%U(i0:i1, j0:j1, STG)
                !forcing_context%forcing_work(:,:,UH:VH) = 0.0_dp

                ! Store the forcing context inside the domain
                domain%forcing_context_cptr = c_loc(forcing_context)
                forcing_context => NULL()
                ! Define the forcing subroutine
                domain%forcing_subroutine => apply_forcing_patch
                ! Re-set the stage deformation to zero, because now we will apply it as a forcing_subroutine
                domain%U(:,:,STG) = 0.0_dp
            end if

        end if


        ! Wall boundaries N/S
        domain%U(:,1:2,ELV) = 100.0_dp
        domain%U(:,(domain%nx(2)-1):domain%nx(2), ELV) = 100.0_dp

        deallocate(x,y)

        ! Wet-dry instabilities fix by flattening send regions with elevation > -10
        !call domain%use_constant_wetdry_send_elevation(elevation_threshold=-10.0_dp)

        if(domain%timestepping_method /= 'linear') then
            domain%manning_squared = 0.02_dp * 0.02_dp !0.0_dp !0.025_dp * 0.025_dp
        end if

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-07_dp)

        !write(log_output_unit,*) 'DEBUG: SETTING STAGE TO ZERO'
        !domain%U(:,:,STG) = max(0.0_dp, domain%U(:,:,ELV) + 1.0e-07_dp)

        write(log_output_unit,*) 'Stage range is: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))
        write(log_output_unit,*) 'Elev range is: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program periodic_multidomain
    !!
    !! Example model using a multidomain with periodic EW boundary conditions.
    !!

    use global_mod, only: ip, dp, minimum_allowed_depth, charlen
    use domain_mod, only: domain_type
    use multidomain_mod, only: multidomain_type, setup_multidomain, test_multidomain_mod
    use boundary_mod, only: flather_boundary, transmissive_boundary
    use local_routines
    use timer_mod
    use logging_mod, only: log_output_unit, send_log_output_to_file
    use stop_mod, only: generic_stop
    use coarray_intrinsic_alternatives, only: swals_mpi_init, swals_mpi_finalize
    use iso_c_binding, only: C_DOUBLE !, C_INT, C_LONG
    implicit none

    ! Type holding all domains 
    type(multidomain_type) :: md

    ! Local timing object
    type(timer_type) :: program_timer

    ! Change this to decrease the cell size by mesh_refine (i.e. for convergence testing)
    integer(ip), parameter :: mesh_refine = 1_ip ! 1_ip --> 4 arc minute

    ! The global (i.e. outer-domain) time-step in the multidomain 
    real(dp) ::  global_dt = 7.82_dp * (1.0_dp/mesh_refine) !5.8_dp * (1.0_dp/mesh_refine) ! 7.82_dp * (1.0_dp/mesh_refine)

    ! Approx timestep between outputs
    real(dp) :: approximate_writeout_frequency = 3600.0_dp
    real(dp) :: final_time = 3600.0_dp * 24.0_dp

    ! Length/width
    real(dp), parameter, dimension(2):: global_lw = [360.0_dp, 135.0_dp]
    ! Lower-left corner coordinate
    !real(dp), parameter, dimension(2):: global_ll = [-180.0_dp, -70.0_dp]
    real(dp), parameter, dimension(2):: global_ll = [-40.0_dp, -70.0_dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = nint(global_lw*15*mesh_refine)

    ! Inner domain 
    integer(ip), parameter :: nest_ratio = 5_ip

    ! Useful misc variables
    integer(ip):: j, i, nd
    real(dp) :: rise_time
    character(len=charlen) :: stage_file, model_name, rise_time_char

    call swals_mpi_init

    ! Assume the stage file was passed to the command line
    call get_command_argument(1, stage_file)
    ! Name for the model output folder, passed to the command line
    call get_command_argument(2, model_name)
    ! The "rise time" for the earthquake in seconds, optionally passed to the commandline
    call get_command_argument(3, rise_time_char)
    if(rise_time_char == '') then
        rise_time = 0.0_dp
    else
        read(rise_time_char, *) rise_time
    end if


    ! Set the model name
    md%output_basedir = './OUTPUTS/' // trim(model_name)

    call program_timer%timer_start('setup')

#ifndef SPHERICAL
    write(log_output_unit,*) 'Code assumes spherical coordinates, but SPHERICAL is not defined'
    call generic_stop
#endif

    ! Set periodic EW boundary condition
    !md%periodic_xs = [-180.0_dp, 180.0_dp]
    md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
    
    
    ! nd domains in this model
    nd = 1  ! 2
    allocate(md%domains(nd))

    !
    ! Setup basic metadata
    !

    ! Linear domain
    md%domains(1)%lw = global_lw
    md%domains(1)%lower_left =global_ll
    md%domains(1)%nx = global_nx
    md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
    md%domains(1)%dx_refinement_factor = 1.0_dp
    md%domains(1)%timestepping_refinement_factor = 1_ip
    md%domains(1)%timestepping_method = 'linear'

    !@ Higher res around region of interest
    !call md%domains(2)%match_geometry_to_parent(&
    !    parent_domain=md%domains(1), &
    !    lower_left  = [139.2_dp, 41.5_dp], &
    !    upper_right = [140.5_dp, 43.5_dp], &
    !    dx_refinement_factor = nest_ratio, &
    !    timestepping_refinement_factor = 1_ip,
    !    rounding_method='nearest')
    !md%domains(2)%timestepping_method = 'midpoint'
    !md%domains(2)%support_elevation_forcing = .true.

    ! Allocate domains and prepare comms
    call md%setup()
    call md%memory_summary()

    ! Give boundary condition to domans that have non-nesting boundaries
    !do j = 1, size(md%domains)
    !    if(any(.not.md%domains(j)%is_nesting_boundary)) md%domains(j)%boundary_subroutine => flather_boundary
    !end do

    ! Set initial conditions
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j), stage_file, rise_time, global_dt)
    end do
    call md%make_initial_conditions_consistent()
    
    ! NOTE: For stability in 'null' regions, we set them to 'high land' that
    ! should be inactive. 
    call md%set_null_regions_to_dry()
   
    write(log_output_unit,*) 'End setup'

    ! Print the gravity-wave CFL limit, to guide timestepping
    do j = 1, size(md%domains)
        write(log_output_unit,*) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
    end do

    call program_timer%timer_end('setup')
    call program_timer%timer_start('evolve')

    flush(log_output_unit)

    !
    ! Evolve the code
    !
    do while (.true.)
        
        ! IO 
        call md%write_outputs_and_print_statistics(&
            approximate_writeout_frequency=approximate_writeout_frequency,&
            timing_tol = 1.0e-06_dp)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

    end do

    call program_timer%timer_end('evolve')
    call md%finalise_and_print_timers

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)

    call swals_mpi_finalize
end program
