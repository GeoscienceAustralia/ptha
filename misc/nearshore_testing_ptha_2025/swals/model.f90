!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module "local_routines" with various helper subroutines.
#include "model_local_routines.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_model
    !!
    !! A global-to-local model with various high-res sites around Australia.
    !! 
    !! Basic usage involves these commandline arguments (designed for a particular study):
    !!     ./model stage_file rise_time model_name run_type load_balance_file offshore_solver_type offshore_manning highres_regions
    !!
    !! where the arguments are:
    !!     stage_file (raster filename with stage perturbation)
    !!     rise_time (time in seconds over which stage perturbation is applied. Can be 0.0 for instantaneous)
    !!     model_name (this will be included within the output_folder name -- use it to help identify outputs)
    !!     run_type (either 'test' or 'test_load_balance' or 'full'):
    !!         'full' is used to run a proper model; 
    !!         the other cases run short models using the load_balance_file ('test_load_balance') or not using any
    !!         load balancing ('test'). The latter case is useful to produce outputs required to make a load_balance_file.
    !!     load_balance_file (file with load balancing metadata. Make it empty '' to use crude defaults)"
    !!     offshore_solver_type ( 'linear_with_manning' or 'leapfrog_nonlinear' or 'linear_with_linear_friction' or
    !!         'linear_with_reduced_linear_friction' or 'linear_with_delayed_linear_friction' or 'linear_with_no_friction') ). 
    !!          Used to test different deep-ocean propagation approaches. 
    !!     offshore_manning (manning coefficient for offshore solver if using 'linear_with_manning' or 'leapfrog_nonlinear')"
    !!     highres_regions (either 'none' [global domain only] or 'australia' [first model using all highres areas, uses "perth" in WA] or 
    !!         'NSW' [only use NSW high-res domains]) or 'perth' [only use perth high-res domains] or 'SWWA' [only use south-west western Australia highres domains]
    !!         or 'australiaSWWA' [like 'australia' but using the more extensive domains in South west western Australia]
    !!         or 'NWWA' or 'australiaWA' or 'WA'
    !!
    !! The geometry of the domains is specified in this program.
    !!

    !
    ! Imports from SWALS 
    !
    use global_mod, only: ip, dp, charlen
        ! Integer/real precision and default character length
    use domain_mod, only: STG, UH, VH, ELV
        ! Indices of stage, depth-integrated velocities 
        ! (UH = east, VH = north) and elevation in the domain%U array.
    use multidomain_mod, only: multidomain_type
        ! The main type that holds all domains and evolves them
    use coarray_intrinsic_alternatives, only: swals_mpi_init, swals_mpi_finalize 
        ! Call at the beginning/end to ensure mpi starts and finishes. 
        ! Does nothing if not compiled for distributed parallel
    use timer_mod, only: timer_type
        ! For local timing of the code
    use logging_mod, only: log_output_unit
        ! Write messages to log_output_unit
    use stop_mod, only: generic_stop
        ! For halting when an error occurs (works in serial or parallel, tries 
        ! to close files, etc).

    !
    ! Case specific helper routines
    !
    use local_routines, only: set_initial_conditions, parse_commandline_args, highres_regions, included_regions

    implicit none

    type(multidomain_type) :: md
        ! md is the main object -- holds all domains, evolves them, etc.

    type(timer_type) :: program_timer
        ! Local code-timing object

    real(dp), parameter:: global_lw(2) = [360.0_dp, 147.0_dp] 
        ! Length/width of multidomain in degrees lon,lat
    real(dp), parameter:: global_ll(2) = [-40.0_dp, -79.0_dp]
        ! Lower-left corner coordinate of multidomain in degrees lon,lat

    integer(ip), parameter :: mesh_refine = 4_ip 
        ! Increase this to decrease the cell side-length by mesh_refine 
        ! (i.e. for convergence testing). 4_ip --> 1 arcmin in the coarse 
        ! global domain

    real(dp), parameter ::  global_dt = 6.0_dp * (1.0_dp/mesh_refine) 
        ! The global time-step in the multidomain. Must satisfy CFL condition 
        ! everywhere (in combination with local timestepping that is specified
        ! when defining the domains below)
    real(dp), parameter :: approximate_writeout_timestep = 30.0_dp
        ! Approx timestep between any outputs (in this case, tide-gauge outputs)
    integer(ip), parameter :: write_grids_every_nth_step = 9999999_ip !3600_ip
        ! Optionally write grids less often than the writeout timestep, to keep file-size down
    integer(ip), parameter :: print_every_nth_step = 10_ip
        ! Optionally print outputs less often than the writeout timestep, to keep the file-size down
    integer(ip), parameter :: write_gauges_every_nth_step = 1_ip
        ! Optionally write gauges less often than the writeout timestep, to keep the file-size down

    real(dp), parameter :: seconds_before_evolve = 0.0_dp 
        ! Some domains might not need to evolve at the start -- use this as
        ! the time before which they evolve. Crude approach (could be domain
        ! specific).

    integer(ip), parameter :: nd_global = 4, &
                              nd_nsw = 8, & ! Lots of domains in NSW
                              nd_victoria = 2, &
                              nd_perth = 2, & ! 
                              nd_SWWA = 5, & ! For models with SWWA, this is in addition to the perth domains
                              nd_geraldton = 2, & 
                              nd_NWWA = 5 

    integer(ip), parameter :: global_eastcoast = 3, global_westcoast = 2
        ! The global domain is split into 'nd_global' pieces initially (not considering parallel partitioning).
        ! Define the index containing the eastcoast and westcoast nested grids. Actually, these could be any integer from 1-nd_global,
        ! because SWALS will figure out the appropriate nesting domain by itself so long as we provide one of the global
        ! domains to inform the initial geometry alignment.

    real(dp) :: final_time
        ! Duration of simulation in seconds (start_time = 0.0_dp)

    integer(ip):: j, nd, nsw_regional, vic_regional, perth_regional, nwwa_regional, &
        geraldton_regional, global_main, last_di
        ! Useful misc local variables
    character(len=charlen) :: stage_file, model_name, run_type, &
        offshore_solver_type
        ! For reading commandline
    real(dp) :: linear_friction_delay_time, linear_friction_delay_value
        ! For case with linear friction after some time

    call swals_mpi_init 
        ! Ensure MPI is initialised

#ifndef SPHERICAL
    write(log_output_unit,*) &
        'Code assumes spherical coordinates, but SPHERICAL is not defined'
    call generic_stop
#endif

    ! This also sets the variable highres_regions 
    call parse_commandline_args(stage_file, run_type, final_time, model_name, &
        md%load_balance_file, offshore_solver_type, md%output_basedir)

    !
    ! Basic definition of multidomain
    !

    call program_timer%timer_start('startup_define_multidomain_geometry')

    !
    ! Figure out how many domains are needed
    !
    if(highres_regions == 'none') then
        ! Do not use any of the regional or high-res domains
        nd = nd_global
        included_regions = [character(len=charlen):: "global"]
    else if(highres_regions == 'australia') then
        ! High res domains everywhere, but with limited coverage in south-west WA
        nd = nd_global + nd_nsw + nd_victoria + nd_perth
        included_regions = [character(len=charlen):: "global", "nsw", "victoria", "perth"]
    else if(highres_regions == 'NSW') then
        ! High res domains in NSW only - plus include regional Victoria domain
        ! so that edge-waves can propagate well around that coast
        nd = nd_global + nd_nsw + 1
        included_regions = [character(len=charlen):: "global", "nsw"]
    else if(highres_regions == 'perth') then
        ! High res domains in Perth only
        nd = nd_global + nd_perth
        included_regions = [character(len=charlen):: "global", "perth"]
    else if(highres_regions == 'SWWA') then
        ! Highres domains in south west WA, more than just Perth region
        nd = nd_global + nd_perth + nd_SWWA + nd_geraldton
        included_regions = [character(len=charlen):: "global", "perth", "swwa", "geraldton"]
    else if(highres_regions == 'NWWA') then
        ! Global domains plus domains in NWWA only
        nd = nd_global + nd_NWWA + nd_geraldton
        included_regions = [character(len=charlen):: "global", "nwwa", "geraldton"]
    else if(highres_regions == 'australiaSWWA') then
        ! High res domains everywhere, including the more extensive southwest WA treatment.
        nd = nd_global + nd_nsw + nd_victoria + nd_perth + nd_SWWA + nd_geraldton
        included_regions = [character(len=charlen):: "global", "nsw", "victoria", "perth", "swwa", "geraldton"]
    else if(highres_regions == 'WA') then
        ! Global domains plus everything in WA
        nd = nd_global + nd_perth + nd_SWWA + nd_geraldton + nd_NWWA
        included_regions = [character(len=charlen):: "global", "perth", "swwa", "geraldton", "nwwa"]
    else if(highres_regions == 'australiaWA') then
        ! High res domains everywhere, including the more extensive WA treatment.
        nd = nd_global + nd_nsw + nd_victoria + nd_perth + nd_SWWA + nd_geraldton + nd_NWWA
        included_regions = [character(len=charlen):: "global", "nsw", "victoria", "perth", "swwa", "geraldton", "nwwa"]
    else
        write(log_output_unit, *) "Unknown highres_regions: ", trim(highres_regions)
        call generic_stop
    end if

    ! nd domains in this model
    allocate(md%domains(nd))

    ! Sanity check: Always including global domains
    if(.not. any(included_regions == "global")) then
        write(log_output_unit, *) included_regions
        write(log_output_unit, *) "error: included_regions must include 'global'"
        call generic_stop
    end if

    !
    ! Enforce periodic EW boundary condition as the values equal
    ! the x-range of the multidomain
    !
    md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
    if(abs(global_lw(1) - 360.0_dp) > 1.0e-06_dp) then
        write(log_output_unit, *) "Model assumes periodic boundary conditions, must have global_lw(1) = 360.0_dp, but it is", &
            global_lw(1)
    end if

    last_di = 0 ! Track the last domain index

    !
    ! Setup domain metadata
    !
    do j = 1, nd_global
        ! Global linear domain, split into (nd_global =4) pieces by longitude
        ! This is used for all model types

        md%domains(j)%lw = [global_lw(1)*1.0_dp/nd_global, global_lw(2)]
        md%domains(j)%lower_left = &
            [global_ll(1) + (j-1)*global_lw(1)*1.0_dp/nd_global, global_ll(2)]
        md%domains(j)%dx = 1/60.0_dp * 4.0_dp * [1.0_dp, 1.0_dp] / mesh_refine 
        md%domains(j)%nx = nint(md%domains(j)%lw/md%domains(j)%dx)
        md%domains(j)%dx_refinement_factor = 1.0_dp
        md%domains(j)%timestepping_refinement_factor = 1_ip
        md%domains(j)%nc_grid_output%spatial_stride = 4_ip 
            ! Reduce output file size by only saving every n'th cell 
        last_di = last_di + 1

        ! A few options for the type of offshore solver
        select case(offshore_solver_type)
        case ("linear_with_manning")
            md%domains(j)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
            md%domains(j)%linear_solver_is_truely_linear = .true.
            ! Try Chezy friction in the offshore domains only -- more energy loss. 
            ! md%domains(j)%friction_type = 'chezy'
        case ("linear_with_linear_friction")
            md%domains(j)%timestepping_method = 'linear' 
            md%domains(j)%linear_friction_coeff = 1.0e-05_dp
            md%domains(j)%linear_solver_is_truely_linear = .true.
        case ("linear_with_reduced_linear_friction")
            md%domains(j)%timestepping_method = 'linear' 
            md%domains(j)%linear_friction_coeff = 1.0_dp/(36.0_dp*3600.0_dp)
            md%domains(j)%linear_solver_is_truely_linear = .true.
        case ("linear_with_delayed_linear_friction")
            md%domains(j)%timestepping_method = 'linear' 
            md%domains(j)%linear_friction_coeff = 0.0_dp ! Later change to 1e-05_dp
            md%domains(j)%linear_solver_is_truely_linear = .true.
            ! Apply linear friction after 12 hours
            linear_friction_delay_time = 12.0_dp * 3600.0_dp 
            linear_friction_delay_value = 1.0e-05_dp
        case ("linear_with_no_friction")
            md%domains(j)%timestepping_method = 'linear' 
            md%domains(j)%linear_friction_coeff = 0.0e-05_dp
            md%domains(j)%linear_solver_is_truely_linear = .true.
        case ("leapfrog_nonlinear")
            md%domains(j)%timestepping_method = 'leapfrog_nonlinear'
        case default
            write(log_output_unit,*) 'Unrecognized offshore_solver_type'
            call generic_stop
        end select

    end do

    ! Below here, the domains are nonlinear. In each region we've got a
    ! regional domain containing one or more high-res nested domains.
    call setup_nested_domains

    ! Sanity check on the number of domains
    if(last_di /= nd) then
        write(log_output_unit, *) included_regions
        write(log_output_unit, *) last_di, nd
        write(log_output_unit, *) "error: last_di /= nd, problem in counting domains"
        call generic_stop
    end if

    ! Minor adjustments to domains
    do j = 1, size(md%domains)

        md%domains(j)%nc_grid_output%time_var_store_flag(STG:VH) = .true.
        md%domains(j)%nc_grid_output%time_var_store_flag(ELV) = .false.
            ! Store stage/UH/VH grids over time, but not ELEVATION (it will 
            ! be stored once anyway).

        if(any(highres_regions == [character(len=charlen) :: 'australiaSWWA', 'SWWA', 'NWWA', 'WA', 'australiaWA'])) then
            ! Reduce storage for the newer models. Not done for the older models to
            ! be backward compatable. With hindsight I should have done this in every case
            md%domains(j)%nc_grid_output%time_var_store_flag = .false.
        end if

        md%domains(j)%theta = 4.0_dp
            ! Use "non-TVD" limiting in nonlinear domains. Less dissipative.
    end do

    call program_timer%timer_end('startup_define_multidomain_geometry')

    !
    ! Setup the multidomain object
    !
    call program_timer%timer_start('startup_md_setup')

    call md%setup()
        ! Allocate domains and prepare comms
    call program_timer%timer_end('startup_md_setup')

    !
    ! Read initial conditions and make them consistent with each other
    !
    call program_timer%timer_start('startup_set_initial_conditions')
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j), stage_file, global_dt)
    end do
    call md%memory_summary()
    call md%make_initial_conditions_consistent()
        ! Perform a parallel halo exchange so initial conditions are consistent.

    call md%set_null_regions_to_dry()
        ! Set 'null' regions (i.e. non-halo areas where other domains have 
        ! priority) to 'high/dry land' that will be inactive. This enhances 
        ! stability. Such regions cannot interact with priority regions
        ! (because the halo update prevents it).

    call program_timer%timer_end('startup_set_initial_conditions')
  
    !
    ! Gauges 
    ! 
    call program_timer%timer_start('startup_set_gauges')

    if(highres_regions == 'SWWA' .or. highres_regions == 'australiaSWWA') then
        ! This includes more gauges in SWWA region -- only for new models to ensure backward compatibility
        call md%set_point_gauges_from_csv("point_gauges_combined_SWWA.csv", &
            skip_header=1_ip)
    else if(highres_regions == 'NWWA' .or. highres_regions == 'WA' .or. highres_regions == 'australiaWA') then
        ! This includes more gauges in NWWA region -- only for new models to ensure backward compatibility
        call md%set_point_gauges_from_csv("point_gauges_combined_NWWA.csv", &
            skip_header=1_ip)
    else
       ! Fewer gauges in WA region
       call md%set_point_gauges_from_csv("point_gauges_combined.csv", &
           skip_header=1_ip)
    end if
    call program_timer%timer_end('startup_set_gauges')
  
    !
    ! Final setup work.
    ! 
    call program_timer%timer_start('startup_end')

    call md%record_initial_volume()
        ! For mass conservation checks

    do j = 1, size(md%domains)
        write(log_output_unit,*) 'domain: ', j, 'ts: ', &
            md%domains(j)%stationary_timestep_max()
            ! Print the gravity-wave CFL limit to guide timestepping
        call md%domains(j)%timer%reset
            ! Reset the domain timers, so that load-balancing only sees the 
            ! evolve info
    end do

    write(log_output_unit,*) 'End setup'
    call program_timer%timer_end('startup_end')
    
    flush(log_output_unit)

    !
    ! Main evolve loop
    !

    do while (.true.)

        call program_timer%timer_start('IO')
        call md%write_outputs_and_print_statistics(&
            ! Print and write outputs
            approximate_writeout_frequency=approximate_writeout_timestep, &
                ! Time between writes is ~= "approximate_writeout_timestep"
            write_grids_less_often = write_grids_every_nth_step, &
                ! Write gridded outputs less often
            write_gauges_less_often = write_gauges_every_nth_step, &
                ! Write gauges every time 
            print_less_often = print_every_nth_step, & !1_ip,&
                ! Print domain statistics less often 
            timing_tol = (global_dt/2.01_dp))
        call program_timer%timer_end('IO')

        if (md%domains(1)%time > final_time) exit

        if(offshore_solver_type == 'linear_with_delayed_linear_friction') then
            ! Set the linear friction coefficient if the delay time has been exceeded
            if(md%domains(1)%time > linear_friction_delay_time) then
                do j = 1, size(md%domains)
                    if(md%domains(j)%timestepping_method == 'linear') then
                        md%domains(j)%linear_friction_coeff = linear_friction_delay_value
                    end if
                end do
            end if
        end if

        call program_timer%timer_start('evolve')
        call md%evolve_one_step(global_dt)
            ! Evolve the model by global_dt seconds
        call program_timer%timer_end('evolve')

    end do

    call md%finalise_and_print_timers
        ! Close files and print timing info (required for load balancing)

    write(log_output_unit,*) ''
    write(log_output_unit, *) 'Program timer'
    write(log_output_unit, *) ''
    call program_timer%print(log_output_unit)
    flush(log_output_unit)

    call swals_mpi_finalize

contains

subroutine setup_nested_domains()
    !
    ! Setup the nested grid geometry
    ! This became complicated so moving here. For future models it is probably better
    ! to use interfaces like for the NSW/WA/Gladstone projects where we provided the geometry in a text file.
    ! 

    if(any(included_regions == "nsw")) then
        !
        ! NSW high-res domains (+ regional victoria, as it might affect NSW in
        ! the east) 
        !
        
        nsw_regional = last_di + 1 
        call md%domains(nsw_regional)%match_geometry_to_parent(&
            ! Better-than-global resolution domain for NSW
            parent_domain=md%domains(global_eastcoast), &
            lower_left  = [149.25_dp, -38.5_dp], &
            upper_right = [152.5_dp, -32.0_dp], &
            dx_refinement_factor = 7_ip, &
            timestepping_refinement_factor = 4_ip, &!3_ip ,& 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional)%timestepping_method = 'rk2' 
        md%domains(nsw_regional)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional)%static_before_time = seconds_before_evolve
        last_di = last_di + 1
        
        call md%domains(nsw_regional+1)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Hawkesbury
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.85_dp, -33.75_dp], &
            upper_right = [151.60_dp, -33.40_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+1)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+1)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+1)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+2)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Sydney
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.85_dp, -34.2_dp], &
            upper_right = [151.40_dp, -33.75_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 5_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+2)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+2)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+2)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+3)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Port Kembla
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.84_dp, -34.6_dp], &
            upper_right = [151.10_dp, -34.35_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 5_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+3)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+3)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+3)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+4)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Jervis Bay
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.62_dp, -35.22_dp], &
            upper_right = [150.90_dp, -34.86_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+4)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+4)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+4)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+5)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Ulladullah
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.44_dp, -35.40_dp], &
            upper_right = [150.60_dp, -35.26_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+5)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+5)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+5)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+6)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Batemans Bay
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [150.08_dp, -35.8_dp], &
            upper_right = [150.35_dp, -35.6_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+6)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+6)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+6)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(nsw_regional+7)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Eden
            parent_domain=md%domains(nsw_regional), &
            lower_left  = [149.85_dp, -37.13_dp], &
            upper_right = [150._dp, -37.0_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 3_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nsw_regional+7)%timestepping_method = 'rk2' 
        md%domains(nsw_regional+7)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(nsw_regional+7)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        !
        ! Regional Victoria domain
        !
        vic_regional = last_di + 1
        call md%domains(vic_regional)%match_geometry_to_parent(&
            ! Better-than-global domain for Victoria
            parent_domain=md%domains(global_eastcoast), &
            lower_left = [141.0_dp, -41.25_dp], &
            upper_right = [149.25_dp, -37.3_dp], &
            dx_refinement_factor = 7_ip, &
            timestepping_refinement_factor = 4_ip, &
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(vic_regional)%timestepping_method = 'rk2'
        md%domains(vic_regional)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(vic_regional)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

    end if

    if(any(included_regions == 'victoria')) then

        ! Sanity check
        if(.not. any(included_regions == 'nsw')) then
            write(log_output_unit, *) included_regions
            write(log_output_unit, *) "error: included_regions must include 'nsw' if it includes 'victoria'"
            call generic_stop
        end if

        !! Include highres Vic if we are doing an 'australia-wide'
        !! model. For many events we only have tidal-gauge data in NSW, and
        !! in those cases this is a waste of compute
        call md%domains(vic_regional+1)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain Portland 
            parent_domain=md%domains(vic_regional), &
            lower_left = [141.45_dp, -38.5_dp], &
            upper_right = [141.75_dp, -38.22_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(vic_regional+1)%timestepping_method = 'rk2'
        md%domains(vic_regional+1)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(vic_regional+1)%static_before_time = seconds_before_evolve
        last_di = last_di + 1
    end if

    if(any(included_regions == "perth")) then
        perth_regional = last_di + 1

        !
        ! Regional Perth domain
        !
        call md%domains(perth_regional)%match_geometry_to_parent(&
            ! Better-than-global res domain around SW WA.
            parent_domain=md%domains(global_westcoast), &
            lower_left = [114.0_dp, -36.0_dp], &
            upper_right = [116.3_dp, -30.0_dp], &
            dx_refinement_factor = 7_ip, &
            timestepping_refinement_factor = 5_ip, &
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional)%timestepping_method = 'rk2'
        md%domains(perth_regional)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(perth_regional+1)%match_geometry_to_parent(&
            ! Higher-res nonlinear domain around Hillarys/Perth
            parent_domain=md%domains(perth_regional), &
            lower_left = [115.6_dp, -32.35_dp], &
            upper_right = [115.9_dp, -31.4_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+1)%timestepping_method = 'rk2'
        md%domains(perth_regional+1)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+1)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

    end if

    if(any(included_regions == "swwa")) then
        !
        ! More extensive SWWA domains
        !

        ! Sanity check
        if(.not. any(included_regions == "perth")) then
            write(log_output_unit, *) included_regions
            write(log_output_unit, *) "error: included_regions must include 'perth' if it includes 'swwa'"
            call generic_stop
        end if

        call md%domains(perth_regional+2)%match_geometry_to_parent(&
            ! Highres-Mandurah
            parent_domain=md%domains(perth_regional), &
            lower_left = [115.45_dp, -32.70_dp], &
            upper_right = [115.8_dp, -32.35_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+2)%timestepping_method = 'rk2'
        md%domains(perth_regional+2)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+2)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(perth_regional+3)%match_geometry_to_parent(&
            ! Highres-Bunbury
            parent_domain=md%domains(perth_regional), &
            lower_left = [115.55_dp, -33.35_dp], &
            upper_right = [115.7_dp, -33.20_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+3)%timestepping_method = 'rk2'
        md%domains(perth_regional+3)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+3)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(perth_regional+4)%match_geometry_to_parent(&
            ! Highres-PortGeographe, more refinement
            parent_domain=md%domains(perth_regional), &
            lower_left = [115.38_dp, -33.64_dp], &
            upper_right = [115.4_dp, -33.61_dp], &
            dx_refinement_factor = 11_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+4)%timestepping_method = 'rk2'
        md%domains(perth_regional+4)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+4)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(perth_regional+5)%match_geometry_to_parent(&
            ! Highres Lancelin
            parent_domain=md%domains(perth_regional), &
            lower_left = [115.31_dp, -31.04_dp], &
            upper_right = [115.34_dp, -30.99_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+5)%timestepping_method = 'rk2'
        md%domains(perth_regional+5)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+5)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(perth_regional+6)%match_geometry_to_parent(&
            ! Highres Jurian Bay
            parent_domain=md%domains(perth_regional), &
            lower_left =  [114.95_dp, -30.36_dp], &
            upper_right = [115.05_dp, -30.21_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(perth_regional+6)%timestepping_method = 'rk2'
        md%domains(perth_regional+6)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(perth_regional+6)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

    end if

    !if(any(highres_regions == [character(len=charlen) :: 'australiaSWWA', 'SWWA'])) then
    if(any(included_regions == "geraldton")) then
        !
        ! Include Geraldton
        ! 

        geraldton_regional = last_di + 1
        call md%domains(geraldton_regional)%match_geometry_to_parent(&
            ! Regional Geraldton
            parent_domain=md%domains(global_westcoast), &
            lower_left = [112.85_dp, -30.0_dp], &
            upper_right = [115.15_dp, -27.0_dp], &
            dx_refinement_factor = 7_ip, &
            timestepping_refinement_factor = 5_ip, &
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(geraldton_regional)%timestepping_method = 'rk2'
        md%domains(geraldton_regional)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(geraldton_regional)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        call md%domains(geraldton_regional+1)%match_geometry_to_parent(&
            ! Highres Geraldton
            parent_domain=md%domains(geraldton_regional), &
            lower_left =  [114.56_dp, -28.81_dp], &
            upper_right = [114.64_dp, -28.58_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(geraldton_regional+1)%timestepping_method = 'rk2'
        md%domains(geraldton_regional+1)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(geraldton_regional+1)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

    end if

    if(any(included_regions == "nwwa")) then
        ! This model only uses the global domains, plus the domains here
        nwwa_regional = last_di + 1

        ! Regional NWWA domain 1
        call md%domains(nwwa_regional)%match_geometry_to_parent(&
            parent_domain=md%domains(global_westcoast), &
            lower_left =  [113.2_dp, -22.7_dp], &
            upper_right = [115.8_dp, -18.5_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nwwa_regional)%timestepping_method = 'rk2'
        md%domains(nwwa_regional)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(nwwa_regional)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        ! Regional NWWA domain 2
        call md%domains(nwwa_regional+1)%match_geometry_to_parent(&
            parent_domain=md%domains(global_westcoast), &
            lower_left =  [115.8_dp, -21.5_dp], &
            upper_right = [122.8_dp, -16.5_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 3_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nwwa_regional+1)%timestepping_method = 'rk2'
        md%domains(nwwa_regional+1)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(nwwa_regional+1)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        ! Nested domain near Exmouth, more refinement
        call md%domains(nwwa_regional+2)%match_geometry_to_parent(&
            parent_domain=md%domains(nwwa_regional), &
            lower_left =  [114.12_dp, -22.0_dp], &
            upper_right = [114.20_dp, -21.9_dp], &
            dx_refinement_factor = 11_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nwwa_regional+2)%timestepping_method = 'rk2'
        md%domains(nwwa_regional+2)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(nwwa_regional+2)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        ! Nested domain near Onslow, more refinement
        call md%domains(nwwa_regional+3)%match_geometry_to_parent(&
            parent_domain=md%domains(nwwa_regional), &
            lower_left =  [115.08_dp, -21.67_dp], &
            upper_right = [115.17_dp, -21.60_dp], &
            dx_refinement_factor = 11_ip, & 
            timestepping_refinement_factor = 4_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nwwa_regional+3)%timestepping_method = 'rk2'
        md%domains(nwwa_regional+3)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(nwwa_regional+3)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

        ! Nested domain near Dampier, regular refinement
        call md%domains(nwwa_regional+4)%match_geometry_to_parent(&
            parent_domain=md%domains(nwwa_regional+1), &
            lower_left =  [116.3_dp, -20.92_dp], &
            upper_right = [117.3_dp, -20.33_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 2_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(nwwa_regional+4)%timestepping_method = 'rk2'
        md%domains(nwwa_regional+4)%nc_grid_output%spatial_stride = 1 
        if(run_type == 'full') &
            md%domains(nwwa_regional+4)%static_before_time = seconds_before_evolve
        last_di = last_di + 1

    end if

end subroutine
end program
