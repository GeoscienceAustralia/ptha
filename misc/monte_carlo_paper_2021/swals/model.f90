!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module "local_routines" with various helper subroutines.
#include "model_local_routines.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_model
    !!
    !! A global-to-local model with various high-res sites around Tonga
    !!
    !! Basic usage involves these commandline arguments (designed for a particular study):
    !!     ./model stage_file rise_time model_name run_type load_balance_file offshore_solver_type offshore_manning highres_regions outer_domain_extent output_style
    !!
    !! where the arguments are:
    !!     stage_file (raster filename with stage perturbation)
    !!     rise_time (time in seconds over which stage perturbation is applied. Can be 0.0 for instantaneous)
    !!     ambient_sea_level (e.g. 0.0 for MSL)
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
    !!     highres_regions 'none' or 'tonga'
    !!     outer_domain_extent 'global' or 'regional' or 'pacific'
    !!     output_style ('animation' or 'few_grids'). If 'animation' we write grids quite often, otherwise rarely.
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
    use boundary_mod, only: flather_boundary
        ! flather_boundary can be used as a passive boundary condition that tries to maintain a given stage, while letting wave
        ! pass.
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
    use local_routines, only: set_initial_conditions, parse_commandline_args, rise_time

    implicit none

    type(multidomain_type) :: md
        ! md is the main object -- holds all domains, evolves them, etc.

    type(timer_type) :: program_timer
        ! Local code-timing object

    real(dp):: global_lw(2) 
        ! Length/width of multidomain in degrees lon,lat
    real(dp):: global_ll(2)
        ! Lower-left corner coordinate of multidomain in degrees lon,lat

    integer(ip), parameter :: mesh_refine = 4_ip 
        ! Increase this to decrease the cell side-length by mesh_refine 
        ! (i.e. for convergence testing). 4_ip --> 1 arcmin in the coarse 
        ! global domain

    real(dp), parameter ::  global_dt = 1.5_dp * (1.0_dp/mesh_refine) 
        ! The global time-step in the multidomain. Must satisfy CFL condition 
        ! everywhere (in combination with local timestepping that is specified
        ! when defining the domains below)
    real(dp), parameter :: approximate_writeout_timestep = 30.0_dp
        ! Approx timestep between any outputs (in this case, tide-gauge outputs)
    integer(ip) :: write_grids_every_nth_step 
        ! Optionally write grids only every nth approximate_writeout_timestep, to keep file-size down
    integer(ip) :: print_every_nth_step = 10_ip
        ! Optionally print outputs only every nth approximate_writeout_timestep, to keep the file-size down
    integer(ip) :: write_gauges_every_nth_step = 1_ip
        ! Optionally write gauges only every nth approximate_writeout_timestep, to keep the file-size down

    real(dp) :: seconds_before_evolve
        ! Some domains might not need to evolve at the start -- use this as
        ! the time before which they evolve. Crude approach (could be domain
        ! specific).

    !integer(ip), parameter :: nd_global = 1, nd_tonga = 4
    integer(ip), parameter :: nd_global = 1, nd_tonga = 6
        ! Number of domains in different regions.
        ! If nd_tonga = 4, then we use two large high-res domains covering Tongatapu -- robust but slow
        ! If nd_tonga = 6, then we use 4 smaller high-res domains covering Tongatapu -- twice as fast.

    real(dp) :: final_time
        ! Duration of simulation in seconds (start_time = 0.0_dp)

    integer(ip):: j, nd, tonga_regional, global_main
        ! Useful misc local variables
    character(len=charlen) :: stage_file, model_name, run_type, &
        offshore_solver_type, highres_regions, outer_domain_extent, output_style
        ! For reading commandline
    real(dp) :: linear_friction_delay_time, linear_friction_delay_value
        ! For case with linear friction after some time
    logical :: energy_is_finite 
        ! Use this to throw an error if the energy becomes infinite or NA/NaN

    call swals_mpi_init 
        ! Ensure MPI is initialised

#ifndef SPHERICAL
    write(log_output_unit,*) &
        'Code assumes spherical coordinates, but SPHERICAL is not defined'
    call generic_stop
#endif

    ! 
    call parse_commandline_args(stage_file, run_type, final_time, model_name, &
        md%load_balance_file, offshore_solver_type, md%output_basedir, &
        highres_regions, outer_domain_extent, output_style)

    if(output_style == 'animation') then
        ! Write grids pretty often (every nth approximate_writeoutput_timestep)
        write_grids_every_nth_step = 2_ip
    else if(output_style == 'few_grids') then
        ! Write grids on the first and last output step.
        ! These can be useful for identifying issues in the initial condition, or
        ! any late-time instabilities. 
        write_grids_every_nth_step = floor(final_time/approximate_writeout_timestep)
    end if

    !
    ! Basic definition of multidomain
    !

    call program_timer%timer_start('startup_define_multidomain_geometry')

    ! Figure out how many domains are needed
    if(highres_regions == 'none') then
        ! Do not use any of the regional or high-res domains
        nd = nd_global
    else if(highres_regions == 'tonga') then
        ! High res domains everywhere
        nd = nd_global + nd_tonga
    end if
    allocate(md%domains(nd))
        ! nd domains in this model

    ! Determine the dimensions of the global domain
    if(outer_domain_extent == 'global') then
        ! This is a global domain with periodic EW boundary conditions
        global_lw = [360.0_dp, 147.0_dp]    
        global_ll = [-40.0_dp, -79.0_dp]

        ! Evolve the local domains right from the start
        seconds_before_evolve = 0.0_dp

    else if(outer_domain_extent == 'regional') then
        ! This would be used for source-models on the Kermadec-Tonga trench

        global_ll = [138.0_dp, -50.0_dp]
        global_lw = [250.0_dp, 0.0_dp] - global_ll    

        ! Evolve the local domains right from the start
        seconds_before_evolve = 0.0_dp

    else if(outer_domain_extent == 'pacific') then
        ! This can be used for Pacific-wide testing

        global_ll = [90.0_dp, -79.0_dp]
        global_lw = [310.0_dp, 65.0_dp] - global_ll    

        ! We do not need to evolve the inner domain for awhile.
        seconds_before_evolve = 30000.0_dp
    end if


    ! Code logic assumes a periodic domain, make sure it holds
    if(outer_domain_extent == 'global') then
        write(log_output_unit,*) 'Assuming periodic EW boundary'
        md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
            ! This will enforce periodic EW boundary condition as the values equal
            ! the x-range of the multidomain
    end if

    ! Code currently uses reflective NS boundaries -- make sure we have an NS range
    ! large enough to justify that. Probably we want to get rid of this for near-field work.

    !
    ! Setup domain metadata
    !

    do j = 1, nd_global
        ! Global linear domain, split into pieces by longitude
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

        ! A few options for the type of offshore solver
        select case(offshore_solver_type)
        case ("linear_with_manning")
            md%domains(j)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
            md%domains(j)%linear_solver_is_truely_linear = .true.
            ! Try Chezy friction in the offshore domains only -- more energy loss. 
            ! md%domains(j)%friction_type = 'chezy'
        case ("leapfrog_nonlinear")
            md%domains(j)%timestepping_method = 'leapfrog_nonlinear'
        case default
            write(log_output_unit,*) 'Unrecognized offshore_solver_type'
            call generic_stop
        end select

    end do

    ! Below here, the domains are nonlinear. In each region we've got a
    ! regional domain containing one or more high-res nested domains.

    if(highres_regions == 'tonga' ) then
        !
        ! Tonga high-res domains 
        !
        
        tonga_regional = nd_global + 1 
        call md%domains(tonga_regional)%match_geometry_to_parent(&
            ! Better-than-global resolution domain
            parent_domain=md%domains(1), & 
            ! This will correctly nest with the global domain, even if the latter is split into pieces.
            lower_left  = [183.0_dp, -23.0_dp], &
            upper_right = [189.0_dp, -17.0_dp], &
            dx_refinement_factor = 7_ip, &
            timestepping_refinement_factor = 1_ip, &!3_ip ,& 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(tonga_regional)%timestepping_method = 'rk2' 
        md%domains(tonga_regional)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(tonga_regional)%static_before_time = seconds_before_evolve

       
        call md%domains(tonga_regional+1)%match_geometry_to_parent(&
            ! Pretty good nonlinear domain Tongatapu, which will nest the 'very-high-res'
            ! regions
            parent_domain=md%domains(tonga_regional), &
            lower_left  = [184.6_dp, -21.30_dp], &
            upper_right = [185.05_dp,-20.95_dp], &
            dx_refinement_factor = 7_ip, & 
            timestepping_refinement_factor = 3_ip, & 
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(tonga_regional+1)%timestepping_method = 'rk2' 
        md%domains(tonga_regional+1)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') &
            md%domains(tonga_regional+1)%static_before_time = seconds_before_evolve

        if(nd_tonga == 4) then
            !
            ! This approach uses only 2 high-res domains that cover all of Tonga.
            ! It is robust but quite slow.
            !
             
            call md%domains(tonga_regional+2)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain eastern Tongatapu
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.75_dp, -21.28_dp], &
                upper_right = [184.98_dp, -21.119_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 9_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+2)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+2)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+2)%static_before_time = seconds_before_evolve

            call md%domains(tonga_regional+3)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain around NE Tongatapu
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.638_dp, -21.214_dp], &
                upper_right = [184.750_dp, -21.053_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 11_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+3)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+3)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+3)%static_before_time = seconds_before_evolve

        else if(nd_tonga == 6) then
            !
            ! This approach uses smaller domains to cover most of Tongatapu, and is faster than the case with nd_tonga=4.
            !

            call md%domains(tonga_regional+2)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain Nuku'alofa
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.75_dp, -21.207_dp], &
                upper_right = [184.8912_dp, -21.119_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 6_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+2)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+2)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+2)%static_before_time = seconds_before_evolve

            call md%domains(tonga_regional+3)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain around NE Tongatapu
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.8912_dp, -21.180_dp], &
                upper_right = [184.9800_dp, -21.119_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 9_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+3)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+3)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+3)%static_before_time = seconds_before_evolve

            call md%domains(tonga_regional+4)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain on plains west of Nuku'alofa
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.65_dp, -21.192_dp], &
                upper_right = [184.75_dp, -21.119_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 9_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+4)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+4)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+4)%static_before_time = seconds_before_evolve

            call md%domains(tonga_regional+5)%match_geometry_to_parent(&
                ! Higher-res nonlinear domain around NW Tongatapu
                parent_domain=md%domains(tonga_regional+1), &
                lower_left  = [184.643_dp, -21.119_dp], &
                upper_right = [184.690_dp, -21.060_dp], &
                dx_refinement_factor = 5_ip, & 
                timestepping_refinement_factor = 9_ip, & 
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(tonga_regional+5)%timestepping_method = 'rk2' 
            md%domains(tonga_regional+5)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(tonga_regional+5)%static_before_time = seconds_before_evolve
        else
            ! Only specified nd_tonga values are supported.
            write(log_output_unit, *) 'ERROR: no high-res geometries have been defined for the specified nd_tonga value'
            call generic_stop
        end if

    end if

    ! Minor adjustments to domains
    do j = 1, size(md%domains)

        md%domains(j)%nc_grid_output%time_var_store_flag(STG:VH) = .true.
        if(run_type == 'full') md%domains(j)%nc_grid_output%time_var_store_flag(ELV) = .false.
            ! Store stage/UH/VH grids over time, but not ELEVATION (it will 
            ! be stored once anyway).

        !md%domains(j)%theta = 4.0_dp 
        if(j <= tonga_regional) md%domains(j)%theta = 4.0_dp
            ! Use "non-TVD" limiting in nonlinear domains. Less dissipative.

        md%domains(j)%local_timestepping_scale = 0.8_dp
            ! Better nesting stability with this < 1.

        if(any(.not. md%domains(j)%is_nesting_boundary)) md%domains(j)%boundary_subroutine => flather_boundary
            ! Assign a flather boundary condition to boundaries on the exterior. Beware
            ! this can still partly reflect -- so keep it far away

        if(rise_time > 0.0_dp ) md%domains(j)%support_elevation_forcing=.true.
            ! To use a rise-time with elevation forcing, some of the finite-volume solvers need this special option.
            ! It makes them evolve the elevation. Note it's probably not to do this [aside from testing, say], because
            ! it complicates the interpretation of max-stage -- in areas with subsidence, the max-stage may end up
            ! being the pre-subsidence topography. That could be worked-around by resetting domain%max_U in the loop,
            ! but for now we don't do anything like that.

        !md%domains(j)%nc_grid_output%flush_every_n_output_steps = 1_ip
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
    call md%set_point_gauges_from_csv("../gauges/point_gauges/point_gauges_tonga.csv", &
        skip_header=1_ip)
    call program_timer%timer_end('startup_set_gauges')
  
    !
    ! Final setup work.
    ! 
    call program_timer%timer_start('startup_end')

    call md%record_initial_volume()
        ! For mass conservation checks

    do j = 1, size(md%domains)
        write(log_output_unit,*) 'domain: ', j, 'ts: ', md%domains(j)%stationary_timestep_max()
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
            timing_tol = (global_dt/2.01_dp), &
            energy_is_finite = energy_is_finite)

        if(.not. energy_is_finite) then
            ! If the energy is infinite or NaN, we should finish the simulation
            write(log_output_unit, *) "ERROR: ENERGY IS NOT FINITE -- HALTING COMPUTATION EARLY"
            exit
        end if

        call program_timer%timer_end('IO')

        if (md%domains(1)%time > final_time) exit

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

end program
