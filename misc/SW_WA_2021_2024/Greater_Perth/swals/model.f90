!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module "local_routines" with various helper subroutines.
#include "model_local_routines.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_model
    !!
    !! A global-to-local model with high-res sites in specified locations.
    !!
    !! Basic usage involves these commandline arguments:
    !!     ./model stage_file model_name run_type coarse_nonlinear_domain_extents_file fine_nonlinear_domain_extents_file highres_nonlinear_domain_extents_file load_balance_file ambient_sea_level
    !!
    !! where the arguments are:
    !!     stage_file (raster filename with stage perturbation)
    !!     model_name (used within the output_folder name) 
    !!     run_type (either 'test' or 'test_load_balance' or 'full') 
    !!     coarse_nonlinear_domain_extents_file (file with extents of coarse nonlinear domains)
    !!     fine_nonlinear_domain_extents_file (file with extents of fine nonlinear domains)
    !!     highres_nonlinear_domain_extents_file (file with extents of very high res nonlinear domains)
    !!     load_balance_file (file with load balancing metadata, or '' to use crude defaults)
    !!     ambient_sea_level (background sea-level in m, e.g. 0.0 or 0.95)
    !!        
    !!

    !
    ! Imports from SWALS 
    !
    use global_mod, only: ip, dp, charlen, gravity, pi
        ! Integer/real precision and default character length
    use domain_mod, only: STG, UH, VH, ELV
        ! Indices of stage, depth-integrated velocities 
        ! (UH = east, VH = north) and elevation in the domain%U array.
    use multidomain_mod, only: multidomain_type
        ! The main type that holds all domains and evolves them
    use boundary_mod, only: flather_boundary
        ! flather_boundary can be used as a passive boundary condition that 
        ! tries to maintain a given stage, while letting the waves pass.
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
    use file_io_mod, only: read_csv_into_array
        ! Read simple csv files
    !
    ! Case specific helper routines
    !
    use local_routines, only: set_initial_conditions, parse_commandline_args, &
        rise_time

    implicit none

    ! md is the main object -- holds all domains, evolves them, etc.
    type(multidomain_type) :: md

    ! Local code-timing object
    type(timer_type) :: program_timer

    ! Lower-left corner coordinate of multidomain in degrees lon,lat
    !real(dp), parameter:: global_ll(2) = [-40.0_dp, -79.0_dp]
    real(dp), parameter:: global_ll(2) = [16.0_dp, -75.0_dp]

    ! Length/width of multidomain in degrees lon,lat
    !real(dp), parameter:: global_lw(2) = [360.0_dp, 147.0_dp] 
    real(dp), parameter:: global_lw(2) = [156.0_dp , 33.0_dp] - global_ll 

    ! Increase this to decrease the cell side-length by mesh_refine 
    ! (i.e. for convergence testing). 4.0_dp --> 1 arcmin in the coarse 
    ! global domain
    real(dp), parameter :: mesh_refine = 4.0_dp

    ! The global time-step in the multidomain. Must satisfy CFL condition 
    ! everywhere (in combination with local timestepping that is specified
    ! when defining the domains below)
    real(dp), parameter ::  global_dt = 4.5_dp / (mesh_refine)

    ! If 'animate_specific_regions', then write grids regularly after time
    ! = grid_animation_start_time, but make them very coarse except on 
    ! selected high-res grids.
    ! Otherwise, write all grids at full resolution, but not very often.
    character(len=charlen), parameter :: grid_output_style = 'regular' !'animate_specific_regions' ! 
    ! Write grids frequently after this time (only when grid_output_style 
    ! = 'animate_specific_regions')
    real(dp), parameter :: grid_animation_start_time = &
        3.0_dp * 3600_dp - global_dt/2.0_dp

    ! Approx timestep between any file outputs
    real(dp), parameter :: approximate_writeout_timestep = 30.0_dp
    ! Optionally write grids less often than the writeout timestep, to keep
    ! file-size down. BEWARE THE FOLLOWING EXCEPTION: If we use 
    ! grid_output_style='animate_specific_regions', then grids will in any 
    ! case be written every time-step after time=grid_animation_start_time
    integer(ip) :: write_grids_every_nth_step = 400_ip 
    ! Optionally print outputs less often than the writeout timestep, to 
    ! keep the file-size down
    integer(ip), parameter :: print_every_nth_step = 10_ip
    ! Optionally write gauges less often than the writeout timestep, to 
    ! keep the file-size down
    integer(ip), parameter :: write_gauges_every_nth_step = 1_ip
    ! Grid-size refinement of each 'nesting level' relative to the one above
    ! The indices in order correspond to the 'global domain' (always 1), the 'coarse nonlinear domains',
    ! the 'fine nonlinear domains', and the 'highres nonlinear domains'
    integer(ip), parameter :: dx_refinement_factors(4) = [ 1_ip,  9_ip, 6_ip,  3_ip]
    ! In this model we can optionally coarsen some domains at each level. Coarsening is applied to
    ! domains that have an entry of 1 in the 'coarsen' column of the nonlinear_domain_extents.
    ! The coarsening_factors correspond to (in order) the 'global domain', the 'coarse nonlinear domains',
    ! the 'fine nonlinear domains', and the 'highres nonlinear domains'. Negative values
    ! denote no coarsening. Positive values mean that the 'coarsend' domains are refined by
    ! (dx_refinement_factors/coarsening_factors), e.g. a value of 2 means the grid-side-length is doubled.
    integer(ip), parameter :: coarsening_factors(4)    = [-1_ip, -1_ip, 2_ip, -1_ip]

    ! Non-global domains might not need to evolve at the start -- this 
    ! specifies the time before which they evolve. Crude approach (could be
    ! domain specific).
    real(dp), parameter :: seconds_before_evolve = 0.0_dp 

    ! Split the global domain into this many different regions.
    integer(ip), parameter :: nd_global = 1

    ! Duration of simulation in seconds (start_time = 0.0_dp)
    real(dp) :: final_time

    ! Useful misc local variables
    integer(ip):: j, nd, p_i, last_setup_domain, &
        num_coarse_nests, num_fine_nests, num_highres_nests, &
        global_domain_index_to_nest_with, coarse_domain_index_to_nest_with, &
        fine_domain_index_to_nest_with
    logical :: is_highres_domain
    real(dp) :: min_ts, stationary_ts, animation_tag(2,6)

    real(dp), parameter :: fourarcmin = 4.0_dp/60.0_dp
    real(dp), parameter :: fourarcmin_dist = 7421.29938622_dp

    ! Some commandline info
    character(len=charlen) :: stage_file, model_name, run_type, &
        coarse_nonlinear_domain_extents_file, fine_nonlinear_domain_extents_file,&
        highres_nonlinear_domain_extents_file

    ! bounding coordinates of domains, read from the 
    ! (coarse/fine)_nonlinear_domain_extents_file
    real(dp), allocatable :: coarse_nonlinear_domain_extents(:,:), &
        fine_nonlinear_domain_extents(:,:), &
        highres_nonlinear_domain_extents(:,:)

    call swals_mpi_init 
        ! Ensure MPI is initialised

#ifndef SPHERICAL
    write(log_output_unit,*) &
        'Code assumes spherical coordinates, but SPHERICAL is not defined'
    call generic_stop
#endif

    call parse_commandline_args(stage_file, run_type, final_time, model_name, &
        coarse_nonlinear_domain_extents_file, &
        fine_nonlinear_domain_extents_file, &
        highres_nonlinear_domain_extents_file, md%load_balance_file, &
        md%output_basedir)

    !
    ! Basic definition of multidomain
    !
    call program_timer%timer_start('startup_define_multidomain_geometry')
    
    if(global_lw(1) == 360.0_dp) then
        ! Periodic EW boundary condition
        md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
    end if

    num_coarse_nests = 0
    if(coarse_nonlinear_domain_extents_file /= '') then
        ! Read lower-left and upper-right extents of 'coarse' nonlinear domains, 
        ! which have a higher level of refinement than the global domain 
        ! [which they are nested inside]
        call read_csv_into_array(coarse_nonlinear_domain_extents, &
            coarse_nonlinear_domain_extents_file, skip_header=1_ip)
        num_coarse_nests = size(coarse_nonlinear_domain_extents, 2)
    end if

    num_fine_nests = 0
    if(fine_nonlinear_domain_extents_file /= '') then
        ! Read lower-left and upper-right extents of 'fine' nonlinear domains -- 
        ! which have a higher level of refinement than the coarse nonlinear domains
        ! [which they are nested inside].
        call read_csv_into_array(fine_nonlinear_domain_extents, &
            fine_nonlinear_domain_extents_file, skip_header=1_ip)
        num_fine_nests = size(fine_nonlinear_domain_extents, 2)
    end if

    num_highres_nests = 0
    if(highres_nonlinear_domain_extents_file /= '') then
        ! Read lower-left and upper-right extents of 'highres' nonlinear domains -- 
        ! which have a higher level of refinement than the fine nonlinear domains
        ! [which they are nested inside].
        call read_csv_into_array(highres_nonlinear_domain_extents, &
            highres_nonlinear_domain_extents_file, skip_header=1_ip)
        num_highres_nests = size(highres_nonlinear_domain_extents, 2)
    end if

    ! Figure out how many domains are needed -- [may be changed by parallel 
    ! partitioning]
    nd = nd_global + num_coarse_nests + num_fine_nests + num_highres_nests
    allocate(md%domains(nd))
    write(log_output_unit,*) "Number of domains initially created: ", nd
    write(log_output_unit,*) "    Parallel partitioning may increase this number"


    !
    ! Setup global domain
    !
    do j = 1, nd_global
        ! Global linear domain, split into "nd_global" pieces by longitude
        ! Note the load balance file can alternatively be used to partition it.

        md%domains(j)%lw = [global_lw(1)*1.0_dp/nd_global, global_lw(2)]
        md%domains(j)%lower_left = &
            [global_ll(1) + (j-1)*global_lw(1)*1.0_dp/nd_global, global_ll(2)]
        md%domains(j)%dx = fourarcmin * [1.0_dp, 1.0_dp] / mesh_refine
        md%domains(j)%nx = nint(md%domains(j)%lw/md%domains(j)%dx)
        md%domains(j)%dx_refinement_factor = dx_refinement_factors(1)
        md%domains(j)%timestepping_refinement_factor = 1_ip

        ! Reduce output file size by only saving every n'th cell 
        md%domains(j)%nc_grid_output%spatial_stride = 1_ip 

        md%domains(j)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
        md%domains(j)%linear_solver_is_truely_linear = .true.

    end do
    last_setup_domain = nd_global

    !
    ! Setup the 'coarse nonlinear domains'
    !
    if(num_coarse_nests > 0) then
        ! All 'coarse nonlinear domains' will use this domain-index as a parent. 
        global_domain_index_to_nest_with = 1_ip 
        call setup_nesting_level_domains(2_ip, global_domain_index_to_nest_with, &
            last_setup_domain, coarse_nonlinear_domain_extents)        
    end if

    !
    ! Setup the fine nonlinear domains
    !
    if(num_fine_nests > 0) then
        ! All fine nonlinear domains will use this domain-index as a parent domain
        coarse_domain_index_to_nest_with = nd_global + 1_ip
        call setup_nesting_level_domains(3_ip, coarse_domain_index_to_nest_with, &
            last_setup_domain, fine_nonlinear_domain_extents)
    end if

    !
    ! Set up "highres nonlinear domains". 
    !
    if(num_highres_nests > 0) then
        ! These are finer than the fine nonlinear domains, and are expected to resolve 
        ! the coast very well. They nest with a "non-coarsened" fine nonlinear domain.
        fine_domain_index_to_nest_with = -HUGE(1_ip)
        do j = 1, size(fine_nonlinear_domain_extents, 2)
            if(fine_nonlinear_domain_extents(7,j) == 0.0_dp) then
                fine_domain_index_to_nest_with = nd_global + num_coarse_nests + j
                exit
            end if
        end do
        call setup_nesting_level_domains(4_ip, fine_domain_index_to_nest_with, &
            last_setup_domain, highres_nonlinear_domain_extents)
    end if

    !
    ! Minor adjustments to domains
    !
    do j = 1, size(md%domains)

        ! Store stage/UH/VH grids over time, but not ELEVATION (it will 
        ! be stored once anyway).
        !md%domains(j)%time_grids_to_store = [character(len=charlen):: 'stage', 'uh', 'vh']

        ! Do not store time-grids
        md%domains(j)%time_grids_to_store = [ '' ]

        ! Useful summary statistics
        md%domains(j)%nontemporal_grids_to_store = [character(len=charlen):: &
            'max_stage', 'max_speed', 'max_flux', 'arrival_time', 'elevation0']

        ! Use "non-TVD" limiting in nonlinear domains. Less dissipative.
        md%domains(j)%theta = 4.0_dp

        ! Assign a flather boundary condition to non-nesting boundaries. 
        ! Beware this can still partly reflect -- so keep it far away
        if(any(.not. md%domains(j)%is_nesting_boundary)) then
            md%domains(j)%boundary_subroutine => flather_boundary
        end if

        ! To use a rise-time with elevation forcing, some of the finite-volume 
        ! solvers need this special option. It makes them evolve the elevation.
        ! Note one should be careful doing this because it complicates the
        ! interpretation of max-stage -- in areas with subsidence, the 
        ! max-stage may end up being the pre-subsidence topography. That could
        ! be worked-around by resetting domain%max_U in the loop, but for now
        ! we don't do anything like that.
        if(rise_time > 0.0_dp ) md%domains(j)%support_elevation_forcing=.true.

        if( grid_output_style == 'animate_specific_regions' ) then
            ! For some grid_output_style values we store grids very often at some 
            ! sites, and coarsely elsewhere

            ! Store stage grids over time, but not UH/VH/ELEVATION
            md%domains(j)%time_grids_to_store = ['stage'] 
 
            write(log_output_unit, *) "Error: Have not yet specified specific regions to animate. Example code commented below"
            call generic_stop

            ! ! By default store grids very coarsely (adjusted for some grids 
            ! ! below)
            ! md%domains(j)%nc_grid_output%spatial_stride = 100_ip 

            ! is_highres_domain = &
            !     (j > nd_global + size(coarse_nonlinear_domain_extents, 2))
            ! 
            ! if(is_highres_domain) then 
            !     ! Newcastle special case -- if domains contain any of these 
            !     ! points, store them at high-res
            !     animation_tag(:,1) = [151.6718_dp, -32.9548_dp]
            !     animation_tag(:,2) = [151.6718_dp, -32.7906_dp]
            !     animation_tag(:,3) = [151.8517_dp, -32.7906_dp]
            !     animation_tag(:,4) = [151.8517_dp, -32.9548_dp]
            !     do p_i = 1, 4
            !         if(all(md%domains(j)%lower_left                    < animation_tag(:,p_i)) .and. &
            !            all(md%domains(j)%lower_left + md%domains(j)%lw > animation_tag(:,p_i))) then
            !            md%domains(j)%nc_grid_output%spatial_stride = 1
            !         end if
            !     end do 
            ! end if
        end if
    end do

    call program_timer%timer_end('startup_define_multidomain_geometry')

    !
    ! Allocate domains and prepare parallel comms data structures
    !
    call program_timer%timer_start('startup_md_setup')
    call md%setup()
    call program_timer%timer_end('startup_md_setup')

    !
    ! Read initial conditions and make them consistent with each other
    !
    call program_timer%timer_start('startup_set_initial_conditions')
    do j = 1, size(md%domains)
        call set_initial_conditions(md%domains(j), stage_file, global_dt)
    end do
    call md%memory_summary()

    ! Perform a parallel halo exchange so initial conditions are consistent.
    call md%make_initial_conditions_consistent()

    ! Prevent flows in null regions (which do not affect priority domains)
    call md%set_null_regions_to_dry()

    call program_timer%timer_end('startup_set_initial_conditions')
  
    !
    ! Gauges 
    ! 
    call program_timer%timer_start('startup_set_gauges')
    ! Make a file
    call md%set_point_gauges_from_csv("../gauges/point_gauges_151221.csv", &
        skip_header=1_ip)
    call program_timer%timer_end('startup_set_gauges')
  
    !
    ! Final setup work.
    ! 
    call program_timer%timer_start('startup_end')

    ! For mass conservation checks we record the initial volume
    call md%record_initial_volume()

    do j = 1, size(md%domains)

        ! Print the gravity-wave CFL limit to guide timestepping
        stationary_ts = md%domains(j)%stationary_timestep_max()
        min_ts = global_dt/md%domains(j)%timestepping_refinement_factor
        write(log_output_unit,*) &
            'domain: ', md%domains(j)%myid, &
            ', ts: ', stationary_ts, &
            ', min_ts:', min_ts, &
            merge(', UNSTABLE ', ',   stable ', (stationary_ts < min_ts)), &
            ', min(elevation):', minval(md%domains(j)%U(:,:,ELV)),&
            ', centre:', md%domains(j)%lower_left + md%domains(j)%lw/2, &
            ', dx:', md%domains(j)%dx

        ! Reset the domain timers, so that load-balancing only sees the 
        ! evolve info
        call md%domains(j)%timer%reset
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
    
        if((grid_output_style == 'animate_specific_regions') .and. &
           (md%domains(1)%time > grid_animation_start_time)) then
            ! Write grids frequently from here on in.
            write_grids_every_nth_step = 1_ip
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

    ! Workhorse subroutine to setup each level of nesting
    subroutine setup_nesting_level_domains(NL, domain_index_to_match, &
        last_setup_domain, nonlinear_domain_extents)

        integer(ip), intent(in) :: NL, & ! Nesting level
            domain_index_to_match ! Index of md%domains(:) used to ensure domain bbox suitable for nesting
        integer(ip), intent(inout) :: last_setup_domain ! Last index of md%domains(:) that was setup
        real(dp), intent(in) :: nonlinear_domain_extents(:,:) ! Metadata for nesting level NL

        integer(ip) :: j, di, local_dx_refinement_factor, local_timestepping_refinement 
        real(dp) :: local_max_depth, local_approx_wave_speed, local_coslat, &
            local_approx_grid_size, local_approx_ts

        do j = 1, size(nonlinear_domain_extents, 2)

            di = last_setup_domain + j 

            local_max_depth = -nonlinear_domain_extents(5,j)
            local_approx_wave_speed = sqrt(gravity * max(&
                local_max_depth*1.25_dp, 100.0_dp))

            local_dx_refinement_factor = merge(dx_refinement_factors(NL), &
                (dx_refinement_factors(NL)/coarsening_factors(NL)), &
                nonlinear_domain_extents(7,j) == 0.0_dp)

            if(local_dx_refinement_factor < 0) then
                write(log_output_unit, *) "Error: negative local_dx_refinement_factor", &
                    NL, j, di, local_dx_refinement_factor, nonlinear_domain_extents(7,j)
                call generic_stop
            end if

            local_coslat = cos(nonlinear_domain_extents(2,j)/180.0_dp * pi)
            local_approx_grid_size = local_coslat * fourarcmin_dist / &
                (product(dx_refinement_factors(1:NL-1)) * &
                 local_dx_refinement_factor * mesh_refine)

            ! The current domain could take (approximately) the following time-step
            local_approx_ts = &
                local_approx_grid_size / (2.0_dp * local_approx_wave_speed)
            local_timestepping_refinement = &
                max(1_ip, ceiling(global_dt / local_approx_ts ))

            call md%domains(di)%match_geometry_to_parent(& 
                parent_domain=md%domains(domain_index_to_match), & 
                lower_left  = nonlinear_domain_extents(1:2,j), &
                upper_right = nonlinear_domain_extents(3:4,j), &
                dx_refinement_factor = local_dx_refinement_factor, &
                timestepping_refinement_factor = local_timestepping_refinement, &
                rounding_method='nearest', &
                recursive_nesting=.false.)
            md%domains(di)%timestepping_method = 'rk2' 
            md%domains(di)%nc_grid_output%spatial_stride = 1
            if(run_type == 'full') &
                md%domains(di)%static_before_time = seconds_before_evolve

            ! Check the impact of a small eddy-viscosity
            !md%domains(di)%use_eddy_viscosity = .true.
            !md%domains(di)%eddy_visc_constants = [0.5_dp, 0.0_dp]
        
        end do
    
        last_setup_domain = di
    
    end subroutine

end program
