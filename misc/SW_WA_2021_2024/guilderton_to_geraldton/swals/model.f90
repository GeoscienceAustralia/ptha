! Definitions of the multidomain geometry
#include "model_multidomain_design_mod.f90" 
! Routines to help setup the elevation, forcing, etc.
#include "model_initial_conditions_mod.f90"

program run_model
!!
!! A global-to-local model with high-res sites in specified locations.
!!
!! Basic usage involves these commandline arguments:
!!     ./model stage_file multidomain_design_namelists model_name run_type ambient_sea_level
!!
!! where the arguments are:
!!     stage_file (raster filename with stage perturbation, or a csv file defining a time-varying source inversion)
!!     multidomain_design_namelists (file with namelists controlling high-level setup)
!!     model_name (string used within the output_folder name)
!!     run_type (either 'test' or 'test_load_balance' or 'full', controls the simulation duration)
!!     ambient_sea_level (background sea-level in m, e.g. 0.0 or 0.95)
!!
!! To change the multidomain design, adjust multidomain_design_namelists or edit variables in "model_multidomain_design_mod.f90"
!! To change finer details of the forcing setup, elevation setup, etc, edit "model_local_routines.f90"
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

! Multidomain design
use multidomain_design_mod, only: global_ll, global_lw, global_dx_arcmin, &
    global_dt, final_time_full_runs, final_time_test_runs, &
    use_periodic_EW_multidomain, swals_point_gauge_file, &
    approximate_writeout_timestep, write_grids_every_nth_writeout_timestep, &
    print_every_nth_writeout_timestep, write_gauges_every_nth_writeout_timestep, & 
    time_grids_to_store, nontemporal_grids_to_store, &
    NNL, nesting_domain_extents_file, nesting_domain_timestepping_method, &
    dx_refinement_factors, coarsening_factors, theta_fv_limiter, &
    load_balance_file, real_table_type, nesting_domain_extents, &
    num_domains_per_nesting_level, domain_index_to_nest_with, &
    write_multidomain_design_variables_to_logfile

! Case specific setup and forcing routines
use initial_conditions_mod, only: set_initial_conditions, rise_time, &
    stage_file_is_raster, ambient_sea_level, parse_inputs

implicit none

! md is the main object -- holds all domains, evolves them, etc.
type(multidomain_type) :: md

! Local code-timing object
type(timer_type) :: program_timer

! If 'animate_specific_regions', then write grids regularly after time =
! grid_animation_start_time, but make them very coarse except on selected
! high-res grids.
! Otherwise, write all grids at full resolution, but not very often.
!character(len=charlen), parameter :: grid_output_style = 'regular' !'animate_specific_regions' !
! Write grids frequently after this time (only when grid_output_style =
! 'animate_specific_regions')
!real(dp), parameter :: grid_animation_start_time = 3.0_dp * 3600_dp - global_dt/2.0_dp

! Duration of simulation in seconds (start_time = 0.0_dp)
real(dp) :: final_time

! Non-global domains might not need to evolve at the start (e.g. if we know the
! tsunami takes XX hours to reach them). We can specify the time they evolve to
! speed up models. Crude approach here (could be domain specific). Beware this
! can mess with load balancing.
real(dp) :: seconds_before_evolve

real(dp), parameter :: arcmin_to_deg = 1.0_dp/60.0_dp
! The following distance is useful for approximate timestepping calcs.
real(dp), parameter :: one_arcmin_dist_approx = (7421.29938622_dp/4.0_dp)

! Useful misc local variables
integer(ip):: j, i, nd, p_i, last_setup_domain
!logical :: is_highres_domain
real(dp) :: min_ts, stationary_ts, animation_tag(2,6)

! For command-line info
character(len=charlen) :: stage_file, model_name, run_type, multidomain_design_namelists

! Ensure MPI is initialised
call swals_mpi_init

!
! Basic definition of multidomain
!
call program_timer%timer_start('startup_define_multidomain_geometry')

call parse_inputs(stage_file, multidomain_design_namelists, &
    model_name, run_type, &
    ambient_sea_level, final_time, md%output_basedir, &
    num_domains_per_nesting_level, nesting_domain_extents_file, &
    nesting_domain_extents)

! Depending on the source, the wave might take time to reach the nested domains,
! so we can wait before timestepping
call set_seconds_before_evolve(stage_file, run_type, seconds_before_evolve)

! Figure out how many domains are needed. Later they may be partitioned,
! depending on md%load_balance_file.
nd = sum(num_domains_per_nesting_level)
allocate(md%domains(nd))
write(log_output_unit,*) "Number of domains initially created: ", nd
write(log_output_unit,*) "    Parallel partitioning may increase this number"

md%load_balance_file = load_balance_file

!
! Preliminary global domain setup
!
md%domains(1)%lw = global_lw
md%domains(1)%lower_left = global_ll
md%domains(1)%nx = nint(md%domains(1)%lw/(global_dx_arcmin * arcmin_to_deg)) !
md%domains(1)%dx_refinement_factor = dx_refinement_factors(1)
md%domains(1)%timestepping_refinement_factor = 1_ip
md%domains(1)%timestepping_method = nesting_domain_timestepping_method(1)
md%domains(1)%linear_solver_is_truely_linear = .true.
md%domains(1)%nc_grid_output%spatial_stride = 1_ip ! Can reduce output res

if(use_periodic_EW_multidomain) then
    ! Tell the multidomain that it has periodic E-W behaviour
    md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
end if

! Preliminary setup for nested domains
domain_index_to_nest_with(1) = 1_ip
last_setup_domain = 1_ip
do j = 2, NNL ! Avoid j==1 (global domain, already done)
    call setup_domains_on_nesting_level_NL(md, NL=j, &
        parent_domain_index=domain_index_to_nest_with(j-1), &
        last_setup_domain=last_setup_domain, &
        nesting_NL_domain_extents=nesting_domain_extents(j)%metadata, &
        domain_index_that_can_parent=domain_index_to_nest_with(j))
end do

!
! Minor adjustments to domains prior to setup
!
do j = 1, size(md%domains)

    ! Store the non-empty time grid variables
    md%domains(j)%time_grids_to_store = pack(time_grids_to_store, &
        mask=time_grids_to_store /= "")

    ! Store the non-empty nontemporal grid variables
    md%domains(j)%nontemporal_grids_to_store = pack(nontemporal_grids_to_store, &
        mask=nontemporal_grids_to_store /= "")

    ! To use a rise-time with elevation forcing, some of the finite-volume
    ! solvers need this special option. It makes them evolve the elevation.
    if(rise_time > 0.0_dp .or. (.not. stage_file_is_raster) ) &
        md%domains(j)%support_elevation_forcing=.true.
    ! One should be careful doing this because it complicates the interpretation
    ! of max-stage -- in areas with subsidence, the max-stage may end up being
    ! the pre-subsidence topography. That could be worked-around by resetting
    ! domain%max_U in the loop, but for now we don't do anything like that.

    !if( grid_output_style == 'animate_specific_regions' ) then
    !    ! For some grid_output_style values we store grids very often at some
    !    ! sites, and coarsely elsewhere
    !    call animate_specific_regions
    !end if
end do

call program_timer%timer_end('startup_define_multidomain_geometry')

!
! Allocate domains and prepare parallel comms data structures
!
call program_timer%timer_start('startup_md_setup')
call md%setup()
call program_timer%timer_end('startup_md_setup')

! Ensure the input namelist variables are logged
call write_multidomain_design_variables_to_logfile

!
! Set boundary conditions
!
do j = 1, size(md%domains)
    ! Assign a flather boundary condition to non-nesting boundaries.
    ! Beware this can still partly reflect -- so keep it far away
    if(any(.not. md%domains(j)%is_nesting_boundary)) then
        md%domains(j)%boundary_subroutine => flather_boundary
    end if
end do

!
! Read initial conditions and make them consistent with each other
!
call program_timer%timer_start('startup_set_initial_conditions')
do j = 1, size(md%domains)
    call set_initial_conditions(md%domains(j), stage_file, global_dt, md%all_dx_md)
end do
! Perform a parallel halo exchange so initial conditions are consistent.
call md%make_initial_conditions_consistent()

! Prevent flows in null regions (which do not affect priority domains)
call md%set_null_regions_to_dry()

! Print some summary statistics on memory usage (NB: does this need updating?)
call md%memory_summary()

call program_timer%timer_end('startup_set_initial_conditions')

!
! Gauges
!
call program_timer%timer_start('startup_set_gauges')
call md%set_point_gauges_from_csv(swals_point_gauge_file, skip_header=1_ip)
call program_timer%timer_end('startup_set_gauges')

!
! A few more minor things before evolving
!
call program_timer%timer_start('startup_end')
! For mass conservation checks we record the initial volume
call md%record_initial_volume()

! Print the gravity-wave CFL for each domain, to guide timestepping
do j = 1, size(md%domains)

    stationary_ts = md%domains(j)%stationary_timestep_max()
    min_ts = global_dt/md%domains(j)%timestepping_refinement_factor

    ! Include a flag "UNSTABLE" or "stable" to record whether the min time-step
    ! is within the CFL-permitted timestep. Also info on max-depth/location 
    ! which can be used to troubleshoot.
    write(log_output_unit,*) &
        'domain: ', md%domains(j)%myid, &
        ', stationary_timestep: ', stationary_ts, &
        ', min_timestep:', min_ts, &
        merge(', UNSTABLE ', ',   stable ', (stationary_ts < min_ts)), &
        ', min(elevation):', minval(md%domains(j)%U(:,:,ELV)),&
        ', centre:', md%domains(j)%lower_left + md%domains(j)%lw/2, &
        ', dx:', md%domains(j)%dx, &
        ', theta: ', md%domains(j)%theta

    ! Reset the domain timers, so that load-balancing only sees the evolve info
    call md%domains(j)%timer%reset
end do

write(log_output_unit,*) 'End setup'
call program_timer%timer_end('startup_end')

flush(log_output_unit)

!
! Main evolve loop
!
do while (.true.)

    ! Print and write outputs
    call program_timer%timer_start('IO')
    call md%write_outputs_and_print_statistics(&
        ! Time between writes is ~= "approximate_writeout_timestep"
        approximate_writeout_frequency=approximate_writeout_timestep, &
        write_grids_less_often = write_grids_every_nth_writeout_timestep, &
        write_gauges_less_often = write_gauges_every_nth_writeout_timestep, &
        print_less_often = print_every_nth_writeout_timestep, & !1_ip,&
        ! If the time is within timing_tol of the desired output time, then
        ! the output will be written even if the time < desired_output_time.
        timing_tol = (global_dt/2.01_dp))
    call program_timer%timer_end('IO')

    if (md%domains(1)%time > final_time) exit

    !if((grid_output_style == 'animate_specific_regions') .and. &
    !   (md%domains(1)%time > grid_animation_start_time)) then
    !    ! Write grids frequently from here on in.
    !    write_grids_every_nth_writeout_timestep = 1_ip
    !end if

    ! Evolve the model by global_dt seconds
    call program_timer%timer_start('evolve')
    call md%evolve_one_step(global_dt)
    call program_timer%timer_end('evolve')

end do

! Close files, and print timing info required for load balancing
call md%finalise_and_print_timers

write(log_output_unit,*) ''
write(log_output_unit, *) 'Program timer'
write(log_output_unit, *) ''
call program_timer%print(log_output_unit)
flush(log_output_unit)

call swals_mpi_finalize

contains

subroutine setup_domains_on_nesting_level_NL(md, NL, parent_domain_index, &
    last_setup_domain, nesting_NL_domain_extents, domain_index_that_can_parent)
    !! Setup nesting level NL for NL = 2, 3, ...

    use global_mod, only: gravity, pi
    use multidomain_design_mod, only: r_ll, r_ur, r_coarse, r_dpth

    type(multidomain_type), intent(inout) :: md !! The multidomain
    integer(ip), intent(in) :: &
        NL, & !! Nesting level
        parent_domain_index 
            !! Index of md%domains(:) in earlier nesting levels to use as a parent domain
    integer(ip), intent(inout) :: last_setup_domain 
        !! Last index of md%domains(:) that was setup (i.e. in a previous nesting level)
    real(dp), intent(in) :: nesting_NL_domain_extents(:,:) 
        !! Metadata on domains in nesting level NL
    integer(ip), intent(out) :: domain_index_that_can_parent
        !! On output this will correspond to a domain on nesting level NL that can serve 
        !! as a parent domain for any subsequent nesting level. Basically any domain will
        !! do, except for coarsened domains.

    integer(ip) :: j, di, local_dx_refinement_factor, local_timestepping_refinement
    real(dp) :: local_max_depth, local_approx_wave_speed, local_coslat, &
        local_approx_grid_size, local_approx_ts
    logical :: have_found_domain_that_can_parent, is_fine_domain

    have_found_domain_that_can_parent = .FALSE.

    do j = 1, size(nesting_NL_domain_extents, 2)

        di = last_setup_domain + j

        ! Record whether domain has regular (fine) resolution, or is flagged for coarsening.
        is_fine_domain = (nesting_NL_domain_extents(r_coarse, j) == 0.0_dp)
        if(coarsening_factors(NL) < 0 .and. (.not. is_fine_domain)) then
            write(log_output_unit, *) 'ERROR: On nesting layer ', NL, ' the input file row ', j, &
                'suggests a domain should be coarsened, but the coarsening_factor(NL) is negative ', &
                '(meaning no coarsening). Maybe you want coarsening_factor(NL) = 1 (???).'
        end if

        ! Find a single domain on this nesting level that
        ! can serve as a parent domain to the next nesting level. Any
        ! domain that is not being coarsened will be OK.
        if((.not. have_found_domain_that_can_parent) .and. is_fine_domain) then
            domain_index_that_can_parent = di
            have_found_domain_that_can_parent = .TRUE.
        end if

        !
        ! Get information for approx timestepping calculation.
        ! At this point we don't know the allowable timestep and cannot compute 
        ! it (even in the still-water case) because we haven't yet read the
        ! elevation data. But we need a useful lower bound to define the
        ! timestepping_refinement_factor and setup the multidomain halos. 
        ! Later, the domain might take larger timesteps (if compiled with 
        ! -DLOCAL_TIMESTEP_PARTITIONED_DOMAINS)
        !

        ! An estimate of the domain's max depth was provided in the input file.
        ! It is usually imprecise (coarse data + rough estimate of nesting halo size).
        local_max_depth = -nesting_NL_domain_extents(r_dpth,j)
        ! Conservative approximate wave speed from sqrt(g * d)
        local_approx_wave_speed = sqrt(gravity * max(local_max_depth*1.25_dp, 100.0_dp))

        ! dx refinement factor for the current domain, which might be coarsened.
        local_dx_refinement_factor = merge(dx_refinement_factors(NL), &
            (dx_refinement_factors(NL)/coarsening_factors(NL)), &
            is_fine_domain)

        if(local_dx_refinement_factor <= 0) then
            write(log_output_unit, *) "Error: non-positive local_dx_refinement_factor", &
                NL, j, di, local_dx_refinement_factor, nesting_NL_domain_extents(7,j)
            call generic_stop
        end if
        
        ! Rough approximation of the local grid size
        local_coslat = cos(nesting_NL_domain_extents(r_ll(2),j)/180.0_dp * pi)
        local_approx_grid_size = local_coslat * &
            one_arcmin_dist_approx * global_dx_arcmin / &
            (product(dx_refinement_factors(1:NL-1)) * local_dx_refinement_factor)

        ! The current domain could take (approximately) the following time-step
        ! This assumes a CFL of 0.5_dp (appropriate for the second order FV schemes)
        local_approx_ts = &
            local_approx_grid_size / (2.0_dp * local_approx_wave_speed)
        local_timestepping_refinement = &
            max(1_ip, ceiling(global_dt / local_approx_ts ))

        !
        ! Set domain parameters and possibly tweak location for nesting
        !
        call md%domains(di)%match_geometry_to_parent(&
            parent_domain=md%domains(parent_domain_index), &
            lower_left  = nesting_NL_domain_extents(r_ll,j), &
            upper_right = nesting_NL_domain_extents(r_ur,j), &
            dx_refinement_factor = local_dx_refinement_factor, &
            timestepping_refinement_factor = local_timestepping_refinement, &
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(di)%timestepping_method = nesting_domain_timestepping_method(NL)
        md%domains(di)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') md%domains(di)%static_before_time = seconds_before_evolve
        md%domains(di)%theta = theta_fv_limiter(NL)

        ! Check the impact of a small eddy-viscosity
        !md%domains(di)%use_eddy_viscosity = .true.
        !md%domains(di)%eddy_visc_constants = [0.5_dp, 0.0_dp]

    end do

    last_setup_domain = di

end subroutine

subroutine set_seconds_before_evolve(stage_file, run_type, seconds_before_evolve)
    ! Depending on the tsunami source, we might be able to delay computations
    ! in the nested domains for a given time period (i.e. before the wave 
    ! reaches them from the global domain).
    character(len=charlen), intent(in) :: stage_file, run_type
    real(dp), intent(inout) :: seconds_before_evolve

    seconds_before_evolve = 0.0_dp

    !if(index(run_type, 'test') /= 0) then
    !    ! Test runs always evolve straight away
    !    seconds_before_evolve = 0.0_dp
    !else if(index(stage_file, '_43731_') /= 0) then
    !    ! Scenario from QFES workshop
    !    ! FIXME: Determine a systematic approach
    !    seconds_before_evolve = 3600 * 3.0_dp
    !end if
end subroutine

!subroutine animate_specific_regions
!    ! Do whatever is needed to setup the animation.
!
!    write(log_output_unit, *) &
!        "Error: Have not yet specified specific regions to animate.", &
!        "Example code commented below"
!    call generic_stop
!
!    ! Store stage grids over time, but not UH/VH/ELEVATION
!    !md%domains(j)%time_grids_to_store = ['stage']
!
!    ! ! By default store grids very coarsely (adjusted for some grids
!    ! ! below)
!    ! md%domains(j)%nc_grid_output%spatial_stride = 100_ip
!
!    ! is_highres_domain = (j > sum(num_domains_per_nesting_level(1:2))
!    !
!    ! if(is_highres_domain) then
!    !     ! Newcastle special case -- if domains contain any of these
!    !     ! points, store them at high-res
!    !     animation_tag(:,1) = [151.6718_dp, -32.9548_dp]
!    !     animation_tag(:,2) = [151.6718_dp, -32.7906_dp]
!    !     animation_tag(:,3) = [151.8517_dp, -32.7906_dp]
!    !     animation_tag(:,4) = [151.8517_dp, -32.9548_dp]
!    !     do p_i = 1, 4
!    !         if(all(md%domains(j)%lower_left                    < animation_tag(:,p_i)) .and. &
!    !            all(md%domains(j)%lower_left + md%domains(j)%lw > animation_tag(:,p_i))) then
!    !            md%domains(j)%nc_grid_output%spatial_stride = 1
!    !         end if
!    !     end do
!    ! end if
!end subroutine

end program
