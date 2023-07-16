!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module "local_routines" with various helper subroutines.
#include "model_local_routines.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_model
!!
!! A global-to-local model with high-res sites in specified locations.
!!
!! Basic usage involves these commandline arguments:
!!     ./model stage_file model_name run_type coarse_nonlinear_domain_extents_file fine_nonlinear_domain_extents_file highres_nonlinear_domain_extents_file extrahighres_nonlinear_domain_extents_file load_balance_file ambient_sea_level
!!
!! where the arguments are:
!!     stage_file (raster filename with stage perturbation, or a csv file defining a time-varying source inversion)
!!     model_name (string used within the output_folder name)
!!     run_type (either 'test' or 'test_load_balance' or 'full', controls the simulation duration)
!!     coarse_nonlinear_domain_extents_file (file with extents of coarse nonlinear domains)
!!     fine_nonlinear_domain_extents_file (file with extents of fine nonlinear domains)
!!     highres_nonlinear_domain_extents_file (file with extents of very high res nonlinear domains)
!!     extrahighres_nonlinear_domain_extents_file (file with extents for yet higher res nonlinear domains)
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
! Case specific utilities
!
use local_routines, only: set_initial_conditions, parse_commandline_args, &
    rise_time, stage_file_is_raster

implicit none

! md is the main object -- holds all domains, evolves them, etc.
type(multidomain_type) :: md

! Local code-timing object
type(timer_type) :: program_timer

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the large-scale domain extent, timestepping and scale for resolution
!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Lower-left corner coordinate of multidomain in degrees lon,lat
real(dp), parameter:: global_ll(2) = [16.0_dp, -75.0_dp]

! Length/width of multidomain in degrees lon,lat
real(dp), parameter:: global_lw(2) = [156.0_dp , 33.0_dp] - global_ll

! Increase this to decrease the cell side-length of ALL domains by
! mesh_refine (i.e. for convergence testing). 4.0_dp --> 1 arcmin
! in the coarse global domain
real(dp), parameter :: mesh_refine = 4.0_dp

! The global time-step in the multidomain. Must satisfy CFL condition
! everywhere (in combination with local timestepping that is specified
! when defining the domains below)
real(dp), parameter ::  global_dt = 4.5_dp / (mesh_refine)

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! High-level nesting information
!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of nesting levels & their refinement (including the global domain as
! one level, and the other domains provided on the commandline)
integer(ip), parameter :: NNL = 5_ip

! Files containing domain definitions (read from command-line arguments)
character(len=charlen) :: coarse_nonlinear_domain_extents_file, &
    fine_nonlinear_domain_extents_file, highres_nonlinear_domain_extents_file, &
    extrahighres_nonlinear_domain_extents_file

! We will make an array with one entry per nest level, holding domain metadata.
! This type helps
type real_table_type
    real(dp), allocatable :: metadata(:,:)
end type
! Row indices in real_table_type%metadata
integer(ip) :: r_ll(2) = [1_ip, 2_ip], & ! Lower left lon,lat
               r_ur(2) = [3_ip, 4_ip], & ! Upper right lon,lat
               r_dpth = 5_ip, &          ! Max depth (estimate) in domain
               r_coarse = 7_ip           ! Coarsen the domain? (0.0 = No, 1.0 = Yes)

! File & data-structure with domain metadata each level of nesting
character(len=charlen) :: nesting_domain_extents_file(NNL)
type(real_table_type) :: nesting_domain_extents(NNL)

! Other useful nesting constants
integer(ip) :: num_nests(NNL), domain_index_to_nest_with(NNL)

! Grid-size refinement of each 'nesting level' relative to the one above
! The indices in order correspond to 1. 'global domain' (refinement == 1_ip),
! 2. 'coarse nonlinear domains', 3. 'fine nonlinear domains',
! 4. 'highres nonlinear domains', 5. 'extrahighres_nonlinear_domains'
integer(ip), parameter :: dx_refinement_factors(NNL) = [1_ip, 9_ip, 6_ip, 3_ip, 3_ip]

! Optionally coarsen some of domains at each nesting_level. This can be useful
! to treat occasional domains in deep water, which can slow models down
! without coarsening, and usually don't require such high resolution.
! Coarsening is applied to the 'j'th domain in nesting-level 'N' if:
!     nesting_domain_extents(N)%metadata(r_coarse,j) = 1.0_dp
! Positive values of coarsening_factor(N) imply these domains are refined by
! (dx_refinement_factors/coarsening_factors)
!     i.e. a value of 2_ip means the cell-side-length is doubled.
! Negative values imply no coarsening will be used.
! Array order is the same as "dx_refinement_factors"
integer(ip), parameter :: coarsening_factors(NNL) = [-1_ip, -1_ip, 2_ip, -1_ip, -1_ip]

! Non-global domains might not need to evolve at the start (e.g. if we know the
! tsunami takes XX hours to reach them). We can specify the time they evolve to
! speed up models. Crude approach here (could be domain specific). Beware this
! can mess with load balancing.
real(dp), parameter :: seconds_before_evolve = 0.0_dp

! Optionally do local smoothing of elevation along coarse-to-fine nesting boundaries
logical, parameter :: smooth_elevation_along_nesting_boundaries = .TRUE. !.FALSE.


!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Frequency of file output for grids, gauges, and printing
!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Approx timestep between any file outputs
real(dp), parameter :: approximate_writeout_timestep = 30.0_dp

! Optionally write grids less often than the writeout timestep, to keep
! file-size down. BEWARE AN EXCEPTION: If we use
! grid_output_style='animate_specific_regions', grids will regardless be written
! every time-step after time=grid_animation_start_time
integer(ip) :: write_grids_every_nth_step = 400_ip

! Optionally print outputs less often than the writeout timestep, to keep the
! file-size down
integer(ip), parameter :: print_every_nth_step = 10_ip

! Optionally write gauges less often than the writeout timestep, to keep the
! file-size down
integer(ip), parameter :: write_gauges_every_nth_step = 1_ip

! If 'animate_specific_regions', then write grids regularly after time =
! grid_animation_start_time, but make them very coarse except on selected
! high-res grids.
! Otherwise, write all grids at full resolution, but not very often.
character(len=charlen), parameter :: grid_output_style = 'regular' !'animate_specific_regions' !
! Write grids frequently after this time (only when grid_output_style =
! 'animate_specific_regions')
real(dp), parameter :: grid_animation_start_time = 3.0_dp * 3600_dp - global_dt/2.0_dp

!@!!!!!!!!!!!!!
! Other details
!@!!!!!!!!!!!!!

! Duration of simulation in seconds (start_time = 0.0_dp)
real(dp) :: final_time

! Useful misc local variables
integer(ip):: j, i, nd, p_i, last_setup_domain
!logical :: is_highres_domain
real(dp) :: min_ts, stationary_ts, animation_tag(2,6)

! Global domain cell-size is fourarcmin if mesh_refine=1.0
real(dp), parameter :: fourarcmin = 4.0_dp/60.0_dp
! The following distance is useful for approximate timestepping calcs.
real(dp), parameter :: fourarcmin_dist = 7421.29938622_dp

! For command-line info
character(len=charlen) :: stage_file, model_name, run_type

!
! Ensure MPI is initialised
!
call swals_mpi_init

#ifndef SPHERICAL
write(log_output_unit,*) &
    'Code assumes spherical coordinates, but SPHERICAL is not defined'
call generic_stop
#endif

!
! Basic definition of multidomain
!
call program_timer%timer_start('startup_define_multidomain_geometry')

call parse_commandline_args(stage_file, &
    run_type, &
    final_time, &
    model_name, &
    coarse_nonlinear_domain_extents_file, &
    fine_nonlinear_domain_extents_file, &
    highres_nonlinear_domain_extents_file, &
    extrahighres_nonlinear_domain_extents_file, &
    md%load_balance_file, &
    md%output_basedir)

! Files with domain metadata for each nesting level (including fake file for
! global domain)
nesting_domain_extents_file = [character(len=charlen) :: &
    '', & ! Global domain --> empty filename
    coarse_nonlinear_domain_extents_file, &
    fine_nonlinear_domain_extents_file, &
    highres_nonlinear_domain_extents_file, &
    extrahighres_nonlinear_domain_extents_file]

! Initialise vars such that we get errors if they are not set later
num_nests(:) = 0_ip
domain_index_to_nest_with(:) = -1_ip

! Number domains in global domain
num_nests(1) = 1_ip
! Read domain metadata for nested domains, and record how many at each nesting
! level
do j = 2, NNL ! Avoid j=1 (global domain, treated separately)
    if(nesting_domain_extents_file(j) /= '') then
        call read_csv_into_array(nesting_domain_extents(j)%metadata, &
            nesting_domain_extents_file(j), skip_header=1_ip)
        num_nests(j) = size(nesting_domain_extents(j)%metadata, 2)
    end if
end do
if(.not. all(num_nests > 0) ) call generic_stop

! Figure out how many domains are needed. Later they may be partitioned,
! depending on md%load_balance_file.
nd = sum(num_nests)
allocate(md%domains(nd))
write(log_output_unit,*) "Number of domains initially created: ", nd
write(log_output_unit,*) "    Parallel partitioning may increase this number"

!
! Setup global domain
!

! Geometry
md%domains(1)%lw = global_lw
md%domains(1)%lower_left = global_ll
md%domains(1)%dx = fourarcmin * [1.0_dp, 1.0_dp] / mesh_refine
md%domains(1)%nx = nint(md%domains(1)%lw/md%domains(1)%dx)
md%domains(1)%dx_refinement_factor = dx_refinement_factors(1)
md%domains(1)%timestepping_refinement_factor = 1_ip

! Timestepping method
md%domains(1)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
md%domains(1)%linear_solver_is_truely_linear = .true.

! Increase to decimate output grids
md%domains(1)%nc_grid_output%spatial_stride = 1_ip

if(global_lw(1) == 360.0_dp) then
    ! Periodic EW boundary condition
    md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
end if

! Used below to define the parent domain for any domain at the next finer
! nesting level.
domain_index_to_nest_with(1) = 1_ip

! Used below to keep track of total number of domains
last_setup_domain = 1_ip

! Read metadata for each level of nesting
do j = 2, NNL ! Avoid j==1 (global domain, treated elsewhere)
    call setup_nesting_level_domains(j, domain_index_to_nest_with(j-1), &
        last_setup_domain, nesting_domain_extents(j)%metadata)

    ! Find a domain from nesting level i that is OK as a parent domain for
    ! nesting level i+1. Avoid domains that have been coarsened
    rowloop: do i = 1, num_nests(j)
        if(nesting_domain_extents(j)%metadata(r_coarse,i) == 0.0_dp) then
            domain_index_to_nest_with(j) = sum(num_nests(1:(j-1))) + i
            exit rowloop
        end if
    end do rowloop
end do

!
! Minor adjustments to domains

!
do j = 1, size(md%domains)

    ! Store stage/UH/VH grids over time, but not ELEVATION (it will be stored
    ! once anyway).
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
    if(rise_time > 0.0_dp .or. (.not. stage_file_is_raster) ) &
        md%domains(j)%support_elevation_forcing=.true.
    ! One should be careful doing this because it complicates the interpretation
    ! of max-stage -- in areas with subsidence, the max-stage may end up being
    ! the pre-subsidence topography. That could be worked-around by resetting
    ! domain%max_U in the loop, but for now we don't do anything like that.

    if( grid_output_style == 'animate_specific_regions' ) then
        ! For some grid_output_style values we store grids very often at some
        ! sites, and coarsely elsewhere
        call animate_specific_regions
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
    if(smooth_elevation_along_nesting_boundaries) then
        call set_initial_conditions(md%domains(j), stage_file, global_dt, md%all_dx_md)
    else
        call set_initial_conditions(md%domains(j), stage_file, global_dt)
    end if
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
call md%set_point_gauges_from_csv("../gauges/point_gauges_2022_12_14.csv", &
    skip_header=1_ip)
call program_timer%timer_end('startup_set_gauges')

!
! Final setup work.
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
        ', ts: ', stationary_ts, &
        ', min_ts:', min_ts, &
        merge(', UNSTABLE ', ',   stable ', (stationary_ts < min_ts)), &
        ', min(elevation):', minval(md%domains(j)%U(:,:,ELV)),&
        ', centre:', md%domains(j)%lower_left + md%domains(j)%lw/2, &
        ', dx:', md%domains(j)%dx

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

    call program_timer%timer_start('IO')

    ! Print and write outputs
    call md%write_outputs_and_print_statistics(&
        ! Time between writes is ~= "approximate_writeout_timestep"
        approximate_writeout_frequency=approximate_writeout_timestep, &
        write_grids_less_often = write_grids_every_nth_step, &
        write_gauges_less_often = write_gauges_every_nth_step, &
        print_less_often = print_every_nth_step, & !1_ip,&
        ! If the time is within timing_tol of the desired output time, then
        ! the output will be written even if the time < desired_output_time.
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

! Close files, and print timing info required for load balancing
call md%finalise_and_print_timers

write(log_output_unit,*) ''
write(log_output_unit, *) 'Program timer'
write(log_output_unit, *) ''
call program_timer%print(log_output_unit)
flush(log_output_unit)

call swals_mpi_finalize


contains

! Workhorse subroutine to setup each level of nesting
subroutine setup_nesting_level_domains(NL, domain_index_to_match, &
    last_setup_domain, nesting_NL_domain_extents)

    integer(ip), intent(in) :: NL, & !! Nesting level
        domain_index_to_match !! Index of md%domains(:) aligned with parent domain
    integer(ip), intent(inout) :: last_setup_domain 
        !! Last index of md%domains(:) that was setup
    real(dp), intent(in) :: nesting_NL_domain_extents(:,:) 
        !! Metadata for nesting level NL

    integer(ip) :: j, di, local_dx_refinement_factor, local_timestepping_refinement
    real(dp) :: local_max_depth, local_approx_wave_speed, local_coslat, &
        local_approx_grid_size, local_approx_ts

    do j = 1, size(nesting_NL_domain_extents, 2)

        di = last_setup_domain + j

        ! Get information on depth and cell-dimensions (for approx timestepping calculation).
        local_max_depth = -nesting_NL_domain_extents(r_dpth,j)
        local_approx_wave_speed = sqrt(gravity * max(local_max_depth*1.25_dp, 100.0_dp))

        local_dx_refinement_factor = merge(dx_refinement_factors(NL), &
            (dx_refinement_factors(NL)/coarsening_factors(NL)), &
            nesting_NL_domain_extents(r_coarse,j) == 0.0_dp)

        if(local_dx_refinement_factor < 0) then
            write(log_output_unit, *) "Error: negative local_dx_refinement_factor", &
                NL, j, di, local_dx_refinement_factor, nesting_NL_domain_extents(7,j)
            call generic_stop
        end if

        local_coslat = cos(nesting_NL_domain_extents(r_ll(2),j)/180.0_dp * pi)
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
            lower_left  = nesting_NL_domain_extents(r_ll,j), &
            upper_right = nesting_NL_domain_extents(r_ur,j), &
            dx_refinement_factor = local_dx_refinement_factor, &
            timestepping_refinement_factor = local_timestepping_refinement, &
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(di)%timestepping_method = 'rk2'
        md%domains(di)%nc_grid_output%spatial_stride = 1
        if(run_type == 'full') md%domains(di)%static_before_time = seconds_before_evolve

        ! Check the impact of a small eddy-viscosity
        !md%domains(di)%use_eddy_viscosity = .true.
        !md%domains(di)%eddy_visc_constants = [0.5_dp, 0.0_dp]

    end do

    last_setup_domain = di

end subroutine

subroutine animate_specific_regions
    ! Do whatever is needed to setup the animation.

    write(log_output_unit, *) &
        "Error: Have not yet specified specific regions to animate.", &
        "Example code commented below"
    call generic_stop

    ! Store stage grids over time, but not UH/VH/ELEVATION
    !md%domains(j)%time_grids_to_store = ['stage']

    ! ! By default store grids very coarsely (adjusted for some grids
    ! ! below)
    ! md%domains(j)%nc_grid_output%spatial_stride = 100_ip

    ! is_highres_domain = (j > sum(num_nests(1:2))
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
end subroutine

end program
