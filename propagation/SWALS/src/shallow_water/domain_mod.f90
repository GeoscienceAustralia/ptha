! Compile with -DTIMER to add timing to the code 
#ifdef TIMER
#   define TIMER_START(tname) call domain%timer%timer_start(tname)
#   define TIMER_STOP(tname)  call domain%timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

module domain_mod
    !
    ! Contains type 'domain_type', which is the main type for solving the
    ! linear/non-linear shallow water equations on a single structured grid.
    !
    ! Takes care of all allocation/file-IO, time-stepping, etc.
    !

    use global_mod, only: dp, ip, charlen, output_precision, &
                          cfl, maximum_timestep, gravity, &
                          advection_beta, &
                          minimum_allowed_depth, &
                          default_timestepping_method, &
                          wall_elevation, &
                          default_output_folder, &
                          send_boundary_flux_data,&
                          force_double, long_long_ip
    use timer_mod, only: timer_type
    use point_gauge_mod, only: point_gauge_type
    use coarray_utilities_mod, only: partitioned_domain_nesw_comms_type
    use nested_grid_comms_mod, only: domain_nesting_type
    use stop_mod, only: generic_stop
    use iso_fortran_env, only: output_unit, int32, int64
    use netcdf_util, only: nc_grid_output_type
    use logging_mod, only: log_output_unit
    use file_io_mod, only: mkdir_p

#ifdef SPHERICAL    
    ! Compile with -DSPHERICAL to get the code to run in spherical coordinates
    use spherical_mod, only: area_lonlat_rectangle, deg2rad, earth_angular_freq
    use global_mod, only: radius_earth
    !
    ! Can only have coriolis if we have spherical. However, the code
    ! will still work if we have spherical and not coriolis. But by default,
    ! we include coriolis with spherical. Add the flag '-DNOCORIOLIS' when compiling
    ! to remove it
    !
#ifndef NOCORIOLIS
#define CORIOLIS
#endif

#endif

#ifndef NOOPENMP
    use omp_lib
#endif

    implicit none

    ! Make everything private, except domain_type, which has its own methods,
    ! and the domain_metadata_type, which is sometimes more convenient to
    ! work with than the domain [since it is 'lightweight']
    private
    public:: domain_type, test_domain_mod

    ! Indices for arrays: Stage, depth-integrated-x-velocity,
    ! depth-integrated-v-velocity, elevation. So e.g. stage
    ! is in domain%U(:,:,STG), and elevation is in domain%U(:,:ELV)
    integer(int32), parameter, public:: STG=1, UH=2, VH=3, ELV=4

    ! Handy constants
    real(dp), parameter :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp
    real(dp), parameter :: QUARTER_dp = HALF_dp * HALF_dp
    real(dp), parameter:: NEG_SEVEN_ON_THREE_dp = -2.0_dp - 1.0_dp/3.0_dp !-7.0_dp/3.0_dp

    ! Use this formatting when converting domain%myid to character
    character(len=charlen), parameter:: domain_myid_char_format = '(I0.20)'

    ! Use this in explicitly vectorized domain routines
    ! A number of methods have been coded once with loops, and once with an
    ! explicitly vectorized alternatives. The latter is sometimes faster, but not always.
    integer, parameter :: vectorization_size = 32

    !
    ! Coefficients controlling the edge extrapolation limiter. 
    !
    ! We control edge extrapolation based on the following quantity:
    !   Q = (min_neighbour_depth - minimum_allowed_depth)/(max_neighbour_depth + limiter_coef3*minimum_allowed_depth)
    ! The idea is that if Q is 'not small', then depth gradients are 'not very
    !   big', and we can use second order extrapolation.
    ! If Q is 'small', then depths are changing rapidly in space, and we
    !   should decay to first-order extrapolation for wet-dry stability.
    ! If the depths are 'very shallow' (locally less than
    !  limiter_coef3*minimum_allowed_depth), then we should also promote decay
    !  to first order to assist stability.
    ! Basically, we 'linearly change' the allowed extrapolation theta, as we
    !   move between 'not small' and 'small' values of Q (see the routines below
    !   where these parameters are used).
    !
    ! We use 'limiter_coef1' for the 'small' threshold, and limiter_coef2 for the 'not small' threshold.
    real(dp), parameter :: limiter_coef1 = 0.1_dp, limiter_coef2 = 0.3_dp, limiter_coef3 = 100.0_dp
    !real(dp), parameter :: limiter_coef1 = 0.00_dp, limiter_coef2 = 0.15_dp, limiter_coef3 = 10.0_dp
    !real(dp), parameter :: limiter_coef1 = 0.01_dp, limiter_coef2 = 0.05_dp, limiter_coef3 = 10.0_dp
    !! This one turns off limiting in practice!
    !real(dp), parameter :: limiter_coef1 = -10000.03_dp, limiter_coef2 = 0.15_dp, limiter_coef3 = 10.0_dp
    ! The following value often occurs in the calculations
    real(dp), parameter :: limiter_coef4 = ONE_dp/(limiter_coef2 - limiter_coef1)

    !
    ! Main type which holds a single grid and associated metadata
    !
    type :: domain_type

        ! [Length,width] of domain in x/y units. 
        real(dp):: lw(2) = -HUGE(1.0_dp)
        ! grid size [number of x-cells, number of y-cells]
        integer(ip):: nx(2) = -1_ip
        ! cell size [dx,dy]
        real(dp):: dx(2) = - HUGE(1.0_dp)
        ! Absolute lower-left coordinate (bottom left corner of cell(1,1))
        real(dp):: lower_left(2) = HUGE(1.0_dp)

        ! Parameter controlling extrapolation for finite volume methods
        ! Defaults depend on timestepping_method
        real(dp) :: theta = -HUGE(1.0_dp) !
        ! CFL number
        real(dp):: cfl = -HUGE(1.0_dp)

        ! This information is useful in the context of nesting, where interior
        ! domains have cell sizes that are an integer divisor of the
        ! coarse-domain dx. Likewise one domain might take integer
        real(dp) :: dx_refinement_factor = ONE_dp
        real(dp) :: dx_refinement_factor_of_parent_domain = ONE_dp
        integer(ip) :: timestepping_refinement_factor = 1_ip

        ! The width of the nesting layer for this domain -- which will depend
        ! on the dx value relative to the neighbouring domains, on the timestepping_method,
        ! and on the timestepping_refinement_factor. See 'get_domain_nesting_layer_thickness'
        integer(ip) :: nest_layer_width = -1_ip

        ! Useful to hold the interior bounding box [ = originally provided bounding box]
        ! since the actual bounding box might be changed to accommodate nesting 
        real(dp) :: interior_bounding_box(4,2) = 0.0_dp
        
        ! Domain ID, which is useful if multiple domains are running
        integer(int64):: myid = 1
        integer(int64):: local_index = 1 ! Useful when we partition a domain in parallel
        !character(len=charlen) :: myid_char = '000001'

        ! Flag to denote boundaries at which nesting occurs: order is N, E, S, W.
        logical :: is_nesting_boundary(4) = .FALSE. 

        ! timestepping_method determines the choice of solver
        character(len=charlen):: timestepping_method = default_timestepping_method

        real(dp) :: max_parent_dx_ratio

        ! Number of quantities in domain%U (stage, uh, vh, elevation)
        integer(ip):: nvar = 4 

        ! Name of ascii file where we output metadata
        character(len=charlen):: metadata_ascii_filename

        ! Subroutine called inside domain%compute_fluxes
        character(len=charlen) :: compute_fluxes_inner_method = 'DE1_low_fr_diffusion' ! 'DE1'

        ! The domain 'interior' is surrounded by 'exterior' cells which are
        ! updated by boundary conditions, or copied from other domains. When
        ! tracking mass conservation, we only want to record inflows/outflows to
        ! interior cells. The interior cells are:
        !    [(1+exterior_cells_width):(domain%nx(1)-exterior_cells_width), &
        !     (1+exterior_cells_width):(domain%nx(2)-exterior_cells_width)]
        ! NOTE: exterior_cells_width should only be used to influence mass conservation tracking calculations.
        ! It is not strictly related to the halo width for parallel computations. When using a multidomain, the
        ! mass conservation tracking is somewhat different (based around subroutines with names matching 'nesting_boundary_flux_*')
        integer(ip):: exterior_cells_width = 2
    
        ! Count number of time steps (useful if we want to do something only
        ! every n'th timestep)
        integer(ip):: nsteps_advanced = 0
       
        ! The evolved time in seconds
        real(dp):: time = ZERO_dp
        ! dt computed from CFL condition (e.g. during flux computation). We might take a smaller timestep (e.g. to match the
        ! timestep on other nested domains)
        real(dp):: max_dt = ZERO_dp
        ! Record dt the last time we ran update_U. Useful if we want to append some other update.
        real(dp) :: dt_last_update = ZERO_dp
        ! Record the time-step evolved by domain%evolve_one_step.
        real(dp):: evolve_step_dt = ZERO_dp
        ! This can set the maximum allowed timestep [used to prevent high timesteps on dry domains]
        real(dp):: maximum_timestep = maximum_timestep
        ! If the time is less than this number, then we will assume the domain flow is stationary, and will not re-compute the flow.
        ! Useful to increase the efficiency of domains where you know nothing happens until some set time (e.g. far-field tsunami ).
        real(dp) :: static_before_time = -HUGE(1.0_dp)
        
        ! Parameter which determines how the 'static depth' is computed in
        ! linear solver [i.e. corresponding to 'h0' in the pressure gradient term 'g h_0 d(free_surface)/dx']
        real(dp):: msl_linear = 0.0_dp
        ! This flag controls whether, for the 'linear' solver, we allow the
        ! pressure-gradient term (g h0 dStage/dx) to have 'h0' varying
        ! over time. If linear_solver_is_truely_linear, then 'h0' is constant.
        ! Otherwise 'h0' varies as the stage varies (which actually makes the equations nonlinear)
        logical :: linear_solver_is_truely_linear = .true.
        ! Useful variable to distinguish staggered-grid and centred-grid numerical methods
        logical :: is_staggered_grid 

        !
        ! Boundary conditions. 
        !
        ! Currently 'boundary_subroutine' is the most important way of imposing
        !    boundary conditions. It must take a domain_type as an INTENT(INOUT) argument
        ! Applications need to define boundary_subroutine, either using a
        !    pre-existing bc, or making a new subroutine.
        ! Some existing boundary conditions rely on 'boundary_function' as well. It takes
        !    the domain_type as well as time, and the location (defined as i, j indices), see the interface below.
        !    Note that the x/y location can easily be obtained from domain%x(i), domain%y(j)
        character(len=charlen):: boundary_type = ''
        procedure(boundary_fun), pointer, nopass:: boundary_function => NULL()
        procedure(boundary_subroutine), pointer, nopass:: boundary_subroutine => NULL() 
        !
        ! Flag whether the boundary is exterior (TRUE) or interior (FALSE). We only need
        ! to apply a boundary condition to the 'exterior boundaries' -- otherwise we have e.g. nesting updates, 
        ! which are different.
        ! Order is North (1), East (2), South (3), West (4) 
        logical :: boundary_exterior(4) = .TRUE. 

        ! 
        ! Mass conservation tracking  -- store as double, even if dp is single prec.
        !
        ! Store the flux through the N, E, S, W boundaries, single domain 
        real(force_double):: boundary_flux_store(4)  = ZERO_dp
        ! Alternative to the above which is appropriate for the nested-grid case
        real(force_double):: boundary_flux_store_exterior(4)  = ZERO_dp
        ! Time integrate the boundary fluxes, single domain
        real(force_double):: boundary_flux_time_integral = ZERO_dp
        ! Alternative to the above which is appropriate for the nested-grid case
        real(force_double):: boundary_flux_time_integral_exterior = ZERO_dp
        ! We need an intermediate variable to take care of time-stepping
        ! This integrates the boundary fluxes within the evolve step only, single domain
        real(force_double):: boundary_flux_evolve_integral = ZERO_dp
        ! Alternative to the above which is appropriate for the nested-grid case
        real(force_double):: boundary_flux_evolve_integral_exterior = ZERO_dp
        
        ! Lower/upper x and y indices over which the SWE computation takes place
        ! These might restrict the SWE update to a fraction of the domain (e.g.
        ! beginning of an earthquake-tsunami run where only a fraction of the domain is
        ! active.
        ! Currently implemented for leap-frog solvers only
        integer(ip):: xL, xU, yL, yU

        ! Output folder units
        ! Unit for time (stored as ascii)
        integer(ip):: output_time_unit_number
        ! Base folder inside which we store domain outputs.
        character(len=charlen):: output_basedir = default_output_folder
        ! Output folder (it will begin with domain%output_basedir, and include another timestamped folder) 
        character(len=charlen):: output_folder_name = ''
        ! Unit number for domain log file
        integer(ip):: logfile_unit = output_unit

        ! Units for home-brew binary output format [better to use netcdf nowadays].
        ! This is only used when the code is compiled with "-DNONETCDF"
        integer(ip), allocatable :: output_variable_unit_number(:)

        ! Count the number of time-steps at which we clip depths.
        ! For a 'well behaved' model this should remain zero -- non-zero values indicate mass conservation issues.
        integer(long_long_ip) :: negative_depth_fix_counter = 0

        ! Type to manage netcdf grid outputs
        type(nc_grid_output_type) :: nc_grid_output

        ! Type to record CPU timings
        type(timer_type):: timer

        ! Type to do single-grid coarray communication. This has been superceeded by multidomain. 
        type(partitioned_domain_nesw_comms_type):: partitioned_comms
        ! Determine whether we call domain%partitioned_comms%communicate. It would be better
        ! to hide this inside the partitioned_comms class, but ifort 2019 segfaults when that is done.
        logical :: use_partitioned_comms = .false. 

        ! Variables controlling the storage of maximum-stage.
        ! We don't have to store the max_U. For some problems (linear solver) that can
        ! take a significant fraction of the total time, or use too much memory.
        logical:: record_max_U = .true.
        integer(ip):: max_U_update_frequency = 1
        
        ! Spatial coordinates, dx/dy distances (useful for spherical coordinates)
        real(dp), allocatable :: x(:), y(:), distance_bottom_edge(:), distance_left_edge(:)
        real(dp), allocatable :: area_cell_y(:)

        ! These variables are only used with spherical coordinates
        real(dp), allocatable :: coslat(:), coslat_bottom_edge(:), tanlat_on_radius_earth(:)
        ! These variables are only used with coriolis
        real(dp), allocatable :: coriolis(:), coriolis_bottom_edge(:)

        ! Type to manage storing of tide gauges
        type(point_gauge_type) :: point_gauges

        !
        ! Type to manage nesting communication. This is required for nesting, and is the standard
        ! approach for distributed-memory parallel runs.
        !
        type(domain_nesting_type) :: nesting
 

        ! Big arrays to hold the domain variables
        !
        ! U holds main flow variables -- STG, UH, VH, ELV 
        ! First 2 dimensions = space, 3rd = number of quantities
        real(dp), allocatable :: U(:,:,:) ! Needed always
        real(dp), allocatable :: max_U(:,:,:) ! Needed if max_U is to be output. FIXME: Consider storing in output_precision
        !
        ! Multi-dimensional arrays below are only required for nonlinear solver. Would be possible
        ! to further reduce memory usage, but not completely trivial. Not all of these are required
        ! for linear type timestepping
        real(dp), allocatable :: flux_NS(:,:,:)  
        real(dp), allocatable :: flux_EW(:,:,:)  
        real(dp), allocatable :: depth(:,:)  
        real(dp), allocatable :: velocity(:,:, :) 
        real(dp), allocatable :: explicit_source(:,:,:) 
        real(dp), allocatable :: explicit_source_VH_j_minus_1(:,:) ! Separate from explicit_source for OPENMP parallel logic
        real(dp), allocatable :: manning_squared(:,:) 
        real(dp), allocatable :: backup_U(:,:,:) ! Needed for some timestepping methods
        real(dp), allocatable :: friction_work(:,:,:) ! Only used for leapfrog_linear_plus_nonlinear_friction
        logical :: friction_work_is_setup = .false. ! Flag used for efficiency with leapfrog_linear_plus_nonlinear_friction solvers
#ifdef DEBUG_ARRAY
        ! For debugging it can be helpful to have this array.  If DEBUG_ARRAY is defined, then it will be allocated with dimensions
        ! (nx, ny), and will be written to the netcdf file at each time.
        real(dp), allocatable :: debug_array(:,:)
#endif

        CONTAINS

        ! Initialisation
        procedure:: allocate_quantities => allocate_quantities

        ! Reporting
        procedure:: print => print_domain_statistics

        ! Core routines that occur within a timestep
        ! (consider making these not type bound -- since the user should not
        !  really call them)
        procedure:: compute_depth_and_velocity => compute_depth_and_velocity
        procedure:: compute_fluxes => compute_fluxes
        !procedure:: compute_fluxes => compute_fluxes_vectorized !! Slower than un-vectorized version on GD home machine
        procedure:: update_U => update_U  ! Slower or faster, depending on wet/dry areas
        !procedure:: update_U => update_U_restructured  ! Slower or faster, depending on wet/dry areas
        !procedure:: update_U => update_U_vectorized !! Slower on an NCI test
        procedure:: backup_quantities => backup_quantities

        ! Setup
        procedure:: precompute_friction_work => precompute_friction_work
        
        ! Timestepping (consider making only 'evolve_one_step' type bound -- since the user should not really call others)
        procedure:: one_euler_step => one_euler_step
        procedure:: one_rk2_step => one_rk2_step 
        procedure:: one_rk2n_step => one_rk2n_step 
        procedure:: one_midpoint_step => one_midpoint_step
        procedure:: one_linear_leapfrog_step => one_linear_leapfrog_step
        procedure:: one_leapfrog_linear_plus_nonlinear_friction_step => one_leapfrog_linear_plus_nonlinear_friction_step
        procedure:: evolve_one_step => evolve_one_step
        procedure:: update_max_quantities => update_max_quantities

        ! Boundary conditions. This just calls whatever domain%boundary_subroutine points to (consider making not type bound)
        procedure:: update_boundary => update_boundary

        ! IO
        procedure:: create_output_files => create_output_files
        procedure:: write_to_output_files => write_to_output_files
        procedure:: write_max_quantities => write_max_quantities
        procedure:: log_outputs => divert_logfile_unit_to_file

        ! Mass conservation reporting
        procedure:: mass_balance_interior => mass_balance_interior
        procedure:: volume_interior => volume_interior

        ! Functions to estimate gravity-wave time-step (often useful to call at the start of a model run). These are not used
        ! in the flow-algorithms, they are just for convenience.
        procedure:: linear_timestep_max => linear_timestep_max
        procedure:: nonlinear_stationary_timestep_max => nonlinear_stationary_timestep_max
        procedure:: stationary_timestep_max => stationary_timestep_max

        ! Setup and store point-gauge outputs
        procedure:: setup_point_gauges => setup_point_gauges
        procedure:: write_gauge_time_series => write_gauge_time_series

        ! finalization (e.g. close netcdf files)
        procedure:: finalise => finalise_domain

        ! Smoothing of elevation. Not recommended in general, but with poor elevation data it might reduce artefacts.
        procedure :: smooth_elevation => smooth_elevation

        ! Nesting.
        ! -- Keep track of fluxes for mass-conservation tracking
        procedure:: nesting_boundary_flux_integral_multiply => nesting_boundary_flux_integral_multiply
        procedure:: nesting_boundary_flux_integral_tstep => nesting_boundary_flux_integral_tstep
        procedure:: nesting_flux_correction_everywhere => nesting_flux_correction_everywhere
        ! -- Determine if a given point is in the 'priority domain'
        procedure:: is_in_priority_domain => is_in_priority_domain
        ! -- If the current grid communicates to a coarser grid, this routine can make elevation constant within each coarse-grid
        ! cell in the send-region. Useful to make elevation consistent in a nesting situation
        procedure:: use_constant_wetdry_send_elevation => use_constant_wetdry_send_elevation
        ! -- When initialising nested domains, set lower-left/upper-right/resolution near the desired locations,
        ! adjusting as required to support nesting (e.g. so that edges align with parent domain, and cell size is an integer divisor
        ! of the parent domain, etc).
        procedure:: match_geometry_to_parent => match_geometry_to_parent


    end type domain_type

    interface

        ! This gives the interface for a generic 'boundary function' which
        ! returns an array of 4 output values (stage/uh/vh/elevation)
        !
        ! It is used in conjunction with a number of different types of
        ! boundary conditions below
        !
        ! @param domain the domain
        ! @param t time 
        ! @param i,j the i/j index we are operating on (presumed on the boundary)
        !
        function boundary_fun(domain, t, i, j) RESULT(stage_uh_vh_elev)
            import dp, ip, domain_type
            implicit none
            type(domain_type), intent(in):: domain 
            real(dp), intent(in):: t
            integer(ip), intent(in) :: i, j
            real(dp):: stage_uh_vh_elev(4)
        end function

        !
        ! The user can provide a boundary subroutine which is supposed to update
        ! the domain boundaries. It is called by domain%update_boundary whenever
        ! a boundary update is required by the timestepping_method. Note that
        ! this may well mean 'boundary_fun' is not required
        !
        subroutine boundary_subroutine(domain)
            import domain_type
            implicit none
            type(domain_type), intent(inout):: domain
        end subroutine
    
    end interface

    contains
  
    ! 
    ! Convenience printing function 
    !
    ! FIXME: For nesting domains, consider modifying this to only use values
    !        where the priority_domain corresponds to the current domain.
    !        In general, I use the multidomain printing routines.
    subroutine print_domain_statistics(domain)
        class(domain_type), intent(inout):: domain
        real(dp):: maxstage, minstage
        integer:: i,j, ecw, nx, ny
        real(dp):: dry_depth_threshold, energy_total, energy_potential, energy_kinetic
        real(dp):: depth, depth_iplus, depth_jplus
        logical, parameter:: report_energy_statistics=.TRUE.

TIMER_START('printing_stats')

        dry_depth_threshold = minimum_allowed_depth

        nx = domain%nx(1)
        ny = domain%nx(2)
            
        ! Min/max stage in wet areas
        maxstage = -huge(ONE_dp)
        minstage = huge(ONE_dp)
        do j = 2, ny-1
            do i = 2, nx-1
                if(domain%U(i,j,STG) > domain%U(i,j,ELV) + dry_depth_threshold) then
                    maxstage = max(maxstage, domain%U(i,j,STG))
                    minstage = min(minstage, domain%U(i,j,STG))
                end if
            end do
        end do

        ! Print main statistics
        write(domain%logfile_unit, *) ''
        write(domain%logfile_unit, *) 'Domain ID: '
        write(domain%logfile_unit, *) '        ', domain%myid
        write(domain%logfile_unit, *) 'Time: '
        write(domain%logfile_unit, *) '        ', domain%time
        write(domain%logfile_unit, *) 'nsteps_advanced:'
        write(domain%logfile_unit, *) '        ', domain%nsteps_advanced
        write(domain%logfile_unit, *) 'max_allowed_dt: '
        write(domain%logfile_unit, *) '        ', domain%max_dt
        write(domain%logfile_unit, *) 'evolve_step_dt: '
        write(domain%logfile_unit, *) '        ', domain%evolve_step_dt
        write(domain%logfile_unit, *) 'Stage: '
        write(domain%logfile_unit, *) '        ', maxstage
        write(domain%logfile_unit, *) '        ', minstage

        if((.not. domain%is_staggered_grid)) then

            ! Typically depth/velocity are not up-to-date with domain%U, so fix
            ! that here.
            call domain%compute_depth_and_velocity()

            write(domain%logfile_unit, *) 'u: '
            write(domain%logfile_unit, *) '        ', maxval(domain%velocity(2:(nx-1), 2:(ny-1),UH))
            write(domain%logfile_unit, *) '        ', minval(domain%velocity(2:(nx-1), 2:(ny-1),UH))
            write(domain%logfile_unit, *) 'v: '
            write(domain%logfile_unit, *) '        ', maxval(domain%velocity(2:(nx-1), 2:(ny-1),VH))
            write(domain%logfile_unit, *) '        ', minval(domain%velocity(2:(nx-1), 2:(ny-1),VH))
            write(domain%logfile_unit, *) '.........'

        end if

        ! Optionally report integrated energy
        if(report_energy_statistics) then

            ecw = domain%exterior_cells_width
            energy_potential = ZERO_dp
            energy_kinetic = ZERO_dp

            if(domain%is_staggered_grid .and. domain%linear_solver_is_truely_linear) then
                !
                ! Linear leap-frog scheme doesn't store velocity
                ! Get total, potential and kinetic energy in domain interior
                !
                ! This does not exactly lead to energy conservation, but with small-enough time-steps
                ! it should be good for the linear solver.
                !

                ! Integrate over the model interior
                do j = 1 + ecw, (domain%nx(2) - ecw)
                    do i = (1+ecw), (domain%nx(1)-ecw)

                        ! Integrate over wet cells
                        if(domain%msl_linear > domain%U(i,j,ELV) + dry_depth_threshold) then
            
                            energy_potential = energy_potential + domain%area_cell_y(j) *&
                                (domain%U(i,j,STG) - domain%msl_linear)**2

                            !
                            ! Kinetic energy integration needs to account for grid staggering -- cell
                            ! areas are constant for constant lat, but changing as lat changes.
                            !
                            depth_iplus = HALF_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + &
                                domain%U(i+1,j,STG) - domain%U(i+1,j,ELV))
                            if(depth_iplus > dry_depth_threshold) then
                                energy_kinetic = energy_kinetic + domain%area_cell_y(j)*&
                                    (domain%U(i,j,UH)**2)/depth_iplus
                            end if

                            depth_jplus = HALF_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + &
                                domain%U(i,j+1,STG) - domain%U(i,j+1,ELV))
                            if(depth_jplus > dry_depth_threshold) then
                                energy_kinetic = energy_kinetic + HALF_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1))*&
                                    (domain%U(i,j,VH)**2)/depth_jplus
                            end if

                        end if
                    end do
                end do
            else if(.not. domain%is_staggered_grid) then
                !
                ! Compute energy statistics using non-staggered grid
                !

                ! Integrate over the model interior
                do j = 1 + ecw, (domain%nx(2) - ecw)
                    do i = (1+ecw), (domain%nx(1)-ecw)
                        depth = domain%U(i,j,STG) - domain%U(i,j,ELV)
                        if(depth > dry_depth_threshold) then
                            energy_potential = energy_potential + domain%area_cell_y(j) * &
                                (domain%U(i,j,STG) - domain%msl_linear)**2
                                !depth * domain%U(i,j,STG) 
                            energy_kinetic = energy_kinetic + domain%area_cell_y(j) * &
                                (domain%U(i,j,UH)**2 + domain%U(i,j,VH)**2)/depth
                        end if 
                    end do
                end do

            end if

            ! In the case we use a "not-truely-linear" variant of the linear solver, the energy
            ! calculations do not really make sense
            if(domain%is_staggered_grid .and. (.not. domain%linear_solver_is_truely_linear)) then
                ! Do nothing
            else
                ! Rescale energy statistics appropriately 
                energy_potential = energy_potential * gravity * HALF_dp
                energy_kinetic = energy_kinetic * HALF_dp
                energy_total = energy_potential + energy_kinetic

                write(domain%logfile_unit, *) 'Energy Total / rho: ', energy_total
                write(domain%logfile_unit, *) 'Energy Potential / rho: ', energy_potential
                write(domain%logfile_unit, *) 'Energy Kinetic / rho: ', energy_kinetic
            end if

        end if

        ! Mass conservation check
        write(domain%logfile_unit, *) 'Mass Balance (domain interior): ', domain%mass_balance_interior()

TIMER_STOP('printing_stats')
    end subroutine

    ! 
    ! Set up the full domain, allocate arrays, etc
    !
    ! @param global_lw real(dp) array size 2. length/width of domain in same units as x,y coordinates
    ! @param global_nx integer array size 2. number of x/y cells in the domain
    ! @param global_ll real(dp) array size 2. lower left x/y coordinate of the domain (at the corner of the lower left cell)
    ! @param create_output_files optional. If .TRUE. or not provided, then make output files
    ! @param co_size_xy optional integer array size 2. Split up domain into sub-tiles of this dimension [coarray only -- note this
    ! is a simple "single domain decomposition", not used with the multidomain approach]
    ! @param ew_periodic Use EW periodic boundaries [coarray only -- note this is a simple "single domain decomposition", not used
    ! with the multidomain approach]
    ! @param ns_periodic Use NS periodic boundaries [coarray only -- note this is a simple "single domain decomposition", not used
    ! with the multidomain approach]
    ! @param verbose Print info about the domain
    !
    subroutine allocate_quantities(domain, global_lw, global_nx, global_ll, create_output_files,&
        co_size_xy, ew_periodic, ns_periodic, verbose)

        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: global_lw(2), global_ll(2)
        integer(ip), intent(in):: global_nx(2)
        logical, optional, intent(in) :: create_output_files
        integer(ip), optional, intent(in):: co_size_xy(2)
        logical, optional, intent(in) :: ew_periodic, ns_periodic, verbose

        integer(ip) :: nx, ny, nvar
        integer(ip) :: i, j, k
        logical :: create_output, use_partitioned_comms, ew_periodic_, ns_periodic_, verbose_
        real(dp):: local_lw(2), local_ll(2)
        integer(ip):: local_nx(2)

        ! Send domain print statements to the default log (over-ridden later if we send the domain log to its own file
        domain%logfile_unit = log_output_unit

        if(present(create_output_files)) then
            create_output = create_output_files
        else
            create_output = .TRUE.
        end if

        if(present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .true.
        end if

        ! Only split the domain if we provided a co_size with > 1 image
        if (present(co_size_xy)) then
            if(maxval(co_size_xy) > 1) then
                use_partitioned_comms = .TRUE.
                domain%use_partitioned_comms = .TRUE.

                ! In parallel, we need to tell the code whether the domain is periodic or not
                if(present(ew_periodic)) then
                    ew_periodic_ = ew_periodic
                else
                    ew_periodic_ = .FALSE.
                end if

                if(present(ns_periodic)) then
                    ns_periodic_ = ns_periodic
                else
                    ns_periodic_ = .FALSE.
                end if

            else
                use_partitioned_comms = .FALSE.
            endif
        else
            use_partitioned_comms = .FALSE.
        endif


        ! Set good defaults for different timestepping methods
        select case(domain%timestepping_method)
        case('rk2')
            if(domain%theta == -HUGE(1.0_dp)) domain%theta = 1.6_dp
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.99_dp
            domain%is_staggered_grid = .false.
        case('rk2n')
            if(domain%theta == -HUGE(1.0_dp)) domain%theta = 1.6_dp
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.9_dp
            domain%is_staggered_grid = .false.
        case('midpoint')
            if(domain%theta == -HUGE(1.0_dp)) domain%theta = 1.6_dp
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.99_dp
            domain%is_staggered_grid = .false.
        case('euler')
            if(domain%theta == -HUGE(1.0_dp)) domain%theta = 0.9_dp
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.9_dp
            domain%is_staggered_grid = .false.
        case('linear')
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.7_dp
            domain%is_staggered_grid = .true.
        case('leapfrog_linear_plus_nonlinear_friction')
            if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = 0.7_dp
            domain%is_staggered_grid = .true.
        case default
            print*, 'domain%timestepping_method = ', trim(domain%timestepping_method), ' not recognized'
            call generic_stop()
        end select


        if(use_partitioned_comms) then
            ! Note that this only supports a single grid, and has largely been replaced by the "multidomain"
            ! parallel infrastructure.

            ! Compute the ll/lw/nx for this sub-domain
            call domain%partitioned_comms%initialise(co_size_xy, global_ll, global_lw, global_nx, &
                local_ll, local_lw, local_nx, &
                ew_periodic=ew_periodic_, ns_periodic=ns_periodic_)
            domain%lower_left = local_ll
            domain%lw = local_lw 
            domain%nx = local_nx

            call domain%partitioned_comms%print()
        
            ! Make sure that communication boundaries are not numerical boundaries
            do i = 1,4
                if(domain%partitioned_comms%neighbour_images(i) > 0) domain%boundary_exterior(i) = .FALSE.
            end do
        else
            ! Not using coarrays (simple case)
            domain%lw = global_lw
            domain%nx = global_nx
            domain%lower_left = global_ll
        end if

        ! For the linear solver, these variables can be modified to evolve domain%U only in the region domain%U(xL:xU, yL:yU). This
        ! can speed up some simulations (but the user must adaptively change the variables).
        domain%xL = 1_ip
        domain%yL = 1_ip
        domain%xU = domain%nx(1)
        domain%yU = domain%nx(2)

        ! Compute dx
        domain%dx = domain%lw/(domain%nx*ONE_dp)

        ! Compute the x/y cell center coordinates. We assume a linear mapping between the x-index and the x-coordinate, and
        ! similarly for y. We support regular cartesian coordinates (m), and lon/lat spherical coordinates (degrees).
        nx = domain%nx(1)
        ny = domain%nx(2)
        nvar = domain%nvar 
        allocate(domain%x(nx), domain%y(ny))
        ! x/y coordinates (only stored along domain edge, assumed constant)
        do i = 1, nx
            domain%x(i) = domain%lower_left(1) + ((i - HALF_dp)/(nx*ONE_dp))*domain%lw(1)
        end do
        do i = 1, ny
            domain%y(i) = domain%lower_left(2) + ((i - HALF_dp)/(ny*ONE_dp))* domain%lw(2) 
        end do
     
#ifdef SPHERICAL
        ! For spherical coordinates it saves computation to have cos(latitude) at cells and edges
        allocate(domain%coslat(ny), domain%coslat_bottom_edge(ny+1), domain%tanlat_on_radius_earth(ny))
        domain%coslat = cos(domain%y * deg2rad)
        domain%coslat_bottom_edge(1:ny) = cos((domain%y - HALF_dp * domain%dx(2))*deg2rad)
        domain%coslat_bottom_edge(ny+1) = cos((domain%y(ny) + HALF_dp*domain%dx(2))*deg2rad)
        ! This term appears in 'extra' spherical shallow water equations terms from Williamson et al., (1992) and other modern
        ! derivations
        domain%tanlat_on_radius_earth = tan(domain%y * deg2rad) / radius_earth
#endif

#ifdef CORIOLIS
        ! Spherical coordinates must be used if coriolis is used.
#ifndef SPHERICAL
        stop 'Cannot define preprocessing flag CORIOLIS without also defining SPHERICAL'
#endif    
        ! Coriolis parameter
        allocate(domain%coriolis_bottom_edge(ny+1), domain%coriolis(ny))
        domain%coriolis = 2.0_dp * sin(domain%y * deg2rad) * earth_angular_freq
        domain%coriolis_bottom_edge(1:ny) = 2.0_dp * earth_angular_freq * &
            sin((domain%y - HALF_dp * domain%dx(2))*deg2rad)
        domain%coriolis_bottom_edge(ny+1) = 2.0_dp * earth_angular_freq * &
            sin((domain%y(ny) + HALF_dp * domain%dx(2))*deg2rad)

#endif


        ! Distances along edges. 
        ! For spherical lon/lat coordinates, distance_bottom_edge changes with y
        ! For cartesian x/y coordinates, they are constants
        ! In general we allow the bottom-edge distance to change with y, and the left-edge distance to change with x.
        allocate(domain%distance_bottom_edge(ny+1), domain%distance_left_edge(nx+1))

        do i = 1, nx + 1
#ifdef SPHERICAL
            ! Spherical coordinates
            domain%distance_left_edge(i) = domain%dx(2) * deg2rad * radius_earth
#else
            ! Cartesian 
            domain%distance_left_edge(i) = domain%dx(2)
#endif
        end do
        do i = 1, ny+1
#ifdef SPHERICAL
            ! Spherical coordinates
            domain%distance_bottom_edge(i) = domain%dx(1) * deg2rad * radius_earth * &
                domain%coslat_bottom_edge(i)
#else
            ! Cartesian
            domain%distance_bottom_edge(i) = domain%dx(1)
#endif
        end do
        
        ! Area of cells. For spherical, this only changes with y
        allocate(domain%area_cell_y(ny))
#ifdef SPHERICAL
        do i = 1, ny
            domain%area_cell_y(i) = area_lonlat_rectangle(&
                domain%x(1) - domain%dx(1) * HALF_dp, &
                domain%y(i) - domain%dx(2) * HALF_dp, &
                domain%dx(1), &
                domain%dx(2), &
                flat=.TRUE.)
        end do
#else
        domain%area_cell_y = product(domain%dx)
#endif

        !
        ! Below we allocate the main arrays, using loops to promote openmp memory affinity (based on the first-touch principle).
        !
        allocate(domain%U(nx, ny, 1:nvar))
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny, nvar)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, ny
            domain%U(:,j,1:nvar) = ZERO_dp
        end do
        !$OMP END DO
        !$OMP END PARALLEL
       
        if(domain%record_max_U) then 
            allocate(domain%max_U(nx, ny, 1))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%max_U(:,j,1) = -huge(ONE_dp)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
        end if

        ! Many other variables are required for the non-staggered-grid solvers (i.e. the nonlinear finite-volume solvers)
        if((.not. domain%is_staggered_grid)) then 

            allocate(domain%depth(nx, ny))
            allocate(domain%velocity(nx, ny, UH:VH))
            
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%depth(:,j) = ZERO_dp
                domain%velocity(:,j,UH:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            
            ! NOTE: We don't need a flux for elevation 
            allocate(domain%flux_NS(nx, ny+1_ip, STG:VH))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny+1
                domain%flux_NS(:,j,STG:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            allocate(domain%flux_EW(nx + 1_ip, ny, STG:VH))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%flux_EW(:, j, STG:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! Pressure gradient applies only to uh and vh (2 and 3 variables in U)
            allocate(domain%explicit_source(nx, ny, UH:VH))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%explicit_source(:,j,UH:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL
           
            ! This stores components of explicit_source related to the index 'j-1', which we must treat separately when using
            ! openmp.
            allocate(domain%explicit_source_VH_j_minus_1(nx, ny+1))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny+1
                domain%explicit_source_VH_j_minus_1(:,j) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! (Manning friction)*(Manning friction)
            allocate(domain%manning_squared(nx, ny))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%manning_squared(:,j) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL
    
            ! If we use complex timestepping we need to back-up U for variables 1-3
            if (domain%timestepping_method /= 'euler') then
                allocate(domain%backup_U(nx, ny, STG:VH))
                !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
                !$OMP DO SCHEDULE(STATIC)
                do j = 1, ny
                    domain%backup_U(:,j,STG:VH) = ZERO_dp
                end do
                !$OMP END DO
                !$OMP END PARALLEL
            end if

        endif

        if(domain%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction') then
            ! The linear_plus_nonlinear_friction needs a manning term, and a work array
            ! which can hold the "constant" part of the friction terms for efficiency

            ! Manning
            allocate(domain%manning_squared(nx, ny))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%manning_squared(:,j) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! A work array that will hold the expensive "constant" part of the friction terms
            allocate(domain%friction_work(nx, ny, UH:VH))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%friction_work(:,j, UH:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        end if

#ifdef DEBUG_ARRAY
        ! Optional for debugging -- this will be written to the netcdf file
        allocate(domain%debug_array(nx, ny))
        domain%debug_array = ZERO_dp
#endif

        if(create_output) then
            CALL domain%create_output_files()
        end if

        if(verbose_) then
            write(domain%logfile_unit, *) ''
            write(domain%logfile_unit, *) 'dx: ', domain%dx
            write(domain%logfile_unit, *) 'nx: ', domain%nx
            write(domain%logfile_unit, *) 'lw: ', domain%lw
            write(domain%logfile_unit, *) 'lower-left: ', domain%lower_left
            write(domain%logfile_unit, *) ''
            write(domain%logfile_unit, *) ''
            write(domain%logfile_unit, *) 'Total area: ', sum(domain%area_cell_y)
            write(domain%logfile_unit, *) 'distance_bottom_edge(1): ', domain%distance_bottom_edge(1)
            write(domain%logfile_unit, *) 'distange_left_edge(1)', domain%distance_left_edge(1)
            write(domain%logfile_unit, *) ''
        end if

    end subroutine

    ! Precompute the friction-work for the solver "leapfrog_linear_plus_nonlinear_friction"
    !
    ! This assumes stage/UH/VH/elev have been set
    !
    ! The friction work term is of the form
    !     g * n^2 / constant_depth^(7/3)
    ! The key point is that when multiplied by ||UH|| * uh, it will be equal to the
    ! standard friction form: g*constant_depth*friction_slope
    !
    ! The pre-computation removes an expensive power-law call from the inner loop
    !
    subroutine precompute_friction_work(domain)
        class(domain_type), intent(inout) :: domain

        integer(ip) :: i, j, jp1, ip1
        real(dp) :: depth_iph, depth_jph, nsq_iph, nsq_jph

        if(domain%timestepping_method /= 'leapfrog_linear_plus_nonlinear_friction') &
            stop 'precompute_friction_work can only be called with timestepping_method=leapfrog_linear_plus_nonlinear_friction'

        if(.not. allocated(domain%friction_work)) &
            stop 'friction_work is not allocated: ensure elevation is set before running this routine'

        if(domain%linear_solver_is_truely_linear) then
            !
            ! Here we evaluate the depth assuming stage=domain%msl_linear
            ! This is the standard logic for the "truely linear" linear shallow water equations.
            !

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)

                    ! UH component
                    ip1 = min(i+1, domain%nx(1))
                    depth_iph = 0.5_dp * (domain%msl_linear - domain%U(i,j,ELV) + domain%msl_linear - domain%U(ip1,j, ELV))
                    depth_iph = max(depth_iph, minimum_allowed_depth)
                    nsq_iph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(ip1,j))))**2
                    domain%friction_work(i,j,UH) = gravity * nsq_iph * (depth_iph**(-7.0_dp/3.0_dp))

                    ! VH component
                    jp1 = min(j+1, domain%nx(2))
                    depth_jph = 0.5_dp * (domain%msl_linear - domain%U(i,j,ELV) + domain%msl_linear - domain%U(i,jp1, ELV))
                    depth_jph = max(depth_jph, minimum_allowed_depth)
                    nsq_jph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(i,jp1))) )**2
                    domain%friction_work(i,j,VH) = gravity * nsq_jph * (depth_jph**(-7.0_dp/3.0_dp))

                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! For a truely-linear solver, this is only required once
            domain%friction_work_is_setup = .true.
        else
            !
            ! This differs from the above, in that we use domain%U(:,:,STG) to compute the depth.
            ! This means we generally need to recompute the friction-work terms at every time-step.
            !
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)

                    ! UH component
                    ip1 = min(i+1, domain%nx(1))
                    depth_iph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(ip1,j,STG) - domain%U(ip1,j, ELV))
                    depth_iph = max(depth_iph, minimum_allowed_depth)
                    nsq_iph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(ip1,j))))**2
                    domain%friction_work(i,j,UH) = gravity * nsq_iph * (depth_iph**(-7.0_dp/3.0_dp))

                    ! VH component
                    jp1 = min(j+1, domain%nx(2))
                    depth_jph = 0.5_dp * (domain%U(i,j,STG) - domain%U(i,j,ELV) + domain%U(i,jp1,STG) - domain%U(i,jp1, ELV))
                    depth_jph = max(depth_jph, minimum_allowed_depth)
                    nsq_jph = (0.5_dp * (sqrt(domain%manning_squared(i,j)) + sqrt(domain%manning_squared(i,jp1))) )**2
                    domain%friction_work(i,j,VH) = gravity * nsq_jph * (depth_jph**(-7.0_dp/3.0_dp))

                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
        end if

    end subroutine

    !
    ! minmod function, which is used in some gradient limiters
    !
    ! @param a,b real numbers
    !
    elemental function minmod(a,b) result(minmod_ab)
        real(dp), intent(in):: a, b
        real(dp):: minmod_ab
        
        minmod_ab = merge(min(abs(a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
        !minmod_ab = sign(one_dp, a) * max(zero_dp, min(abs(a), sign(one_dp, a)*b))

    end function

    !
    ! minmod subroutine, which is used in some gradient limiters
    !
    ! @param a,b real numbers
    !
    elemental subroutine minmod_sub(a,b, minmod_ab)
        real(dp), intent(in):: a, b
        real(dp), intent(out):: minmod_ab

        if((a>0.0_dp .and. b>0.0_dp).or.(a<0.0_dp .and. b<0.0_dp)) then
            minmod_ab = min(abs(a), abs(b))*sign(ONE_dp, a)
        else
            minmod_ab = ZERO_dp
        end if

    end subroutine

    !
    ! Get the "gradient times dx" around U_local
    ! @param U_local -- variable at the cell of interest, say x
    ! @param U_lower -- variable at (x - dx)
    ! @param U_upper -- variable at (x + dx)
    ! @param theta -- limiter parameter
    ! @param gradient_dx -- hold output
    ! @param n -- length of each vector
    !
    subroutine limited_gradient_dx_vectorized(U_local, U_lower, U_upper, theta, gradient_dx, n)
        integer(ip), intent(in):: n
        real(dp), intent(in):: U_local(n), U_lower(n), U_upper(n), theta(n) 
        real(dp), intent(out) :: gradient_dx(n)

        character(len=charlen), parameter :: limiter_type = 'MC' !'Minmod2' !'Superbee_variant' ! 'MC'

        integer(ip) :: i
        real(dp):: a, b, c, d, e, th, sa, sb, half_sasb

        if(limiter_type == 'MC') then

            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                th = theta(i)
                sa = sign(ONE_dp,a)
                sb = sign(ONE_dp,b)
                half_sasb = HALF_dp * (sa + sb)
                d = min(abs(a), abs(b)) * half_sasb * th ! Limit on local gradient
                e = HALF_dp * (a + b)
                c = merge(ZERO_dp, e, d == ZERO_dp) 
                ! NOTE: IF d /= 0, then clearly d, c have the same sign
                ! We exploit this to avoid a further minmod call (which seems
                ! expensive)
                gradient_dx(i) = merge(min(c, d), max(c, d), d > ZERO_dp)
            end do

        else if(limiter_type == "Superbee_variant") then

            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                ! Divide by 1.6 which is the default 'max theta' in the rk2 algorithms
                th = theta(i) * 2.0_dp/1.6_dp
                !d = minmod(a, th*b)
                d = merge(min(abs(a), abs(th*b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
                !e = minmod(th*a, b)
                e = merge(min(abs(th*a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
                if(abs(e) > abs(d)) then
                    b = e
                else
                    b = d
                endif
                gradient_dx(i) = b 
            end do

        else if(limiter_type == "Minmod2") then
            ! Same as 'MC' (!!!)
            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                th = theta(i)
                e = HALF_dp * (a + b)
                a = a * th
                b = b * th
                if(b > ZERO_dp .and. a > ZERO_dp) then
                    ! Positive slopes
                    if(b < a) a = b
                    d = min(e, a)
                else if(b < ZERO_dp .and. a < ZERO_dp) then
                    ! Negative slopes
                    if(b > a) a = b
                    d = max(e, a)
                else
                    d = ZERO_dp
                endif
                gradient_dx(i) = d
            end do

        else 
            gradient_dx = ZERO_dp
        end if
            

    end subroutine

    !
    ! Get the NS gradients for stage, depth, u-vel, v-vel, at row j
    !
    subroutine get_NS_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_NS, dstage_NS, ddepth_NS, du_NS, dv_NS)

        class(domain_type), intent(in) :: domain
        integer(ip), intent(in) :: j, nx, ny
        real(dp), intent(inout) :: theta_wd_NS(nx), dstage_NS(nx), ddepth_NS(nx), du_NS(nx), dv_NS(nx)

        integer(ip) :: i
        real(dp) :: mindep, maxdep, theta_local

        if(j > 1 .and. j < ny) then
            ! Typical case
            
            ! limiter coefficient
            !$OMP SIMD
            do i = 1, nx
                mindep = min(domain%depth(i, j-1), domain%depth(i,j), domain%depth(i,j+1)) - minimum_allowed_depth
                maxdep = max(domain%depth(i, j-1), domain%depth(i,j), domain%depth(i,j+1)) + &
                    limiter_coef3*minimum_allowed_depth
                theta_local = limiter_coef4 * (mindep/maxdep - limiter_coef1)
                theta_wd_NS(i) = max(domain%theta * min(ONE_dp, theta_local), ZERO_dp)
            end do

            ! stage 
            call limited_gradient_dx_vectorized(domain%U(:,j,STG), domain%U(:,j-1,STG), domain%U(:,j+1,STG), &
                theta_wd_NS, dstage_NS, nx)
            ! depth 
            call limited_gradient_dx_vectorized(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1), &
                theta_wd_NS, ddepth_NS, nx)
            ! u velocity
            call limited_gradient_dx_vectorized(domain%velocity(:,j,UH), domain%velocity(:,j-1, UH), &
                domain%velocity(:,j+1, UH), theta_wd_NS, du_NS, nx)
            ! v velocity
            call limited_gradient_dx_vectorized(domain%velocity(:,j,VH), domain%velocity(:,j-1, VH), &
                domain%velocity(:,j+1, VH), theta_wd_NS, dv_NS, nx)

        else
            ! Border case -- all gradients = zero
            theta_wd_NS = ZERO_dp
            dstage_NS = ZERO_dp
            ddepth_NS = ZERO_dp
            du_NS = ZERO_dp
            dv_NS = ZERO_dp
        end if

    end subroutine

    !
    ! Get the EW gradients for stage, depth, u-vel, v-vel, at row j
    !
    subroutine get_EW_limited_gradient_dx(domain, j, nx, ny, &
            theta_wd_EW, dstage_EW, ddepth_EW, du_EW, dv_EW)

        class(domain_type), intent(in) :: domain
        integer(ip), intent(in) :: j, nx, ny
        real(dp), intent(inout) :: theta_wd_EW(nx), dstage_EW(nx), ddepth_EW(nx), du_EW(nx), dv_EW(nx)

        real(dp) :: mindep, maxdep, theta_local
        integer(ip) :: i
        
        ! limiter coefficient
        !$OMP SIMD
        do i = 2, nx-1
            mindep = min(domain%depth(i-1, j), domain%depth(i,j), domain%depth(i+1,j)) - minimum_allowed_depth
            maxdep = max(domain%depth(i-1, j), domain%depth(i,j), domain%depth(i+1,j)) + &
                limiter_coef3*minimum_allowed_depth
            theta_local = limiter_coef4 * (mindep/maxdep - limiter_coef1)
            theta_wd_EW(i) = max(domain%theta * min(ONE_dp, theta_local), ZERO_dp)
        end do

        ! stage 
        call limited_gradient_dx_vectorized(domain%U(2:(nx-1),j,STG), domain%U(1:(nx-2),j,STG), domain%U(3:nx,j,STG), &
            theta_wd_EW(2:(nx-1)), dstage_EW(2:(nx-1)), nx)
        dstage_EW(1) = ZERO_dp
        dstage_EW(nx) = ZERO_dp

        ! depth 
        call limited_gradient_dx_vectorized(domain%depth(2:(nx-1),j), domain%depth(1:(nx-2),j), domain%depth(3:nx,j), &
            theta_wd_EW(2:(nx-1)), ddepth_EW(2:(nx-1)), nx)
        ddepth_EW(1) = ZERO_dp
        ddepth_EW(nx) = ZERO_dp

        ! u velocity
        call limited_gradient_dx_vectorized(domain%velocity(2:(nx-1),j,UH), domain%velocity(1:(nx-2),j,UH), &
            domain%velocity(3:nx,j, UH), theta_wd_EW(2:(nx-1)), du_EW(2:(nx-1)), nx)
        du_EW(1) = ZERO_dp
        du_EW(nx) = ZERO_dp

        ! u velocity
        call limited_gradient_dx_vectorized(domain%velocity(2:(nx-1),j,VH), domain%velocity(1:(nx-2),j,VH), &
            domain%velocity(3:nx,j,VH), theta_wd_EW(2:(nx-1)), dv_EW(2:(nx-1)), nx)
        dv_EW(1) = ZERO_dp
        dv_EW(nx) = ZERO_dp

    end subroutine

    !
    ! Compute domain%depth and domain%velocity, so they are consistent with domain%U
    !
    ! We also implement a mass conservation check in this routine, and zero velocities
    ! at sites with depth < minimum_allowed_depth
    !
    subroutine compute_depth_and_velocity(domain)

        class(domain_type), intent(inout) :: domain

        real(dp) :: depth_inv
        integer(ip) :: i, j, masscon_error

        ! Recompute depth and velocity
        ! Must be updated before the main loop when done in parallel

        masscon_error = 0_ip
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain) REDUCTION(+:masscon_error)
        masscon_error = 0_ip
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)

                domain%depth(i,j) = domain%U(i,j,STG) - domain%U(i,j,ELV)

                if(domain%depth(i,j) > minimum_allowed_depth) then
                    ! Typical case
                    depth_inv = ONE_dp/domain%depth(i,j)
                    domain%velocity(i,j,UH) = domain%U(i,j,UH) * depth_inv
                    domain%velocity(i,j,VH) = domain%U(i,j,VH) * depth_inv
                else
                    ! Nearly dry or dry case
                    if(domain%depth(i,j) < ZERO_dp) then
                        ! Clip 'round-off' type errors.
                        domain%depth(i,j) = ZERO_dp
                        domain%U(i,j, STG) = domain%U(i,j,ELV)
                        ! Record that clipping occurred
                        masscon_error = masscon_error + 1_ip
                    end if
                    ! Zero velocities when the depth <= minimum-allowed_depth
                    domain%velocity(i,j,UH) = ZERO_dp
                    domain%velocity(i,j,VH) = ZERO_dp
                    domain%U(i,j,UH) = ZERO_dp
                    domain%U(i,j,VH) = ZERO_dp
                end if

            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        domain%negative_depth_fix_counter = domain%negative_depth_fix_counter + masscon_error

    end subroutine

! Experimental energy-conservative flux computation.
#include "domain_compute_fluxes_EEC_include.f90"        

    
    ! Regular DE1 flux
    subroutine compute_fluxes_DE1(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .false.
        logical, parameter :: upwind_transverse_momentum = .false.

#include "domain_compute_fluxes_DE1_include.f90"        

    end subroutine

    ! Low-diffusion DE1 flux
    subroutine compute_fluxes_DE1_low_fr_diffusion(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .true.
        logical, parameter :: upwind_transverse_momentum = .false.

#include "domain_compute_fluxes_DE1_include.f90"        

    end subroutine

    ! Regular DE1 flux with upwind transverse momentum flux
    subroutine compute_fluxes_DE1_upwind_transverse(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .false.
        logical, parameter :: upwind_transverse_momentum = .true.

#include "domain_compute_fluxes_DE1_include.f90"        

    end subroutine

    ! Low-diffusion DE1 flux with upwind transverse momentum flux
    subroutine compute_fluxes_DE1_low_fr_diffusion_upwind_transverse(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
        ! Providing this at compile time leads to substantial optimization. 
        logical, parameter :: reduced_momentum_diffusion = .true.
        logical, parameter :: upwind_transverse_momentum = .true.

#include "domain_compute_fluxes_DE1_include.f90"        

    end subroutine

    !
    ! Compute the fluxes, and other things, in preparation for an update of U 
    !
    ! Update values of:
    ! domain%flux_NS, domain%flux_EW, domain%max_dt, domain%explicit_source, 
    ! domain%explicit_source_VH_j_minus_1, domain%boundary_flux_store
    !
    ! @param domain the model domain type for which fluxes etc will be computed
    ! @param max_dt_out optional real scalar, if provided is set to the max_dt
    !     allowed based on CFL limit.
    subroutine compute_fluxes(domain, max_dt_out)

        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out

        real(dp):: max_dt_out_local
        integer(ip) :: nx, ny, n_ext

        ! Need to have depth/velocity up-to-date for flux computation
        call domain%compute_depth_and_velocity()

        !
        ! Flux computation -- optionally use a few different methods
        !
        select case(domain%compute_fluxes_inner_method)
        case('DE1_low_fr_diffusion')
            ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers
            call compute_fluxes_DE1_low_fr_diffusion(domain, max_dt_out_local)
        case('DE1')
            ! Something like ANUGA (sort of..)
            call compute_fluxes_DE1(domain, max_dt_out_local)
        case('DE1_low_fr_diffusion_upwind_transverse')
            ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers
            ! Also using an upwind treatment of the transverse flux term
            call compute_fluxes_DE1_low_fr_diffusion_upwind_transverse(domain, max_dt_out_local)
        case('DE1_upwind_transverse')
            ! Something like ANUGA (sort of..)
            ! Also using an upwind treatment of the transverse flux term
            call compute_fluxes_DE1_upwind_transverse(domain, max_dt_out_local)
        case('EEC')
            ! Experimental energy conservative method, no wetting/drying
            call compute_fluxes_EEC(domain, max_dt_out_local)
        case default
            write(domain%logfile_unit, *) 'ERROR: Cannot recognize domain%compute_fluxes_inner_method=',&
                TRIM(domain%compute_fluxes_inner_method)
            call generic_stop()
        end select

        ! Update the time-step
        if(present(max_dt_out)) max_dt_out = max_dt_out_local

        !
        ! Spatially integrate boundary fluxes. Order is N, E, S, W. 
        ! Note we integrate in the INTERIOR (i.e. ignoring the outer n_ext rows/columns which can be updated by boundary-condition
        ! routines, then integrate the boundary of what remains). 
        !
        n_ext = domain%exterior_cells_width
        nx = domain%nx(1)
        ny = domain%nx(2) 

        ! Outward boundary flux over the north
        domain%boundary_flux_store(1) = sum(domain%flux_NS( (1+n_ext):(nx-n_ext), ny+1-n_ext, STG))
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG))
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG))
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG))

        !
        ! Compute fluxes relevant to the multidomain case. The key difference to above is that we ignore fluxes associated with
        ! other "priority domain index/image" cells
        !
        if(domain%nesting%my_index > 0) then
            ! We are in a multidomain -- carefully compute fluxes through
            ! exterior boundaries
            domain%boundary_flux_store_exterior = ZERO_dp
       
            ! Here we implement masked versions of the boundary flux sums above, only counting cells where the priority domain is
            ! receiving/sending the fluxes on actual physical boundaries 

            ! North boundary
            if(domain%boundary_exterior(1)) then
                domain%boundary_flux_store_exterior(1) = sum(&
                    domain%flux_NS( (1+n_ext):(nx-n_ext), ny+1-n_ext, STG),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), ny-n_ext) == domain%nesting%my_image ))
            end if

            ! East boundary
            if(domain%boundary_exterior(2)) then
                domain%boundary_flux_store_exterior(2) = sum(&
                    domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (&
                        domain%nesting%priority_domain_index(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(nx-n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if

            ! South boundary
            if(domain%boundary_exterior(3)) then
                domain%boundary_flux_store_exterior(3) = -sum(&
                    domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG),&
                    mask = (&
                        domain%nesting%priority_domain_index((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image((1+n_ext):(nx-n_ext), 1+n_ext) == domain%nesting%my_image ))
            end if

            ! West boundary
            if(domain%boundary_exterior(4)) then
                domain%boundary_flux_store_exterior(4) = -sum(&
                    domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (&
                        domain%nesting%priority_domain_index(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_index .and. &
                        domain%nesting%priority_domain_image(1+n_ext, (1+n_ext):(ny-n_ext)) == domain%nesting%my_image ))
            end if 
        else
            ! We are not doing nesting
            domain%boundary_flux_store_exterior = domain%boundary_flux_store
        end if

    end subroutine

    !
    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
    !
    ! @param domain the domain to be updated
    ! @param dt timestep to advance
    !
    subroutine update_U(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        real(dp):: inv_cell_area_dt, depth, implicit_factor, dt_gravity, fs
        integer(ip):: j, i, kk

        domain%dt_last_update = dt


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !$OMP SIMD
            do i = 1, domain%nx(1)
                !! Fluxes
                do kk = 1, 3
                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
                end do
            end do

            !$OMP SIMD
            do i = 1, domain%nx(1)

                ! Velocity clipping
                depth = domain%depth(i,j)

                if (domain%U(i,j,STG) <= (domain%U(i,j,ELV) + minimum_allowed_depth)) then
                    domain%U(i,j,UH) = ZERO_dp 
                    domain%U(i,j,VH) = ZERO_dp 
                else
#ifndef NOFRICTION
                    ! Implicit friction slope update
                    ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new

                    ! If we multiply this by UH or VH, we get the associated friction slope term
                    fs = domain%manning_squared(i,j) * &
                        sqrt(domain%velocity(i,j,UH) * domain%velocity(i,j,UH) + &
                             domain%velocity(i,j,VH) * domain%velocity(i,j,VH)) * &
                        (max(depth, minimum_allowed_depth)**(NEG_SEVEN_ON_THREE_dp))

                    implicit_factor = ONE_dp/(ONE_dp + dt_gravity*max(depth, ZERO_dp)*fs)
                    !if(domain%U(i,j,STG) <= (domain%U(i,j,ELV) + minimum_allowed_depth)) implicit_factor = ZERO_dp
#else
                    implicit_factor = ONE_dp
#endif
                    ! Pressure gradients 
                    domain%U(i,j,UH) = domain%U(i,j,UH) + inv_cell_area_dt * domain%explicit_source(i, j, UH)
                    ! Friction
                    domain%U(i,j,UH) = domain%U(i,j,UH) * implicit_factor

                    ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
                    ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
                    domain%U(i,j,VH) = domain%U(i,j,VH) + inv_cell_area_dt * &
                        (domain%explicit_source(i, j, VH) + domain%explicit_source_VH_j_minus_1(i, j+1))
                    ! Friction
                    domain%U(i,j,VH) = domain%U(i,j,VH) * implicit_factor
                end if
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + dt

        call domain%update_boundary()
    
    end subroutine

    !
    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
    ! Same as update_U but slightly restructured (faster in some cases)
    !
    ! @param domain the domain to be updated
    ! @param dt timestep to advance
    !
    subroutine update_U_restructured(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        real(dp):: inv_cell_area_dt, depth, implicit_factor(domain%nx(1)), dt_gravity, fs, depth_neg7on3
        integer(ip):: j, i, kk


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !$OMP SIMD
            do i = 1, domain%nx(1)
        
                !! Fluxes
                do kk = 1, 3
                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
                end do
            end do

            !$OMP SIMD
            do i = 1, domain%nx(1)

                depth = domain%depth(i,j)
                depth_neg7on3 = max(depth, minimum_allowed_depth)**(NEG_SEVEN_ON_THREE_dp)

#ifndef NOFRICTION
                    ! Implicit friction slope update
                    ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new

                    ! If we multiply this by UH or VH, we get the associated friction slope term
                    fs = domain%manning_squared(i,j) * &
                        sqrt(domain%velocity(i,j,UH) * domain%velocity(i,j,UH) + &
                             domain%velocity(i,j,VH) * domain%velocity(i,j,VH)) * &
                        depth_neg7on3

                    ! Velocity clipping
                    implicit_factor(i) = merge( ONE_dp/(ONE_dp + dt_gravity*max(depth, ZERO_dp)*fs), ZERO_dp, &
                        (domain%U(i,j,STG) > (domain%U(i,j,ELV) + minimum_allowed_depth)) )
#else
                    implicit_factor(i) = ONE_dp
#endif
            end do

            !$OMP SIMD
            do i = 1, domain%nx(1)
                    ! Pressure gradients 
                    domain%U(i,j,UH) = domain%U(i,j,UH) + inv_cell_area_dt * domain%explicit_source(i, j, UH)
                    ! Friction
                    domain%U(i,j,UH) = domain%U(i,j,UH) * implicit_factor(i)

                    ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
                    ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
                    domain%U(i,j,VH) = domain%U(i,j,VH) + inv_cell_area_dt * &
                        (domain%explicit_source(i, j, VH) + domain%explicit_source_VH_j_minus_1(i, j+1))
                    ! Friction
                    domain%U(i,j,VH) = domain%U(i,j,VH) * implicit_factor(i)
                !end if
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + dt

        call domain%update_boundary()
    
    end subroutine
    
    !
    ! Vectorized version of update_U (i.e. using arrays to avoid the 'i' loop)
    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
    !
    ! @param domain the domain to be updated
    ! @param dt timestep to advance
    !
    subroutine update_U_vectorized(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        real(dp):: inv_cell_area_dt, implicit_factor(domain%nx(1)), dt_gravity, fs(domain%nx(1))
        integer(ip):: j, kk


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(GUIDED)
        do j = 1, domain%nx(2)
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !! Fluxes
            do kk = 1, 3
                domain%U(:,j,kk) = domain%U(:,j,kk) - inv_cell_area_dt * ( & 
                    (domain%flux_NS(:, j+1, kk) - domain%flux_NS(:, j, kk)) + &
                    (domain%flux_EW(2:(domain%nx(1)+1), j, kk) - domain%flux_EW(1:domain%nx(1), j, kk) ))
            end do

#ifndef NOFRICTION
            ! Implicit friction slope update
            ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new

            !! If we multiply this by UH or VH, we get the associated friction slope term
            fs = domain%manning_squared(:,j) * &
                sqrt(domain%velocity(:,j,UH) * domain%velocity(:,j,UH) + &
                     domain%velocity(:,j,VH) * domain%velocity(:,j,VH)) * &
                !norm2(domain%velocity(i,j,UH:VH)) * &
                (max(domain%depth(:,j), minimum_allowed_depth)**(NEG_SEVEN_ON_THREE_dp))

            implicit_factor = ONE_dp/(ONE_dp + dt_gravity*max(domain%depth(:,j), ZERO_dp)*fs)
#else
            implicit_factor = ONE_dp
#endif
            ! Velocity clipping here
            implicit_factor = merge(implicit_factor, ZERO_dp, domain%U(:,j,STG) > (domain%U(:,j,ELV) + minimum_allowed_depth))
            ! Pressure gradients 
            domain%U(:,j,UH) = domain%U(:,j,UH) + inv_cell_area_dt * domain%explicit_source(:, j, UH)
            ! Friction
            domain%U(:,j,UH) = domain%U(:,j,UH) * implicit_factor

            ! Here we add an extra source to take care of the j-1 pressure gradient addition, which
            ! could not be updated in parallel (since j-1 might be affected by another OMP thread) 
            domain%U(:,j,VH) = domain%U(:,j,VH) + inv_cell_area_dt * &
                (domain%explicit_source(:, j, VH) + domain%explicit_source_VH_j_minus_1(:, j+1))
            ! Friction
            domain%U(:,j,VH) = domain%U(:,j,VH) * implicit_factor
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + dt

        call domain%update_boundary()
    
    end subroutine


    ! 
    ! Standard forward euler 1st order timestepping scheme 
    !
    ! This is also used as a component of more advanced timestepping schemes
    ! Argument 'timestep' is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced. If not
    !    provided, advance by the CFL permitted timestep, computed by domain%compute_fluxes
    ! @param update_nesting_boundary_flux logical. If TRUE (default), then call 
    !    domain%nesting_boundary_flux_integral_tstep inside the routine. Sometimes it is 
    !    more straightforward to switch this off, and do the computations separately (e.g in rk2).
    subroutine one_euler_step(domain, timestep, update_nesting_boundary_flux_integral)
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep
        logical, optional, intent(in) :: update_nesting_boundary_flux_integral

        character(len=charlen):: timer_name
        real(dp):: ts
        logical :: nesting_bf

        ! By default, update the nesting_boundary_flux_integral tracker
        if(present(update_nesting_boundary_flux_integral)) then
            nesting_bf = update_nesting_boundary_flux_integral
        else
            nesting_bf = .TRUE.
        end if
            

        !TIMER_START('flux')
        call domain%compute_fluxes(ts)
        !TIMER_STOP('flux')

        !TIMER_START('update')

        ! Store the internally computed max timestep
        domain%max_dt = ts

        if(present(timestep)) then
            ts = timestep
        end if

        call domain%update_U(ts)

        ! Track flux through boundaries
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral + &
            ts * sum(domain%boundary_flux_store)
        ! Track flux through 'exterior' boundaries (i.e. not nesting)
        domain%boundary_flux_evolve_integral_exterior = domain%boundary_flux_evolve_integral_exterior + &
            ts * sum(domain%boundary_flux_store_exterior, mask=domain%boundary_exterior)

        !TIMER_STOP('update')

        if(nesting_bf) then 
            ! Update the nesting boundary flux
            call domain%nesting_boundary_flux_integral_tstep(&
                ts,&
                flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
                flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
                var_indices=[STG, VH],&
                flux_already_multiplied_by_dx=.TRUE.)
        end if

        ! Coarray communication, if required (this has been superceeded by the multidomain approach)
        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

    end subroutine

    !
    ! Standard 2-step second order timestepping runge-kutta scheme
    ! Argument timestep is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced
    !
    !
    subroutine one_rk2_step(domain, timestep)
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep

        real(dp):: backup_time, dt_first_step, max_dt_store
        integer(ip):: j
        character(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()

        !
        ! rk2 involves taking 2 euler steps, then halving the change.
        !

        ! First euler step -- 
        if(present(timestep)) then
            call domain%one_euler_step(timestep, &
                update_nesting_boundary_flux_integral = .FALSE.)
            dt_first_step = timestep
        else
            call domain%one_euler_step(update_nesting_boundary_flux_integral = .FALSE.)
            dt_first_step = domain%max_dt
        end if

        ! Update the nesting boundary flux, by only 1/2 tstep. By keeping it outside of
        ! the euler-step, we can accumulate the result over multiple time-steps without
        ! needing a temporary
        call domain%nesting_boundary_flux_integral_tstep(&
            dt=(dt_first_step*HALF_dp),&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=[STG, VH],&
            flux_already_multiplied_by_dx=.TRUE.)


        max_dt_store = domain%max_dt
       
        ! Second euler step with the same timestep 
        call domain%one_euler_step(dt_first_step, &
            update_nesting_boundary_flux_integral=.FALSE.)

        ! Update the nesting boundary flux, by the other 1/2 tstep
        call domain%nesting_boundary_flux_integral_tstep(&
            dt=(dt_first_step*HALF_dp),&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=[STG, VH],&
            flux_already_multiplied_by_dx=.TRUE.)


        !TIMER_START('average')

        ! Take average (but allow for openmp)
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            domain%U(:, j, STG) = HALF_dp * (domain%U(:, j, STG) +&
                domain%backup_U(:, j, STG))
            domain%U(:, j, UH) = HALF_dp * (domain%U(:, j, UH) +&
                domain%backup_U(:, j, UH))
            domain%U(:, j, VH) = HALF_dp * (domain%U(:, j, VH) +&
                domain%backup_U(:, j, VH))
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time (since we updated twice) and boundary flux integral
        domain%time = backup_time + dt_first_step
        domain%boundary_flux_evolve_integral = HALF_dp * domain%boundary_flux_evolve_integral
        domain%boundary_flux_evolve_integral_exterior = HALF_dp * domain%boundary_flux_evolve_integral_exterior

        ! We want the CFL timestep that is reported to always be based on the first step -- so force that here
        domain%max_dt = max_dt_store

        !TIMER_STOP('average')

    end subroutine
    

    !
    ! An n-step 2nd order runge-kutta timestepping scheme
    !
    ! This involves less flux calls per timestep advance than rk2.
    !
    ! BEWARE: It advances (n-1) timesteps, so a call one_rk2n_step(domain, 1.0_dp)
    ! will advance 4.0 seconds. 
    !
    ! Argument timestep is optional, but if provided should satisfy the CFL condition
    !
    ! FIXME: Still need to implement nesting boundary flux integral timestepping, if we
    ! want to use this inside a multidomain
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep, by which the solution is advanced (n-1) times
    !
    subroutine one_rk2n_step(domain, timestep)
        ! Advance (n-1) * timesteps in this routine, with 2nd order in time accuracy
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep

        integer(ip), parameter:: n = 5 ! number of substeps to take, must be > 2
        real(dp), parameter:: n_inverse = ONE_dp / (ONE_dp * n)
        real(dp):: backup_time, dt_first_step, reduced_dt, max_dt_store
        real(dp):: backup_flux_integral_exterior, backup_flux_integral
        integer(ip):: j,k
        character(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
     
        if(present(timestep)) then 
            ! first step 
            call domain%one_euler_step(timestep,&
                update_nesting_boundary_flux_integral=.FALSE.)
            ! store timestep
            dt_first_step = timestep
        else
            ! first step 
            call domain%one_euler_step(update_nesting_boundary_flux_integral=.FALSE.)
            ! store timestep
            dt_first_step = domain%max_dt
        end if

        max_dt_store = domain%max_dt

        ! Steps 2, n-1
        do j = 2, n-1        
            call domain%one_euler_step(dt_first_step,&
                update_nesting_boundary_flux_integral=.FALSE.)
        end do

        ! Store 1/n * (original_u - u)

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        do k = 1, 3
            do j = 1, domain%nx(2)
                domain%backup_U(:,j,k) = (domain%backup_U(:,j, k) - domain%U(:,j,k)) * n_inverse
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Store this for boundary flux tracking
        backup_flux_integral = domain%boundary_flux_evolve_integral*n_inverse
        backup_flux_integral_exterior = domain%boundary_flux_evolve_integral_exterior*n_inverse

        ! Now take one step of duration (n-1)/n * dt
        reduced_dt = (ONE_dp * n - ONE_dp) * n_inverse * dt_first_step
        call domain%one_euler_step(reduced_dt, &
            update_nesting_boundary_flux_integral=.FALSE.)
        
        !TIMER_START('final_update')

        ! Final update
        ! domain%U = domain%U + domain%backup_U
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        do k = 1, 3
            do j = 1, domain%nx(2)
                domain%U(:, j, k) = domain%U(:, j, k) + domain%backup_U(:, j, k)
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Get final boundary flux integral
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral - &
            backup_flux_integral
        domain%boundary_flux_evolve_integral_exterior = domain%boundary_flux_evolve_integral_exterior - &
            backup_flux_integral_exterior

        ! Fix time and timestep (since we updated (n-1)*dt regular timesteps)
        domain%time = backup_time + dt_first_step * (n-1) 

        ! We want the max_dt reported to always be based on the first-sub-step CFL (since that is
        ! what is required for stability, even though at later sub-steps the max_dt might reduce)
        domain%max_dt = max_dt_store

        !TIMER_STOP('final_update')

    end subroutine


    !
    ! Another 2nd order timestepping scheme (like the trapezoidal rule) 
    ! Argument timestep is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced
    !
    subroutine one_midpoint_step(domain, timestep)
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in) :: timestep

        real(dp):: backup_time, dt_first_step
        integer(ip):: j

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
        
        !TIMER_START('flux')
        call domain%compute_fluxes(dt_first_step)
        !TIMER_STOP('flux')

        domain%max_dt = dt_first_step

        !TIMER_START('update')
        ! First euler sub-step
        if(present(timestep)) then

            dt_first_step = timestep
            call domain%update_U(dt_first_step*HALF_dp)
        else
            ! First euler sub-step
            call domain%update_U(dt_first_step*HALF_dp)
        end if
        !TIMER_STOP('update')

        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

        ! Compute fluxes 
        !TIMER_START('flux')
        call domain%compute_fluxes()
        !TIMER_STOP('flux')

        ! Set U back to backup_U
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            domain%U(:, j, STG) = domain%backup_U(:, j, STG)
            domain%U(:, j,  UH) = domain%backup_U(:, j, UH)
            domain%U(:, j,  VH) = domain%backup_U(:, j, VH)
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time
        domain%time = backup_time

        ! Update U
        !TIMER_START('update')
        call domain%update_U(dt_first_step)
        !TIMER_STOP('update')

        ! Update the nesting boundary flux
        call domain%nesting_boundary_flux_integral_tstep(&
            dt_first_step,&
            flux_NS=domain%flux_NS, flux_NS_lower_index=1_ip, &
            flux_EW=domain%flux_EW, flux_EW_lower_index=1_ip, &
            var_indices=[STG, VH],&
            flux_already_multiplied_by_dx=.TRUE.)


        !TIMER_START('partitioned_comms')
        if(domain%use_partitioned_comms) call domain%partitioned_comms%communicate(domain%U)
        !TIMER_STOP('partitioned_comms')

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store)*&
            dt_first_step
        domain%boundary_flux_evolve_integral_exterior = sum(domain%boundary_flux_store_exterior, mask=domain%boundary_exterior)*&
            dt_first_step


    end subroutine

    ! 
    ! Linear shallow water equations leap-frog update
    !
    ! Update domain%U by timestep dt, using the linear shallow water equations.
    ! Note that unlike the other timestepping routines, this does not require a
    ! prior call to domain%compute_fluxes 
    !
    ! @param domain the domain to advance
    ! @param dt the timestep to advance. Should remain constant in between repeated calls
    !     to the function (since the numerical method assumes constant timestep)
    !
    subroutine one_linear_leapfrog_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        if(domain%linear_solver_is_truely_linear) then
            call one_truely_linear_leapfrog_step(domain, dt)
        else
            call one_not_truely_linear_leapfrog_step(domain, dt)
        end if

    end subroutine

    ! Truely-linear leap-frog solver
    subroutine one_truely_linear_leapfrog_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .TRUE.

        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 
#include "domain_linear_solver_include.f90"

    end subroutine

    ! Not-truely-linear leap-frog solver
    subroutine one_not_truely_linear_leapfrog_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .FALSE.

        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 
#include "domain_linear_solver_include.f90"

    end subroutine

    ! Linear combined with nonlinear friction. This might be useful to allow some dissipation
    ! in very large scale tsunami models.
    subroutine one_leapfrog_linear_plus_nonlinear_friction_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt

        ! Compute some expensive parts of the friction term
        ! If domain%linear_solver_is_truely_linear, this will only happen
        ! once (on the first timestep). Otherwise it should happen every timestep
        if(.not. domain%friction_work_is_setup) call domain%precompute_friction_work

        if(domain%linear_solver_is_truely_linear) then
            call one_leapfrog_truely_linear_plus_nonlinear_friction_step(domain, dt)
        else
            call one_leapfrog_not_truely_linear_plus_nonlinear_friction_step(domain, dt)
        end if

    end subroutine

    ! Truely-linear leap-frog solver PLUS NONLINEAR FRICTION
    subroutine one_leapfrog_truely_linear_plus_nonlinear_friction_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .true.
        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 

#define LINEAR_PLUS_NONLINEAR_FRICTION
#include "domain_linear_solver_include.f90"
#undef LINEAR_PLUS_NONLINEAR_FRICTION

    end subroutine

    ! Not-truely-linear leap-frog solver PLUS NONLINEAR FRICTION
    subroutine one_leapfrog_not_truely_linear_plus_nonlinear_friction_step(domain, dt)
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt
        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .false.
        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 

#define LINEAR_PLUS_NONLINEAR_FRICTION
#include "domain_linear_solver_include.f90"
#undef LINEAR_PLUS_NONLINEAR_FRICTION

    end subroutine
    
    !
    ! Routine to run all boundary conditions
    !
    subroutine update_boundary(domain)
        class(domain_type), intent(inout):: domain 

        !TIMER_START('boundary_update')

        if(associated(domain%boundary_subroutine)) then
            CALL domain%boundary_subroutine(domain)
        end if

        !TIMER_STOP('boundary_update')

    end subroutine

    !
    ! Copy domain%U to domain%backup_U
    !
    subroutine backup_quantities(domain)
        
        class(domain_type), intent(inout):: domain
        integer(ip):: j, k

        !TIMER_START('backup')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        do k = STG, VH
            do j = 1, domain%nx(2)
                domain%backup_U(:, j, k) = domain%U(:, j, k)
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        !TIMER_STOP('backup')

    end subroutine

    !
    ! Main 'high-level' evolve routine
    !
    ! 'timestep' is an optional argument. If provided it must satisfy the CFL condition.
    ! Otherwise the maximum timestep satisfying the CFL condition is used (considering
    ! also the cfl number stored in domain%cfl)
    !
    ! The actual timestepping method used is determined by the value of
    ! domain%timestepping_method ('euler', 'rk2', 'rk2n', 'midpoint')
    !
    subroutine evolve_one_step(domain, timestep)
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep

        real(dp):: time0
        character(len=charlen) :: timestepping_method
        real(dp):: static_before_time

TIMER_START('evolve_one_step')

        timestepping_method = domain%timestepping_method
        time0 = domain%time
        static_before_time = domain%static_before_time

        ! Check whether we should 'skip' the evolve
        if(time0 < static_before_time) then
            timestepping_method = 'static'
        end if

        ! Reset the boundary fluxes integrated over the evolve step to zero
        domain%boundary_flux_evolve_integral = ZERO_dp
        domain%boundary_flux_evolve_integral_exterior = ZERO_dp
   
        select case (timestepping_method) 

        case ('static')
            ! Do nothing -- except update the time if the timestep was provided
            if(present(timestep)) then
                domain%time = domain%time + timestep 
            end if
        case ('euler')
            if(present(timestep)) then
                call domain%one_euler_step(timestep)
            else
                call domain%one_euler_step()
            end if
        case ('rk2')
            if(present(timestep)) then
                call domain%one_rk2_step(timestep)
            else
                call domain%one_rk2_step()
            end if
        case('rk2n')
            if(present(timestep)) then
                call domain%one_rk2n_step(timestep)
            else
                call domain%one_rk2n_step()
            end if
        case ('midpoint')
            if(present(timestep)) then
                call domain%one_midpoint_step(timestep)
            else
                call domain%one_midpoint_step()
            end if
        case ('linear')
            if(present(timestep)) then
                call domain%one_linear_leapfrog_step(timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for linear evolve_one_step'
                call generic_stop()
            end if
        case ('leapfrog_linear_plus_nonlinear_friction')
            if(present(timestep)) then
                call domain%one_leapfrog_linear_plus_nonlinear_friction_step(timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for ', &
                    'leapfrog_linear_plus_nonlinear_friction evolve_one_step'
                call generic_stop()
            end if

        case default
            write(domain%logfile_unit,*) 'ERROR: domain%timestepping_method not recognized'
            call generic_stop()
        end select

        ! For some problems updating max U can take a significant fraction of the time,
        ! so we allow it to only be done occasionally
        if(timestepping_method /= 'static') then
            if(mod(domain%nsteps_advanced, domain%max_U_update_frequency) == 0) then
                call domain%update_max_quantities()
            end if
        end if

        ! Record the timestep here
        domain%evolve_step_dt = domain%time - time0

        ! Update the boundary flux time integral
        domain%boundary_flux_time_integral = domain%boundary_flux_time_integral + &
            domain%boundary_flux_evolve_integral
        domain%boundary_flux_time_integral_exterior = domain%boundary_flux_time_integral_exterior + &
            domain%boundary_flux_evolve_integral_exterior
        

        domain%nsteps_advanced = domain%nsteps_advanced + 1

TIMER_STOP('evolve_one_step')
    end subroutine

    !
    ! Convenience function to compute the volume of water in the 'interior' of the domain
    ! This involves all parts of the domain that are more than domain%exterior_cells_width
    ! from the edge. For multidomains, there are other purpose-built routines -- this one
    ! is appropriate for the single-domain case.
    !
    function volume_interior(domain) result(domain_volume)
        class(domain_type), intent(in):: domain
        real(force_double) :: domain_volume

        integer(ip) :: j, n_ext
        real(force_double) :: local_sum
        
        ! Volume on the interior. At the moment the interior is all
        ! but the outer cells of the domain, but that could change.

        !!TIMER_START('volume_interior')
        n_ext = domain%exterior_cells_width
        domain_volume = ZERO_dp
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, n_ext) REDUCTION(+:domain_volume)
        domain_volume = ZERO_dp
        !$OMP DO SCHEDULE(STATIC)
        do j = (1+n_ext), (domain%nx(2) - n_ext) 
            ! Be careful about precision
            local_sum = sum(&
                real(domain%U((1+n_ext):(domain%nx(1)-n_ext),j,STG), force_double) - &
                real(domain%U((1+n_ext):(domain%nx(1)-n_ext),j,ELV), force_double))

            domain_volume = domain_volume + real(domain%area_cell_y(j) * local_sum, force_double)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        !!TIMER_STOP('volume_interior')

    end function

    !
    ! Compute the volume on interior cells and add to the integrated boundary flux. 
    !
    ! This should sum to a constant in the absence of mass sources, for single domains. Note this relies on the time-stepping
    ! routines correctly computing the boundary_flux_time_integral. 
    !
    ! For multidomains, there are other purpose-built mass conservation routines (which take into account the "priority domain"
    ! information)
    !
    function mass_balance_interior(domain) result(mass_balance)
        class(domain_type), intent(in):: domain
        real(dp) :: mass_balance

        mass_balance = domain%volume_interior() + domain%boundary_flux_time_integral 

    end function

    !
    ! Function to compute the max timestep allowed for the linear shallow water equations, using the provided CFL, assuming the
    ! leap-frog time-stepping scheme
    !
    function linear_timestep_max(domain) result(timestep)
        class(domain_type), intent(in):: domain
        real(dp) :: timestep

        real(dp) :: ts_max
        integer(ip) :: i, j

        ts_max = huge(ONE_dp)
    
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                ! max timestep = Distance along latitude / wave-speed <= dt
                ts_max = min(ts_max, &
                    0.5_dp * min(&
                        (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)),& 
                        (domain%distance_left_edge(i+1)   + domain%distance_left_edge(i)  ) ) / &
                    sqrt(gravity * max(domain%msl_linear-domain%U(i,j,ELV), minimum_allowed_depth)) )
            end do
        end do
        timestep = ts_max * domain%cfl

    end function

    !
    ! Function to compute the max timestep allowed for the nonlinear shallow water equations, for a stationary flow (i.e. gravity
    ! wave only), using the provided CFL. 
    !
    ! Results are not exactly the same as the gravity-wave timestep limit for nonlinear-FV solvers, because the latter compute wave
    ! speeds at edges rather than cell centres. But results are very close in general for a stationary flow.
    !
    function nonlinear_stationary_timestep_max(domain) result(timestep)
        class(domain_type), intent(in):: domain
        real(dp) :: timestep

        real(dp) :: ts_max
        integer(ip) :: i, j

        ts_max = huge(ONE_dp)
    
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                ! max timestep = Distance along latitude / wave-speed <= dt
                ! Actually the non-linear schemes can often only take half this step, but we correct that elsewhere
                ts_max = min(ts_max, &
                    0.5_dp * min(&
                        (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)),& 
                        (domain%distance_left_edge(i+1)   + domain%distance_left_edge(i)  ) ) / &
                        sqrt(gravity * max(domain%U(i,j,STG)-domain%U(i,j,ELV), minimum_allowed_depth)) )
            end do
        end do
        timestep = ts_max * domain%cfl

    end function

    ! Convenience wrapper which uses 'linear_timestep_max' or 'nonlinear_stationary_timestep_max' depending on the numerical method.
    !
    ! This is mainly useful when setting up models, as a guide to the timestep that one might be able to take.
    !
    ! Beware the time-step calculations are not exactly the same as the results from the nonlinear-FV solvers (because the latter use
    ! the edge-based wave speed calculations). But they should be very similar. 
    function stationary_timestep_max(domain) result(timestep)
        class(domain_type), intent(in):: domain
        real(dp) :: timestep
        ! Leapfrog type numerical methods
        character(len=charlen) :: leapfrog_type_solvers(2) = [ character(len=charlen) :: &
            'linear', 'leapfrog_linear_plus_nonlinear_friction' ]
        ! Typical DE1-type finite volume methods
        character(len=charlen) :: nonlinear_FV_solvers_1(3) = [ character(len=charlen) :: &
            'euler', 'rk2', 'midpoint' ]
        ! Less typical finite-volume methods
        character(len=charlen) :: nonlinear_FV_solvers_2(1) = [ character(len=charlen) :: &
            'rk2n' ]

        if( any(domain%timestepping_method == leapfrog_type_solvers) ) then
            !
            ! Leap-frog type solvers
            !
            if(domain%linear_solver_is_truely_linear) then
                timestep = domain%linear_timestep_max()
            else
                timestep = domain%nonlinear_stationary_timestep_max()
            end if

        else if(any(domain%timestepping_method == nonlinear_FV_solvers_1)) then
            !
            ! Standard nonlinear solvers in SWALS
            !
            timestep = domain%nonlinear_stationary_timestep_max() * 0.5_dp
        else if(any(domain%timestepping_method == nonlinear_FV_solvers_2)) then
            !
            ! The unusual 'rk2n' timestepping method. 
            !
            ! Note 'rk2n' will advance 4 timesteps with one call to the timestepping routine -- consider changing that interface, in
            ! which case the '1.0' here would change to '4.0'. Currently 'rk2n' is not supported for multidomain, and is not
            ! commonly used, so leave as-is for now.
            !
            timestep = domain%nonlinear_stationary_timestep_max() * 0.5_dp * 1.0_dp
        else
            write(log_output_unit, *) 'Error in stationary_timestep_max: Unrecognized timestepping_method'
            write(log_output_unit, *) trim(domain%timestepping_method)
            call generic_stop
        end if

    end function

    !
    ! Initialise output files
    !
    subroutine create_output_files(domain)
        class(domain_type), intent(inout):: domain

        character(len=charlen):: mkdir_command, cp_command, t1, t2, t3, t4, &
                                 output_folder_name
        integer(ip):: i, metadata_unit, natt
        character(len=charlen), allocatable :: attribute_names(:), attribute_values(:)

        ! Create output directory
        call date_and_time(t1, t2, t3)
        ! Get domain id as a character
        write(t3, domain_myid_char_format) domain%myid
        write(t4, '(I0.5)') domain%local_index

        output_folder_name = trim(domain%output_basedir) // '/RUN_ID' // trim(t3) // &
            '_' // trim(t4) // '_' // trim(t1) // '_' // trim(t2)
        call mkdir_p(output_folder_name)

        ! Make a filename to hold domain metadata, and write the metadata
        t1 = trim(output_folder_name) // '/' // 'Domain_info_ID' // trim(t3) // '.txt'
        domain%metadata_ascii_filename = t1
        open(newunit = metadata_unit, file=domain%metadata_ascii_filename)
        write(metadata_unit, *) 'id :', domain%myid
        write(metadata_unit, *) 'nx :', domain%nx
        write(metadata_unit, *) 'dx :', domain%dx
        write(metadata_unit, *) 'lower_left_corner: ', domain%lower_left
        write(metadata_unit, *) 'dp_precision: ', dp
        write(metadata_unit, *) 'ip_precision: ', ip
        write(metadata_unit, *) 'output_precision: ', output_precision
        close(metadata_unit)

#ifdef NONETCDF
        ! Write to a home-brew binary format
        allocate(domain%output_variable_unit_number(domain%nvar))
        do i=1, domain%nvar
            ! Make output_file_name
            write(t1,'(I1)') i
            write(t2, *) 'Var_', trim(t1), '_ID', trim(t3)
            t1 = trim(output_folder_name) // '/' // adjustl(trim(t2))

            ! Binary. Using 'stream' access makes it easy to read in R
            open(newunit = domain%output_variable_unit_number(i), file = t1, &
                access='stream', form='unformatted')
        end do

#else
        !
        ! Write to netcdf
        !

        t1 = trim(output_folder_name) // '/' // 'Grid_output_ID' // trim(t3) // '.nc'

        ! Character attributes
        natt = 7!8
        allocate(attribute_names(natt), attribute_values(natt))
        ! 
        attribute_names(1) = 'timestepping_method'
        attribute_values(1) = trim(domain%timestepping_method)
        !
        attribute_names(2) = 'myid'
        write(attribute_values(2), domain_myid_char_format) domain%myid
        !
        attribute_names(3) = 'dx_refinement_factor'
        write(attribute_values(3), *) domain%dx_refinement_factor
        !
        attribute_names(4) = 'max_parent_dx_ratio'
        write(attribute_values(4), *) domain%max_parent_dx_ratio
        !
        attribute_names(5) = 'timestepping_refinement_factor'
        write(attribute_values(5), *) domain%timestepping_refinement_factor
        !
        attribute_names(6) = 'interior_bounding_box'
        write(attribute_values(6), *) domain%interior_bounding_box
        ! 
        attribute_names(7) = 'is_nesting_boundary_N_E_S_W' 
        write(attribute_values(7), *) domain%is_nesting_boundary

        call domain%nc_grid_output%initialise(filename=t1,&
            output_precision = output_precision, &
            record_max_U = domain%record_max_U, &
            xs = domain%x, ys = domain%y, &
            attribute_names=attribute_names, attribute_values=attribute_values)

        deallocate(attribute_names, attribute_values)

        ! If we have setup domain%nesting, then store locations where the current
        ! domain is the priority domain
        if(allocated(domain%nesting%priority_domain_index)) then
            call domain%nc_grid_output%store_priority_domain_cells(&
                domain%nesting%priority_domain_index,&
                domain%nesting%priority_domain_image,&
                domain%nesting%my_index,&
                domain%nesting%my_image)
        end if
#endif

        ! Make a time file_name. Store as ascii
        t1 = trim(output_folder_name) // '/' // 'Time_ID' // trim(t3) // '.txt'
        open(newunit = domain%output_time_unit_number, file = t1)

        ! Copy code to the output directory

        mkdir_command = 'mkdir -p ' // trim(output_folder_name) // '/Code'
        !call execute_command_line(trim(mkdir_command))
        call system(trim(mkdir_command))
        
        cp_command = 'cp *.f* make* ' // trim(output_folder_name) // '/Code'
        !call execute_command_line(trim(cp_command))
        call system(trim(cp_command))

        domain%output_folder_name = output_folder_name

    end subroutine create_output_files

    !
    ! Write output to files (must call create_output_files first). 
    !
    ! If time_only=.TRUE., only write the time. This is used to avoid writing the main model grids. Typically useful when values at
    ! point-gauges are being recorded, and we want to store the time too, but it would take too much disk to store the model grids
    !
    subroutine write_to_output_files(domain, time_only)
        class(domain_type), intent(inout):: domain
        logical, optional, intent(in) :: time_only
        integer(ip):: i, j
        logical:: to

TIMER_START('fileIO')
        if(present(time_only)) then
            to = time_only
        else
            to = .FALSE.
        end if

       
        if(to .eqv. .FALSE.) then 
#ifdef NONETCDF
            ! Write to home brew binary format
            do i = 1, domain%nvar
                do j = 1, domain%nx(2)
                    ! Binary
                    write(domain%output_variable_unit_number(i)) real(domain%U(:,j,i), output_precision)
                end do
            end do
#else

#ifdef DEBUG_ARRAY
            ! Write to netcdf, with debug_array
            call domain%nc_grid_output%write_grids(domain%time, domain%U, domain%debug_array)    
#else
            ! Write to netcdf (regular case)
            call domain%nc_grid_output%write_grids(domain%time, domain%U)    
#endif

#endif
        end if

        ! Time too, as ascii
        write(domain%output_time_unit_number, *) domain%time

TIMER_STOP('fileIO')

    end subroutine write_to_output_files

    !
    ! Set the domain%logfile_unit to point to an actual file
    !
    subroutine divert_logfile_unit_to_file(domain)
        class(domain_type), intent(inout):: domain
        character(len=charlen):: logfile_name, domain_ID

        write(domain_ID, domain_myid_char_format) domain%myid

        if(domain%output_folder_name == '') then
            logfile_name = './logfile_' // TRIM(domain_ID) // '.log'
        else
            logfile_name = TRIM(domain%output_folder_name) // '/logfile_' // TRIM(domain_ID) // '.log'
        end if

        open(newunit=domain%logfile_unit, file=logfile_name)

    end subroutine

    !
    ! Keep track of the maxima of stage (only) 
    !
    subroutine update_max_quantities(domain)
        class(domain_type), intent(inout):: domain
        integer(ip):: j, k, i

        if(domain%record_max_U) then
            !TIMER_START('update_max_quantities')

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)

            !! Previously we stored all quantities in max_U. However, generally
            !! max uh and max vh don't have much meaning by themselves. So now
            !! we just store max stage
            !!$OMP DO SCHEDULE(GUIDED), COLLAPSE(2)
            !DO k = 1, domain%nvar

            !! UPDATE: Only record max stage
            !$OMP DO SCHEDULE(STATIC)
                do j = domain%yL, domain%yU !1, domain%nx(2)
                    do i = domain%xL, domain%xU !1, domain%nx(1)
                        domain%max_U(i,j,STG) = max(domain%max_U(i,j,STG), domain%U(i,j,STG))
                    end do
                end do
            !END DO
            !$OMP END DO
            !$OMP END PARALLEL
            
            !TIMER_STOP('update_max_quantities')
        end if

    end subroutine update_max_quantities

    !
    ! Write max quantities to a file (usually just called once at the end of a simulation). 
    !
    ! Currently we only write stage, followed by the elevation. Although the latter usually doesn't evolve, we generally want both
    ! the max stage and the elevation for plotting purposes, so it is saved here too.
    !
    subroutine write_max_quantities(domain)
        class(domain_type), intent(in):: domain
        character(len=charlen):: max_quantities_filename
        integer:: i, j, k

        if(domain%record_max_U) then
#ifdef NONETCDF
            ! Use home-brew binary format to output max stage

            max_quantities_filename = trim(domain%output_folder_name) // '/Max_quantities'

            open(newunit=i, file=max_quantities_filename, access='stream', form='unformatted')
        
            !DO k = 1, domain%nvar
            !! Update: Only record max stage
            k = STG
                do j = 1, domain%nx(2)
                    write(i) real(domain%max_U(:,j,k), output_precision)
                end do
            !END DO

            ! Also store elevation since it is typically useful, and we might not store it otherwise
            do j = 1, domain%nx(2)
                write(i) real(domain%U(:,j,ELV), output_precision)
            end do
            
            close(i)
#else
            if(allocated(domain%manning_squared)) then
                call domain%nc_grid_output%store_max_stage(max_stage=domain%max_U(:,:,STG), &
                    elevation=domain%U(:,:,ELV), manning_squared = domain%manning_squared)

            else
                call domain%nc_grid_output%store_max_stage(max_stage=domain%max_U(:,:,STG), &
                    elevation=domain%U(:,:,ELV))
            end if

            
#endif
        end if
 
    end subroutine

    ! Determine whether points x,y are in the priority domain
    !
    !@param domain
    !@param x, y -- real arrays with x, y coordinates
    !@param is_in_priority -- logical array with the same length as x and y
    !
    subroutine is_in_priority_domain(domain, x, y, is_in_priority)
        class(domain_type), intent(in) :: domain 
        real(dp), intent(in) :: x(:), y(:)
        logical, intent(inout) :: is_in_priority(:)

        integer(ip) :: i, i0, j0

        if(domain%nesting%my_index == 0) then
            write(log_output_unit, *) ' Error in is_in_priority_domain: Nesting is not active'
            call generic_stop
        end if

        if(size(x) /= size(y) .or. size(x) /= size(is_in_priority)) then
            write(log_output_unit, *) ' Error in is_in_priority_domain: All inputs must have same length'
            call generic_stop
        end if

        is_in_priority = .false.
        do i = 1, size(x)

            i0 = ceiling( (x(i) - domain%lower_left(1))/domain%dx(1) )
            j0 = ceiling( (y(i) - domain%lower_left(2))/domain%dx(2) )

            if(i0 > 0 .and. i0 <= domain%nx(1) .and. j0 > 0 .and. j0 <= domain%nx(2)) then
                if( domain%nesting%priority_domain_index(i0, j0) == domain%nesting%my_index .and. &
                    domain%nesting%priority_domain_image(i0, j0) == domain%nesting%my_image) then
                   is_in_priority(i) = .true.
                end if
            end if
            
        end do

    end subroutine


    !
    ! Set up point gauges, which record values of variables at cells nearest the given xy points.
    ! 
    ! @param domain The domain within which we will record outputs (in domain%U)
    ! @param xy_coords numeric array of x,y coordinates with 2 rows and as many columns as points. These must all be inside the
    ! extents of the domain
    ! @param time_series_var (optional) array with the indices of U to store at gauges each timestep, Default is [STG, UH, VH] to store
    ! stage, uh, vh
    ! @param static_var (optional) array with the indices of U to store at gauges only once. Default is [ELV] to store elevation
    ! @param gauge_ids (optional) an REAL ID for each gauge. Default gives (1:size(xy_coords(1,:))) * 1.0. Even though integers are
    ! natural for IDS, we store as REAL, just to avoid precision loss issues if someone decides to use large numbers.
    ! @param attribute_names optional character vector with the names of global attributes for the netdf file
    ! @param attribute_values optional character vector with the values for global attributes in the netcdf file (corresponding to
    ! attribute_names)
    !
    subroutine setup_point_gauges(domain, xy_coords, time_series_var, static_var, gauge_ids, &
        attribute_names, attribute_values)

        class(domain_type), intent(inout):: domain
        real(dp), intent(in) :: xy_coords(:,:)
        integer(ip), optional, intent(in):: time_series_var(:), static_var(:)
        real(dp), optional, intent(in):: gauge_ids(:)
        character(charlen), optional, intent(in):: attribute_names(:), attribute_values(:)

        character(charlen) :: netcdf_gauge_output_file, t3
        integer(ip), allocatable:: tsv(:), sv(:)
        real(dp), allocatable:: gauge_ids_local(:)
        real(dp) :: bounding_box(2,2)
        integer(ip):: i, i0, j0
        logical, allocatable :: priority_gauges(:)

        ! Ensure domain was already allocated
        if(.not. allocated(domain%U)) then
            write(domain%logfile_unit,*) 'Error in domain%setup_point_gauges -- the domain must be allocated and have' 
            write(domain%logfile_unit,*) 'elevation etc initialised, before trying to setup the point_gauges' 
            flush(domain%logfile_unit)
            call generic_stop
        end if

        ! Set variables in domain%U(:,:,k) that are stored at gauges every output timestep
        if(present(time_series_var)) then
            allocate(tsv(size(time_series_var)))
            tsv = time_series_var
        else
            ! Default case -- store stage/uh/vh every output step
            allocate(tsv(3))
            tsv = [STG, UH, VH]
        end if

        ! Set variables in domain%U(:,:,k) that are stored only once (at the start)
        if(present(static_var)) then
            allocate(sv(size(static_var)))
            sv = static_var
        else 
            ! Default case -- store elevation once
            allocate(sv(1))
            sv = [ELV]
        end if

        ! Setup or create some 'IDs' for each gauge. These come in handy.
        if(present(gauge_ids)) then
            if(size(gauge_ids) /= size(xy_coords(1,:))) then
                write(domain%logfile_unit,*) 'Number of gauge ids does not equal number of coordinates' 
                flush(domain%logfile_unit)
                call generic_stop
            end if
            allocate(gauge_ids_local(size(gauge_ids)))
            gauge_ids_local = gauge_ids
        else
            ! Default case -- give sequential integer ids
            allocate(gauge_ids_local(size(xy_coords(1,:))))
            do i = 1, size(gauge_ids_local)
                gauge_ids_local(i) = i*ONE_dp
            end do
        end if
    
        ! Get domain id as a character and put it in the output file name
        write(t3, domain_myid_char_format) domain%myid
        netcdf_gauge_output_file = trim(domain%output_folder_name) // '/' // &
            'Gauges_data_ID' // trim(t3) // '.nc'

        ! Pass the bounding box
        bounding_box(1,1:2) = domain%lower_left
        bounding_box(2,1:2) = domain%lower_left + domain%lw

        if(domain%nesting%my_index > 0) then
            ! If we are nesting, tell the model which gauges are in the priority domain
            allocate(priority_gauges(size(xy_coords(1,:))))
            call domain%is_in_priority_domain(xy_coords(1,:), xy_coords(2,:), priority_gauges)
            ! Allocate the gauges
            call domain%point_gauges%allocate_gauges(xy_coords, tsv, sv, gauge_ids_local, &
                bounding_box=bounding_box, priority_gauges = priority_gauges)
            deallocate(priority_gauges)
        else
            ! Not nesting
            ! Allocate the gauges
            call domain%point_gauges%allocate_gauges(xy_coords, tsv, sv, gauge_ids_local, &
                bounding_box=bounding_box)
        end if

        ! Add attributes and initialise the point_gauges object
        if((present(attribute_names)).AND.(present(attribute_values))) then
            call domain%point_gauges%initialise_gauges(domain%lower_left, domain%dx, &
                domain%nx, domain%U, netcdf_gauge_output_file, &
                attribute_names=attribute_names, attribute_values=attribute_values)
        else
            call domain%point_gauges%initialise_gauges(domain%lower_left, domain%dx, &
                domain%nx, domain%U, netcdf_gauge_output_file)
        end if

    end subroutine

    !
    ! Write the values of point gauges to a file
    !
    subroutine write_gauge_time_series(domain)

        class(domain_type), intent(inout):: domain
TIMER_START('write_gauge_time_series')
        if(allocated(domain%point_gauges%time_series_values)) then
            call domain%point_gauges%write_current_time_series(domain%U, domain%time)
        end if
TIMER_STOP('write_gauge_time_series')
    end subroutine

    !
    ! Routine to call once we no longer need the domain. 
    !
    ! One case where this is important is when using netcdf output -- since if the files are not closed, then they may not be
    ! completely written out.
    !
    subroutine finalise_domain(domain)

        class(domain_type), intent(inout):: domain

        integer(ip) :: i
        logical :: is_open

        ! Close the gauges netcdf file -- since otherwise it might not finish
        ! writing.
        call domain%point_gauges%finalise()
    
        ! Flush all open file units -- FIXME here I assume we have only a fixed number. Likely true, but beware if it is not
        ! satisfied.
        do i = 1, 10000
            inquire(i, opened = is_open) 
            if(is_open) call flush(i)
        end do

#ifndef NONETCDF
        ! Close the grids netcdf file
        call domain%nc_grid_output%finalise()
#endif

    end subroutine


    ! If doing nesting, we need to track boundary flux integrals through send/recv regions
    !
    ! This routine multiplies all boundary_flux_integrals in send/recv regions by a constant 'c', which
    ! is required for the flux tracking.
    !
    subroutine nesting_boundary_flux_integral_multiply(domain, c)
    
        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: c

        integer(ip) :: i
    
TIMER_START('nesting_boundary_flux_integral_multiply')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, c)

        ! The 'send communicator' fluxes
        if(allocated(domain%nesting%send_comms)) then
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%send_comms)
                call domain%nesting%send_comms(i)%boundary_flux_integral_multiply(c)
            end do 
            !$OMP END DO
        end if


        ! The 'receive communicator' fluxes
        if(allocated(domain%nesting%recv_comms)) then
            !$OMP DO SCHEDULE (DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)
                call domain%nesting%recv_comms(i)%boundary_flux_integral_multiply(c)
            end do 
            !$OMP END DO
        end if

        !$OMP END PARALLEL

TIMER_STOP('nesting_boundary_flux_integral_multiply')
    end subroutine
    
    ! If doing nesting, we probably want to track boundary flux integrals through send/recv regions -- e.g. to allow for flux
    ! correction
    !
    ! This routine adds "dt * current_value_of_fluxes * dx" to all boundary_flux_integrals in send/recv regions, as required for
    ! flux tracking.
    !
    ! @param domain
    ! @param dt real (timestep)
    ! @param flux_NS rank 3 real array with north-south fluxes
    ! @param flux_NS_lower_index integer. Assume flux_NS(:,1,:) contains the bottom edge flux for cells with j index =
    !   flux_NS_lower_index. For example, "flux_NS_lower_index=1" implies that flux_NS includes fluxes along the models' southern
    !   boundary.
    !   For some of our solvers this is not true [e.g. linear leapfrog, because the 'mass flux' terms are effectively stored in the
    !   domain%U variable].  Hence, it is useful to sometimes have flux_NS_lower_index != 1. In that case, it is assumed we never
    !   try to access the flux_NS value for a location where it doesn't exist.
    ! @param flux_EW rank 3 real array with east-west fluxes
    ! @param flux_EW_lower_index integer. Assume flux_EW(1,:,:) contains the left-edge flux for cells with i index =
    !   flux_EW_lower_index. For example, "flux_EW_lower_index=1" implies that flux_EW includes fluxes right to the boundary. 
    !   For some of our solvers this is not true [e.g. linear leapfrog, because the 'mass flux terms are effectively stored in the
    !   domain%U variable]. 
    !   Hence, it is useful to sometimes have flux_EW_lower_index != 1. In that case, it is assumed we never try to access the
    !   flux_EW value for a location where it doesn't exist.
    ! @param var_indices integer rank1 array of length 2. Gives [lower, upper] indices of the 3rd rank of flux_NS/flux_EW that we
    !   use in the time-stepping
    ! @param flux_already_multiplied_by_dx logical. If TRUE, assume that flux_NS/EW have already been multiplied by their
    !   distance_bottom_edge/ distance_left_edge terms.
    !   Because some solvers store 'flux' and some store 'flux . dx', this allows us to work with whatever fluxes are available,
    !   without manual transformation.
    !
    subroutine nesting_boundary_flux_integral_tstep(domain, dt, &
        flux_NS, flux_NS_lower_index,&
        flux_EW, flux_EW_lower_index,&
        var_indices, flux_already_multiplied_by_dx)
    
        class(domain_type), intent(inout) :: domain
        integer(ip), intent(in) :: flux_NS_lower_index, flux_EW_lower_index, var_indices(2)
        real(dp), intent(in) :: dt, flux_NS(:,:,:), flux_EW(:,:,:)
        logical, intent(in) :: flux_already_multiplied_by_dx

        integer(ip) :: i

!TIMER_START('nesting_boundary_flux_integral_tstep')

        ! Apply to both the send and recv comms. This means that after
        ! communication, we can compare fluxes that were computed with
        ! different numerical methods, and potentially apply flux correction.

        ! Deal with 'send comminicators'
        if(allocated(domain%nesting%send_comms)) then
            do i = 1, size(domain%nesting%send_comms)
                call domain%nesting%send_comms(i)%boundary_flux_integral_tstep( dt,&
                    flux_NS, flux_NS_lower_index, domain%distance_bottom_edge, &
                    flux_EW, flux_EW_lower_index, domain%distance_left_edge, & 
                    var_indices, flux_already_multiplied_by_dx)
            end do 
        end if

        ! Deal with 'receive comminicators'
        if(allocated(domain%nesting%recv_comms)) then
            do i = 1, size(domain%nesting%recv_comms)
                call domain%nesting%recv_comms(i)%boundary_flux_integral_tstep( dt,&
                    flux_NS, flux_NS_lower_index, domain%distance_bottom_edge, &
                    flux_EW, flux_EW_lower_index, domain%distance_left_edge, &
                    var_indices, flux_already_multiplied_by_dx)
            end do
        end if
!TIMER_STOP('nesting_boundary_flux_integral_tstep')

    end subroutine


    ! Apply nesting flux correction all throughout the domain, even at non-priority-domain sites.
    !
    ! This is a key routine for mainintaining mass conservation in the multidomain case.
    !
    ! The general idea is that after we have received from other domains, we can apply flux correction "just like the other domains
    ! would do". This means we avoid having to do multiple nesting communications to make the fluxes consistent.
    !
    ! @param domain the domain
    ! @param all_dx_md rank 3 array with dx(1:2) for all domains in the multidomain
    ! @param all_is_staggered_grid_md rank 2 array with integer recording whether grid is staggered, for all domains in multidomain
    ! @param fraction_of real in [0.0,1.0], apply some fraction of the flux correction. By default apply completely (i.e. 1.0)
    !
    subroutine nesting_flux_correction_everywhere(domain, all_dx_md, all_is_staggered_grid_md, fraction_of)
        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: all_dx_md(:,:,:)
        integer(ip), intent(in) :: all_is_staggered_grid_md(:,:)
        real(dp), optional, intent(in) :: fraction_of

        integer(ip) :: i, j, k, n0, n1, m0, m1, dm, dn, dir, ni, mi
        integer(ip) :: dx_ratio, p0, p1, ii, jj, ic, jc, jL, iL, i1, j1, ind
        integer(ip) :: my_index, my_image, nbr_index, nbr_image
        integer(ip) :: out_index, out_image, sg, ic_last, jc_last
        integer(ip) :: var1, varN, dm_outside, dn_outside, inv_cell_ratios_ip(2), nip, nim, njp, njm
        real(dp) :: fraction_of_local, du_di_p, du_di_m, du_dj_p, du_dj_m, dmin, dmax
        real(dp) :: gradient_scale_x, gradient_scale_y, d_i_jp, d_i_jm, d_ip_j, d_im_j, d_i_j
        real(dp) :: du_di, du_dj, area_scale
        integer, parameter :: NORTH = 1, SOUTH = 2, EAST = 3, WEST = 4
        logical :: equal_sizes

        if(present(fraction_of)) then
            fraction_of_local = fraction_of
        else
            fraction_of_local = 1.0_dp
        end if

        if(send_boundary_flux_data .and. allocated(domain%nesting%recv_comms)) then

TIMER_START('nesting_flux_correction')

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, all_dx_md, all_is_staggered_grid_md, fraction_of_local)

            !
            ! NORTH BOUNDARIES. 
            !
            ! By doing each box direction separately (i.e. all north, then all south, ...), we can ensure that multiple openmp
            ! threads do not try to update the same domain%U cells at once.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)

                my_index = domain%nesting%recv_comms(i)%my_domain_index
                my_image = domain%nesting%recv_comms(i)%my_domain_image_index
                ! Note that nbr_index and nbr_image correspond to the priority domain, because
                ! we always receive from the priority domain
                nbr_index = domain%nesting%recv_comms(i)%neighbour_domain_index
                nbr_image = domain%nesting%recv_comms(i)%neighbour_domain_image_index

                !
                ! North boundary
                !
                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)
                dm_outside = 1
                dir = NORTH

                if(m1 < domain%nx(2) ) then
                    do ni = n0, n1
                        ! Look at the priority domain just to the north
                        out_index = domain%nesting%priority_domain_index(ni, m1+dm_outside)
                        out_image = domain%nesting%priority_domain_image(ni, m1+dm_outside)

                        if(out_index == nbr_index .and. out_image == nbr_image) then
                            ! Do nothing if the send/recv areas are from the same domain
                        else
                            ! Compute offset 'dm' such that 'm1 + dm' is outside or inside the recv box, as appropriate.
                            ! Also compute other quantities required to do the update
                            ! Note either 'dm = dm_outside' or 'dm = 0'. 
                            call compute_offset_inside_or_out(dm, dm_outside, out_index, out_image, nbr_index, nbr_image, &
                                my_index, my_image, is_ew = .false., dx_ratio=dx_ratio, equal_sizes=equal_sizes, &
                                var1 = var1, varN = varN)

                            ! If the area to be 'flux-corrected' has a resolution higher than
                            ! the current domain, then the flux-correction would anyway not have been
                            ! applied to the central cell. So no changes should occur in that case (dx_ratio = 0)
                            if(dx_ratio > 0) then

                                ! Indices to help us 'spread' the correction of >= 1 cell
                                if(dm > 0) then
                                    p0 = m1 + dm
                                    p1 = m1 + dm + dx_ratio - 1
                                    sg = 1
                                else
                                    p0 = m1 + dm - dx_ratio + 1
                                    p1 = m1 + dm
                                    ! If we update 'inside', must flip the sign of the flux correction
                                    sg = -1
                                end if

                                if( .not. (out_index == my_index .and. out_image == my_image)) then
                                    !! If we are correcting a boundary which does not touch the current domain,
                                    !! then there will be a matching north/south edge from 2 different recv_boxes.
                                    !! The corrections will be:
                                    !!     A)   my_domain_flux - nbr_domain_flux
                                    !!    and
                                    !!     B)   my_domain_flux - out_domain_flux
                                    !!
                                    !! We want to correct with:
                                    !!        (nbr_domain_flux - out_domain_flux) = B - A
                                    !! In this case we need to flip the sign of flux_correction 'A' to get the desired
                                    !! correction. The box associated with 'A' occurs when (dm /= 0).
                                    if(dm /= 0) sg = -1 * sg
                                end if

                                if((my_index /= domain%nesting%priority_domain_index(ni, p0)) .or. &
                                   (my_image /= domain%nesting%priority_domain_image(ni, p0)) ) then
                                   sg = -sg
                                end if

                                area_scale = sum(domain%area_cell_y(p0:p1))
                                do k = var1, varN
                                    domain%U(ni, p0:p1, k) = domain%U(ni, p0:p1, k) - &
                                        sg * &
                                        real(domain%nesting%recv_comms(i)%recv_box_flux_error(dir)%x(ni - n0 + 1, k)/&
                                            area_scale, dp) * fraction_of_local
                                end do

                            end if
                        end if
                    end do

                end if

            end do
            !$OMP END DO

            !
            ! SOUTH BOUNDARIES.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)

                my_index = domain%nesting%recv_comms(i)%my_domain_index
                my_image = domain%nesting%recv_comms(i)%my_domain_image_index
                nbr_index = domain%nesting%recv_comms(i)%neighbour_domain_index
                nbr_image = domain%nesting%recv_comms(i)%neighbour_domain_image_index

                !
                ! South boundary
                !
                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)
                dm_outside = -1
                dir = SOUTH

                if(m0 > 1) then
                    do ni = n0, n1
                        ! Look at the priority domain just to the south
                        out_index = domain%nesting%priority_domain_index(ni, m0+dm_outside)
                        out_image = domain%nesting%priority_domain_image(ni, m0+dm_outside)

                        if(out_index == nbr_index .and. out_image == nbr_image) then
                            ! Do nothing if the send/recv areas are from the same domain,
                        else
                            !
                            ! For comments, see the NORTH case above
                            !

                            call compute_offset_inside_or_out(dm, dm_outside, out_index, out_image, nbr_index, nbr_image, &
                                my_index, my_image, is_ew = .false., dx_ratio=dx_ratio, equal_sizes=equal_sizes, &
                                var1=var1, varN=varN)

                            if(dx_ratio > 0) then
   
                                if(dm < 0) then
                                    p0 = m0 + dm - dx_ratio + 1
                                    p1 = m0 + dm
                                    sg = 1
                                else
                                    p0 = m0 + dm
                                    p1 = m0 + dm + dx_ratio - 1
                                    sg = -1
                                end if

                                if( .not. (out_index == my_index .and. out_image == my_image)) then
                                    if(dm /= 0) sg = -1 * sg
                                end if

                                if((my_index /= domain%nesting%priority_domain_index(ni, p0)) .or. &
                                   (my_image /= domain%nesting%priority_domain_image(ni, p0)) ) then
                                   sg = -sg
                                end if

                                area_scale = sum(domain%area_cell_y(p0:p1))
                                do k = var1, varN
                                    domain%U(ni, p0:p1, k) = domain%U(ni, p0:p1, k) + &
                                        sg * &
                                        real(domain%nesting%recv_comms(i)%recv_box_flux_error(dir)%x(ni - n0 + 1, k)/&
                                            area_scale, dp) * fraction_of_local
                                end do
                            end if
                        end if
                    end do

                end if
            end do
            !$OMP END DO

            !
            ! EAST BOUNDARIES.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)
                my_index = domain%nesting%recv_comms(i)%my_domain_index
                my_image = domain%nesting%recv_comms(i)%my_domain_image_index
                nbr_index = domain%nesting%recv_comms(i)%neighbour_domain_index
                nbr_image = domain%nesting%recv_comms(i)%neighbour_domain_image_index

                !
                ! East boundary
                !
                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)
                dn_outside = 1
                dir = EAST

                if(n1 < domain%nx(1)) then
                    do mi = m0, m1
                        ! Look at the priority domain just to the east
                        out_index = domain%nesting%priority_domain_index(n1 + dn_outside, mi)
                        out_image = domain%nesting%priority_domain_image(n1 + dn_outside, mi)

                        if(out_index == nbr_index .and. out_image == nbr_image) then
                            ! Do nothing if the send/recv areas are from the same domain
                        else
                            ! Compute offset 'dn' such that 'n1 + dn' is outside or inside the recv box, as appropriate
                            ! Also compute other quantities required to do the update
                            ! Note either 'dn = dn_outside' or 'dn = 0'.
                            call compute_offset_inside_or_out(dn, dn_outside, out_index, out_image, nbr_index, nbr_image, &
                                my_index, my_image, is_ew = .true., dx_ratio=dx_ratio, equal_sizes=equal_sizes, &
                                var1=var1, varN=varN)

                            ! If the area to be 'flux-corrected' has a finer resolution than
                            ! the current domain, then the flux-correction would anyway not have been
                            ! applied to the central cell. So no changes should occur in that case (dx_ratio = 0)
                            if(dx_ratio > 0) then
                                ! Indices to help us 'spread' the correction of >= 1 cell
                                if(dn > 0) then
                                    p0 = n1 + dn
                                    p1 = n1 + dn + dx_ratio - 1
                                    sg = 1
                                else
                                    p0 = n1 + dn - dx_ratio + 1
                                    p1 = n1 + dn
                                    ! If we update 'inside', must flip the sign of the flux correction
                                    sg = -1
                                end if

                                if( .not. (out_index == my_index .and. out_image == my_image)) then
                                    !! If we are correcting a boundary which does not touch the current domain,
                                    !! then there will be a matching east/west edge from 2 different recv_boxes.
                                    !! The corrections will be:
                                    !!     A)   my_domain_flux - nbr_domain_flux
                                    !!    and
                                    !!     B)   my_domain_flux - out_domain_flux
                                    !!
                                    !! We want to correct with:
                                    !!        (nbr_domain_flux - out_domain_flux) = B - A
                                    !! In this case we need to flip the sign of flux_correction 'A' to get the desired
                                    !! correction. The box associated with 'A' occurs when (dn /= 0).
                                    if(dn /= 0) sg = -1 * sg
                                end if

                                if((my_index /= domain%nesting%priority_domain_index(p0, mi)) .or. &
                                   (my_image /= domain%nesting%priority_domain_image(p0, mi)) ) then
                                   sg = -sg
                                end if

                                area_scale = (domain%area_cell_y(mi)*dx_ratio)
                                do k = var1, varN
                                    domain%U(p0:p1, mi, k) = domain%U(p0:p1, mi, k) - &
                                        sg * &
                                        real(domain%nesting%recv_comms(i)%recv_box_flux_error(dir)%x(mi - m0 + 1, k)/&
                                             area_scale, dp) * fraction_of_local
                                end do
                            end if
                        end if
                    end do

                end if
            end do
            !$OMP END DO

            !
            ! WEST BOUNDARIES.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)
                my_index = domain%nesting%recv_comms(i)%my_domain_index
                my_image = domain%nesting%recv_comms(i)%my_domain_image_index
                nbr_index = domain%nesting%recv_comms(i)%neighbour_domain_index
                nbr_image = domain%nesting%recv_comms(i)%neighbour_domain_image_index
                !
                !
                ! West boundary
                !
                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)
                dn_outside = -1
                dir = WEST

                if(n0 > 1) then
                    do mi = m0, m1
                        ! Look at the priority domain just to the west
                        out_index = domain%nesting%priority_domain_index(n0 + dn_outside, mi)
                        out_image = domain%nesting%priority_domain_image(n0 + dn_outside, mi)

                        if(out_index == nbr_index .and. out_image == nbr_image) then
                            ! Do nothing if the send/recv areas are from the same domain
                        else
                            !
                            ! For comments, see 'EAST' case above.
                            !
                            call compute_offset_inside_or_out(dn, dn_outside, out_index, out_image, nbr_index, nbr_image, &
                                my_index, my_image, is_ew=.true., dx_ratio=dx_ratio, equal_sizes=equal_sizes, &
                                var1=var1, varN=varN)

                            if(dx_ratio > 0) then

                                if(dn < 0) then
                                    p0 = n0 + dn - dx_ratio + 1
                                    p1 = n0 + dn
                                    sg = 1
                                else
                                    p0 = n0 + dn
                                    p1 = n0 + dn + dx_ratio - 1
                                    sg = -1
                                end if

                                if( .not. (out_index == my_index .and. out_image == my_image)) then
                                    if(dn /= 0) sg = -1 * sg
                                end if

                                if((my_index /= domain%nesting%priority_domain_index(p0, mi)) .or. &
                                   (my_image /= domain%nesting%priority_domain_image(p0, mi)) ) then
                                   sg = -sg
                                end if

                                area_scale = (domain%area_cell_y(mi)*dx_ratio)
                                do k = var1, varN
                                    domain%U(p0:p1, mi, k) = domain%U(p0:p1, mi, k) + &
                                        sg * &
                                        real(domain%nesting%recv_comms(i)%recv_box_flux_error(dir)%x(mi - m0 + 1, k)/&
                                            area_scale, dp) * fraction_of_local
                                end do
                            end if
                        end if
                    end do

                end if
            end do
            !$OMP END DO

            !
            ! Later we do some interpolation inside some recv boxes.
            ! This may require the use of cell values outside the receive box, which hypothetically might also be being updated.
            ! To avoid any issues, we copy such data here, using 'recv_box_flux_error' as a scratch space.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)

                ! If my domain is not finer, we do not need to do anything
                if(.not. domain%nesting%recv_comms(i)%my_domain_is_finer) cycle

                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)

                inv_cell_ratios_ip = nint(1.0_dp/domain%nesting%recv_comms(i)%cell_ratios)

                !
                ! Copy values of U north of the box, using recv_box_flux_error as a scratch space
                ! 
                ind = min(m1 + 1 + inv_cell_ratios_ip(2)/2, domain%nx(2))
                domain%nesting%recv_comms(i)%recv_box_flux_error(NORTH)%x(1:(n1-n0+1),1:4) = domain%U(n0:n1, ind, 1:4)
                
                !
                ! Copy values of U south of the box, using recv_box_flux_error as a scratch space
                ! 
                ind = max(m0 - 1 - inv_cell_ratios_ip(2)/2, 1)
                domain%nesting%recv_comms(i)%recv_box_flux_error(SOUTH)%x(1:(n1-n0+1),1:4) = domain%U(n0:n1, ind, 1:4)

                !
                ! Copy values of U east of the box, using recv_box_flux_error as a scratch space
                ! 
                ind = min(n1 + 1 + inv_cell_ratios_ip(1)/2, domain%nx(1))
                domain%nesting%recv_comms(i)%recv_box_flux_error(EAST)%x(1:(m1-m0+1),1:4) = domain%U(ind, m0:m1, 1:4)

                !
                ! Copy values of U west of the box, using recv_box_flux_error as a scratch space
                !
                ind = max(n0 - 1 - inv_cell_ratios_ip(1)/2, 1)
                domain%nesting%recv_comms(i)%recv_box_flux_error(WEST)%x(1:(m1-m0+1),1:4) = domain%U(ind, m0:m1, 1:4)

            end do
            !$OMP END DO


            !
            ! At this point we do the interpolation
            !

            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms)

                ! If my domain is not finer, we do not need to do anything
                if(.not. domain%nesting%recv_comms(i)%my_domain_is_finer) cycle


                n0 = domain%nesting%recv_comms(i)%recv_inds(1,1)
                n1 = domain%nesting%recv_comms(i)%recv_inds(2,1)
                m0 = domain%nesting%recv_comms(i)%recv_inds(1,2)
                m1 = domain%nesting%recv_comms(i)%recv_inds(2,2)

                inv_cell_ratios_ip = nint(1.0_dp/domain%nesting%recv_comms(i)%cell_ratios)

                ! Interpolation
                ! Loop over 'coarse parent grid' cells by taking steps of inv_cell_ratios_ip
                do jL = m0, (m1 - inv_cell_ratios_ip(2) + 1), inv_cell_ratios_ip(2)
                    do iL = n0, (n1 - inv_cell_ratios_ip(1) + 1), inv_cell_ratios_ip(1)

                        ! Get the 'central cell col index'. Deliberate integer divisions
                        jc = jL + inv_cell_ratios_ip(2)/2
                        ! Get the 'central cell row index'. Deliberate integer divisions
                        ic = iL + inv_cell_ratios_ip(1)/2

                        ! Positive/negative derivative indices, with out-of-bounds protection
                        nip = min(ic + inv_cell_ratios_ip(1), domain%nx(1))
                        nim = max(ic - inv_cell_ratios_ip(1), 1)
                        njp = min(jc + inv_cell_ratios_ip(2), domain%nx(2))
                        njm = max(jc - inv_cell_ratios_ip(2), 1)

                        !
                        ! Below, compute depth change on the 'coarse grid' we receive from,
                        ! and suppress gradients where depths are changing rapidly
                        !

                        d_i_j = domain%U(ic, jc, STG) - domain%U(ic, jc, ELV)
                        ! depth at i+, j
                        if(nip < n1) then
                            d_ip_j = domain%U(nip, jc, STG) - domain%U(nip, jc, ELV)
                        else
                            ! Use the copied version to ensure we don't access cells with in-process updates
                            d_ip_j = domain%nesting%recv_comms(i)%recv_box_flux_error(EAST)%x(jc - m0 + 1, STG) - &
                                     domain%nesting%recv_comms(i)%recv_box_flux_error(EAST)%x(jc - m0 + 1, ELV)
                        end if

                        ! depth at i-, j
                        if(nim > n0) then
                            d_im_j = domain%U(nim, jc, STG) - domain%U(nim, jc, ELV)
                        else
                            ! Use the copied version to ensure we don't access cells with in-process updates
                            d_im_j = domain%nesting%recv_comms(i)%recv_box_flux_error(WEST)%x(jc - m0 + 1, STG) - &
                                     domain%nesting%recv_comms(i)%recv_box_flux_error(WEST)%x(jc - m0 + 1, ELV)
                        end if

                        ! Limit on x-gradient
                        dmax = max(d_i_j, max(d_ip_j, d_im_j))
                        dmin = min(d_i_j, min(d_ip_j, d_im_j))
                        if( dmax <= minimum_allowed_depth * 100.0_dp) then
                            gradient_scale_x = 0.0_dp
                        else
                            ! Smooth transition from 1.0 to 0.0 as dmin/dmax drops from 0.25 to 0.1
                            gradient_scale_x = &
                                max(0.0_dp, (dmin/dmax - 0.1_dp))/(0.25_dp - 0.1_dp)
                            gradient_scale_x = min(1.0_dp, gradient_scale_x)
                        end if

                        ! depth at i, j+
                        if(njp < m1) then
                            d_i_jp = domain%U(ic, njp, STG) - domain%U(ic, njp, ELV)
                        else
                            ! Use the copied version to ensure we don't access cells with in-process updates
                            d_i_jp = domain%nesting%recv_comms(i)%recv_box_flux_error(NORTH)%x(ic - n0 + 1, STG) - &
                                     domain%nesting%recv_comms(i)%recv_box_flux_error(NORTH)%x(ic - n0 + 1, ELV)
                        end if

                        ! depth at i, j-
                        if(njm > m0) then
                            d_i_jm = domain%U(ic, njm, STG) - domain%U(ic, njm, ELV)
                        else
                            ! Use the copied version to ensure we don't access cells with in-process updates
                            d_i_jm = domain%nesting%recv_comms(i)%recv_box_flux_error(SOUTH)%x(ic - n0 + 1, STG) - &
                                     domain%nesting%recv_comms(i)%recv_box_flux_error(SOUTH)%x(ic - n0 + 1, ELV)
                        end if

                        ! Limit on y-gradient
                        dmax = max(d_i_j, max(d_i_jp, d_i_jm))
                        dmin = min(d_i_j, min(d_i_jp, d_i_jm))
                        if( dmax <= minimum_allowed_depth * 100.0_dp) then
                            gradient_scale_y = 0.0_dp
                        else
                            ! Smooth transition from 1.0 to 0.0 as dmin/dmax drops from 0.25 to 0.1
                            gradient_scale_y = &
                                max(0.0_dp, (dmin/dmax - 0.1_dp))/(0.25_dp - 0.1_dp)
                            gradient_scale_y = min(1.0_dp, gradient_scale_y)
                        end if

                        !
                        ! Interpolation
                        !
                        do k = STG, ELV

                            ! Firstly compute interpolation gradients.
                            ! The gradients should only use U values at the center of the coarser cells we receive from
                            ! (i.e. U(ic, jc, ...) and index offsets by inv_cell_ratio_ip). 
                            ! Those values should not be changed by the operations below.

                            ! Gradient to the west
                            if(nim > n0) then
                                du_di_m = (domain%U(ic, jc, k) - domain%U(nim, jc, k))/max((ic - nim), 1)
                            else
                                ! Use the copied version to ensure we don't access cells with in-process updates
                                du_di_m = (domain%U(ic, jc, k) - &
                                     domain%nesting%recv_comms(i)%recv_box_flux_error(WEST)%x(jc - m0 + 1, k)) / &
                                     max((ic - nim), 1)
                            end if

                            ! Gradient to the east
                            if(nip < n1) then
                                du_di_p = (domain%U(nip, jc, k) - domain%U(ic, jc, k))/max((nip - ic), 1)
                            else
                                ! Use the copied version to ensure we don't access cells with in-process updates
                                du_di_p = (domain%nesting%recv_comms(i)%recv_box_flux_error(EAST)%x(jc - m0 + 1, k) - &
                                    domain%U(ic, jc, k))/max((nip - ic), 1)
                            end if

                            ! Gradient to the south
                            if(njm > m0) then
                                du_dj_m = (domain%U(ic, jc, k) - domain%U(ic, njm, k))/max((jc - njm), 1)
                            else
                                ! Use the copied version to ensure we don't access cells with in-process updates
                                du_dj_m = (domain%U(ic, jc, k) - &
                                    domain%nesting%recv_comms(i)%recv_box_flux_error(SOUTH)%x(ic - n0 + 1, k)) / &
                                    max((jc - njm), 1)
                            end if

                            ! Gradient to the north
                            if(njp < m1) then
                                du_dj_p = (domain%U(ic, njp, k) - domain%U(ic, jc, k))/max((njp - jc), 1)
                            else
                                ! Use the copied version to ensure we don't access cells with in-process updates
                                du_dj_p = (domain%nesting%recv_comms(i)%recv_box_flux_error(NORTH)%x(ic - n0 + 1, k) - &
                                    domain%U(ic, jc, k))/max((njp - jc), 1)
                            end if

                            !! Loop over the 'fine grid' cells (ii,jj) inside each coarse grid cell (ic,jc)
                            do j1 = 1, inv_cell_ratios_ip(2)

                                jj = jL - 1 + j1

                                do i1 = 1, inv_cell_ratios_ip(1)

                                    ii = iL - 1 + i1

                                    ! x gradient
                                    if(ii < ic) then
                                        du_di = du_di_m
                                    else
                                        du_di = du_di_p
                                    end if

                                    ! y gradient
                                    if(jj < jc) then
                                        du_dj = du_dj_m
                                    else
                                        du_dj = du_dj_p
                                    end if

                                    domain%U(ii,jj,k) = domain%U(ic, jc, k) + &
                                        (ii-ic) * du_di * gradient_scale_x  + &
                                        (jj-jc) * du_dj * gradient_scale_y

                                end do
                            end do
                        end do

                    end do
                end do
                
            end do
            !$OMP END DO

            !$OMP END PARALLEL

TIMER_STOP('nesting_flux_correction')

        end if

        contains
            
            ! The flux correction should either be applied inside or outside of the recv box boundary, depending
            ! on whether the priority domain in the recv box is coarser (inside) or finer (outside) than the priority
            ! domain just outside. 
            !
            ! This can be represented by changing the 'index perturbation' to 0 (inside) or +1 (outside, north/east boundary)
            ! or -1 (outside, west/south boundary)
            !
            ! Supposing the index-offset for "outside" is dm_outside, this routine makes 'dm' have either the
            ! value of 'dm_outside' (for outside) or 0 (for inside), by checking the resolutions
            !
            ! @param dm resultant perturbation
            ! @param dm_outside +1 if north/east boundary, -1 if west/south boundary
            ! @param out_index, out_image the index/image of the priority domain just outside the boundary
            ! @param nbr_index, nbr_image the index/image of the priority domain inside the boundary (i.e. the one we receive from)
            ! @param is_ew logical, .true. if comparing EW directions, .false. if comparing NS directions
            ! @param dx_ratio integer ratio of 'my cell size' with the 'to-be-corrected cell size' in either the EW or NS directions
            ! (depending on value of is_ew)
            ! @param equal_sizes report whether the nbr/out cells are of equal size. Note this is not the same as dx=1,
            ! because that refers to the 'my/out' cells, not 'nbr/out' cells.
            ! @param var1 index of first variable to update (normally STG=1)
            ! @param varN index of last variable to update (normally VH=3, except for staggered grid, where it is STG=1)
            subroutine compute_offset_inside_or_out(dm, dm_outside, out_index, out_image, nbr_index, nbr_image, &
                    my_index, my_image, is_ew, dx_ratio, equal_sizes, var1, varN)

                integer(ip), intent(out) :: dm
                integer(ip), intent(in) :: dm_outside, out_index, out_image, nbr_index, nbr_image, my_index, my_image
                logical, intent(in) :: is_ew
                integer(ip), intent(out) :: dx_ratio
                logical, intent(out) :: equal_sizes
                integer(ip), intent(out) :: var1, varN

                integer(ip) :: i1, cor_index, cor_image
                real(dp) :: area_out, area_nbr

                ! If is_ew, then compare dx in the x direction. Otherwise compare dx in the y direction
                if(is_ew) then
                    i1 = 1
                else
                    i1 = 2
                end if

                equal_sizes = .false.

                area_out = all_dx_md(1, out_index, out_image) * all_dx_md(2, out_index, out_image)
                area_nbr = all_dx_md(1, nbr_index, nbr_image) * all_dx_md(2, nbr_index, nbr_image)

                ! Figure out if the priority domain in the recv box is finer/coarser/same
                ! than the priority domain out of the recv box
                if(area_out > 1.25_dp*area_nbr) then
                    ! outside priority domain is coarser. Note the '1.25' factor is a kludge -- the dx
                    ! values should always differ by at least a factor of 2, unless they are identical. 
                    ! The 1.25 factor simply provides protection against round-off
                    dm = dm_outside
                else if (area_out < 0.8_dp*area_nbr) then
                    ! outside priority domain is finer. The '0.8' factor is a kludge to protect
                    ! against round-off. See above.
                    dm = 0
                else
                    ! Both 'nbr' and 'out' have the same cell size. 
                    equal_sizes = .true.
                    ! Need a rule to decide which one to correct
                    ! Select dm based on the order of the index, with tie-breaking by image
                    if(out_index > nbr_index) then
                        dm = dm_outside
                    else
                        if(out_index < nbr_index) then
                            dm = 0
                        else if(out_index == nbr_index) then
                            ! "Corner" case, decide based on the image 
                            if(out_image > nbr_image) then
                                dm = dm_outside
                            else
                                dm = 0
                            end if
                        end if
                    end if
                end if

                ! Get indices for 'the domain to be corrected'  = cor_index, cor_image
                if(dm == dm_outside) then
                    cor_image = out_image
                    cor_index = out_index
                else
                    cor_image = nbr_image
                    cor_index = nbr_index
                end if
                
                ! Get dx ratio of 'domain to be corrected' vs 'my domain', with round-off protection as above
                if(all_dx_md(i1, cor_index, cor_image) > 0.8_dp * all_dx_md(i1, my_index, my_image)) then
                    ! Note the factor 0.8 is a kludge to protect against round-off. The cell-size ratio should either
                    ! be an integer, or 1/integer
                    dx_ratio = nint(all_dx_md(i1, cor_index, cor_image)/all_dx_md(i1, my_index, my_image))
                else
                    dx_ratio = 0
                end if

                !! Define variables to adjust
                var1 = STG ! Always adjust stage, for mass conservation.
                if(all_is_staggered_grid_md(cor_index, cor_image) == 1) then
                    ! When correcting linear domains (i.e. on the staggerd grid), we do not correct momentum fluxes.
                    varN = STG
                else
                    ! When correcting properly nonlinear domains (i.e. without staggered grid), we correct momentum fluxes.
                    varN = VH
                end if

            end subroutine 

    end subroutine

    ! Make elevation constant in nesting send_regions that go to a single coarser cell, if the maximum elevation is above
    ! elevation_threshold
    ! 
    ! This was done (in the past) to avoid wet-dry instabilities, caused by aggregating over wet-and-dry cells on a finer domain,
    ! which is then sent to a coarser domain.  Such an operation will break the hydrostatic balance, unless the elevation in the
    ! fine cells is constant
    !
    ! NOTE: Instead of using this, a better approach is to send the data from the centre cell to the coarse grid. That way we do not
    ! need to hack the elevation data to avoid these instabilities. With the latter approach, this routine seems defunct.
    !
    subroutine use_constant_wetdry_send_elevation(domain, elevation_threshold)

        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: elevation_threshold

        integer(ip) :: i, ic, jc
        integer(ip) :: send_inds(2,3), cr(2)
        real(dp) :: mean_elev, max_elev

        ! Move on if we do not send data
        if(.not. allocated(domain%nesting%send_comms) ) return

        ! Loop over all 'send' regions
        do i = 1, size(domain%nesting%send_comms)

            ! Only operate on finer domains
            if(.not. domain%nesting%send_comms(i)%my_domain_is_finer) cycle
        
            ! [num_x, num_y] cells per coarse domain cell
            cr = nint(ONE_dp/domain%nesting%send_comms(i)%cell_ratios)

            send_inds = domain%nesting%send_comms(i)%send_inds

            ! Loop over j cells that are sent, in steps corresponding to
            ! the coarse j cells
            do jc = send_inds(1,2) - 1, send_inds(2,2) - cr(2), cr(2)

                ! Loop over i cells that are sent, in steps corresponding
                ! to the coarse i cells
                do ic = send_inds(1,1) - 1, send_inds(2,1) - cr(1), cr(1)

                    ! Maximum elevation of cells which will be aggregated
                    ! and sent to the coarser domain 
                    max_elev = maxval( &
                        domain%U((ic+1):(ic+cr(1)), (jc+1):(jc+cr(2)), ELV) )
                    ! Mean elevation of cells which will be aggregated and
                    ! sent to the coarser domain
                    mean_elev = sum( &
                        domain%U((ic+1):(ic+cr(1)), (jc+1):(jc+cr(2)), ELV) )/&
                        (product(cr))

                    ! Replace all elevation values with well behaved solution, if
                    ! (max_elevation > threshold)
                    if(max_elev > elevation_threshold) then
                        domain%U((ic+1):(ic+cr(1)), (jc+1):(jc+cr(2)), ELV) = mean_elev
                    end if

                end do
            end do

        end do


    end subroutine

    !
    ! Convenience routine for setting up nested domain geometries.
    !
    ! Set the domain "lower-left, lw, dx" etc, in a way that ensures the domain can nest with its parent domain, while having an
    ! extent and resolution close to the desired values.
    !
    ! @param domain A domain_type which is unallocated, and does not have lw, dx, nx set
    ! @param parent_domain Another domain_type which DOES have lw, dx, nx set. We want to
    !   give the new domain boundaries that are consistent with the parent domain (for nesting purposes), 
    !   i.e. the corners of each parent domain cell correspond to corners of child domain cells
    ! @param lower_left The desired lower left of the new domain. We will change this for consistency with nesting
    ! @param upper_right The desired upper right of the new domain. We will change this for consistency with nesting
    ! @param dx_refinement_factor Integer such that parent_domain%dx/dx_refinement_factor = new_domain%dx
    ! @param timestepping refinement factor How many time-steps should the new domain take, for each global time-step
    ! @param rounding method optional character controlling how we adjust lower-left/upper-right. 
    !   If rounding_method = 'expand' (DEFAULT), then we adjust the new domain lower-left/upper-right so that the provided
    !   lower-left/upper-right are definitely contained in the new domain. If rounding_method = 'nearest', we move
    !   lower-left/upper-right onto the nearest cell corner of the parent domain. This can be preferable if we want to have multiple
    !   child domains which share boundaries with each other -- but does not ensure the provided lower-left/upper-right are within
    !   the new domain
    ! @param recursive_nesting optional logical. If TRUE(default), the domain's dx_refinement_factor will be multiplied
    !   by its parent domain's dx_refinement_factor before storing in domain%dx_refinement_facor. This will not affect 
    !   domain%dx. But it should should allow the domain to communicate 'cleanly' with "parents of its parent" if it is
    !   partitioned, SO LONG AS the lower_left is also on a corner of the earlier generation's domain. If .FALSE., then just
    !   use the provided dx_refinement_factor, and tell the domain that it's parent-domain's dx_refinement_factor=1.
    !   This is less likely to allow communicating with grandparent domains, but might allow for a more efficient split-up
    !   of the domain.
    !
    subroutine match_geometry_to_parent(domain, parent_domain, lower_left, &
        upper_right, dx_refinement_factor, timestepping_refinement_factor, &
        rounding_method, recursive_nesting)
    
        class(domain_type), intent(inout) :: domain
        class(domain_type), intent(in) :: parent_domain  
        real(dp), intent(in) :: lower_left(2), upper_right(2)
        integer(ip), intent(in) :: dx_refinement_factor, timestepping_refinement_factor
        character(*), intent(in), optional :: rounding_method
        logical, intent(in), optional :: recursive_nesting

        real(dp) :: ur(2), parent_domain_dx(2)
        character(len=charlen) :: rounding
        logical :: recursive_nest

        if(present(rounding_method)) then
            rounding = rounding_method
        else
            rounding = 'expand'
        end if

        if(present(recursive_nesting)) then
            recursive_nest = recursive_nesting
        else
            recursive_nest = .true.
        end if

        ! Check the parent_domain was setup ok
        if(any(parent_domain%nx <= 0) .or. any(parent_domain%lw <= 0) .or. &
           any(parent_domain%lower_left == HUGE(1.0_dp))) then
            write(log_output_unit, *) 'Parent domain must have lw, lower_left and nx already set' 
            call generic_stop
        end if

        ! Parent domain's dx might not have been defined
        parent_domain_dx = parent_domain%lw/parent_domain%nx

        select case (rounding)
        case ('expand')
            ! Ensure that after rounding, the originally requested domain is contained in the final domain

            !  domain%lower_left is on a cell boundary of the parent domain -- and is further 'west/south' than 'lower_left'
            domain%lower_left = parent_domain%lower_left + &
                floor((lower_left - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx 

            ! upper_right = (domain%lower_left + domain%lw) is on a cell boundary of the parent domain
            ur = parent_domain%lower_left + &
                ceiling((upper_right - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

        case('nearest')
            ! Find the 'nearest' match in parent domain. This might mean we reduce the requested size of the domain
            domain%lower_left = parent_domain%lower_left + &
                nint((lower_left - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx 

            ur = parent_domain%lower_left + &
                nint((upper_right - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

        case default

            write(domain%logfile_unit, *) ' rounding_method ', TRIM(rounding), ' not recognized '
            call generic_stop()

        end select

        ! Now we can set the child domain's properties
        domain%lw =  ur - domain%lower_left
        domain%dx = parent_domain_dx/dx_refinement_factor
        domain%nx = nint(domain%lw / domain%dx) ! Is a multiple of dx_refinement_factor
        domain%timestepping_refinement_factor = timestepping_refinement_factor

        if(recursive_nest) then
            domain%dx_refinement_factor = real(dx_refinement_factor * parent_domain%dx_refinement_factor, dp)
            domain%dx_refinement_factor_of_parent_domain = real(parent_domain%dx_refinement_factor, dp)
        else
            domain%dx_refinement_factor = real(dx_refinement_factor, dp)
            domain%dx_refinement_factor_of_parent_domain = ONE_dp
        end if

    end subroutine

    !
    ! Basic smooth of elevation, applied along x/y axes separately
    !
    ! @param method is an optional character controlling the smoothing type.
    !        values are '9pt_average'
    subroutine smooth_elevation(domain, smooth_method)
        class(domain_type), intent(inout) :: domain
        character(*), optional, intent(in) :: smooth_method
    
        integer(ip) :: i, j
        real(dp), allocatable:: elev_block(:,:)
        character(len=charlen) :: method

        if(present(smooth_method)) then
            method = smooth_method
        else
            method = '9pt_average'
        end if

        select case(method)
        case('9pt_average')
            ! 3-point average along y, followed by 3-point average along x

            ! Smooth bathymetry along 'i' axis
            do j = 2, domain%nx(2) - 1
                domain%U(2:(domain%nx(1)-1),j,ELV) = (1.0_dp/3.0_dp)*(&
                    domain%U(2:(domain%nx(1)-1),j,ELV) + &
                    domain%U(1:(domain%nx(1)-2),j,ELV) + &
                    domain%U(3:(domain%nx(1)-0),j,ELV) )
            end do
            ! Smooth bathymetry along 'j' axis
            do i = 2, domain%nx(1) - 1
                domain%U(i, 2:(domain%nx(2)-1),ELV) = (1.0_dp/3.0_dp)*(&
                    domain%U(i, 2:(domain%nx(2)-1),ELV) + &
                    domain%U(i, 1:(domain%nx(2)-2),ELV) + &
                    domain%U(i, 3:(domain%nx(2)-0),ELV) )
            end do
        case default
            stop 'Smoothing method not recognized'
        end select

    end subroutine

    !
    ! Very limited testing of the domain routines.
    ! More important are the validation tests.
    !
    subroutine test_domain_mod
        type(domain_type) :: domain
    
        integer(ip), parameter :: N = 10
        real(dp) :: U(N), U_lower(N), U_upper(N)
        real(dp) :: theta(N), extrapolation_sign(N)
        real(dp) :: desired_answer(N), answer(N)
        integer(ip) :: i

        theta = 1.0_dp
        U = (/(i*1.0_dp - 5.5_dp, i = 1, N)/)
        U_lower = U - 1.0_dp
        U_upper = U + 1.0_dp

        ! Basic extrapolation tests, positive side
        extrapolation_sign=1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * (U + U_upper)
        call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side, reduced theta
        theta = 0.5_dp
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * U + 0.5_dp*0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side, zero theta
        theta = 0.0_dp
        extrapolation_sign = 1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 1.0_dp * U + 0.0_dp*0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Local min.
        theta = 4.0_dp
        extrapolation_sign = 1.0_dp
        U_lower = U + 2.0_dp
        U_upper = U + 3.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U 
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Local max.
        theta = 4.0_dp
        U_lower = U - 2.0_dp
        U_upper = U - 3.0_dp
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U 
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Limited gradients
        ! 
        ! Central gradient = 2, lower-gradient=1, upper-gradient=3
        U_lower = U - 1.0_dp
        U_upper = U + 3.0_dp
        extrapolation_sign = 1.0_dp
        theta = 1.0_dp ! Limited-gradient = 1
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U + 0.5_dp * 1.0_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        theta = 1.5_dp ! Limited-gradient = 1.5
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U + 0.5_dp * 1.5_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        
        extrapolation_sign = -1.0_dp ! As above, negative side
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U - 0.5_dp * 1.5_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        contains

            subroutine assert_equal_within_tol(desired_answer, answer, line)
                real(dp) :: desired_answer(N), answer(N)
                integer :: line

                real(dp) :: tol

                tol = spacing(1.0_dp) * 100

                if(all(abs(desired_answer - answer) < tol)) then
                    print*, 'PASS'
                else
                    print*, 'FAIL', line, desired_answer - answer
                end if

            end subroutine

    end subroutine

end module domain_mod
