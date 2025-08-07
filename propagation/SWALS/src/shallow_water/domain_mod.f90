! Compile with -DTIMER to add timing to the code
#ifdef TIMER
#   define TIMER_START(tname) call domain%timer%timer_start(tname)
#   define TIMER_STOP(tname)  call domain%timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

! Compile with -DEVOLVE_TIMER to add detailed timing of the evolve loop
#ifdef EVOLVE_TIMER
#   define EVOLVE_TIMER_START(tname) call domain%evolve_timer%timer_start(tname)
#   define EVOLVE_TIMER_STOP(tname)  call domain%evolve_timer%timer_end(tname)
#else
#   define EVOLVE_TIMER_START(tname)
#   define EVOLVE_TIMER_STOP(tname)
#endif

module domain_mod
    !!
    !! Contains type 'domain_type', which is the main type for solving the
    !! linear/non-linear shallow water equations on a single structured grid.
    !!

    use global_mod, only: dp, ip, charlen, output_precision, &
        force_double, force_long_double, long_long_ip, &
        maximum_timestep, gravity, &
        advection_beta, &
        minimum_allowed_depth, &
        default_nonlinear_timestepping_method, &
        wall_elevation, &
        default_output_folder, &
        send_boundary_flux_data
    use timer_mod, only: timer_type
    use timestepping_metadata_mod, only: timestepping_metadata, &
        setup_timestepping_metadata, timestepping_method_index
    use point_gauge_mod, only: point_gauge_type
    use nested_grid_comms_mod, only: domain_nesting_type, process_received_data
    use coarray_point2point_comms_mod, only: p2p_comms_type
    use stop_mod, only: generic_stop
    use iso_fortran_env, only: output_unit, int32, int64
    use iso_c_binding, only: c_ptr, c_null_ptr
    use netcdf_util, only: nc_grid_output_type, &
        gridded_output_variables_time, & ! Names of time-varying grids that SWALS can output
        gridded_output_variables_max_U, & ! Names of grids that SWALS can output once, which require tracking the flow
        gridded_output_variables_other ! Names of other grids that SWALS can output once
    use logging_mod, only: log_output_unit
    use file_io_mod, only: mkdir_p
    use cliffs_tolkova_mod, only: cliffs, setSSLim
    use extrapolation_limiting_mod, only: limited_gradient_dx_vectorized
    use linear_dispersive_solver_mod, only: dispersive_solver_type

#ifdef SPHERICAL
    ! Compile with -DSPHERICAL to get the code to run in spherical coordinates
    use spherical_mod, only: area_lonlat_rectangle, deg2rad, earth_angular_freq, distance_haversine
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

    ! Make everything private except domain_type (which has its own methods) and
    ! the indices of key variables in the 3rd dimension of domain%U(:,:,:)
    private
    public:: domain_type, STG, UH, VH, ELV, test_domain_mod

    integer(int32), parameter :: STG=1, UH=2, VH=3, ELV=4
        ! Indices for arrays: Stage, depth-integrated-x-velocity, depth-integrated-y-velocity, elevation. So e.g. stage is in
        ! domain%U(:,:,STG), and elevation is in domain%U(:,:,ELV)

    ! Handy constants
    real(dp), parameter :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp
    real(dp), parameter :: QUARTER_dp = HALF_dp * HALF_dp
    real(dp), parameter :: NEG_SEVEN_ON_THREE_dp = -7.0_dp/3.0_dp

    ! Make 'rk2n' involve this many flux calls -- it will evolve (rk2n_n_value -1) steps,
    ! with one extra step used for a second-order correction
    integer(ip), parameter :: rk2n_n_value = 5_ip

    ! Use this formatting when converting domain%myid to character
    character(len=charlen), parameter:: domain_myid_char_format = '(I0.20)'

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
    ! We use 'limiter_coef1' for the 'small' threshold, and limiter_coef2 for the 'not small' threshold.
    real(dp), parameter :: limiter_coef1 = 0.1_dp, limiter_coef2 = 0.3_dp, limiter_coef3 = 100.0_dp
    !real(dp), parameter :: limiter_coef1 = 0.00_dp, limiter_coef2 = 0.15_dp, limiter_coef3 = 10.0_dp
    !real(dp), parameter :: limiter_coef1 = 0.01_dp, limiter_coef2 = 0.05_dp, limiter_coef3 = 10.0_dp
    !! This one turns off limiting in practice!
    !real(dp), parameter :: limiter_coef1 = -10000.03_dp, limiter_coef2 = 0.15_dp, limiter_coef3 = 10.0_dp
    ! The following value often occurs in the calculations
    real(dp), parameter :: limiter_coef4 = ONE_dp/(limiter_coef2 - limiter_coef1)


    type :: forcing_term_type
        !!
        !! In SWALS, domains can have user-specified forcing terms defined by setting
        !!    domain%forcing_subroutine => subroutine_matching_forcing_subroutine_interface
        !!    domain%forcing_context_cptr = c_loc(data_required_by_subroutine_above)
        !! But what happens if we want to use multiple forcing terms? In that case the current
        !! type is used to store an array of forcing subroutine/context pairs, which we loop over.
        !!
        procedure(forcing_subroutine), pointer, nopass:: forcing_subroutine => NULL()
            !! The user can use this to create source-terms. If associated, is called inside domain%apply_forcing
        type(c_ptr) :: forcing_context_cptr = c_null_ptr
            !! This can be used to pass arbitrary data to the user-specified forcing_subroutine.
            !! Before we call the forcing_subroutine, the code will set domain%forcing_context_cptr to be this value
            !! The forcing subroutine takes the domain as an argument, so this allows the forcing_context to be accessed
            !! FIXME: Consider whether this can be made class(*), so that the memory cleanup is automatic.
    end type

    type :: domain_type
        !!
        !! Main type for solving the linear/nonlinear shallow water equations on a single
        !! structured grid.
        !!

        !
        ! Geometry
        !
        real(dp):: lw(2) = -HUGE(1.0_dp)
            !! [Length,width] of domain in x/y units.
        integer(ip):: nx(2) = -1_ip
            !! grid size [number of x-cells, number of y-cells]
        real(dp):: dx(2) = - HUGE(1.0_dp)
            !! cell size [dx,dy] (= lw/nx)
        real(dp):: lower_left(2) = HUGE(1.0_dp)
            !! Absolute lower-left coordinate (bottom left corner of cell(1,1))

        !
        ! Solver related
        !
        integer(ip):: nvar = 4
            !! Number of quantities in domain%U (stage, uh, vh, elevation)
        character(len=charlen):: timestepping_method = default_nonlinear_timestepping_method
            !! timestepping_method determines the choice of solver
        integer(ip) :: timestepping_method_nesting_thickness_for_one_sw_timestep = -1
            !! Holds "nesting_thickness_for_one_sw_timestep" for the shallow-water timestepping
            !! method, as defined in timestepping_metadata_mod
        character(len=charlen) :: compute_fluxes_inner_method = 'DE1_low_fr_diffusion' ! 'DE1'
            !! Subroutine called inside domain%compute_fluxes, to allow variations on flux computation.
        real(dp) :: theta = -HUGE(1.0_dp) !
            !! Parameter controlling extrapolation for finite volume methods
        real(dp):: cfl = -HUGE(1.0_dp)
            !! CFL number
        integer(ip):: nsteps_advanced = 0
            !! Count number of time steps (useful if we want to do something only
            !! every n'th timestep)
        real(dp):: time = ZERO_dp
            !! The evolved time in seconds
        real(dp):: max_dt = ZERO_dp
            !! dt computed from CFL condition (e.g. during flux computation). We might take a smaller timestep (e.g. to match the
            !! timestep on other nested domains)
        real(dp) :: dt_last_update = ZERO_dp
            !! Record dt the last time we ran update_U. Useful if we want to append some other update.
        real(dp):: evolve_step_dt = ZERO_dp
            !! Record the time-step evolved by domain%evolve_one_step.
        real(dp):: maximum_timestep = maximum_timestep
            !! This can set the maximum allowed timestep [used to prevent high timesteps on dry domains]
        real(dp) :: static_before_time = -HUGE(1.0_dp)
            !! If the time is less than this number, then we will assume the domain flow is stationary, and will not re-compute the
            !! flow. Useful to increase the efficiency of domains where you know nothing happens until some set time (e.g. far-field
            !! tsunami).
        real(dp):: msl_linear = 0.0_dp
            !! Parameter which determines how the 'static depth' is computed in
            !! linear solver [i.e. corresponding to 'h0' in the pressure gradient term 'g h_0 d(free_surface)/dx'],
            !! AND the reference level for potential energy computation (for all solvers)
        logical :: linear_solver_is_truely_linear = .true.
            !! This flag controls whether, for the 'linear' solver, we allow the
            !! pressure-gradient term (g h0 dStage/dx) to have 'h0' varying
            !! over time. If linear_solver_is_truely_linear, then 'h0' = (msl_linear - elevation).
            !! Otherwise 'h0' varies as the stage varies (which actually makes the equations nonlinear)
        logical :: is_staggered_grid
            !! Useful variable to distinguish staggered-grid and centred-grid numerical methods
        logical :: adaptive_timestepping
            !! Can the time-step vary over time?
        real(dp) :: cliffs_minimum_allowed_depth = 0.001_dp
            !! The CLIFFS solver seems to require tuning of the minimum allowed depth (and often it should
            !! be larger than for SWALS finite-volume solvers). For field-scale problems, Tolkova (2014) mentions
            !! values of 0.1 m on an inundation grid, and 0.5m on a larger-scale grid. For experimental type
            !! problems, Tolkova (2014) mentions values of say 1 - 2 mm.
        real(dp) :: cliffs_bathymetry_smoothing_alpha = 2.0_dp
            !! For bathymetry smoothing with the cliffs method, Tolkova (in CLIFFS manual) suggests this is typically a reasonable
            !! value.
        real(dp) :: leapfrog_froude_limit = 10.0_dp
            !! Froude-limiter for nonlinear leapfrog scheme. Velocity will be supressed at higher
            !! froude numbers, so as not to exceed this.
        character(len=charlen) :: friction_type = 'manning'
            !! If friction_type = 'manning' then interpret domain%manning_squared as (manning's n)**2
            !! If friction_type = 'chezy' then interpret domain%manning_squared as (1/chezy_friction)**2, and use the Chezy friction
            !! model
        real(dp) :: linear_friction_coeff = 0.0_dp
            !! Linear friction coefficient similar to Fine et al., (2012),  Kulikov et al., (2014), and others.
            !! For the UH (and VH) equations we add a term '- linear_friction_coeff * UH (or VH)' to the right-hand-side.
            !! This kind of linear decay seems heuristically consistent with the observed exponential decay of tsunami energy once
            !! it has spread globally. In the aforementioned papers which do not use nonlinear friction, they use
            !! linear_friction_coeff = 1.0e-05_dp.
            !! Note this parameterization is different to that used often in the global tidal modelling literature for
            !! linear friction -- where instead the term is like "drag_coefficient * U (or V)" when the momentum equation
            !! is written in terms of the depth-integrated velocity (e.g. Egbert and Erofeeva 2002)
        integer(ip):: xL, xU, yL, yU
            !! Lower/upper x and y indices over which the SWE computation takes place
            !! These might restrict the SWE update to a fraction of the domain (e.g.
            !! beginning of an earthquake-tsunami run where only a fraction of the domain is
            !! active. Currently implemented for linear leap-frog solvers only

        !
        ! Eddy viscosity terms -- only for nonlinear finite-volume solvers
        !
        logical :: eddy_viscosity_is_supported
            !! Can eddy-viscosity be used with this solver? Depends on the time-stepping method
        logical :: use_eddy_viscosity = .false.
            !! Use a flux routine that includes a contribution from eddy-viscosity?
        real(dp) :: eddy_visc_constants(2) = [0.0_dp, 0.5_dp]
            !! Constants used to control the eddy-viscosity model. All eddy-visc-models are of the form
            !! "constant + another_model". The first coefficient is the constant eddy viscosity,
            !! while the second relates to either the smagorinksy eddy-viscocity model, or the
            !! shiono/knight eddy-viscosity model.

        !
        ! Dispersion
        !
        logical :: use_dispersion = .false. ! By default do not use dispersion
        type(dispersive_solver_type) :: ds

        !
        ! Nesting-related
        !
        type(domain_nesting_type) :: nesting
            !! Type to manage nesting communication. This is required for nesting, and is the standard
            !! approach for distributed-memory parallel runs.
        real(dp) :: dx_refinement_factor = ONE_dp
        real(dp) :: dx_refinement_factor_of_parent_domain = ONE_dp
            !! This information is useful in the context of nesting, where interior
            !! domains have cell sizes that are an integer divisor of the
            !! coarse-domain dx.
        integer(ip) :: timestepping_refinement_factor = 1_ip
            !! Number of timesteps taken by this domain inside each multidomain timestep.
            !! Only used in conjunction with the multidomain class
        integer(ip) :: nest_layer_width = -1_ip
            !! The width of the nesting layer for this domain in the multidomain. This will depend
            !! on the dx value relative to the neighbouring domains, on the timestepping_method,
            !! and on the timestepping_refinement_factor. See 'get_domain_nesting_layer_thickness'
        real(dp) :: interior_bounding_box(4,2) = 0.0_dp
            !! Useful to hold the interior bounding box [ = originally provided bounding box]
            !! since the actual bounding box might be changed to accommodate nesting
        integer(int64):: myid = 1
            !! Domain ID, which is useful if multiple domains are running
        integer(int64):: local_index = 1
            !! Useful when we partition a domain in parallel
        logical :: is_nesting_boundary(4) = .FALSE.
            !! Flag to denote boundaries at which nesting occurs: order is N, E, S, W.
        real(dp) :: max_parent_dx_ratio = -1.0_dp
            !! Maximum value of (dx-from-domains-we-receive-from)/my-dx. Useful for nesting.

        !
        ! Boundary conditions.
        !
        character(len=charlen):: boundary_type = '' !! FIXME: DEFUNCT (??)
        procedure(boundary_fun), pointer, nopass:: boundary_function => NULL()
            !! Some existing boundary conditions rely on 'boundary_function' as well. It takes
            !! the domain_type as well as time, and the location (defined as i, j indices), see the interface below.
            !! Note that the x/y location can easily be obtained from domain%x(i), domain%y(j)
        procedure(boundary_subroutine), pointer, nopass:: boundary_subroutine => NULL()
            !! Currently 'boundary_subroutine' is the most important way of imposing
            !! boundary conditions. It must take a domain_type as an INTENT(INOUT) argument.
            !! Applications need to define boundary_subroutine, either using a
            !! pre-existing bc, or making a new subroutine.
        type(c_ptr) :: boundary_context_cptr = c_null_ptr
            !! This can be used to pass arbitrary data to the user-specified boundary subroutine.
            !! The latter takes the domain as an argument, so one can access the boundary_context
            !! via that. FIXME: Consider making this class(*) if the compiler support is OK, so that
            !! the memory cleanup is automatic.
        logical :: boundary_exterior(4) = .TRUE.
            !! Flag whether the boundary is exterior (TRUE) or interior (FALSE). We only need
            !! to apply a boundary condition to the 'exterior boundaries' -- otherwise we have e.g. nesting updates,
            !! which are different.
            !! Order is North (1), East (2), South (3), West (4)

        !
        ! Forcing terms
        !
        procedure(forcing_subroutine), pointer, nopass:: forcing_subroutine => NULL()
            !! The user can use this to create source-terms. If associated, is called inside domain%apply_forcing
            !! If you need to add multiple forcing terms, then for each forcing term, set this (and forcing_context_cptr if
            !! required) and then "call domain%store_forcing()"
        type(c_ptr) :: forcing_context_cptr = c_null_ptr
            !! This can be used to pass arbitrary data to the user-specified forcing_subroutine.
            !! The latter takes the domain as an argument, so one can access the forcing_context
            !! via that. FIXME: Consider making this class(*) if compiler support is OK, so the memory
            !! cleanup is automatic
        type(forcing_term_type), allocatable :: forcing_terms_storage(:)
            !! This is used internally in the code to store multiple forcing terms. A call to domain%store_forcing
            !! will move domain%forcing_subroutine and domain%forcing_context_cptr into an element of this array. This
            !! can be repeated to add in mutliple forcing terms.
        logical :: support_elevation_forcing = .FALSE.
            !! If the forcing_patch_type tries to evolve the elevation, but support_elevation_forcing is FALSE, then
            !! the code will throw an error. The value of this variable will be overridden automatically for schemes
            !! where forcing the elevation is always OK. For other schemes one can manually set this to TRUE, and
            !! the code will do extra calculations to enable elevation forcing.

        !
        ! Mass conservation tracking  -- store as double, even if dp is single prec.
        !
        real(force_double):: boundary_flux_store(4)  = ZERO_dp
            !! Store the flux through the N, E, S, W boundaries, single domain
        real(force_double):: boundary_flux_time_integral = ZERO_dp
            !! Time integrate the boundary fluxes, single domain
        real(force_double):: boundary_flux_evolve_integral = ZERO_dp
            !! We need an intermediate variable to take care of boundary_flux time-stepping.
            !! This integrates the boundary fluxes within the evolve step only, single domain
        real(force_double):: boundary_flux_store_exterior(4)  = ZERO_dp
            !! Store the flux through priority_domain parts of N, E, S, W boundaries, as required for the nested-grid case
        real(force_double):: boundary_flux_time_integral_exterior = ZERO_dp
            !! Time integrated boundary fluxes in nested-grid case
        real(force_double):: boundary_flux_evolve_integral_exterior = ZERO_dp
            !! Alternative to the above which is appropriate for the nested-grid case
        integer(long_long_ip) :: negative_depth_fix_counter = 0
            !! Count the number of time-steps at which we clip depths.
            !! For a 'well behaved' model this should remain zero -- non-zero values indicate mass conservation issues.
        integer(ip):: exterior_cells_width = 2
            !! The domain 'interior' is surrounded by 'exterior' cells which are
            !! updated by boundary conditions, or copied from other domains. When
            !! tracking mass conservation, we only want to record inflows/outflows to
            !! interior cells. The interior cells are:
            !!    [(1+exterior_cells_width):(domain%nx(1)-exterior_cells_width), &
            !!     (1+exterior_cells_width):(domain%nx(2)-exterior_cells_width)]
            !! NOTE: exterior_cells_width should only be used to influence mass conservation tracking calculations.
            !! It is not strictly related to the halo width for parallel computations. When using a multidomain, the
            !! mass conservation tracking is somewhat different (based around subroutines with names matching
            !! 'nesting_boundary_flux_*')
        integer(ip) :: minimum_nesting_layer_thickness = 0_ip
            !! There are situations where we might want to force the nesting layer thickness for parallel computations.
            !! Set this to the desired thickness. The actual value will not be less than this.

        !
        ! Output file content -- bespoke binary format (that is no longer properly supported, consider removing)
        !
        character(len=charlen):: metadata_ascii_filename
            !! Name of ascii file where we output metadata
        integer(ip):: output_time_unit_number
            !! Output file unit for time (stored as ascii)
        character(len=charlen):: output_basedir = default_output_folder
            !! Base folder inside which we store domain outputs.
        character(len=charlen):: output_folder_name = ''
            !! Output folder (it will begin with domain%output_basedir, and include another timestamped folder)
        logical :: output_folders_were_created = .FALSE.
            !! To avoid trying to repeatedly create output folders
        integer(ip):: logfile_unit = output_unit
            !! Unit number for domain log file
        integer(ip), allocatable :: output_variable_unit_number(:)
            !! Units for home-brew binary output format [better to use netcdf nowadays].
            !! This is only used when the code is compiled with "-DNONETCDF"
        !
        ! netcdf-output
        !
        type(nc_grid_output_type) :: nc_grid_output
            !! Type to manage netcdf grid outputs
        character(len=charlen), allocatable :: time_grids_to_store(:)
            !! Specify which gridded variables to store over time. For example:
            !!   [character(len=charlen):: 'stage', 'uh', 'vh', 'elev'] to store everything;
            !!   [''] to store nothing;
            !!   ['stage'] to just store the stage.
            !! If unallocated then it will be set in domain%nc_grid_output%initialise
        character(len=charlen), allocatable :: nontemporal_grids_to_store(:)
            !! Specify which 'nontemporal' variables to store. For example:
            !!   [character(len=charlen):: 'max_stage', 'max_flux', 'max_speed', 'arrival_time', 'manning_squared', 'elevation0', 'elevation_source_file_index', 'time_of_max_stage', 'min_stage']
            !! to store everything;
            !!   [''] to store nothing;
            !!   ['max_stage'] to just store the maximum stage.
        logical:: record_max_U = .true.
            !! (Deprecated) approach to specifying whether we store maximum-stage (.true.) or not (.false.).
            !! For new code please specify 'nontemporal_grids_to_store' instead, because the latter is more flexible
        character(len=charlen), allocatable :: max_U_variables(:)
            !! Names of variables contained in max_U. Currently if it has size > 0, then the first variable is
            !! ALWAYS 'max_stage' (for backward compatability - consider revising)
        integer(ip):: max_U_update_frequency = 1
            !! Only update domain%max_U every n time-steps. Values other than 1 introduce
            !! some approximation error, but the option might be useful for speed in some cases.
        real(dp) :: arrival_stage_threshold_above_msl_linear = 0.01_dp
            !! The "arrival time" is defined as the time at which the stage exceeds
            !! (msl_linear + arrival_stage_threshold_above_msl_linear) AND the cell is wet
        type(point_gauge_type) :: point_gauges
            !! Type to manage storing of tide gauges

        !
        ! Timing
        !
        type(timer_type):: timer, evolve_timer
            !! Types to record CPU timings

        !
        ! 1D arrays
        !
        real(dp), allocatable :: x(:), y(:), distance_bottom_edge(:), distance_left_edge(:)
            !! Spatial coordinates, dx/dy distances (useful for spherical coordinates)
        real(dp), allocatable :: area_cell_y(:)
            !! Area of cells (varies with 'y' only)
        real(dp), allocatable :: coslat(:), coslat_bottom_edge(:), tanlat_on_radius_earth(:)
            !! These variables are only used with spherical coordinates
        real(dp), allocatable :: coriolis(:), coriolis_bottom_edge(:)
            !! These variables are only used with coriolis
        real(force_double), allocatable :: volume_work(:), energy_potential_on_rho_work(:), energy_kinetic_on_rho_work(:)
            !! Work arrays which permit an openmp-reproducible summation
            !! (since for standard REDUCE(+: variable), the summation order may change).

        !
        ! Big arrays to hold the domain variables
        !
        real(dp), allocatable :: U(:,:,:)
            !! U holds main flow variables -- STG, UH, VH, ELV
            !! First 2 dimensions = space, 3rd = number of quantities
        real(dp), allocatable :: max_U(:,:,:)
            !! Store various flow maxima. FIXME: Consider storing in reduced precision.
        !
        ! Multi-dimensional arrays below are only required for nonlinear solver. Would be possible
        ! to further reduce memory usage, but not completely trivial. Not all of these are required
        ! for linear type timestepping
        real(dp), allocatable :: flux_NS(:,:,:)  !! NS fluxes
        real(dp), allocatable :: flux_EW(:,:,:)  !! EW fluxes
        real(dp), allocatable :: depth(:,:) !! Depths at cell centre
        real(dp), allocatable :: velocity(:,:, :) !! Velocity
        real(force_double), allocatable :: explicit_source(:,:,:) !! Used for pressure gradient terms
        real(force_double), allocatable :: explicit_source_VH_j_minus_1(:,:) !! Separate from explicit_source for OPENMP parallel logic
        real(dp), allocatable :: manning_squared(:,:) !! Friction
        real(dp), allocatable :: backup_U(:,:,:) !! Needed for some timestepping methods
        real(dp), allocatable :: friction_work(:,:,:) !! Friction for leapfrog_nonlinear and leapfrog_linear_plus_nonlinear_friction
        logical :: friction_work_is_setup = .false. !! Flag used for efficiency with leapfrog_linear_plus_nonlinear_friction solvers
        real(dp), allocatable :: advection_work(:,:,:) !! Work array
        integer(ip) :: advection_work_dim3_size = 5 !! Size of 3rd dimension of advection_work
#ifdef DEBUG_ARRAY
        real(dp), allocatable :: debug_array(:,:)
            !! For debugging it can be helpful to have this array.  If DEBUG_ARRAY is defined, then it will be allocated with
            !! dimensions (nx, ny), and will be written to the netcdf file at each time.
#endif
        logical :: friction_with_ambient_fluxes = .false.
            !! Optionally include other currents (e.g. tides) in the friction term. This is an experiment, currently only supported for
            !! leapfrog_linear_with_nonlinear_friction
        real(dp), allocatable :: ambient_flux(:,:,:)
            !! Ambient (e.g. tidal) depth-integrated velocities with easting/northing, used if friction_with_ambient_fluxes=.true.
            !! FIXME: Consider removing this.
        real(dp), allocatable :: elevation_source_file_index(:,:)
            !! Store the file index used to set the elevation data at each point. 
            !! This can be recorded when using the multi_raster_type
            !! to set the elevation. Although conceptually an integer, it is real 
            !! because less code modification was required.

        CONTAINS

        ! Initialisation
        procedure, non_overridable :: allocate_quantities => allocate_quantities

        ! Reporting
        procedure, non_overridable :: print => print_domain_statistics

        ! Core routines that occur within a timestep
        ! (consider making these not type bound -- since the user should not really call them)
        procedure, non_overridable :: compute_depth_and_velocity => compute_depth_and_velocity
        procedure, non_overridable :: compute_fluxes => compute_fluxes
        !procedure, non_overridable :: compute_fluxes => compute_fluxes_vectorized !! Slower than un-vectorized version on GD home machine
        procedure, non_overridable :: update_U => update_U  ! Slower or faster, depending on wet/dry areas
        !procedure, non_overridable :: update_U => update_U_restructured  ! Slower or faster, depending on wet/dry areas
        !procedure, non_overridable :: update_U => update_U_vectorized !! Slower on an NCI test
        procedure, non_overridable :: backup_quantities => backup_quantities

        ! Timestepping (consider making only 'evolve_one_step' type bound -- since the user should not really call others)
        procedure, non_overridable :: evolve_one_step => evolve_one_step
        procedure, non_overridable :: update_max_quantities => update_max_quantities

        ! Boundary conditions. This just calls whatever domain%boundary_subroutine points to (consider making not type bound)
        procedure, non_overridable :: update_boundary => update_boundary

        ! Forcing term
        procedure, non_overridable :: apply_forcing => apply_forcing
        procedure, non_overridable :: store_forcing => store_forcing

        ! Dispersive terms
        procedure, non_overridable :: copy_U_to_dispersive_last_timestep => copy_U_to_dispersive_last_timestep
        procedure, non_overridable :: solve_dispersive_terms => solve_dispersive_terms

        ! IO
        procedure, non_overridable :: create_output_folders => create_output_folders
        procedure, non_overridable :: create_output_files => create_output_files
        procedure, non_overridable :: write_to_output_files => write_to_output_files
        procedure, non_overridable :: write_max_quantities => write_max_quantities
        procedure, non_overridable :: log_outputs => divert_logfile_unit_to_file

        ! Mass conservation and other reporting
        procedure, non_overridable :: mass_balance_interior => mass_balance_interior
        procedure, non_overridable :: volume_interior => volume_interior
        procedure, non_overridable :: volume_in_priority_domain => volume_in_priority_domain
        procedure, non_overridable :: compute_domain_statistics => compute_domain_statistics_new

        ! Functions to estimate gravity-wave time-step (often useful to call at the start of a model run). These are not used
        ! in the flow-algorithms, they are just for convenience.
        procedure, non_overridable :: linear_timestep_max => linear_timestep_max
        procedure, non_overridable :: nonlinear_stationary_timestep_max => nonlinear_stationary_timestep_max
        procedure, non_overridable :: stationary_timestep_max => stationary_timestep_max

        ! Setup and store point-gauge outputs
        procedure, non_overridable :: setup_point_gauges => setup_point_gauges
        procedure, non_overridable :: write_gauge_time_series => write_gauge_time_series

        ! finalization (e.g. close netcdf files)
        procedure, non_overridable :: finalise => finalise_domain

        ! Smoothing of elevation. Not recommended in general, but with poor elevation data it might reduce artefacts.
        procedure, non_overridable :: smooth_elevation => smooth_elevation
        ! Sometimes local smoothing near nesting boundaries can aid stability
        procedure, non_overridable :: smooth_elevation_near_point => smooth_elevation_near_point
        procedure, non_overridable :: smooth_elevation_near_nesting_fine2coarse_boundaries => &
            smooth_elevation_near_nesting_fine2coarse_boundaries

        ! Sending and receiving halos
        procedure, non_overridable :: send_halos => send_halos
        procedure, non_overridable :: receive_halos => receive_halos

        ! Nesting.
        ! -- Keep track of fluxes for mass-conservation tracking
        procedure, non_overridable :: nesting_boundary_flux_integral_multiply => nesting_boundary_flux_integral_multiply
        procedure, non_overridable :: nesting_boundary_flux_integral_tstep => nesting_boundary_flux_integral_tstep
        procedure, non_overridable :: nesting_flux_correction_everywhere => nesting_flux_correction_everywhere
        ! -- Determine if a given point is in the 'priority domain'
        procedure, non_overridable :: is_in_priority_domain => is_in_priority_domain
        ! -- When initialising nested domains, set lower-left/upper-right/resolution near the desired locations,
        ! adjusting as required to support nesting (e.g. so that edges align with parent domain, and cell size is an integer divisor
        ! of the parent domain, etc).
        procedure, non_overridable :: match_geometry_to_parent => match_geometry_to_parent


    end type domain_type

    interface

        function boundary_fun(domain, t, i, j) RESULT(stage_uh_vh_elev)
            !! This gives the interface for a generic 'boundary function' which
            !! returns an array of 4 output values (stage/uh/vh/elevation)
            !!
            !! It is used in conjunction with a number of different types of
            !! boundary conditions.
            !!
            import dp, ip, domain_type
            implicit none
            type(domain_type), intent(in):: domain
            real(dp), intent(in):: t ! Time
            integer(ip), intent(in) :: i, j
                !! i,j index of the domain at which we compute the boundary condition. Generally these values
                !! would be at the domain boundary (e.g. i == 1 or i == domain%nx(1) or j == 1 or j == domain%nx(2))
            real(dp):: stage_uh_vh_elev(4) !! Stage, uh, vh, elevation according to the boundary condition.
        end function

        subroutine boundary_subroutine(domain)
            !! The user can provide a boundary subroutine which is supposed to update
            !! the domain boundaries. It is called by domain%update_boundary whenever
            !! a boundary update is required by the timestepping_method. Note that
            !! this may well mean 'boundary_fun' is not required
            import domain_type
            implicit none
            type(domain_type), intent(inout):: domain
        end subroutine

        subroutine forcing_subroutine(domain, dt)
            !! The user can provide a boundary subroutine which is supposed to update
            !! the domain boundaries. It is called by domain%update_boundary whenever
            !! a boundary update is required by the timestepping_method. Note that
            !! this may well mean 'boundary_fun' is not required
            import domain_type, dp
            implicit none
            type(domain_type), intent(inout):: domain
            real(dp), intent(in) :: dt
        end subroutine
    end interface

    contains

    subroutine print_domain_statistics(domain)
        !!
        !! Convenience printing function.
        !!
        class(domain_type), intent(inout):: domain

        real(dp):: maxstage, minstage, minspeed, maxspeed
        integer:: ecw
        real(dp) :: energy_total_on_rho, energy_potential_on_rho, energy_kinetic_on_rho

TIMER_START('printing_stats')

        call domain%compute_domain_statistics(maxstage, maxspeed, minstage, minspeed, energy_potential_on_rho,&
            energy_kinetic_on_rho, energy_total_on_rho)

        ! Print main statistics
        write(domain%logfile_unit, *) ''
        write(domain%logfile_unit, *) 'Domain ID: '
        write(domain%logfile_unit, *) '        ', domain%myid
        write(domain%logfile_unit, *) 'Time: '
        write(domain%logfile_unit, *) '        ', domain%time
        write(domain%logfile_unit, *) 'nsteps_advanced:'
        write(domain%logfile_unit, *) '        ', domain%nsteps_advanced
        write(domain%logfile_unit, *) 'negative_depth_fix_counter:'
        write(domain%logfile_unit, *) '        ', domain%negative_depth_fix_counter
        write(domain%logfile_unit, *) 'max_allowed_dt: '
        write(domain%logfile_unit, *) '        ', domain%max_dt
        write(domain%logfile_unit, *) 'evolve_step_dt: '
        write(domain%logfile_unit, *) '        ', domain%evolve_step_dt
        write(domain%logfile_unit, *) 'Stage: '
        write(domain%logfile_unit, *) '        ', maxstage
        write(domain%logfile_unit, *) '        ', minstage
        write(domain%logfile_unit, *) 'Speed: '
        write(domain%logfile_unit, *) '        ', maxspeed
        write(domain%logfile_unit, *) '        ', minspeed
        write(domain%logfile_unit, "(A, ES25.12E3)") 'Energy Potential / rho (zero when stage=domain%msl_linear): ', &
                                      energy_potential_on_rho
        write(domain%logfile_unit, "(A, ES25.12E3)") 'Energy Kinetic / rho: ', energy_kinetic_on_rho
        write(domain%logfile_unit, "(A, ES25.12E3)") 'Energy Total / rho: ', energy_total_on_rho
        ! Mass conservation check
        write(domain%logfile_unit, *) 'Mass Balance (domain interior): ', domain%mass_balance_interior()

TIMER_STOP('printing_stats')
    end subroutine


    subroutine compute_domain_statistics_new(domain, maxstage, maxspeed, minstage, minspeed, &
        energy_potential_on_rho, energy_kinetic_on_rho, energy_total_on_rho)
        !!
        !! Compute various statistics in the domain. If nesting is occurring, then only consider
        !! priority_domain cells. Never include cells within 'exterior_cells_width' from the boundary.
        !!
        class(domain_type), intent(inout) :: domain
            !! The domain for which we compute statistics
        real(dp), intent(out) :: maxstage, maxspeed, minstage, minspeed
            !! Stage and speed statistics in wet areas
        real(dp), intent(out) :: energy_potential_on_rho
            !! Potential-energy/water-density. Potential energy includes an arbitrary offset, and herein that
            !! is defined so that the potential energy is zero if stage=msl_linear (or dry) everywhere.
        real(dp), intent(out) :: energy_kinetic_on_rho
            !! Kinetic energy/water-density.
        real(dp), intent(out) :: energy_total_on_rho


        logical :: is_nesting, is_wet, use_truely_linear_method
        integer(ip) :: ecw, i, j
        real(dp) :: depth_C, depth_N, depth_E, depth_N_inv, depth_E_inv, speed_sq
        ! Ensure the energy accumulation uses high precision, easy to lose precision here.
        real(force_double) :: e_flow, e_constant, w0, d1, w2, d3

TIMER_START("compute_statistics")
        !print*, "force_double: ", force_double
        !print*, "domain%msl_linear: ", domain%msl_linear

        ! If nesting is occurring, the above stats are only computed
        ! in the priority domain area
        is_nesting = (domain%nesting%my_index > 0)

        ! Reset stats of interest
        maxstage = -HUGE(1.0_dp)
        maxspeed = 0.0_dp ! speed is always > 0
        minstage = HUGE(1.0_dp)
        minspeed = HUGE(1.0_dp) ! speed is always > 0

        ecw = domain%exterior_cells_width

        ! Use these arrays to store the row-wise sums of potential and kinetic energy.
        domain%energy_potential_on_rho_work = 0.0_force_double
        domain%energy_kinetic_on_rho_work = 0.0_force_double

        ! Compute stats of interest in parallel
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(is_nesting, domain, ecw), &
        !$OMP REDUCTION(max: maxspeed) REDUCTION(min: minspeed), &
        !$OMP REDUCTION(max: maxstage) REDUCTION(min: minstage)

        ! Reset stats of interest
        maxstage = -HUGE(1.0_dp)
        maxspeed = 0.0_dp ! speed is always > 0
        minstage = HUGE(1.0_dp)
        minspeed = HUGE(1.0_dp)

        ! If we are using a 'truely-linear' solver then the depth is always recorded from MSL for certain
        ! calculations (pressure gradient term, and wetting/drying)
        use_truely_linear_method = domain%is_staggered_grid .and. domain%linear_solver_is_truely_linear .and. &
            any(domain%timestepping_method == [character(len=charlen) :: 'linear', 'leapfrog_linear_plus_nonlinear_friction'])

        !$OMP DO SCHEDULE(STATIC)
        do j = (1 + ecw), (domain%nx(2) - ecw )
            do i = (1 + ecw), (domain%nx(1) - ecw )

                if(is_nesting) then
                    ! Cycle if this cell is not on the priority domain
                    if( domain%nesting%is_priority_domain_not_periodic(i, j) == 0_ip ) cycle
                end if

                if(use_truely_linear_method) then
                    is_wet = (domain%msl_linear > domain%U(i,j,ELV) + minimum_allowed_depth)
                        ! For the 'truely linear' equations, the "depth" is relative to MSL, irrespective
                        ! of the stage. Note that for 'truely-linear' equations it is valid to have
                        ! the stage below the bed (we can create such a model -- just multiply the
                        ! solution by a sufficiently large number -- which must be a solution according to linearity).
                else
                    is_wet = (domain%U(i,j,STG) > domain%U(i,j,ELV) + minimum_allowed_depth)
                end if

                ! The potential energy is defined as:
                !     Integral of [ g * depth*(elev + 1/2 depth) ]
                ! This can be rearranged noting (depth = stage - elev) to give
                !     Integral of [ (g/2 stage^2 - g/2 elev^2) ]
                ! ('difference of 2 squares')
                w0 = real(domain%U(i,j,STG), force_double) + real(domain%U(i,j,ELV), force_double)
                d1 = real(domain%U(i,j,STG), force_double) - real(domain%U(i,j,ELV), force_double)
                e_flow = merge( w0*d1, 0.0_force_double, is_wet)
                ! Once integrated the g/2 elev^2 term will be a constant if the bathymetry is fixed and there is no wetting and
                ! drying. We are only interested in changes in integrated energy, and subtracting that constant recovers the
                ! (g/2 stage^2) formulation of potential energy, which is typically used for the linear shallow water equations in
                ! deep ocean tsunami type cases.
                w2 = real(domain%msl_linear, force_double) + real(domain%U(i,j,ELV), force_double)
                d3 = real(domain%msl_linear, force_double) - real(domain%U(i,j,ELV), force_double)
                e_constant = merge( w2*d3, 0.0_force_double, domain%msl_linear > domain%U(i,j,ELV) + minimum_allowed_depth)
                    ! e_constant must be COMPLETELY UNAFFECTED by changes in the flow - it is a constant offset.

                ! ! Here is the classical form of available potential energy (which works without wetting/drying)
                !w0 = (real(domain%U(i,j,STG), force_double) - real(domain%msl_linear, force_double))**2
                !e_flow = merge(w0, 0.0_force_double, is_wet)
                !e_constant = 0.0_force_double

                domain%energy_potential_on_rho_work(j) = domain%energy_potential_on_rho_work(j) + &
                    real(domain%area_cell_y(j) * gravity * 0.5_dp, force_double) * (e_flow - e_constant)

                if(is_wet) then

                    maxstage = max(maxstage, domain%U(i,j,STG))
                    minstage = min(minstage, domain%U(i,j,STG))

                    if(domain%is_staggered_grid) then

                        ! Get depths at UH point (E) and VH point (N).
                        ! The interpretation of 'depth' varies depending on whether we are using the 'truely-linear'
                        ! treatment of the pressure gradient term in the linear shallow water equations.
                        if(use_truely_linear_method) then
                            ! Depth is always relative to MSL.
                            depth_E = 0.5_dp * sum( (domain%msl_linear - domain%U(i:i+1,j,ELV)))
                            depth_N = 0.5_dp * sum( (domain%msl_linear - domain%U(i,j:j+1,ELV)))
                        else
                            ! Spatially averaged depth on staggered grid
                            depth_E = 0.5_dp * sum( (domain%U(i:i+1,j,STG) - domain%U(i:i+1,j,ELV)))
                            depth_N = 0.5_dp * sum( (domain%U(i,j:j+1,STG) - domain%U(i,j:j+1,ELV)))
                        end if

                        depth_E_inv = merge(1.0_dp/depth_E, 0.0_dp, depth_E > minimum_allowed_depth)
                        depth_N_inv = merge(1.0_dp/depth_N, 0.0_dp, depth_N > minimum_allowed_depth)

                        speed_sq = domain%U(i,j,UH) * domain%U(i,j,UH) * depth_E_inv**2  + &
                                   domain%U(i,j,VH) * domain%U(i,j,VH) * depth_N_inv**2

                        ! The kinetic energy is integral of ( 1/2 depth speed^2)
                        domain%energy_kinetic_on_rho_work(j) = domain%energy_kinetic_on_rho_work(j) + real(0.5_dp * &
                            (domain%area_cell_y(j) * domain%U(i,j,UH) * domain%U(i,j,UH) * depth_E_inv + &
                            0.5_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1)) * &
                                domain%U(i,j,VH) * domain%U(i,j,VH) * depth_N_inv  ), &
                            force_double)

                    else
                        ! Co-loated grid
                        depth_C = domain%U(i,j,STG) - domain%U(i,j,ELV)
                        speed_sq = (domain%U(i,j,UH) * domain%U(i,j,UH) + domain%U(i,j,VH) * domain%U(i,j,VH)   )/ &
                                    (depth_C*depth_C)

                        ! The kinetic energy is integral of ( 1/2 depth speed^2)
                        domain%energy_kinetic_on_rho_work(j) = domain%energy_kinetic_on_rho_work(j) + &
                            real(domain%area_cell_y(j) * 0.5_dp * depth_C * speed_sq, force_double)
                    end if

                    maxspeed = max(maxspeed, speed_sq) ! Will undo the sqrt later.
                    minspeed = min(minspeed, speed_sq) ! Will undo the sqrt later

                end if
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        ! Convert 'speed**2' to 'speed'
        maxspeed = sqrt(maxspeed)
        minspeed = sqrt(minspeed)
        ! Copy energies (in force_double) precision back to possibly lower precision input var
        energy_kinetic_on_rho = sum(domain%energy_kinetic_on_rho_work)
        energy_potential_on_rho = sum(domain%energy_potential_on_rho_work)
        energy_total_on_rho = energy_kinetic_on_rho + energy_potential_on_rho

TIMER_STOP("compute_statistics")
    end subroutine


    subroutine allocate_quantities(domain, global_lw, global_nx, global_ll, create_output_files, verbose)
        !!
        !! Set up the full domain, allocate arrays, etc
        !!
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: global_lw(2) !! length/width of domain in same units as x,y coordinates
        real(dp), intent(in):: global_ll(2) !! lower-left x/y coordinates of domain (at corner of the lower-left cell)
        integer(ip), intent(in):: global_nx(2) !! Number of cells in x/y directions
        logical, optional, intent(in) :: create_output_files !! If .TRUE. or not provided, then make output files
        logical, optional, intent(in) :: verbose !! Print info about the domain

        integer(ip) :: nx, ny, nvar
        integer(ip) :: i, j, k, tsi
        logical :: create_output, verbose_
        real(dp):: local_lw(2), local_ll(2)
        integer(ip):: local_nx(2)

        ! Send domain print statements to the default log. Over-ridden later if we send the domain log to its own file.
        domain%logfile_unit = log_output_unit

        create_output = .true.
        if(present(create_output_files)) create_output = create_output_files

        verbose_ = .true.
        if(present(verbose)) verbose_ = verbose

        ! Set default parameters for different timestepping methods
        !
        ! Find index corresponding to domain%timestepping_method in the timestepping_metadata
        tsi = timestepping_method_index(domain%timestepping_method)
        ! Set slope-limiter theta parameter
        if(domain%theta == -HUGE(1.0_dp)) domain%theta = timestepping_metadata(tsi)%default_theta
        ! Set CFL number
        if(domain%cfl == -HUGE(1.0_dp)) domain%cfl = timestepping_metadata(tsi)%default_cfl
        ! Flag to note if the grid should be interpreted as staggered or cell-centred
        domain%is_staggered_grid = (timestepping_metadata(tsi)%is_staggered_grid == 1)
        ! Does the timestepping method allow for adaptive timestepping?
        domain%adaptive_timestepping = timestepping_metadata(tsi)%adaptive_timestepping
        ! Does the timestepping method allow for eddy viscosity?
        domain%eddy_viscosity_is_supported = timestepping_metadata(tsi)%eddy_viscosity_is_supported
        ! ! Uncomment the following line to use eddy viscosity by default if supported (but beware this
        ! ! will overwrite what the user specified)
        !if(domain%eddy_viscosity_is_supported) domain%use_eddy_viscosity = .true.
        if(domain%use_eddy_viscosity .and. (.not. domain%eddy_viscosity_is_supported)) then
            write(log_output_unit, *) 'Error: Eddy viscosity is not supported for the numerical scheme ', &
                trim(domain%timestepping_method)
            call generic_stop
        end if

        ! Do we allow domain%forcing_subroutine to update the elevation? 
        if(timestepping_metadata(tsi)%forcing_elevation_is_allowed == 'always') then
            ! For these schemes, elevation forcing is always possible
            ! (without special treatment)
            domain%support_elevation_forcing=.TRUE.
        end if
        ! Note: if (timestepping_metadata(tsi)%forcing_elevation_is_allowed == 'optional') then
        ! one can manually set domain%support_elevation_forcing=TRUE, prior to this function.
        ! However, one cannot do that if (timestepping_metadata(tsi)%forcing_elevation_is_allowed == 'never')
        if(timestepping_metadata(tsi)%forcing_elevation_is_allowed == 'never' .and. &
            domain%support_elevation_forcing) then
            write(log_output_unit, *) 'Error: The numerical scheme ', trim(domain%timestepping_method)
            write(log_output_unit, *) 'does not allow forcing of elevation.'
            call generic_stop
        end if

        ! Grid geometry
        domain%lw = global_lw ! Same units as coordinate system
        domain%nx = global_nx ! Number of cells
        domain%lower_left = global_ll ! Same units as coordinate system

        ! Cell size (same units as coordinate system)
        domain%dx = domain%lw/(domain%nx*ONE_dp)

        ! For the linear solver, these variables can be modified to evolve domain%U only in 
        ! the region domain%U(xL:xU, yL:yU, ...). This has been used to speed up some single-grid 
        ! simulations, but the user must adaptively change the variables.
        domain%xL = 1_ip
        domain%yL = 1_ip
        domain%xU = domain%nx(1)
        domain%yU = domain%nx(2)

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

        ! Distances along cell edges.
        ! For spherical lon/lat coordinates, distance_bottom_edge changes with y
        ! For cartesian x/y coordinates, the cell edge lengths are constant
        ! In general we allow the bottom-edge distance to change with y, and the left-edge distance to change with x,
        ! although to date I never used the latter case, and beware that other parts of the code may assume
        ! distance_left_edge is constant in space.
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

        ! Work-space to do volume calculations
        allocate(domain%volume_work(ny))
        allocate(domain%energy_potential_on_rho_work(ny))
        allocate(domain%energy_kinetic_on_rho_work(ny))

        !
        ! If we specified the time_grids_to_store, check for spelling errors
        !
        if(allocated(domain%time_grids_to_store)) then
            do i = 1, size(domain%time_grids_to_store)
                if(.not. ((domain%time_grids_to_store(i) == '') .or. &
                    any(domain%time_grids_to_store(i) == gridded_output_variables_time))) then

                    write(log_output_unit, *) 'Unsupported variable in time_grids_to_store: ', &
                        trim(domain%time_grids_to_store(i))
                    call generic_stop
                end if
            end do
        end if
        !
        ! Determine which flow statistics we store that are not time-series.
        !
        if(.not. allocated(domain%nontemporal_grids_to_store)) then
            ! Mimic 'old behaviour' before we allowed specification of domain%nontempral_grids_to_store.
            ! In new code please avoid setting domain%record_max_U -- instead set domain%nontemporal_grids_to_store
            if(domain%record_max_U) then
                ! Old default where we store maxima
                domain%nontemporal_grids_to_store = [character(len=charlen):: 'max_stage', 'elevation0', 'manning_squared']
            else
                ! Old default where we don't store maxima or elevation or manning
                domain%nontemporal_grids_to_store = ['']
            end if
        end if
        ! Check for spelling errors in nontemporal_grids_to_store
        do i = 1, size(domain%nontemporal_grids_to_store)
            if(.not. ((domain%nontemporal_grids_to_store(i) == '') .or. &
                any(domain%nontemporal_grids_to_store(i) == gridded_output_variables_max_U) .or. &
                any(domain%nontemporal_grids_to_store(i) == gridded_output_variables_other) )) then

                write(log_output_unit, *) 'Unsupported variable in nontemporal_grids_to_store: ', &
                    trim(domain%nontemporal_grids_to_store(i))
                call generic_stop
            end if
        end do
        ! Define the names of nontemporal variables to store that require tracking flow statistics
        allocate(domain%max_U_variables(0))
        do i = 1, size(gridded_output_variables_max_U)
            if(any(domain%nontemporal_grids_to_store == gridded_output_variables_max_U(i))) &
                domain%max_U_variables = [domain%max_U_variables, gridded_output_variables_max_U(i)]
        end do
        if(size(domain%max_U_variables) > 0) then
            ! Ensure the first max_U variable is 'max_stage'
            if(domain%max_U_variables(1) /= 'max_stage') then
                write(log_output_unit, *) 'Error: "max_stage" MUST be added to domain%nontemporal_grids_to_store.'
                write(log_output_unit, *) 'While not logically essential, some older code assumes that if'
                write(log_output_unit, *) 'domain%max_U is allocated, then domain%max_U(:,:,1) has max-stage'
                call generic_stop
            end if

            ! Ensure that if we store "time_of_max_stage", then we also store "max_stage". This is 
            ! required with the current calculation method (although in principle isn't necessary).
            if(any(domain%max_U_variables == 'time_of_max_stage') .and. &
               all(domain%max_U_variables /= 'max_stage')) then
                write(log_output_unit,*) 'Error: if "time_of_max_stage" is being stored then "max_stage"'
                write(log_output_unit,*) 'must be stored too (limitation of current implementation).'
            end if
        end if

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

        !
        ! Storage for various flow maxima statistics
        !
        if(size(domain%max_U_variables) > 0) then
            allocate(domain%max_U(nx, ny, size(domain%max_U_variables)))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                ! The default value can be stored even at reduced precision.
                ! Although I expected -huge(1.0) could be stored, I ran into problems with netcdf
                ! which did not occur for slightly smaller magnitude numbers -- hence the factor 0.99
                ! in the following.
                domain%max_U(:,j,:) = -0.99_output_precision * huge(real(1.0, kind=output_precision))
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! Set the initial values for the min_stage to a large number
            do k = 1, size(domain%max_U_variables)
                if (domain%max_U_variables(k) == 'min_stage') then
                    ! No need to do this in parallel since already first touched
                    domain%max_U(:,:,k) = 0.99_output_precision * huge(real(1.0, kind=output_precision))
                end if
            end do
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

            ! (Manning friction)*(Manning friction) {Or (1/chezy_friction)^2 if domain%friction_type == 'chezy'}
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
                ! For these schemes, if domain%forcing_subroutine modifies the elevation, then
                ! we need to include elevation in domain%backup_U, and time-step the elevation like
                ! other variables. This adds some storage and computation, so we do not enable it by
                ! default.
                if(domain%support_elevation_forcing) then
                    allocate(domain%backup_U(nx, ny, STG:ELV))
                else
                    allocate(domain%backup_U(nx, ny, STG:VH))
                end if
                !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
                do k = 1, size(domain%backup_U, 3)
                    !$OMP DO SCHEDULE(STATIC)
                    do j = 1, ny
                        domain%backup_U(:,j,k) = ZERO_dp
                    end do
                    !$OMP END DO NOWAIT
                end do
                !$OMP END PARALLEL
            end if

        endif

        if(domain%timestepping_method == 'leapfrog_linear_plus_nonlinear_friction' .or. &
           domain%timestepping_method == 'leapfrog_nonlinear') then
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

        if(domain%timestepping_method == 'leapfrog_nonlinear') then

            ! A work array for nonlinear advection terms
            allocate(domain%advection_work(nx, ny, 1:domain%advection_work_dim3_size))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%advection_work(:,j, 1:domain%advection_work_dim3_size) = ZERO_dp
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

        end if

        !
        ! Setup for Chezy or Manning friction
        !
        select case(domain%friction_type)
        case('manning')
            ! Nothing to do here, but this avoids hitting "case default"
        case('chezy')
            if(domain%timestepping_method == 'cliffs') then
                write(log_output_unit, *) "ERROR: chezy_friction has not yet been implemented for our cliffs solver"
                call generic_stop
            end if
        case default
            write(log_output_unit, *) ' Error: Unrecognized value of domain%friction_type '
            call generic_stop
        end select

#ifdef DEBUG_ARRAY
        ! Optional for debugging -- this will be written to the netcdf file
        allocate(domain%debug_array(nx, ny))
        domain%debug_array = ZERO_dp
#endif

        if(domain%friction_with_ambient_fluxes) then
            ! Include ambient fluxes (e.g. from some external tide model) when computing friction terms for leapfrog solver.

            if(.not. any(domain%timestepping_method == ['leapfrog_linear_plus_nonlinear_friction'])) then
                write(log_output_unit, *) 'Error: domain%friction_with_ambient_fluxes = TRUE is only supported when using'
                write(log_output_unit, *) '       domain%timestepping_method = "leapfrog_linear_plus_nonlinear_friction" '
                call generic_stop
            end if

            allocate(domain%ambient_flux(nx, ny, UH:VH))
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, ny
                domain%ambient_flux(:, j, UH:VH) = ZERO_dp
            end do
            !$OMP END DO
            !$OMP END PARALLEL
        end if

        ! Setup the dispersive solver type's data
        if(domain%use_dispersion) then
            ! Dispersive terms are not supported for all solvers
            if(.not. any(domain%timestepping_method == [character(len=charlen):: &
                'linear', 'leapfrog_linear_plus_nonlinear_friction', 'leapfrog_nonlinear', &
                'rk2', 'euler', 'midpoint', 'cliffs']) ) then
                write(log_output_unit, *) 'Dispersion not supported for timestepping_method: ', trim(domain%timestepping_method)
                call generic_stop
            end if
            if(any(domain%timestepping_method == [character(len=charlen):: 'rk2', 'euler'])) then
                write(log_output_unit, *) 'Using finite volume solver ', trim(domain%timestepping_method), &
                    ' with dispersion. But the finite volume solver "midpoint" may be more accurate,', &
                    ' as the other options involve some first-order (rather than second order) approximations ', &
                    ' when including dispersion in the timestepping.'
            end if
            call domain%ds%setup(nx, ny, domain%is_staggered_grid)
        end if

        if(any(domain%nontemporal_grids_to_store == 'elevation_source_file_index')) then
            ! Include space to record which file was used to set the elevation at each site.
            ! The resulting array can be passed as an optional argument to 
            ! read_raster_mod::multi_raster_type%get_values. 
            allocate(domain%elevation_source_file_index(nx, ny))
            ! Typically only used once in serial, so no need for "first touch" openmp.
            domain%elevation_source_file_index = -1.0_dp
        end if


        if(create_output) call domain%create_output_files()

        if(verbose_) then
            write(domain%logfile_unit, *) ''
            write(domain%logfile_unit, *) 'dx: ', domain%dx
            write(domain%logfile_unit, *) 'nx: ', domain%nx
            write(domain%logfile_unit, *) 'lw: ', domain%lw
            write(domain%logfile_unit, *) 'lower-left: ', domain%lower_left
            write(domain%logfile_unit, *) 'upper-right:', domain%lower_left + domain%lw
            write(domain%logfile_unit, *) 'Total area: ', sum(domain%area_cell_y)
            write(domain%logfile_unit, *) 'distance_bottom_edge(1): ', domain%distance_bottom_edge(1)
            write(domain%logfile_unit, *) 'distance_left_edge(1): ', domain%distance_left_edge(1)
            write(domain%logfile_unit, *) ''
        end if

    end subroutine

    subroutine compute_depth_and_velocity(domain, min_depth)
        !! Compute domain%depth and domain%velocity, so they are consistent with domain%U.
        !! This only applies to the finite-volume solvers on non-staggered grids
        !! We also implement a mass conservation check in this routine, and zero velocities
        !! at sites with depth < minimum_allowed_depth.

        class(domain_type), intent(inout) :: domain
        real(dp), optional, intent(in) :: min_depth

        real(dp) :: depth_inv, min_depth_local
        integer(ip) :: i, j, masscon_error
EVOLVE_TIMER_START('compute_depth_and_velocity')
        ! Recompute depth and velocity
        ! Must be updated before the main loop when done in parallel

        if(present(min_depth)) then
            min_depth_local = min_depth
        else
            min_depth_local = minimum_allowed_depth
        end if

        masscon_error = 0_ip
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, min_depth_local) REDUCTION(+:masscon_error)
        masscon_error = 0_ip
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)

                domain%depth(i,j) = domain%U(i,j,STG) - domain%U(i,j,ELV)

                if(domain%depth(i,j) > min_depth_local) then
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

EVOLVE_TIMER_STOP('compute_depth_and_velocity')
    end subroutine

! Get details of the compute_fluxes routines here
#include "domain_mod_compute_fluxes_alternatives_include.f90"


    subroutine compute_fluxes(domain, max_dt_out)
        !!
        !! Compute the fluxes, and other things, in preparation for an update of domain%U.
        !! Update values of:
        !! domain%flux_NS, domain%flux_EW, domain%max_dt, domain%explicit_source,
        !! domain%explicit_source_VH_j_minus_1, domain%boundary_flux_store
        !!
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(inout) :: max_dt_out
            !! If provided, it is used to store the computed max_dt based on the CFL limit. Note this is
            !! not setting the timestep, it is simply storing the computed timestep.

        real(dp):: max_dt_out_local
        integer(ip) :: nx, ny, n_ext

        ! Need to have depth/velocity up-to-date for flux computation
        call domain%compute_depth_and_velocity()

EVOLVE_TIMER_START('compute_fluxes')
        !
        ! Flux computation -- optionally use a few different methods
        !
        select case(domain%compute_fluxes_inner_method)
        case('DE1_low_fr_diffusion')
            if(.not. domain%use_eddy_viscosity) then
                ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers
                call compute_fluxes_DE1_low_fr_diffusion(domain, max_dt_out_local)
            else
                ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers, and eddy-viscosity
                call compute_fluxes_DE1_low_fr_diffusion_eddyvisc(domain, max_dt_out_local)
            end if
        case('DE1')
            if(.not. domain%use_eddy_viscosity) then
                ! Something like ANUGA (sort of..)
                call compute_fluxes_DE1(domain, max_dt_out_local)
            else
                ! Something like ANUGA (sort of..), and eddy-viscosity
                call compute_fluxes_DE1_eddyvisc(domain, max_dt_out_local)
            end if
        case('DE1_low_fr_diffusion_upwind_transverse')
            if(.not. domain%use_eddy_viscosity) then
                ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers
                ! Also using an upwind treatment of the transverse flux term
                call compute_fluxes_DE1_low_fr_diffusion_upwind_transverse(domain, max_dt_out_local)
            else
                ! Something like ANUGA (sort of..) but with diffusion scaled for low-froude numbers
                ! Also using an upwind treatment of the transverse flux term, and eddy-viscosity
                call compute_fluxes_DE1_low_fr_diffusion_upwind_transverse_eddyvisc(domain, max_dt_out_local)
            end if
        case('DE1_upwind_transverse')
            if(.not. domain%use_eddy_viscosity) then
                ! Something like ANUGA (sort of..)
                ! Also using an upwind treatment of the transverse flux term
                call compute_fluxes_DE1_upwind_transverse(domain, max_dt_out_local)
            else
                ! Something like ANUGA (sort of..)
                ! Also using an upwind treatment of the transverse flux term, and eddy-viscosity
                call compute_fluxes_DE1_upwind_transverse_eddyvisc(domain, max_dt_out_local)
            end if
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
                    mask = (domain%nesting%is_priority_domain_not_periodic((1+n_ext):(nx-n_ext), ny-n_ext) == 1_ip) )
            end if

            ! East boundary
            if(domain%boundary_exterior(2)) then
                domain%boundary_flux_store_exterior(2) = sum(&
                    domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (domain%nesting%is_priority_domain_not_periodic(nx-n_ext, (1+n_ext):(ny-n_ext)) == 1_ip) )
            end if

            ! South boundary
            if(domain%boundary_exterior(3)) then
                domain%boundary_flux_store_exterior(3) = -sum(&
                    domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG),&
                    mask = (domain%nesting%is_priority_domain_not_periodic((1+n_ext):(nx-n_ext), 1+n_ext) == 1_ip) )
            end if

            ! West boundary
            if(domain%boundary_exterior(4)) then
                domain%boundary_flux_store_exterior(4) = -sum(&
                    domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG),&
                    mask = (domain%nesting%is_priority_domain_not_periodic(1+n_ext, (1+n_ext):(ny-n_ext)) == 1_ip) )
            end if
        else
            ! We are not doing nesting
            domain%boundary_flux_store_exterior = domain%boundary_flux_store
        end if

EVOLVE_TIMER_STOP('compute_fluxes')
    end subroutine

! Get details of the update_U routine here
#include "domain_mod_update_U_DE1_alternatives_include.f90"

    subroutine update_U(domain, dt)
        !!
        !! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
        !!
        class(domain_type), intent(inout):: domain
        real(dp), intent(in):: dt !! Timestep to advance

EVOLVE_TIMER_START('update_U')
        if(domain%friction_type == 'manning') then
            call update_U_manning(domain, dt)
        else if(domain%friction_type == 'chezy') then
            call update_U_chezy(domain, dt)
        end if
EVOLVE_TIMER_STOP('update_U')

        call domain%update_boundary()
        call domain%apply_forcing(dt)

    end subroutine

    subroutine update_boundary(domain)
        !!
        !! Routine to update all boundary conditions
        !!
        class(domain_type), intent(inout):: domain

EVOLVE_TIMER_START('boundary_update')

        if(associated(domain%boundary_subroutine)) then
            CALL domain%boundary_subroutine(domain)
        end if

EVOLVE_TIMER_STOP('boundary_update')

    end subroutine

    subroutine apply_forcing(domain, dt)
        !!
        !! Routine to call a forcing term (i.e. user-specified source term)
        !! In general this will require knowledge of the timestep
        !!
        class(domain_type), intent(inout):: domain
        real(dp), intent(in) :: dt
        integer(ip) :: i

EVOLVE_TIMER_START('apply_forcing')

        if(associated(domain%forcing_subroutine)) then
            ! In case domain%forcing_subroutine is associated, add it to the forcing array
            ! and disassociate it. This is allowed where only one forcing_subroutine exists,
            ! to support code that existed before domain%store_forcing() was written

            if(allocated(domain%forcing_terms_storage)) then
                ! This should only happen once
                write(domain%logfile_unit, *) &
                    'Error in forcing term setup: There are multiple forcing terms, but at least ', &
                    'one was defined without a subsequent call to domain%store_forcing(). If multiple ', &
                    'forcing terms are used, then domain%store_forcing() must be called after defining each.'
                call generic_stop
            end if

            call domain%store_forcing()
        end if

        if(allocated(domain%forcing_terms_storage)) then

            ! Apply the forcing terms in order of storage
            do i = 1, size(domain%forcing_terms_storage)
                ! Use the forcing data for the ith forcing term
                domain%forcing_subroutine => domain%forcing_terms_storage(i)%forcing_subroutine
                domain%forcing_context_cptr = domain%forcing_terms_storage(i)%forcing_context_cptr

                if(associated(domain%forcing_subroutine)) then
                    call domain%forcing_subroutine(domain, dt)
                end if

                ! Disassociate the forcing data
                domain%forcing_subroutine => NULL()
                domain%forcing_context_cptr = c_null_ptr
            end do
        end if

EVOLVE_TIMER_STOP('apply_forcing')

    end subroutine

    subroutine backup_quantities(domain)
        !!
        !! Copy domain%U to domain%backup_U
        !!

        class(domain_type), intent(inout):: domain
        integer(ip):: j, k

EVOLVE_TIMER_START('backup')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        do k = STG, size(domain%backup_U, 3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, domain%nx(2)
                domain%backup_U(:, j, k) = domain%U(:, j, k)
            end do
            !$OMP END DO NOWAIT
        end do
        !$OMP END PARALLEL

EVOLVE_TIMER_STOP('backup')

    end subroutine


! Get the timestepping routine details for evolve_one_step
#include "domain_mod_timestepping_alternatives_include.f90"


    subroutine evolve_one_step(domain, timestep)
        !
        ! Main 'high-level' evolve routine.
        ! The actual timestepping method used is determined by the value of
        ! domain%timestepping_method ('euler', 'rk2', 'rk2n', 'midpoint', ...).
        !
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep
            !! Optional for some domain%timestepping_method's, necessary for others.
            !! If provided it must satisfy the CFL condition.

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

        ! update_max_quantities usually occurs after the timestep, but that would miss the 
        ! initial condition (prior to any timestepping) if we didn't include the following
        if(domain%nsteps_advanced == 0) call domain%update_max_quantities()

        select case (timestepping_method)

        case ('static')
            ! Do nothing -- except update the time if the timestep was provided
            if(present(timestep)) then
                domain%time = domain%time + timestep
            end if
        case ('euler')
            if(present(timestep)) then
                call one_euler_step(domain, timestep)
            else
                call one_euler_step(domain)
            end if
        case ('rk2')
            if(present(timestep)) then
                call one_rk2_step(domain, timestep)
            else
                call one_rk2_step(domain)
            end if
        case('rk2n')
            if(present(timestep)) then
                call one_rk2n_step(domain, timestep)
            else
                call one_rk2n_step(domain)
            end if
        case ('midpoint')
            if(present(timestep)) then
                call one_midpoint_step(domain, timestep)
            else
                call one_midpoint_step(domain)
            end if
        case ('linear')
            if(present(timestep)) then
                call one_linear_leapfrog_step(domain, timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for linear evolve_one_step'
                call generic_stop()
            end if
        case ('leapfrog_linear_plus_nonlinear_friction')
            if(present(timestep)) then
                call one_leapfrog_linear_plus_nonlinear_friction_step(domain, timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for ', &
                    'leapfrog_linear_plus_nonlinear_friction evolve_one_step'
                call generic_stop()
            end if

        case('leapfrog_nonlinear')

            if(present(timestep)) then
                call one_leapfrog_nonlinear_step(domain, timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for ', &
                    'leapfrog_nonlinear evolve_one_step'
                call generic_stop()
            end if

        case ('cliffs')
            if(present(timestep)) then
                call one_cliffs_step(domain, timestep)
            else
                write(domain%logfile_unit,*) 'ERROR: timestep must be provided for ', &
                    'cliffs evolve_one_step'
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

    function volume_interior(domain) result(domain_volume)
        !!
        !! Convenience function to compute the volume of water in the 'interior' of the domain
        !! This involves all parts of the domain that are more than domain%exterior_cells_width
        !! from the edge. For multidomains, there are other purpose-built routines -- this one
        !! is appropriate for the single-domain case.
        !!
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

    subroutine volume_in_priority_domain(domain, domain_volume)
        !!
        !! Compute the volume of water in "priority domain" cells in the domain, ignoring
        !! parts of the domain that are within the domain%exterior_cells_width boundary (in
        !! general those regions are boundary condition sites).
        !! Control the summation order in a way that is reproducible even with openmp.
        !!
        class(domain_type), intent(inout):: domain
        real(force_double), intent(inout) :: domain_volume
        integer(ip) :: n_ext, ny, nx, j
        real(force_double) :: local_sum

TIMER_START('compute_volume_in_priority_domain')

            domain%volume_work = 0.0_force_long_double

            n_ext = domain%exterior_cells_width
            ny = domain%nx(2)
            nx = domain%nx(1)
            ! Integrate over the area with this domain = priority domain, which is not in periodic regions.
            ! Be careful about precision ( see calls to real(...., force_double) ).
            ! Use domain%volume_work (rather than a reduction) for reproducibility of the sum in parallel.
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, n_ext, nx, ny)
            !$OMP DO SCHEDULE(STATIC)
            do j = (n_ext+1), ny - n_ext
                local_sum = sum(&
                        real(domain%U( (n_ext+1):(nx-n_ext) ,j ,STG), force_long_double) - &
                        real(domain%U( (n_ext+1):(nx-n_ext) ,j ,ELV), force_long_double ), &
                    mask=(domain%nesting%is_priority_domain_not_periodic( (n_ext+1):(nx-n_ext) ,j) == 1_ip) )
                domain%volume_work(j) = real(domain%area_cell_y(j) * local_sum, force_long_double)
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            domain_volume = sum(domain%volume_work)

TIMER_STOP('compute_volume_in_priority_domain')

    end subroutine

    function mass_balance_interior(domain) result(mass_balance)
        !!
        !! Compute the volume on interior cells and add to the integrated boundary flux.
        !! This should sum to a constant in the absence of mass sources, for single domains. Note this relies on the time-stepping
        !! routines correctly computing the boundary_flux_time_integral.
        !! For multidomains, there are other purpose-built mass conservation routines (which take into account the "priority domain"
        !! information).
        !!
        class(domain_type), intent(in):: domain
        real(dp) :: mass_balance

        mass_balance = domain%volume_interior() + domain%boundary_flux_time_integral

    end function

    function linear_timestep_max(domain) result(timestep)
        !!
        !! Function to compute the max timestep allowed for the linear shallow water equations, using the provided CFL, assuming the
        !! leap-frog time-stepping scheme
        !!
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
                        sum(domain%distance_bottom_edge(j:j+1)),&
                        sum(domain%distance_left_edge(i:i+1)) ) / &
                    sqrt(gravity * max(domain%msl_linear-domain%U(i,j,ELV), minimum_allowed_depth)) )
            end do
        end do
        timestep = ts_max * domain%cfl

    end function

    function nonlinear_stationary_timestep_max(domain) result(timestep)
        !!
        !! Function to compute the max timestep allowed for the nonlinear shallow water equations, for a stationary flow (i.e. gravity
        !! wave only), using the provided CFL, similar to [cfl * min(cell_dx, cell_dy)/sqrt(gravity * max(depth)) ]. The numerical
        !! methods can often only take half this step, but we correct for that later.
        !! Results are not exactly the same as the gravity-wave timestep limit for nonlinear-FV solvers, because the latter compute
        !! wave speeds at edges rather than cell centres. But results are very close in general for a stationary flow.
        !!
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

    function stationary_timestep_max(domain) result(timestep)
        !! Convenience wrapper which uses 'linear_timestep_max' or 'nonlinear_stationary_timestep_max', depending on the numerical
        !! method, to compute the max timestep allowed for a stationary flow (i.e. just gravity waves, ignoring flow velocities for
        !! nonlinear shallow water equations).
        !! This is mainly useful when setting up models, as a guide to the timestep that one might be able to take.
        !! Beware the time-step calculations are not exactly the same as the results from the nonlinear-FV solvers (because the
        !! latter use the edge-based wave speed calculations). But they should be very similar.
        class(domain_type), intent(in):: domain
        real(dp) :: timestep

        ! Leapfrog type numerical methods
        character(len=charlen) :: leapfrog_type_solvers(3) = [ character(len=charlen) :: &
            'linear', 'leapfrog_linear_plus_nonlinear_friction', 'leapfrog_nonlinear']
        ! Typical DE1-type finite volume methods
        character(len=charlen) :: nonlinear_FV_solvers_1(3) = [ character(len=charlen) :: &
            'euler', 'rk2', 'midpoint' ]
        ! Less typical finite-volume methods
        character(len=charlen) :: nonlinear_FV_solvers_2(1) = [ character(len=charlen) :: &
            'rk2n' ]
        ! Cliffs
        character(len=charlen) :: cliffs_solver(1) = [ character(len=charlen) :: &
            'cliffs' ]

        if( any(domain%timestepping_method == leapfrog_type_solvers) ) then
            !
            ! Leap-frog type solvers
            !
            if(domain%timestepping_method == 'leapfrog_nonlinear' .or. &
                (.not. domain%linear_solver_is_truely_linear)) then
                ! Nonlinear variant
                timestep = domain%nonlinear_stationary_timestep_max()
            else
                ! Truely-linear variant
                timestep = domain%linear_timestep_max()
            end if

        else if(any(domain%timestepping_method == nonlinear_FV_solvers_1)) then
            !
            ! Standard nonlinear solvers in SWALS
            !
            timestep = domain%nonlinear_stationary_timestep_max() * 0.5_dp
        else if(any(domain%timestepping_method == nonlinear_FV_solvers_2)) then
            !
            ! The unusual 'rk2n' timestepping method -- we need (timestep/(rk2n_n_value-1)) to satisfy
            ! the typical FV CFL condition, because we take a number of substeps
            !
            timestep = domain%nonlinear_stationary_timestep_max() * 0.5_dp * (rk2n_n_value-1.0_dp)
        else if(any(domain%timestepping_method == cliffs_solver)) then
            !
            ! Cliffs (Tolkova)
            !
            timestep = domain%nonlinear_stationary_timestep_max()
        else
            write(log_output_unit, *) 'Error in stationary_timestep_max: Unrecognized timestepping_method'
            write(log_output_unit, *) trim(domain%timestepping_method)
            call generic_stop
        end if

    end function

    subroutine create_output_folders(domain, copy_code)
        !!
        !! Initialise output folders for the domain
        !!
        class(domain_type), intent(inout):: domain
        logical, optional, intent(in) :: copy_code

        character(len=charlen):: mkdir_command, cp_command, t1, t2, t3, t4, &
                                 output_folder_name, code_folder
        logical :: copy_code_local
        integer :: local_status

        ! Quick exit if folders already exist
        if(domain%output_folders_were_created) return

        if(present(copy_code)) then
            copy_code_local = copy_code
        else
            copy_code_local = .TRUE.
        end if

        ! Create output directory
        call date_and_time(t1, t2, t3)
        ! Get domain id as a character
        write(t3, domain_myid_char_format) domain%myid
        write(t4, '(I0.5)') domain%local_index

        output_folder_name = trim(domain%output_basedir) // '/RUN_ID' // trim(t3) // &
            '_' // trim(t4) // '_' // trim(t1) // '_' // trim(t2)
        domain%output_folder_name = output_folder_name
        call mkdir_p(domain%output_folder_name)

        if(copy_code_local) then
            ! Copy code to the output directory
            code_folder = trim(domain%output_folder_name) // '/Code'
            call mkdir_p(code_folder)

            cp_command = 'cp *.f* make* ' // trim(domain%output_folder_name) // '/Code'
            call execute_command_line(trim(cp_command), exitstat=local_status)
            !call system(trim(cp_command))
        end if
        domain%output_folders_were_created = .TRUE.

    end subroutine

    subroutine create_output_files(domain)
        !!
        !! Initialise output files for the domain
        !!
        class(domain_type), intent(inout):: domain

        character(len=charlen):: mkdir_command, cp_command, t1, t2, t3, t4, &
                                 output_folder_name
        integer(ip):: i, metadata_unit, natt
        character(len=charlen), allocatable :: attribute_names(:), attribute_values(:)


        call domain%create_output_folders()

        ! Get domain id as a character
        write(t3, domain_myid_char_format) domain%myid

        ! Make a filename to hold domain metadata, and write the metadata
        t1 = trim(domain%output_folder_name) // '/' // 'Domain_info_ID' // trim(t3) // '.txt'
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

        ! Make a time file_name. Store as ascii
        t1 = trim(domain%output_folder_name) // '/' // 'Time_ID' // trim(t3) // '.txt'
        open(newunit = domain%output_time_unit_number, file = t1)

        allocate(domain%output_variable_unit_number(domain%nvar))
        do i=1, domain%nvar
            ! Make output_file_name
            write(t1,'(I1)') i
            write(t2, *) 'Var_', trim(t1), '_ID', trim(t3)
            t1 = trim(domain%output_folder_name) // '/' // adjustl(trim(t2))

            ! Binary. Using 'stream' access makes it easy to read in R
            open(newunit = domain%output_variable_unit_number(i), file = t1, &
                access='stream', form='unformatted')
        end do

#else
        !
        ! Write to netcdf
        !

        t1 = trim(domain%output_folder_name) // '/' // 'Grid_output_ID' // trim(t3) // '.nc'
        ! Character attributes
        natt = 7
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
            xs = domain%x, ys = domain%y, &
            time_grids_to_store = domain%time_grids_to_store, &
            nontemporal_grids_to_store = domain%nontemporal_grids_to_store, &
            attribute_names=attribute_names, attribute_values=attribute_values)

        deallocate(attribute_names, attribute_values)

        ! If we have setup domain%nesting, then store locations where the current
        ! domain is the priority domain
        if(allocated(domain%nesting%priority_domain_index)) then
            call domain%nc_grid_output%store_priority_domain_cells(&
                domain%nesting%priority_domain_index,&
                domain%nesting%priority_domain_image,&
                domain%nesting%is_priority_domain_not_periodic,&
                domain%nesting%my_index,&
                domain%nesting%my_image)
        end if
#endif


    end subroutine create_output_files

    subroutine write_to_output_files(domain, time_only)
        !!
        !! Write domain data to output files (must call create_output_files first).
        !!
        class(domain_type), intent(inout):: domain
        logical, optional, intent(in) :: time_only
            !! If time_only=.TRUE., only write the time. This is used to avoid writing the main model grids. Typically useful when
            !! values at point-gauges are being recorded, and we want to store the time too, but it would take too much disk to
            !! store the model grids

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

        ! Time too, as ascii
        write(domain%output_time_unit_number, *) domain%time
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

TIMER_STOP('fileIO')

    end subroutine write_to_output_files

    subroutine divert_logfile_unit_to_file(domain)
        !!
        !! Set the domain%logfile_unit to point to an actual file, rather than stdout
        !!
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

    subroutine update_max_quantities(domain)
        !!
        !! Keep track of the maxima of stage, i.e. domain%U(:,:,STG)
        !!
        class(domain_type), intent(inout):: domain
        integer(ip):: j, k, i, ip1, jp1, toms

        real(dp) :: local_depth, local_depth_inv, arrival_stage, &
            depth_E, depth_N, depth_E_inv, depth_N_inv
        logical :: use_truely_linear_method

        if(domain%record_max_U) then
EVOLVE_TIMER_START('update_max_quantities')

            toms = findloc(domain%max_U_variables, value='time_of_max_stage', dim=1) ! 0 if not present

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, toms)
            do k = 1, size(domain%max_U_variables)

                select case(domain%max_U_variables(k))

                case('max_stage')

                    if(toms > 0) then
                        ! Store max-stage and the time of max-stage
                        !$OMP DO SCHEDULE(STATIC)
                        do j = domain%yL, domain%yU !1, domain%nx(2)
                            do i = domain%xL, domain%xU !1, domain%nx(1)
                                if(domain%U(i,j,STG) > domain%max_U(i,j,k)) then
                                    domain%max_U(i,j,k) = domain%U(i,j,STG)
                                    domain%max_U(i,j,toms) = domain%time
                                end if
                            end do
                        end do
                        !$OMP END DO NOWAIT
                    else
                        ! Only store max stage
                        !$OMP DO SCHEDULE(STATIC)
                        do j = domain%yL, domain%yU !1, domain%nx(2)
                            do i = domain%xL, domain%xU !1, domain%nx(1)
                                domain%max_U(i,j,k) = max(domain%max_U(i,j,k), domain%U(i,j,STG))
                            end do
                        end do
                        !$OMP END DO NOWAIT
                    end if

                case('time_of_max_stage')
                    ! Treated in case('max_stage'), do nothing else.

                case('max_speed')

                    if(domain%is_staggered_grid) then
                        ! Speed on staggered grid
                        use_truely_linear_method = (domain%linear_solver_is_truely_linear .and. &
                            any(domain%timestepping_method == &
                            [character(len=charlen) :: 'linear', 'leapfrog_linear_plus_nonlinear_friction']))

                        !$OMP DO SCHEDULE(STATIC)
                        do j = domain%yL, domain%yU !1, domain%nx(2)
                            do i = domain%xL, domain%xU !1, domain%nx(1)
                                ! Get depths at UH point (E) and VH point (N).
                                ! The interpretation of 'depth' varies depending on whether we are using the 'truely-linear'
                                ! treatment of the pressure gradient term in the linear shallow water equations.
                                ip1 = min(i+1, domain%nx(1))
                                jp1 = min(j+1, domain%nx(2))
                                if(use_truely_linear_method) then
                                    ! Depth is always relative to MSL.
                                    depth_E = 0.5_dp * sum( (domain%msl_linear - domain%U(i:ip1,j,ELV)))
                                    depth_N = 0.5_dp * sum( (domain%msl_linear - domain%U(i,j:jp1,ELV)))
                                else
                                    ! Spatially averaged depth on staggered grid
                                    depth_E = 0.5_dp * sum( (domain%U(i:ip1,j,STG) - domain%U(i:ip1,j,ELV)))
                                    depth_N = 0.5_dp * sum( (domain%U(i,j:jp1,STG) - domain%U(i,j:jp1,ELV)))
                                end if

                                depth_E_inv = merge(1.0_dp/depth_E, 0.0_dp, depth_E > minimum_allowed_depth)
                                depth_N_inv = merge(1.0_dp/depth_N, 0.0_dp, depth_N > minimum_allowed_depth)

                                domain%max_U(i,j,k) = max(domain%max_U(i,j,k), sqrt(&
                                    domain%U(i,j,UH) * domain%U(i,j,UH) * depth_E_inv**2  + &
                                    domain%U(i,j,VH) * domain%U(i,j,VH) * depth_N_inv**2))
                            end do
                        end do
                        !$OMP END DO NOWAIT
                    else
                        ! Speed on co-located grid
                        !$OMP DO SCHEDULE(STATIC)
                        do j = domain%yL, domain%yU !1, domain%nx(2)
                            do i = domain%xL, domain%xU !1, domain%nx(1)
                                local_depth = domain%U(i,j,STG) - domain%U(i,j,ELV)
                                local_depth_inv = merge(1.0_dp/local_depth, 0.0_dp, local_depth > minimum_allowed_depth)
                                domain%max_U(i,j,k) = max(&
                                    domain%max_U(i,j,k), &
                                    sqrt((domain%U(i,j,UH)**2 + domain%U(i,j,VH)**2)*local_depth_inv**2))
                            end do
                        end do
                        !$OMP END DO NOWAIT
                    end if

                case('max_flux')

                    ! Track max flux
                    !$OMP DO SCHEDULE(STATIC)
                    do j = domain%yL, domain%yU !1, domain%nx(2)
                        do i = domain%xL, domain%xU !1, domain%nx(1)
                            domain%max_U(i,j,k) = max(&
                                domain%max_U(i,j,k), &
                                sqrt(domain%U(i,j,UH)**2 + domain%U(i,j,VH)**2))
                        end do
                    end do
                    !$OMP END DO NOWAIT

                case('arrival_time')

                    !! Track the arrival time
                    !$OMP DO SCHEDULE(STATIC)
                    do j = domain%yL, domain%yU !1, domain%nx(2)
                        do i = domain%xL, domain%xU !1, domain%nx(1)
                            local_depth = domain%U(i,j,STG) - domain%U(i,j,ELV)
                            arrival_stage = domain%msl_linear + domain%arrival_stage_threshold_above_msl_linear
                            if( (local_depth > minimum_allowed_depth) .and. &
                                (domain%U(i,j,STG) > arrival_stage  ) .and. &
                                ! The arrival time was initialised to a negative number
                                ( domain%max_U(i,j,k) < 0.0_dp ) ) then
                                domain%max_U(i,j,k) = domain%time
                            end if
                        end do
                    end do
                    !$OMP END DO NOWAIT

                case('min_stage')

                    !! Track the minimum stage
                    !$OMP DO SCHEDULE(STATIC)
                    do j = domain%yL, domain%yU !1, domain%nx(2)
                        do i = domain%xL, domain%xU !1, domain%nx(1)
                            domain%max_U(i,j,k) = min(domain%max_U(i,j,k), domain%U(i,j,STG))
                        end do
                    end do
                    !$OMP END DO NOWAIT

                case default
                    write(log_output_unit, *) ' Unknown value in domain%max_U_variables '
                    call generic_stop
                end select

            end do

            !$OMP END PARALLEL

EVOLVE_TIMER_STOP('update_max_quantities')
        end if

    end subroutine update_max_quantities

    subroutine write_max_quantities(domain)
        !!
        !! Write max quantities to a file (usually just called once at the end of a simulation).
        !! Currently we only write stage, followed by the elevation. Although the latter usually doesn't evolve, we generally want
        !! both the max stage and the elevation for plotting purposes, so it is saved here too.
        !!
        class(domain_type), intent(in):: domain
        character(len=charlen):: max_quantities_filename
        integer:: i, j, k

#ifdef NONETCDF
        if(domain%record_max_U) then
            ! Use home-brew binary format to output max stage (only)
            ! This is not up-to-date (consider removing support for the home-brew file format)

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
        end if
#else

            ! Save the nontemporal variables that require tracking the flow
            do j = 1, size(domain%max_U_variables)
                call domain%nc_grid_output%store_static_variable(domain%max_U_variables(j), domain%max_U(:,:,j))
            end do

            ! Save manning^2
            if(allocated(domain%manning_squared) .and. &
                any(domain%nontemporal_grids_to_store == 'manning_squared')) then
                call domain%nc_grid_output%store_static_variable('manning_squared', domain%manning_squared)
            end if

            ! Save elevation
            if(any(domain%nontemporal_grids_to_store == 'elevation0')) then
                call domain%nc_grid_output%store_static_variable('elevation0', domain%U(:,:,ELV))
            end if

            ! Save the elevation source file index
            if(allocated(domain%elevation_source_file_index) .and. &
               any(domain%nontemporal_grids_to_store == 'elevation_source_file_index')) then
                call domain%nc_grid_output%store_static_variable('elevation_source_file_index', &
                    domain%elevation_source_file_index)
            end if

#endif

    end subroutine

    subroutine is_in_priority_domain(domain, x, y, is_in_priority)
        !! Determine whether points x,y are in the priority domain of "domain".
        class(domain_type), intent(in) :: domain
        real(dp), intent(in) :: x(:), y(:) !! x/y coordinates
        logical, intent(inout) :: is_in_priority(:) !! Has the same length as x and y

        integer(ip) :: i, i0, j0

        if(domain%nesting%my_index == 0) then
            write(log_output_unit, *) ' Error in is_in_priority_domain: Nesting is not active'
            call generic_stop
        end if

        if(size(x, kind=ip) /= size(y, kind=ip) .or. size(x, kind=ip) /= size(is_in_priority, kind=ip)) then
            write(log_output_unit, *) ' Error in is_in_priority_domain: All inputs must have same length'
            call generic_stop
        end if

        is_in_priority = .false.
        do i = 1, size(x, kind=ip)

            i0 = ceiling( (x(i) - domain%lower_left(1))/domain%dx(1) )
            j0 = ceiling( (y(i) - domain%lower_left(2))/domain%dx(2) )

            if(i0 > 0 .and. i0 <= domain%nx(1) .and. j0 > 0 .and. j0 <= domain%nx(2)) then
                ! Here we INCLUDE points in any periodic regions that receive from their own domain.
                is_in_priority(i) = ( &
                    domain%nesting%priority_domain_index(i0, j0) == domain%nesting%my_index .and. &
                    domain%nesting%priority_domain_image(i0, j0) == domain%nesting%my_image )
            end if

        end do

    end subroutine


    subroutine setup_point_gauges(domain, xy_coords, time_series_var, static_var, gauge_ids, &
        attribute_names, attribute_values)
        !!
        !! Set up point gauges, which record values of variables at cells nearest the given xy points.
        !!

        class(domain_type), intent(inout):: domain !! The domain within which we will record outputs (in domain%U)
        real(dp), intent(in) :: xy_coords(:,:)
            !! numeric array of x,y coordinates with 2 rows and as many columns as points. These must all be inside the
            !! extents of the domain
        integer(ip), optional, intent(in):: time_series_var(:)
            !! (optional) array with the indices of U to store at gauges each timestep, Default is [STG, UH, VH] to store
            !! stage, uh, vh
        integer(ip), optional, intent(in):: static_var(:)
            !! (optional) array with the indices of U to store at gauges only once. Default is [ELV] to store elevation
        real(dp), optional, intent(in):: gauge_ids(:)
            !! (optional) an REAL ID for each gauge. Default gives (1:size(xy_coords(1,:))) * 1.0. Even though integers are
            !! natural for IDS, we store as REAL, just to avoid precision loss issues if someone decides to use large numbers.
        character(charlen), optional, intent(in):: attribute_names(:), attribute_values(:)
            !! character vectors with the names and values of of global attributes for the netdf file

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
            allocate(tsv(size(time_series_var, kind=ip)))
            tsv = time_series_var
        else
            ! Default case -- store stage/uh/vh every output step
            allocate(tsv(3))
            tsv = [STG, UH, VH]
        end if

        ! Set variables in domain%U(:,:,k) that are stored only once (at the start)
        if(present(static_var)) then
            allocate(sv(size(static_var, kind=ip)))
            sv = static_var
        else
            ! Default case -- store elevation once
            allocate(sv(1))
            sv = [ELV]
        end if

        ! Setup or create some 'IDs' for each gauge. These come in handy.
        if(present(gauge_ids)) then
            if(size(gauge_ids, kind=ip) /= size(xy_coords(1,:), kind=ip)) then
                write(domain%logfile_unit,*) 'Number of gauge ids does not equal number of coordinates'
                flush(domain%logfile_unit)
                call generic_stop
            end if
            allocate(gauge_ids_local(size(gauge_ids, kind=ip)))
            gauge_ids_local = gauge_ids
        else
            ! Default case -- give sequential integer ids
            allocate(gauge_ids_local(size(xy_coords(1,:), kind=ip)))
            do i = 1, size(gauge_ids_local, kind=ip)
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
            allocate(priority_gauges(size(xy_coords(1,:), kind=ip)))
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

    subroutine write_gauge_time_series(domain)
        !!
        !! Write the values of point gauges to a file
        !!

        class(domain_type), intent(inout):: domain
TIMER_START('write_gauge_time_series')
        if(allocated(domain%point_gauges%time_series_values)) then
            call domain%point_gauges%write_current_time_series(domain%U, domain%time)
        end if
TIMER_STOP('write_gauge_time_series')
    end subroutine

    subroutine finalise_domain(domain)
        !!
        !! Routine to call once we no longer need the domain.
        !! One case where this is important is when using netcdf output -- since if the files are not closed, then they may not be
        !! completely written out.
        !!
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
            if(is_open) flush(i)
        end do

#ifndef NONETCDF
        ! Close the grids netcdf file
        call domain%nc_grid_output%finalise()
#endif

    end subroutine


    subroutine nesting_boundary_flux_integral_multiply(domain, c)
        !! If doing nesting, we need to track boundary flux integrals through send/recv regions.
        !! This routine multiplies all boundary_flux_integrals in send/recv regions by a constant 'c', which
        !! is required for the flux tracking.

        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: c

        integer(ip) :: i

TIMER_START('nesting_boundary_flux_integral_multiply')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, c)

        ! The 'send communicator' fluxes
        if(allocated(domain%nesting%send_comms)) then
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%send_comms, kind=ip)
                call domain%nesting%send_comms(i)%boundary_flux_integral_multiply(c)
            end do
            !$OMP END DO NOWAIT
        end if


        ! The 'receive communicator' fluxes
        if(allocated(domain%nesting%recv_comms)) then
            !$OMP DO SCHEDULE (DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms, kind=ip)
                call domain%nesting%recv_comms(i)%boundary_flux_integral_multiply(c)
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

TIMER_STOP('nesting_boundary_flux_integral_multiply')
    end subroutine


    subroutine nesting_boundary_flux_integral_tstep(domain, dt, &
        flux_NS, flux_NS_lower_index,&
        flux_EW, flux_EW_lower_index,&
        var_indices, flux_already_multiplied_by_dx)
        !! If doing nesting, we probably want to track boundary flux integrals through send/recv regions -- e.g. to allow for flux
        !! correction.
        !! This routine adds "dt * current_value_of_fluxes * dx" to all boundary_flux_integrals in send/recv regions, as required for
        !! flux tracking.
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

        class(domain_type), intent(inout) :: domain
        integer(ip), intent(in) :: flux_NS_lower_index, flux_EW_lower_index, var_indices(2)
        real(dp), intent(in) :: dt, flux_NS(:,:,:), flux_EW(:,:,:)
        logical, intent(in) :: flux_already_multiplied_by_dx

        integer(ip) :: i

EVOLVE_TIMER_START('nesting_boundary_flux_integral_tstep')

        ! Apply to both the send and recv comms. This means that after
        ! communication, we can compare fluxes that were computed with
        ! different numerical methods, and potentially apply flux correction.

        ! Deal with 'send comminicators'
        if(allocated(domain%nesting%send_comms)) then
            do i = 1, size(domain%nesting%send_comms, kind=ip)
                call domain%nesting%send_comms(i)%boundary_flux_integral_tstep( dt,&
                    flux_NS, flux_NS_lower_index, domain%distance_bottom_edge, &
                    flux_EW, flux_EW_lower_index, domain%distance_left_edge, &
                    var_indices, flux_already_multiplied_by_dx)
            end do
        end if

        ! Deal with 'receive comminicators'
        if(allocated(domain%nesting%recv_comms)) then
            do i = 1, size(domain%nesting%recv_comms, kind=ip)
                call domain%nesting%recv_comms(i)%boundary_flux_integral_tstep( dt,&
                    flux_NS, flux_NS_lower_index, domain%distance_bottom_edge, &
                    flux_EW, flux_EW_lower_index, domain%distance_left_edge, &
                    var_indices, flux_already_multiplied_by_dx)
            end do
        end if
EVOLVE_TIMER_STOP('nesting_boundary_flux_integral_tstep')

    end subroutine


    subroutine nesting_flux_correction_everywhere(domain, all_dx_md, all_timestepping_methods_md, fraction_of)
        !!
        !! Apply nesting flux correction all throughout the domain, even at non-priority-domain sites.
        !! This is a key routine for mainintaining mass conservation in the multidomain case.
        !! The general idea is that after we have received from other domains, we can apply flux correction "just like the other
        !! domains would do". This means we avoid having to do multiple nesting communications to make the fluxes consistent.
        !!
        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: all_dx_md(:,:,:)
            !! contains dx for all domains in the multidomain
        integer(ip), intent(in) :: all_timestepping_methods_md(:,:)
            !! records the timestepping_method index for all domains in the multidomain
        real(dp), optional, intent(in) :: fraction_of
            !! Apply some fraction of the flux correction. By default apply completely (i.e. 1.0)

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

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, all_dx_md, all_timestepping_methods_md, fraction_of_local)

            !
            ! NORTH BOUNDARIES.
            !
            ! By doing each box direction separately (i.e. all north, then all south, ...), we can ensure that multiple openmp
            ! threads do not try to update the same domain%U cells at once.
            !
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%recv_comms, kind=ip)

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
                            !! FIXME: What about the (unusual) case with a nesting_ratio = 2?
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
            do i = 1, size(domain%nesting%recv_comms, kind=ip)

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
            do i = 1, size(domain%nesting%recv_comms, kind=ip)
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
            do i = 1, size(domain%nesting%recv_comms, kind=ip)
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
            do i = 1, size(domain%nesting%recv_comms, kind=ip)

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
            do i = 1, size(domain%nesting%recv_comms, kind=ip)

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
            ! @param var1 index of first variable to update (normally STG=1) -- negative for no update
            ! @param varN index of last variable to update (normally VH=3, except for staggered grid, where it is STG=1, or where we
            ! don't do any updates, in which case it is negative and less than var1)
            subroutine compute_offset_inside_or_out(dm, dm_outside, out_index, out_image, nbr_index, nbr_image, &
                    my_index, my_image, is_ew, dx_ratio, equal_sizes, var1, varN)

                integer(ip), intent(out) :: dm
                integer(ip), intent(in) :: dm_outside, out_index, out_image, nbr_index, nbr_image, my_index, my_image
                logical, intent(in) :: is_ew
                integer(ip), intent(out) :: dx_ratio
                logical, intent(out) :: equal_sizes
                integer(ip), intent(out) :: var1, varN

                integer(ip) :: i1, cor_index, cor_image, cor_tsi, notcor_tsi
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
                    ! Need a rule to decide which one to correct. All images should see the same rule.
#ifdef OLD_PROCESS_DATA_TO_SEND_B4FEB22
                    ! Solution: Select dm based on the order of the index, with tie-breaking by image
                    ! OLD METHOD -- it has a weakness:
                    !     With this rule, if we run two models that are identical except image indices are reordered,
                    !     then this could behave differently. While there is no 'correct' choice, seems better to have
                    !     the same solution irrespective (perhaps based on geometric criteria)? Note that if the equal
                    !     sized domains are using local time-stepping, so have different time-steps, then presumably
                    !     flux correction will be needed.
                    !     Alternative: using the value of `dm_outside` [e.g. dm = merge(0, dm_outside, dm_outside > 0)].
                    !
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
#else
                    ! This approach will not depend on the ordering of md%domains, or their distribution among images
                    dm = merge(0, dm_outside, dm_outside > 0)
#endif
                end if

                ! Get indices for 'the domain to be corrected'  = cor_index, cor_image
                if(dm == dm_outside) then
                    cor_image = out_image
                    cor_index = out_index
                    ! Get the timestepping-method-index for both the 'to be corrected' and 'not to be corrected' domains
                    cor_tsi = all_timestepping_methods_md(cor_index, cor_image)
                    notcor_tsi = all_timestepping_methods_md(nbr_index, nbr_image)
                else
                    cor_image = nbr_image
                    cor_index = nbr_index
                    ! Get the timestepping-method-index for both the 'to be corrected' and 'not to be corrected' domains
                    cor_tsi = all_timestepping_methods_md(cor_index, cor_image)
                    notcor_tsi = all_timestepping_methods_md(out_index, out_image)
                end if

                ! Get dx ratio of 'domain to be corrected' vs 'my domain', with round-off protection as above
                if(all_dx_md(i1, cor_index, cor_image) > 0.8_dp * all_dx_md(i1, my_index, my_image)) then
                    ! Note the factor 0.8 is a kludge to protect against round-off. The cell-size ratio should either
                    ! be an integer, or 1/integer
                    dx_ratio = nint(all_dx_md(i1, cor_index, cor_image)/all_dx_md(i1, my_index, my_image))
                else
                    dx_ratio = 0
                end if


                if(timestepping_metadata(cor_tsi)%flux_correction_is_unsupported .or. &
                   timestepping_metadata(notcor_tsi)%flux_correction_is_unsupported) then
                   ! One or other solver cannot use flux correction. This is typically the case
                   ! when the solver doesn't track mass fluxes, so it stores fluxes as zero, and there
                   ! would be a large spurious correction if it neighbours another solver that does track fluxes.
                   ! Prevent looping using a reversed loop range, i.e.
                   !    do i = var1, varN
                   ! never enters the loop if varN < var1
                   var1 = -1
                   varN = -2
                else
                    ! Some flux correction should be applied
                    if(timestepping_metadata(cor_tsi)%flux_correction_of_mass_only) then
                        ! Treat solvers that can only do mass-flux correction, not advective momentum fluxes.
                        ! (e.g. linear solvers where the momentum advection terms are ignored)
                        var1 = STG
                        varN = STG
                    else
                        ! Typical case
                        var1 = STG
                        varN = VH
                    end if

                end if

            end subroutine

    end subroutine

    subroutine match_geometry_to_parent(domain, parent_domain, lower_left, &
        upper_right, dx_refinement_factor, timestepping_refinement_factor, &
        rounding_method, recursive_nesting)
        !!
        !! Convenience routine for setting up nested domain geometries.
        !! Set the domain "lower-left, lw, dx" etc, in a way that ensures the domain can nest with its parent domain, while having
        !! an extent and resolution close to the desired values.
        !!
        class(domain_type), intent(inout) :: domain
            !! A domain_type which is unallocated, and does not have lw, dx, nx set
        class(domain_type), intent(in) :: parent_domain
            !! Another domain_type which DOES have lw, dx, nx set. We want to
            !! give the new domain boundaries that are consistent with the parent domain (for nesting purposes),
            !! i.e. the corners of each parent domain cell correspond to corners of child domain cells
        real(dp), intent(in) :: lower_left(2), upper_right(2)
            !! The desired lower-left and upper-right of the new domain. We will change this for consistency with nesting
        integer(ip), intent(in) :: dx_refinement_factor
            !! Integer such that parent_domain%dx/dx_refinement_factor = new_domain%dx
        integer(ip), intent(in) :: timestepping_refinement_factor
            !! How many time-steps should the new domain take, for each global time-step in the multidomain.
        character(*), intent(in), optional :: rounding_method
            !! optional character controlling how we adjust lower-left/upper-right.
            !! If rounding_method = 'expand', then we adjust the new domain lower-left/upper-right so that the provided
            !! lower-left/upper-right are definitely contained in the new domain. If rounding_method = 'nearest' (DEFAULT), we move
            !! lower-left/upper-right onto the nearest cell corner of the parent domain. This can be preferable if we want to have
            !! multiple child domains which share boundaries with each other -- but does not ensure the provided
            !! lower-left/upper-right are within the new domain
        logical, intent(in), optional :: recursive_nesting
            !! If TRUE(default), the domain's dx_refinement_factor will be multiplied
            !! by its parent domain's dx_refinement_factor before storing in domain%dx_refinement_facor. This will not affect
            !! domain%dx. But it should should allow the domain to communicate 'cleanly' with "parents of its parent" if it is
            !! partitioned, SO LONG AS the lower_left is also on a corner of the earlier generation's domain. If .FALSE., then just
            !! use the provided dx_refinement_factor, and tell the domain that it's parent-domain's dx_refinement_factor=1.
            !! This is less likely to allow communicating with grandparent domains, but might allow for a more efficient split-up
            !! of the domain.

        real(force_double) :: ur(2), parent_domain_dx(2), ll(2)
        character(len=charlen) :: rounding
        logical :: recursive_nest
    
        rounding = 'nearest'
        if(present(rounding_method)) rounding = rounding_method

        recursive_nest = .true.
        if(present(recursive_nesting)) recursive_nest = recursive_nesting

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
            ll = real(parent_domain%lower_left, force_double) + &
                floor((lower_left - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

            ! upper_right = (domain%lower_left + domain%lw) is on a cell boundary of the parent domain
            ur = real(parent_domain%lower_left, force_double) + &
                ceiling((upper_right - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

        case('nearest')
            ! Find the 'nearest' match in parent domain. This might mean we reduce the requested size of the domain
            ll = real(parent_domain%lower_left, force_double) + &
                nint((lower_left - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

            ur = real(parent_domain%lower_left, force_double) + &
                nint((upper_right - parent_domain%lower_left)/parent_domain_dx)*parent_domain_dx

        case default

            write(domain%logfile_unit, *) ' rounding_method ', TRIM(rounding), ' not recognized '
            call generic_stop()

        end select

        ! Now we can set the child domain's properties
        domain%lower_left = ll
        domain%lw =  ur - ll
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

    subroutine receive_halos(domain, p2p, skip_if_time_before_this_time, skip_if_time_after_this_time, &
        skip_if_unequal_cell_ratios)
        !!
        !! Loop over domain%nesting%receive_comms, and receive the halo data
        !!
        class(domain_type), intent(inout) :: domain
        type(p2p_comms_type), intent(in) :: p2p
        real(dp), optional, intent(in) :: skip_if_time_before_this_time
        real(dp), optional, intent(in) :: skip_if_time_after_this_time
        logical, optional, intent(in) :: skip_if_unequal_cell_ratios

        integer(ip) :: i, j, ii, jj, iL, iU, jL, jU, ip1, jp1, nbr_j, nbr_ti
        real(dp) :: elev_lim

TIMER_START('receive_halos')

        !$OMP PARALLEL DEFAULT(PRIVATE) &
        !$OMP SHARED(domain, p2p, skip_if_time_before_this_time, skip_if_time_after_this_time, skip_if_unequal_cell_ratios)
        !$OMP DO SCHEDULE(DYNAMIC)
        do i = 1, size(domain%nesting%recv_comms, kind=ip)

            !! Avoid use of type-bound-procedure in openmp region
            !call domain%nesting%recv_comms(i)%process_received_data(domain%U)
            call process_received_data(domain%nesting%recv_comms(i), p2p, domain%U, &
                skip_if_time_before_this_time = skip_if_time_before_this_time, &
                skip_if_time_after_this_time  = skip_if_time_after_this_time, &
                skip_if_unequal_cell_ratios = skip_if_unequal_cell_ratios)
        end do
        !$OMP END DO

        !$OMP DO SCHEDULE(DYNAMIC)
        do i = 1, size(domain%nesting%recv_comms, kind=ip)
            !
            ! Need some more processing of staggered-grid receive regions to deal with potential wetting/drying issues.
            ! This cannot be included in the loop above because for some solvers, we refer to stage values outside
            ! the receive region (extreme values of indices ip1, jp1), which could lead to a race condition.
            ! However one could restructure the code to put much of this content inside the loop above (truely linear solvers,
            ! updates of non-boundary cells).
            !
            if( (.not. domain%is_staggered_grid) .or. (.not. domain%nesting%recv_comms(i)%recv_active) ) cycle
            ! Below here we are definitely using a leapfrog-type solver.

            ! Need to avoid receiving non-zero UH/VH at wet/dry boundaries, in a way that would violate
            ! the solver logic

            ! Indices of the region that received data
            iL = domain%nesting%recv_comms(i)%recv_inds(1,1)
            iU = domain%nesting%recv_comms(i)%recv_inds(2,1)
            jL = domain%nesting%recv_comms(i)%recv_inds(1,2)
            jU = domain%nesting%recv_comms(i)%recv_inds(2,2)

            ! Depending on the numerical method, we need to treat the wet/dry issues differently
            if( any( domain%timestepping_method == [character(len=charlen) :: &
                    'linear', 'leapfrog_linear_plus_nonlinear_friction'] ) .and. &
                domain%linear_solver_is_truely_linear) then
                ! Here we treat the 'truely linear' solvers where all cells above
                !     (domain%msl_linear - minimum_allowed_depth)
                ! are dry. Flux is zero if either neighbouring elevation is dry

                ! Elevation above which flux = 0 for truely-linear domain
                elev_lim = (domain%msl_linear - minimum_allowed_depth)

                ! Loop over the received region, and ensure UH/VH on wet-dry boundaries are 0
                do jj = jL, jU

                    ! jj+1, avoiding out-of-bounds
                    jp1 = min(jj+1, size(domain%U, 2, kind=ip))

                    do ii = iL, iU

                        ! ii+1, avoiding out-of-bounds
                        ip1 = min(ii+1, size(domain%U, 1, kind=ip))

                        ! Condition for zero EW flux on staggered grid
                        if( (domain%U(ii,jj,ELV) >= elev_lim) .or. &
                            (domain%U(ip1,jj,ELV) >= elev_lim) ) then
                            domain%U(ii,jj,UH) = 0.0_dp
                        end if

                        ! Condition for zero NS flux on staggered grid
                        if( (domain%U(ii,jj,ELV) >= elev_lim) .or. &
                            (domain%U(ii,jp1,ELV) >= elev_lim) ) then
                            domain%U(ii,jj,VH) = 0.0_dp
                        end if

                    end do
                end do

            else if(domain%timestepping_method /= 'leapfrog_nonlinear') then
                ! This treats 'not-truely-linear' staggered grid solvers, which are not fully nonlinear.
                ! In this case the UH/VH terms should be zero if either neighbouring cell has
                ! ( stage <= (elev + minimum_allowed_depth) )

                ! Loop over the received region, and ensure UH/VH on
                ! wet-dry boundaries are 0
                do jj = jL, jU

                    ! jj+1, avoiding out-of-bounds
                    jp1 = min(jj+1, size(domain%U, 2, kind=ip))

                    do ii = iL, iU

                        ! ii+1, avoiding out-of-bounds
                        ip1 = min(ii+1, size(domain%U, 1, kind=ip))

                        ! No EW flux if either neighbouring cell is dry
                        if( (domain%U(ii,jj,STG)  - domain%U(ii,jj,ELV)  <= minimum_allowed_depth) .or. &
                            (domain%U(ip1,jj,STG) - domain%U(ip1,jj,ELV) <= minimum_allowed_depth) ) then
                            domain%U(ii,jj,UH) = 0.0_dp
                        end if

                        ! No NS flux if either neighbouring cell is dry
                        if( (domain%U(ii,jj,STG)  - domain%U(ii,jj,ELV)  <= minimum_allowed_depth) .or. &
                            (domain%U(ii,jp1,STG) - domain%U(ii,jp1,ELV) <= minimum_allowed_depth) ) then
                            domain%U(ii,jj,VH) = 0.0_dp
                        end if

                    end do
                end do

            else if(domain%timestepping_method == 'leapfrog_nonlinear') then
                ! For the fully nonlinear staggered-grid solver, we just need to ensure there is no outflow
                ! from dry cells
                do jj = jL, jU

                    ! jj+1, avoiding out-of-bounds
                    jp1 = min(jj+1, size(domain%U, 2, kind=ip))

                    do ii = iL, iU

                        ! ii+1, avoiding out-of-bounds
                        ip1 = min(ii+1, size(domain%U, 1, kind=ip))

                        ! No easterly outflow from dry cell
                        if( (domain%U(ii,jj,STG)  - domain%U(ii,jj,ELV)  <= minimum_allowed_depth) .and. &
                            (domain%U(ii,jj, UH) > 0.0_dp) ) domain%U(ii,jj,UH) = 0.0_dp

                        ! No westerly outflow from dry cell -- not applicable if ii+1 extends outside eastern boundary
                        if( (ip1 == ii + 1) .and. &
                            (domain%U(ip1,jj,STG) - domain%U(ip1,jj,ELV) <= minimum_allowed_depth) .and. &
                            (domain%U(ii,jj, UH) < 0.0_dp) ) domain%U(ii,jj,UH) = 0.0_dp

                        ! No northerly outflow from dry cell
                        if( (domain%U(ii,jj,STG)  - domain%U(ii,jj,ELV)  <= minimum_allowed_depth) .and. &
                            (domain%U(ii,jj,VH) > 0.0_dp ) ) domain%U(ii,jj,VH) = 0.0_dp

                        ! No southerly outflow from dry cell -- not applicable if jj+1 extends outside northern boundary
                        if( (jp1 == jj + 1) .and. &
                            (domain%U(ii,jp1,STG) - domain%U(ii,jp1,ELV) <= minimum_allowed_depth) .and. &
                            (domain%U(ii,jj,VH) < 0.0_dp ) ) domain%U(ii,jj,VH) = 0.0_dp

                    end do
                end do
            end if
        end do ! End of loop over nesting receive comms
        !$OMP END DO
        !$OMP END PARALLEL

TIMER_STOP('receive_halos')

    end subroutine

    subroutine send_halos(domain, p2p, send_to_recv_buffer, time)
        !!
        !! Loop over domain%nesting%send_comms, and send the halo data
        !!
        class(domain_type), intent(inout) :: domain
        type(p2p_comms_type), intent(inout) :: p2p
        logical, intent(in) :: send_to_recv_buffer
            !! If .TRUE., then do parallel communication. If .FALSE., then only copy data
            !! to the send buffer (in that case, one can later call "communicate_p2p"
            !! to do all communications at once, which can be faster as it enables
            !! multiple sends between the same images to be collapsed into a single send).
        real(dp), intent(in), optional :: time

        integer(ip) :: i

TIMER_START('send_halos')
            !! Seems faster to do the openmp inside "process_data_to_send"?
            !! So comment out openmp here. Also, ifort19 seg-faulted when I used openmp here.
            !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, send_to_recv_buffer_local, time)
            !!$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(domain%nesting%send_comms, kind=ip)

                call domain%nesting%send_comms(i)%process_data_to_send(domain%U, time=time)
                call domain%nesting%send_comms(i)%send_data(p2p, send_to_recv_buffer=send_to_recv_buffer)
            end do
            !!$OMP END DO
            !!$OMP END PARALLEL
TIMER_STOP('send_halos')

    end subroutine

    subroutine smooth_elevation(domain, smooth_method)
        !!
        !! Basic smooth of elevation
        !!
        class(domain_type), intent(inout) :: domain
        character(*), optional, intent(in) :: smooth_method
            !! Control the smoothing method. Values are '9pt_average', or 'cliffs'

        integer(ip) :: i, j, nx, ny
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

            nx = domain%nx(1)
            ny = domain%nx(2)
            ! Smooth bathymetry along 'i' axis
            !do j = 2, domain%nx(2) - 1
            domain%U(2:(nx-1),1:ny,ELV) = (1.0_dp/3.0_dp)*(&
                domain%U(2:(nx-1),1:ny,ELV) + &
                domain%U(1:(nx-2),1:ny,ELV) + &
                domain%U(3:(nx-0),1:ny,ELV) )
            !end do
            ! Smooth bathymetry along 'j' axis
            !do i = 2, domain%nx(1) - 1
            domain%U(1:nx, 2:(ny-1),ELV) = (1.0_dp/3.0_dp)*(&
                domain%U(1:nx, 2:(ny-1),ELV) + &
                domain%U(1:nx, 1:(ny-2),ELV) + &
                domain%U(1:nx, 3:(ny-0),ELV) )
            !end do

        case('cliffs')
            ! Cliffs bathymetry smoother. Tolkova has shown (e.g. user manual) that this is important for stability.
            ! This is written with a different elevation-sign-convention than what I use, so reverse the sign
            domain%U(:,:,ELV) = -1*domain%U(:,:,ELV)
            call setSSLim(domain%nx(1), domain%nx(2), domain%U(:,:,ELV), &
                domain%cliffs_bathymetry_smoothing_alpha, domain%cliffs_minimum_allowed_depth)
            ! Move back to our elevation sign convention
            domain%U(:,:,ELV) = -1*domain%U(:,:,ELV)

        case default
            stop 'Smoothing method not recognized'
        end select

    end subroutine

    subroutine smooth_elevation_near_point(domain, number_of_9pt_smooths, pt, &
            smooth_region_radius_meters, transition_region_radius_meters)
        !!
        !! Smooth the elevation grid near a specified point.
        !!
        !! We apply the 9pt_smoother to the elevation a specified number of times. At sites within a distance
        !! "smooth_region_radius_meters" of the point, we use the smoothed elevation value. For points with distance between
        !! smooth_region_radius_meters and transition_region_radius_meters, we linearly blend the smooth DEM value and the original
        !! DEM value by a factor that decreases from 1 to 0 with distance.
        !!
        class(domain_type), intent(inout) :: domain
            !! Input domain
        integer(ip), intent(in) :: number_of_9pt_smooths
            !! Number of calls to domain%smooth_elevation().
        real(dp), intent(in) :: pt(2)
            !! The point in the middle of the region to be smoothed
        real(dp), intent(in) :: smooth_region_radius_meters
            !! Use the smooth DEM within this distance from pt
        real(dp), intent(in) :: transition_region_radius_meters
            !! Transition between the smooth DEM and the original DEM at distances between smooth_region_radius_meters and
            !! transition_region_radius_meters, so that beyond transition_region_radius_meters, we are just using the original DEM.

        real(dp), allocatable :: initial_elevation(:,:)
        integer(ip) :: i, j, nx, ny
        real(dp) :: distance_to_pt, wt, dlon, dlat, cosfac
        logical :: is_nearby

        if(number_of_9pt_smooths < 1) return

        if(smooth_region_radius_meters > transition_region_radius_meters) then
            write(log_output_unit, *) 'smooth_region_radius_meters must be <= transition_region_radius_meters'
            call generic_stop
        end if

        if(smooth_region_radius_meters < 0.0_dp) then
            write(log_output_unit, *) 'smooth_region_radius_meters must be non-negative'
            call generic_stop
        end if

        nx = domain%nx(1)
        ny = domain%nx(2)

#ifdef SPHERICAL
        ! Check whether pt is 'likely' close enough to the domain for smoothing to matter,
        ! by searching in an 'extended' box.

        if(abs(pt(1) - domain%x(1)) >= 360.0_dp .or. abs(pt(1) - domain%x(nx)) >= 360.0_dp) then
            write(log_output_unit, *) 'Code does not yet treat longitude wrapping in spherical coordinates'
            call generic_stop
        end if

        ! Approximate lat change corresponding to the transition_region_radius_meters
        dlat = transition_region_radius_meters/radius_earth * 1.0_dp/DEG2RAD
        ! For dlon, this depends on latitude
        cosfac = min(cos((domain%y(1)-dlat)*DEG2RAD), cos((domain%y(ny) + dlat)*DEG2RAD)) ! Conservative
        if(cosfac <= 0.0_dp) then
            dlon = 360.0_dp
        else
            dlon = min( dlat/cosfac, 360.0_dp)
        end if

        is_nearby = ( (pt(1) >= domain%x(1 ) - dlon) .and. &
                      (pt(1) <= domain%x(nx) + dlon) .and. &
                      (pt(2) >= domain%y(1 ) - dlat) .and. &
                      (pt(2) <= domain%y(ny) + dlat) )

#else
        ! Check whether pt is 'likely' close enough to the domain for smoothing to matter,
        ! by searching in an 'extended' box
        is_nearby = ( (pt(1) >= domain%x(1 ) - transition_region_radius_meters) .and. &
                      (pt(1) <= domain%x(nx) + transition_region_radius_meters) .and. &
                      (pt(2) >= domain%y(1 ) - transition_region_radius_meters) .and. &
                      (pt(2) <= domain%y(ny) + transition_region_radius_meters) )
#endif

        if(.not. is_nearby) return

        !
        ! Apply smoothing
        !

        ! Backup
        initial_elevation = domain%U(:,:, ELV)

        ! Smooth the entire elevation grid. This could be made more efficient if only a small patch of domain
        ! needs to be smoothed.
        do i = 1, number_of_9pt_smooths
            call domain%smooth_elevation(smooth_method='9pt_average')
        end do

        !$OMP PARALLEL DO PRIVATE(distance_to_pt, wt) &
        !$OMP SHARED(initial_elevation, nx, ny, domain, pt, transition_region_radius_meters, smooth_region_radius_meters)
        do j = 1, ny
            do i = 1, nx
                ! Weight the smoothed version, based on the distance to pt
                distance_to_pt = pt_distance(domain%x(i), domain%y(j))
                wt = (transition_region_radius_meters - distance_to_pt) / &
                     (transition_region_radius_meters - smooth_region_radius_meters)
                wt = min(1.0_dp, max(wt, 0.0_dp))
                domain%U(i,j,ELV) = wt * domain%U(i,j,ELV) + (1.0_dp - wt) * initial_elevation(i,j)
            end do
        end do
        !$OMP END PARALLEL DO

        deallocate(initial_elevation)

        contains

            function pt_distance(x, y) result(distance_to_pt)
            real(dp), intent(in) :: x, y
            real(dp) :: distance_to_pt
#ifdef SPHERICAL
                distance_to_pt = distance_haversine(x, y, pt(1), pt(2))
#else
                distance_to_pt = sqrt((x - pt(1))**2 + (y - pt(2))**2)
#endif
            end function

    end subroutine

    subroutine smooth_elevation_near_nesting_fine2coarse_boundaries(domain, all_dx_md, &
            coarse_side_ncells, fine_side_ncells, number_of_9pt_smooths)
        !! Smooth the domain elevation locally near fine2coarse nesting boundaries (for domains that
        !! are part of a multidomain).
        !!
        !! Elevation smoothing at specific sites can improve stability in some cases. In practice
        !! problems most often arise at sites near nesting boundaries (where the grid size changes).
        !! This routine applies "local" smoothing near all coarse to fine boundaries in the domain.
        !! The number of cells involved around the nesting boundary differs for the fine domain vs
        !! the coarse domain.
        !!
        !! See unit-test in the multidomain module. 
        !! See domain%smooth_elevation_near_point for smoothing near a specified point.
        !! See domain%smooth_elevation for smoothing the entire domain.
        !!
        class(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: all_dx_md(:,:,:)
            !! Array with rank (2, max_number_domains_on_any_image, number_images), which in practice
            !! is stored in the multidomain object md%all_dx_md.
            !! It contains the cell size (dx & dy) for each domain on every image in a multidomain.
            !! The first argument "domain" is also part of the multidomain.
            !! At any cell i,j, in our domain, the priority_domain_cellsize(1:2) can be computed as:
            !!     n = domain%nesting%priority_domain_index(i,j)
            !!     m = domain%nesting%priority_domain_image(i,j)
            !!     priority_domain_cellsize = all_dx_md(1:2, n, m)
        integer(ip), optional, intent(in) :: coarse_side_ncells
            !! Smooth within "coarse_side_ncells" of a coarse-2-fine nesting boundary on the coarse domain.
        integer(ip), optional, intent(in) :: fine_side_ncells
            !! Smooth within "fine_side_ncells" of a coarse-2-fine nesting boundary on the fine domain.
        integer(ip), optional, intent(in) :: number_of_9pt_smooths
            !! Smoothing corresponds to this many applications of a 9pt smooth

        real(dp), allocatable :: temp_elev(:,:)
        integer(ip) :: i, j, n, m, i0, j0
        real(dp) :: smooth_elev, nonsmooth_elev
        integer(ip) :: coarse_window_size, fine_window_size, nsmooth, smooth_footprint
        logical:: smooth_footprint_exceeded_boundary

        ! Smooth within this many cells on the coarse side of a nesting boundary
        coarse_window_size = 2_ip
        if(present(coarse_side_ncells)) coarse_window_size = coarse_side_ncells

        ! Smooth within this many cells on the fine side of a nesting boundary
        fine_window_size = max(nint(domain%max_parent_dx_ratio), 2_ip)
        if(present(fine_side_ncells)) fine_window_size = fine_side_ncells

        ! Smooth using this many 9pt smooths
        nsmooth = 1_ip
        if(present(number_of_9pt_smooths)) nsmooth = number_of_9pt_smooths

        ! At cell (i,j), smoothing leads to an influence from cells
        !  [i-smooth_footprint, ..., i+smooth_footprint] .outer_product. &
        !  [j-smooth_footprint, ..., j+smooth_footprint]
        ! We can warn if any of these cells touched the domain boundary (as
        ! changes to the domain partition might interact with the smoothing, potentially
        ! making models more sensitive to the partition than they should be).
        smooth_footprint = 3 + (nsmooth-1)

        ! Store the non-smooth elevation
        allocate(temp_elev(domain%nx(1), domain%nx(2)))
        temp_elev = domain%U(:,:,ELV)

        ! Smooth elevation everywhere.
        ! Later we ignore most of these values -- could be made more efficient.
        do j = 1, nsmooth
            call domain%smooth_elevation('9pt_average')
        end do

        smooth_footprint_exceeded_boundary = .false. 

        !$OMP PARALLEL DEFAULT(PRIVATE),&
        !$OMP SHARED(domain, temp_elev, all_dx_md, coarse_window_size, fine_window_size, nsmooth, smooth_footprint), &
        !$OMP REDUCTION(.or.: smooth_footprint_exceeded_boundary)

        ! Swap temp_elev and domain%U(:,:,ELV) -- so the domain elevation is
        ! non-smooth, and the temp_elev is smooth
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                smooth_elev = domain%U(i,j,ELV)
                nonsmooth_elev = temp_elev(i,j)

                domain%U(i,j,ELV) = nonsmooth_elev
                temp_elev(i,j) = smooth_elev
            end do
        end do
        !$OMP END DO

        ! These may be redefined in loop below
        smooth_footprint_exceeded_boundary = .false.

        ! Find cells that need smoothing, and set their elevations to the smooth value
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            cell_loop: do i = 1, domain%nx(1)

                ! Skip non-priority-domain cells
                if(domain%nesting%is_priority_domain_not_periodic(i,j) == 0_ip) cycle cell_loop

                ! Detect cells on the coarse side of a coarse-2-fine boundary
                do j0 = max(j-coarse_window_size, 1), min(j+coarse_window_size, domain%nx(2))
                    do i0 = max(i-coarse_window_size, 1), min(i+coarse_window_size, domain%nx(1))
                        n = domain%nesting%priority_domain_index(i0, j0)
                        m = domain%nesting%priority_domain_image(i0, j0)
                        ! Finer cells will be an integer divisor finer (1/2, 1/3, ...) than cells in this domain.
                        ! Below we use a factor slightly larger than 1/2 to protect against
                        ! small floating point differences in dx.
                        if(any(all_dx_md(1:2, n, m) < 0.55_dp * domain%dx(1:2))) then
                            domain%U(i,j,ELV) = temp_elev(i,j)
                            if(i - smooth_footprint < 1 .or. i + smooth_footprint > domain%nx(1) .or. &
                               j - smooth_footprint < 1 .or. j + smooth_footprint > domain%nx(2)) then
                               ! Trigger a warning about boundary interactions
                                smooth_footprint_exceeded_boundary = .true.
                            end if
                            cycle cell_loop
                        end if
                    end do
                end do

                ! Detect cells on the fine side of a coarse-2-fine boundary
                do j0 = max(j-fine_window_size, 1), min(j+fine_window_size, domain%nx(2))
                    do i0 = max(i-fine_window_size, 1), min(i+fine_window_size, domain%nx(1))
                        n = domain%nesting%priority_domain_index(i0, j0)
                        m = domain%nesting%priority_domain_image(i0, j0)
                        ! Coarser cells will be an integer factor coarser (2, 3, ...) than cells in this domain.
                        ! Below we use a factor slightly smaller than 2 to protect against
                        ! small floating point differences in dx.
                        if(any(all_dx_md(1:2, n, m) > 1.9_dp * domain%dx(1:2))) then
                            domain%U(i,j,ELV) = temp_elev(i,j)
                            if(i - smooth_footprint < 1 .or. i + smooth_footprint > domain%nx(1) .or. &
                               j - smooth_footprint < 1 .or. j + smooth_footprint > domain%nx(2)) then
                               ! Trigger a warning about boundary interactions
                                smooth_footprint_exceeded_boundary = .true.
                            end if
                            cycle cell_loop
                        end if
                    end do
                end do

            end do cell_loop
        end do
        !$OMP END DO

        !$OMP END PARALLEL

        deallocate(temp_elev) ! Should happen automatically but have seen buggy compilers (years ago).

        ! Check whether the smoothing would have been contained in the nesting layer. Warn if not. 
        if( smooth_footprint_exceeded_boundary ) then  
            write(log_output_unit, *) 'WARNING: ', &
            'Elevation smoothing should have been affected by regions outside the domain. ', &
            'Thus the smooth elevation may differ using a different domain partition. ', &
            'This may not matter - but could be significant if comparing models with different domain partitions.'
        end if


    end subroutine


    subroutine store_forcing(domain)
        !! Adds new forcing terms to domain%forcing_terms_storage.
        !!
        !! Recall that a forcing term can be used by setting domain%forcing_subroutine => my_subroutine,
        !! and optionally populating domain%forcing_context_cptr with a c_ptr to required data. If only
        !! one forcing term is applied then nothing else needs to be done. But what about if we want to
        !! add multiple forcing terms? In that case we add the first term as above. Then
        !! "call domain%store_forcing()" will append those terms to the array
        !! "domain%forcing_terms_storage(:)" and also clear the forcing data (i.e. it will set
        !! domain%forcing_context_cptr=c_null_pointer and domain%foring_subroutine => NULL() ).
        !! We can then optionally add another forcing term by again setting domain%forcing_subroutine and
        !! domain%forcing_context_cptr and calling domain%store_forcing(). And repeat for multiple forcing terms.
        !!
        class(domain_type), intent(inout) :: domain

        type(forcing_term_type) :: temp_forcing
        type(forcing_term_type), allocatable :: swap_array(:)

        if(.not. associated(domain%forcing_subroutine)) return

        ! Store the forcing
        temp_forcing%forcing_subroutine => domain%forcing_subroutine
        temp_forcing%forcing_context_cptr = domain%forcing_context_cptr

        ! Append to the forcing terms storage
        if(allocated(domain%forcing_terms_storage)) then
            ! Re-allocate so we can extend the array-size by 1 element
            allocate(swap_array(size(domain%forcing_terms_storage) + 1))
            swap_array(1:size(domain%forcing_terms_storage)) = domain%forcing_terms_storage
            swap_array(size(domain%forcing_terms_storage) + 1) = temp_forcing
            call move_alloc(swap_array, domain%forcing_terms_storage)
        else
            allocate(domain%forcing_terms_storage(1))
            domain%forcing_terms_storage(1) = temp_forcing
        end if

        ! Remove the forcing from the regular domain terms
        domain%forcing_subroutine => NULL()
        domain%forcing_context_cptr = c_null_ptr

    end subroutine

    subroutine copy_U_to_dispersive_last_timestep(domain)
        !! Copy domain%U to the variable domain%ds%last_U, and record timing info.
        !!
        !! See documentation blurb in domain%solve_dispersive_terms for more info. 
        !!
        class(domain_type), intent(inout) :: domain

        if(domain%use_dispersion) then
TIMER_START('dispersive_store')
            call domain%ds%store_last_U(domain%U)
TIMER_STOP('dispersive_store')
        end if
    end subroutine

    subroutine solve_dispersive_terms(domain, rhs_is_up_to_date, estimate_solution_forward_in_time, forward_time)
        !! Solve the dispersive terms.
        !!
        !! In practice we evolve models with dispersive terms by:
        !!    1. Store domain%U in domain%ds%last_U (via domain%copy_U_to_dispersive_last_timestep). Say this time is tLast
        !!    2. Evolve one or more shallow water steps (finishing at time tNew)
        !!    3. Solve the dispersive terms with domain%solve_dispersive_terms (which this routine does). 
        !!       This updates the solution with dispersive terms, with the time discretization centred from 
        !!       tLast to tNew (assuming U at tLast is stored in domain%ds%last_U).
        !!
        class(domain_type), intent(inout) :: domain
        logical, intent(in) :: rhs_is_up_to_date
            !! Are the RHS terms in the dispersive solve already up to date? This is useful in the multidomain
            !! context when using iteration, if the RHS does not change between iterations.
        logical, intent(in) :: estimate_solution_forward_in_time
            !! Should we use extrapolation in time to guess the solution before iteration
        real(dp), intent(in) :: forward_time
            !! Time for the guessed solution (only used if estimate_solution_forward_in_time is TRUE)

        if(domain%use_dispersion .and. (domain%time >= domain%static_before_time)) then
TIMER_START('dispersive_solve')
            call domain%ds%solve_tridiag(domain%is_staggered_grid, &
                domain%U, domain%dx(1), domain%dx(2), &
                domain%distance_bottom_edge, domain%distance_left_edge, &
                domain%msl_linear, &
                rhs_is_up_to_date = rhs_is_up_to_date, &
                estimate_solution_forward_in_time = estimate_solution_forward_in_time, &
                forward_time = forward_time)
TIMER_STOP('dispersive_solve')
        end if
    end subroutine

    subroutine test_domain_mod
        ! Unit-tests (basic stuff only -- the solver quality is checked by the validation tests)
        type(domain_type) :: domain
        integer(ip), parameter :: nx = 20, ny = 10
        real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp]
        real(dp), parameter :: global_lw(2) = [100.0_dp, 50.0_dp]
        real(dp), parameter :: pt(2) = [50.0_dp, 20.0_dp]
        real(dp), parameter :: r1 = 10.0_dp, r2 = 15.0_dp

        integer(ip) :: i, j
        real(dp) :: err, tol, dist2, wt
        logical :: has_failed

        call domain%allocate_quantities(global_lw, [nx, ny], global_ll, create_output_files = .FALSE., verbose=.FALSE.)

        !
        ! Check the domain was allocated with zero U
        !
        if(all(domain%U == 0.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        !
        ! Check smooth_elevation
        !
        do j = 1, ny
            do i = 1, nx
                ! Ensure that in every 9x9 block there is a single elevation value of 1.0 (the rest are 0.0)
                if(mod(i, 3) == 2 .and. mod(j, 3) == 2) domain%U(i,j,ELV) = 1.0_dp
            end do
        end do
        !do j = 1, ny
        !    print*, j, ':', domain%U(:,j,ELV)
        !end do

        ! With this above configuration, 9pt average should lead to the interior having an elevation of 1.0/9.0
        call domain%smooth_elevation(smooth_method='9pt_average')
        err = maxval(abs(domain%U(2:nx-1, 2:ny-1, ELV) - 1.0_dp/9.0_dp))
        tol = 100.0_dp * spacing(1.0_dp)
        if( err < tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, tol
        end if
        !do j = 1, ny
        !    print*, j, ':', domain%U(:,j,ELV)
        !end do


        !
        ! Check smooth_elevation_near_pt
        !
        domain%U(:,:,ELV) = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                ! Ensure that in every 9x9 block there is a single elevation value of 1.0 (the rest are 0.0)
                if(mod(i, 3) == 2 .and. mod(j, 3) == 2) domain%U(i,j,ELV) = 1.0_dp
            end do
        end do
        !do j = 1, ny
        !    print*, j, ':', domain%U(:,j,ELV)
        !end do

        call smooth_elevation_near_point(domain, number_of_9pt_smooths = 1_ip, pt = pt, &
            smooth_region_radius_meters = r1, transition_region_radius_meters = r2)

        has_failed = .false.
        do j = 2, ny-1
            do i = 2, nx-1

                dist2 = (domain%x(i) - pt(1))**2 + (domain%y(j) - pt(2))**2

                if( ( dist2 <= r1**2 ) ) then
                    ! Inner radius
                    if(abs(domain%U(i,j,ELV) - 1.0_dp/9.0_dp) > tol) has_failed = .true.
                else if(dist2 >= r2**2) then
                    ! Outside smoothing zone
                    if(.not. any(domain%U(i,j,ELV) == [0.0_dp, 1.0_dp])) has_failed = .true.
                else
                    ! Transition zone
                    wt = (r2 - sqrt(dist2)) / (r2 - r1)
                    ! Value should be a weighted average of 1.0/9.0, and either 0, or 1.0
                    if(.not. any(abs( domain%U(i,j,ELV) - (wt * 1.0_dp/9.0_dp + (1.0_dp - wt) * [0.0_dp, 1.0_dp])) < tol)) then
                        has_failed = .true.
                        !print*, domain%U(i,j,ELV) - (wt * 1.0_dp/9.0_dp + (1.0_dp - wt) * [0.0_dp, 1.0_dp]), tol
                    end if
                end if
            end do
        end do

        if(.not. has_failed) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if
        !do j = 1, ny
        !    print*, j, ':', domain%U(:,j,ELV)
        !end do

    end subroutine

end module domain_mod
