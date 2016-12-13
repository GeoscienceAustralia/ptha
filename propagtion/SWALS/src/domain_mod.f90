! Compile with -DTIMER to add timing to the code 
#ifdef TIMER
#   define TIMER_START(tname) call domain%timer%timer_start(tname)
#   define TIMER_STOP(tname)  call domain%timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

MODULE domain_mod

    USE global_mod, only: dp, ip, charlen, output_precision, &
                          cfl, maximum_timestep, gravity, &
                          advection_beta, &
                          minimum_allowed_depth, &
                          minimum_allowed_depth, &
                          default_timestepping_method, extrapolation_theta,&
                          wall_elevation, &
                          default_output_folder
    USE timer_mod, only: timer_type
    USE point_gauge_mod, only: point_gauge_type
    !USE coarray_utilities_mod, only: partitioned_domain_NESW_comms_type
    USE stop_mod, only: generic_stop
    USE iso_fortran_env, only: output_unit

! Compile with -DSPHERICAL to get the code to run in spherical coordinates
#ifdef SPHERICAL    
    USE spherical_mod, only: area_lonlat_rectangle, deg2rad
    USE global_mod, only: radius_earth
#endif

    IMPLICIT NONE

    ! Make everything private, except domain_type, which has its own methods
    PRIVATE
    PUBLIC:: domain_type

    ! Indices for arrays: Stage, depth-integrated-x-velocity,
    ! depth-integrated-v-velocity, elevation. So e.g. stage
    ! is in domain%U(:,:,STG), and elevation is in domain%U(:,:ELV)
    INTEGER(4), PARAMETER, PUBLIC:: STG=1, UH=2, VH=3, ELV=4

    REAL(dp), PARAMETER :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp
    REAL(dp), PARAMETER:: NEG_SEVEN_ON_THREE_dp = -7.0_dp/3.0_dp

    !
    ! Main type used in the program. It holds the model arrays, information on the domain,
    ! etc
    !
    TYPE:: domain_type

        ! [Length,width] of domain in x/y units
        REAL(dp):: lw(2) 
        ! grid size [number of x-cells, number of y-cells]
        INTEGER(ip):: nx(2) 
        ! cell size [dx,dy]
        REAL(dp):: dx(2)
        ! Absolute lower-left coordinate (bottom left corner of cell(1,1))
        REAL(dp):: lower_left(2)

        ! Number of quantities (stage, uh, vh, elevation)
        INTEGER(ip):: nvar = 4 !global_nvar

        ! Domain ID, which is useful if multiple domains are running
        INTEGER(4):: myid = 1

        ! Name of ascii file where we output metadata
        CHARACTER(charlen):: metadata_ascii_filename

        ! The domain 'interior' is surrounded by 'exterior' cells which are
        ! updated by boundary conditions, or copied from other domains. When
        ! tracking mass conservation, we only want to record inflows/outflows to
        ! interior cells. The interior cells are:
        !    [(1+exterior_cells_width):(domain%nx(1)-exterior_cells_width), &
        !     (1+exterior_cells_width):(domain%nx(2)-exterior_cells_width)]
        ! NOTE: exterior_cells_width should only be used to influence mass conservation tracking calculations.
        ! It is not strictly related to the halo width for parallel computations, although we may well want 
        ! the that to be the same.
        INTEGER(ip):: exterior_cells_width = 2
    
        ! Count number of time steps (useful if we want to do something only
        ! every n'th timestep)
        INTEGER(ip):: nsteps_advanced = 0
        
        REAL(dp):: time = ZERO_dp
        ! dt computed from CFL condition (e.g. during flux computation)
        REAL(dp):: max_dt = ZERO_dp
        ! the time-step evolved by domain%evolve_one_step.
        REAL(dp):: evolve_step_dt = ZERO_dp
        ! the maximum allowed timestep [used to prevent high timesteps on dry domains]
        REAL(dp):: maximum_timestep = maximum_timestep
        ! CFL number
        REAL(dp):: cfl = cfl
        ! timestepping_method determines the choice of solver
        CHARACTER(len=charlen):: timestepping_method = default_timestepping_method
        ! Parameter controlling extrapolation for finite volume methods
        REAL(dp):: theta = extrapolation_theta
        ! Parameter which determines how the 'static depth' is computed in
        ! linear solver [corresponding to 'h0' in 'g h_0 d(free_surface)/dx']
        REAL(dp):: MSL_linear = 0.0_dp

        !
        ! Boundary conditions. (FIXME: Consider making a 'boundary_data_type' which
        !    could hold the variables prefixed with 'boundary_'. Not essential but
        !    might reduce proliferation of variables in domain.)
        !
        ! Currently 'boundary_subroutine' is the most important way of imposing
        !    boundary conditions. It must take a domain_type as an INTENT(INOUT) argument
        ! Applications need to define boundary_subroutine, either using a
        !    pre-existing bc, or making a new subroutine.
        ! Some existing boundary conditions rely on a function as well, hence
        !    'boundary_function', but its use should probably be discouraged. It takes
        !    the domain_type as well as time, x, y, see the interface below
        CHARACTER(len=charlen):: boundary_type = ''
        PROCEDURE(boundary_fun), pointer, nopass:: boundary_function => NULL()
        PROCEDURE(boundary_subroutine), pointer, nopass:: boundary_subroutine => default_boundary_subroutine
        !
        ! Flag whether the boundary is exterior (TRUE) or interior (FALSE). N, E, S, W order
        ! The interior boundaries are 'overlaps between domains', and are dealt with by communication. 
        ! The numerical boundary conditions (e.g. reflective, etc) are applied to the exterior boundaries
        LOGICAL :: boundary_exterior(4) = [.TRUE., .TRUE., .TRUE., .TRUE.]

        ! 
        ! Mass conservation tracking 
        !
        ! Store the flux through the N, E, S, W boundaries 
        REAL(dp):: boundary_flux_store(4)  = ZERO_dp
        ! Time integrate the boundary fluxes
        REAL(dp):: boundary_flux_time_integral = ZERO_dp
        ! We need an intermediate variable to take care of time-stepping
        ! This integrates the boundary fluxes within the evolve step only
        REAL(dp):: boundary_flux_evolve_integral = ZERO_dp
        
        ! Lower/upper x and y indices over which the SWE computation takes place
        ! These might restrict the SWE update to a fraction of the domain (e.g.
        ! beginning of an earthquake-tsunami run where only a fraction of the domain is
        ! active.
        ! Currently implemented for linear SWE only
        INTEGER(ip):: xL, xU, yL, yU

        ! Output folder units
        ! (FIXME: Consider making a 'output_data_type' which
        !    could hold the variables prefixed with 'output_'. Not essential but
        !    might reduce proliferation of variables in domain.)

        INTEGER(ip), ALLOCATABLE :: output_variable_unit_number(:)
        INTEGER(ip):: output_time_unit_number
        CHARACTER(len=charlen):: output_basedir = default_output_folder
        CHARACTER(len=charlen):: output_folder_name = ''
        INTEGER(ip):: logfile_unit = output_unit

        ! Type to record CPU timings
        TYPE(timer_type):: timer

        ! Type to do single-grid coarray communication
        !TYPE(partitioned_domain_NESW_comms_type):: comms
        
        ! We don't have to store the max_U. For some problems (linear solver) that can
        ! take a significant fraction of the total time, or use too much memory.
        LOGICAL:: record_max_U = .true.
        INTEGER(ip):: max_U_update_frequency = 1
        
        ! Spatial coordinates, dx/dy distances (useful for spherical coordinates)
        REAL(dp), ALLOCATABLE :: x(:), y(:), distance_bottom_edge(:), distance_left_edge(:)
        REAL(dp), ALLOCATABLE :: area_cell_y(:)
#ifdef SPHERICAL
        REAL(dp), ALLOCATABLE :: coslat(:), coslat_bottom_edge(:)
#endif

        ! Tide gauges
        TYPE(point_gauge_type) :: point_gauges
 
        ! Big arrays
        !
        ! U holds quantities
        ! First 2 dimensions = space, 3rd = number of quantities
        REAL(dp), ALLOCATABLE :: U(:,:,:) ! Needed always
        REAL(dp), ALLOCATABLE :: max_U(:,:,:) ! Needed if max_U is to be output, but can reduce precision
        !
        ! Multi-dimensional arrays below are only required for nonlinear. Would be possible
        ! to further reduce memory usage, but not completely trivial
        REAL(dp), ALLOCATABLE :: flux_NS(:,:,:) ! Could avoid
        REAL(dp), ALLOCATABLE :: flux_EW(:,:,:) ! Could avoid
        REAL(dp), ALLOCATABLE :: depth(:,:) ! Could avoid
        REAL(dp), ALLOCATABLE :: explicit_source(:,:,:) ! Could avoid
        REAL(dp), ALLOCATABLE :: explicit_source_VH_j_minus_1(:,:) ! Separate from explicit_source for OPENMP parallel logic
        REAL(dp), ALLOCATABLE :: velocity(:,:, :) ! Could avoid
        REAL(dp), ALLOCATABLE :: manning_squared(:,:) ! Needed for variable manning
        REAL(dp), ALLOCATABLE :: backup_U(:,:,:) ! Needed

        CONTAINS

        ! Initialisation
        PROCEDURE:: allocate_quantities => allocate_quantities

        ! Reporting
        PROCEDURE:: print => print_domain_statistics

        ! Core routines that occur within a timestep
        PROCEDURE:: get_bottom_edge_values => get_bottom_edge_values
        PROCEDURE:: get_left_edge_values => get_left_edge_values
        PROCEDURE:: compute_fluxes => compute_fluxes
        PROCEDURE:: update_U => update_U 
        PROCEDURE:: backup_quantities => backup_quantities

        ! Timestepping
        PROCEDURE:: one_euler_step => one_euler_step
        PROCEDURE:: one_rk2_step => one_rk2_step 
        PROCEDURE:: one_rk2n_step => one_rk2n_step 
        PROCEDURE:: one_midpoint_step => one_midpoint_step
        PROCEDURE:: one_linear_leapfrog_step => one_linear_leapfrog_step
        PROCEDURE:: evolve_one_step => evolve_one_step

        ! Boundary conditions. This just calls whatever domain%boundary_subroutine points to
        PROCEDURE:: update_boundary => update_boundary

        ! IO
        PROCEDURE:: create_output_files => create_output_files
        PROCEDURE:: write_to_output_files => write_to_output_files
        PROCEDURE:: update_max_quantities => update_max_quantities
        PROCEDURE:: write_max_quantities => write_max_quantities
        PROCEDURE:: log_outputs => divert_logfile_unit_to_file

        ! Mass conservation tracking
        PROCEDURE:: mass_balance_interior => mass_balance_interior
        PROCEDURE:: volume_interior => volume_interior

        ! Time-step for linear solver (a constant)
        PROCEDURE:: linear_timestep_max => linear_timestep_max

        ! Gauges
        PROCEDURE:: setup_point_gauges => setup_point_gauges
        PROCEDURE:: write_gauge_time_series => write_gauge_time_series

        ! Finalization (e.g. close netcdf files)
        PROCEDURE:: finalise => finalise_domain

    END TYPE domain_type

    INTERFACE

        ! This gives the interface for a generic 'boundary function' which
        ! can return 4 output values (stage/uh/vh/elevation)
        !
        ! It is used in conjunction with a number of different types of
        ! boundary conditions below
        !
        FUNCTION boundary_fun(domain, t, x, y) RESULT(stage_uh_vh_elev)
            IMPORT dp, domain_type
            IMPLICIT NONE
            TYPE(domain_type), INTENT(IN):: domain 
            REAL(dp), INTENT(IN):: t, x, y
            REAL(dp):: stage_uh_vh_elev(4)
        END FUNCTION

        !
        ! The user can provide a boundary subroutine which is supposed to update
        ! the domain boundaries. It is called by domain%update_boundary whenever
        ! a boundary update is required by the timestepping_method. Note that
        ! this may well mean 'boundary_fun' is not required
        !
        SUBROUTINE boundary_subroutine(domain)
            IMPORT domain_type
            IMPLICIT NONE
            TYPE(domain_type), INTENT(INOUT):: domain
        END SUBROUTINE
    
    END INTERFACE

    contains
  
    ! 
    ! Convenience printing function 
    !
    SUBROUTINE print_domain_statistics(domain)
        CLASS(domain_type), INTENT(IN):: domain
        REAL(dp):: maxstage, minstage
        INTEGER:: i,j, ecw
        REAL(dp):: dry_depth_threshold, energy_total, energy_potential, energy_kinetic
        REAL(dp):: depth, depth_iplus, depth_jplus
        LOGICAL, PARAMETER:: report_energy_statistics=.FALSE.

        dry_depth_threshold = minimum_allowed_depth
            

        maxstage = -HUGE(1.0_dp)
        minstage = HUGE(1.0_dp)
        DO j = 1, domain%nx(2)
            DO i = 1, domain%nx(1)
                if(domain%U(i,j,STG) > domain%U(i,j,ELV) + dry_depth_threshold) then
                    maxstage = max(maxstage, domain%U(i,j,STG))
                    minstage = min(minstage, domain%U(i,j,STG))
                end if
            END DO
        END DO


        write(domain%logfile_unit, *) ''
        write(domain%logfile_unit, *) 'Domain ID: '
        write(domain%logfile_unit, *) '        ', domain%myid
        write(domain%logfile_unit, *) 'Time: '
        write(domain%logfile_unit, *) '        ', domain%time
        write(domain%logfile_unit, *) 'nsteps_advanced:'
        write(domain%logfile_unit, *) '        ', domain%nsteps_advanced
        write(domain%logfile_unit, *) 'dt: '
        write(domain%logfile_unit, *) '        ', domain%max_dt
        write(domain%logfile_unit, *) 'evolve_step_dt: '
        write(domain%logfile_unit, *) '        ', domain%evolve_step_dt
        write(domain%logfile_unit, *) 'Stage: '
        write(domain%logfile_unit, *) '        ', maxstage
        write(domain%logfile_unit, *) '        ', minstage

        if(domain%timestepping_method /= 'linear') then

            write(domain%logfile_unit, *) 'u: '
            write(domain%logfile_unit, *) '        ', maxval(domain%velocity(:,:,UH))
            write(domain%logfile_unit, *) '        ', minval(domain%velocity(:,:,UH))
            write(domain%logfile_unit, *) 'v: '
            write(domain%logfile_unit, *) '        ', maxval(domain%velocity(:,:,VH))
            write(domain%logfile_unit, *) '        ', minval(domain%velocity(:,:,VH))
            write(domain%logfile_unit, *) '.........'
        end if

        if(report_energy_statistics) then

            ecw = domain%exterior_cells_width
            energy_potential = ZERO_dp
            energy_kinetic = ZERO_dp

            if(domain%timestepping_method == 'linear') then
                !
                ! Linear leap-frog scheme doesn't store velocity
                ! Get total, potential and kinetic energy in domain interior
                !
                ! FIXME: This does not exactly lead to energy conservation,
                ! probably because of the details of the numerics.
                !

                ! Integrate over the model interior
                DO j = 1 + ecw, (domain%nx(2) - ecw)
                    DO i = (1+ecw), (domain%nx(1)-ecw)

                        ! Integrate over wet cells
                        if(domain%MSL_linear > domain%U(i,j,ELV) + dry_depth_threshold) then
            
                            energy_potential = energy_potential + domain%area_cell_y(j) *&
                                (domain%U(i,j,STG) - domain%MSL_linear)**2

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
                    END DO
                END DO
            else
                !
                ! Compute energy statistics using non-staggered grid
                !

                ! Integrate over the model interior
                DO j = 1 + ecw, (domain%nx(2) - ecw)
                    DO i = (1+ecw), (domain%nx(1)-ecw)
                        depth = domain%U(i,j,STG) - domain%U(i,j,ELV)
                        if(depth > dry_depth_threshold) then
                            energy_potential = energy_potential + domain%area_cell_y(j) * &
                                (domain%U(i,j,STG) - domain%MSL_linear)**2
                            energy_kinetic = energy_kinetic + domain%area_cell_y(j) * &
                                (domain%U(i,j,UH)**2 + domain%U(i,j,VH)**2)/depth
                        end if 
                    END DO
                END DO

            end if

            ! Rescale energy statistics appropriately 
            energy_potential = energy_potential * gravity * HALF_dp
            energy_kinetic = energy_kinetic * HALF_dp
            energy_total = energy_potential + energy_kinetic

            write(domain%logfile_unit, *) 'Energy Total / rho: ', energy_total
            write(domain%logfile_unit, *) 'Energy Potential / rho: ', energy_potential
            write(domain%logfile_unit, *) 'Energy Kinetic / rho: ', energy_kinetic

        end if

        write(domain%logfile_unit, *) 'Boundary flux time integral: ', domain%boundary_flux_time_integral

    END SUBROUTINE

    ! 
    ! Set up the full domain, allocate arrays, etc
    !
    ! @param global_lw real(dp) array size 2. length/width of domain in same
    !  units as x,y coordinates
    ! @param global_nx integer array size 2. number of x/y cells in the domain
    ! @param global_ll real(dp) array size 2. lower left x/y coordinate of the
    !  domain (at the corner of the lower left cell)
    ! @param create_output_files optional. If .TRUE. or not provided, then make
    !  output files.
    !
    SUBROUTINE allocate_quantities(domain, global_lw, global_nx, global_ll, create_output_files,&
        co_size_xy)
        CLASS(domain_type), TARGET, INTENT(INOUT):: domain
        REAL(dp), INTENT(IN):: global_lw(2), global_ll(2)
        INTEGER(ip), INTENT(IN):: global_nx(2)
        LOGICAL, OPTIONAL, INTENT(IN) :: create_output_files
        INTEGER(ip), OPTIONAL, INTENT(IN):: co_size_xy(2)

        INTEGER(ip), POINTER:: nx, ny, nvar
        INTEGER(ip) :: i
        LOGICAL :: create_output
        REAL(dp):: local_lw(2), local_ll(2)
        INTEGER(ip):: local_nx(2)

        if(present(create_output_files)) then
            create_output = create_output_files
        else
            create_output = .TRUE.
        end if

        ! Use coarrays if co_size_xy is provided
        if(present(co_size_xy)) then
            !! Compute the ll/lw/nx for this sub-domain
            !CALL domain%comms%initialise(co_size_xy, global_ll, global_lw, global_nx, &
            !    local_ll, local_lw, local_nx)
            !domain%lower_left = local_ll
            !domain%lw = local_lw 
            !domain%nx = local_nx

            !CALL domain%comms%print()
        
            !! Make sure that communication boundaries are not
            !! numerical boundaries
            !DO i = 1,4
            !    if(domain%comms%neighbour_images(i) > 0) domain%boundary_exterior(i) = .FALSE.
            !END DO
            print*, 'ERROR: Coarrays not supported in this version of the code'
            call generic_stop()
        else
            ! Not using coarrays (simple case)
            domain%lw = global_lw
            domain%nx = global_nx
            domain%lower_left = global_ll
        end if

        domain%xL = 1_ip
        domain%yL = 1_ip
        domain%xU = domain%nx(1)
        domain%yU = domain%nx(2)

        ! Compute dx
        domain%dx = domain%lw/(domain%nx*ONE_dp)
        write(domain%logfile_unit, *) ''
        write(domain%logfile_unit, *) 'dx: ', domain%dx
        write(domain%logfile_unit, *) 'nx: ', domain%nx
        write(domain%logfile_unit, *) 'lw: ', domain%lw
        write(domain%logfile_unit, *) 'lower-left: ', domain%lower_left
        write(domain%logfile_unit, *) ''

        nx => domain%nx(1)
        ny => domain%nx(2)
        nvar => domain%nvar 

        ! x/y coordinates (only stored along domain edge, assumed constant)
        ALLOCATE(domain%x(nx), domain%y(ny))
        DO i = 1, nx
            domain%x(i) = domain%lower_left(1) + (i - HALF_dp)/(nx*ONE_dp)*domain%lw(1)
        END DO
        DO i = 1, ny
            domain%y(i) = domain%lower_left(2) + (i - HALF_dp)/(ny*ONE_dp)* domain%lw(2) 
        END DO
     
#ifdef SPHERICAL
        ! For spherical coordinates it saves computation to have cos(latitude)
        ! at cells and edges
        ALLOCATE(domain%coslat(ny), domain%coslat_bottom_edge(ny+1))
        domain%coslat = cos(domain%y * deg2rad)
        domain%coslat_bottom_edge(1:ny) = cos((domain%y - HALF_dp * domain%dx(2))*deg2rad)
        domain%coslat_bottom_edge(ny+1) = cos((domain%y(ny) + HALF_dp*domain%dx(2))*deg2rad)
#endif

        ! Distances along edges 
        ! For spherical coordinates, distance_xedge changes with y 
        ALLOCATE(domain%distance_bottom_edge(ny+1), domain%distance_left_edge(nx+1))

        DO i = 1, nx + 1
#ifdef SPHERICAL
            ! Spherical coordinates
            domain%distance_left_edge(i) = domain%dx(2) * deg2rad * radius_earth
#else
            ! Cartesian 
            domain%distance_left_edge(i) = domain%dx(2)
#endif
        END DO
        DO i = 1, ny+1
#ifdef SPHERICAL
            ! Spherical coordinates
            domain%distance_bottom_edge(i) = domain%dx(1) * deg2rad * radius_earth * &
                domain%coslat_bottom_edge(i)
#else
            ! Cartesian
            domain%distance_bottom_edge(i) = domain%dx(1)
#endif
        END DO
        
        ! Area of cells. For spherical, this only changes with y
        ALLOCATE(domain%area_cell_y(ny))
#ifdef SPHERICAL
        DO i = 1, ny
            domain%area_cell_y(i) = area_lonlat_rectangle(&
                domain%x(1) - domain%dx(1) * HALF_dp, &
                domain%y(i) - domain%dx(2) * HALF_dp, &
                domain%dx(1), &
                domain%dx(2), &
                flat=.true.)
        END DO
#else
        domain%area_cell_y = product(domain%dx)
#endif

        write(domain%logfile_unit, *) ''
        write(domain%logfile_unit, *) 'Total area: ', sum(domain%area_cell_y)
        write(domain%logfile_unit, *) 'distance_bottom_edge(1): ', domain%distance_bottom_edge(1)
        write(domain%logfile_unit, *) 'distange_left_edge(1)', domain%distance_left_edge(1)
        write(domain%logfile_unit, *) ''

        ALLOCATE(domain%U(nx, ny, 1:nvar))
        domain%U = ZERO_dp
       
        if(domain%record_max_U) then 
            ALLOCATE(domain%max_U(nx, ny, 1))
            domain%max_U = -HUGE(1.0_dp)
        end if

        ! Many other variables are required for non-linear FV, but not for
        ! linear leap-frog
        if(domain%timestepping_method /= 'linear') then 

            ALLOCATE(domain%velocity(nx, ny, UH:VH))
            domain%velocity = ZERO_dp
            ALLOCATE(domain%depth(nx, ny))
            domain%depth = ZERO_dp
            
            ! NOTE: We don't need a flux for elevation 
            ALLOCATE(domain%flux_NS(nx, ny+1_ip, 1:3))
            domain%flux_NS = ZERO_dp
            ALLOCATE(domain%flux_EW(nx + 1_ip, ny, 1:3))
            domain%flux_EW = ZERO_dp


            ! Pressure gradient applies only to uh and vh (2 and 3 variables in U)
            ALLOCATE(domain%explicit_source(nx, ny, UH:VH))
            domain%explicit_source = ZERO_dp
            
            ALLOCATE(domain%explicit_source_VH_j_minus_1(nx, ny+1))
            domain%explicit_source_VH_j_minus_1 = ZERO_dp

            ALLOCATE(domain%manning_squared(nx, ny))
            domain%manning_squared = ZERO_dp

            ! If we use complex timestepping we need to back-up U for variables 1-3
            IF (domain%timestepping_method /= 'euler') THEN
                ALLOCATE(domain%backup_U(nx, ny, 1:3))
                domain%backup_U = ZERO_dp
            END IF

        endif

        
        IF(create_output) THEN
            CALL domain%create_output_files()
        END IF

    END SUBROUTINE

    !
    ! minmod function, which is used in some gradient limiters
    !
    ! @param a,b real numbers
    !
    ELEMENTAL FUNCTION minmod(a,b) result(minmod_ab)
        REAL(dp), INTENT(IN):: a, b
        REAL(dp):: minmod_ab, sa, sb
        
        minmod_ab = merge(min(abs(a), abs(b))*SIGN(ONE_dp,a), ZERO_dp, SIGN(ONE_dp,a) == SIGN(ONE_dp,b))

    END FUNCTION

    !
    ! Extrapolation to cell edges on a regular mesh
    !
    ! edge = U_local + limited_gradient * (edge-to-centroid-distance) * extrapolation_sign
    !
    ! @param U_local value of quantity at cell centre (e.g. index i,j)
    ! @param U_lower value of quantity at 'lower' neighbouring cell centre (e.g. index i-1,j)
    ! @param U_upper value of quantity at 'upper' neighbour cell centre (e.g. index i+1, j)
    ! @param theta Parameter controlling the way that 'limited_gradient' is defined
    ! @param extrapolation_sign Integer (+1 or -1) indicating whether to extrapolate forward (+1) or backward(-1)
    !
    ELEMENTAL FUNCTION extrapolate_edge_second_order(U_local, U_lower, U_upper, theta, &
        extrapolation_sign) RESULT(edge_value)
        REAL(dp), INTENT(IN):: U_local, U_lower, U_upper, theta 
        INTEGER(ip), INTENT(IN):: extrapolation_sign
        REAL(dp) :: edge_value
        !LOGICAL, PARAMETER :: min_limit = .FALSE.
        CHARACTER(len=charlen), PARAMETER:: gradient_type = 'standard'
        REAL(dp), PARAMETER:: LIMITER_PAR_dp = 1.9_dp
        REAL(dp):: a, b , c, d, sa, sb

        !! This is much faster
        !FIXME: deliberate bug
        !edge_value = U_local
        !return

        !SELECT CASE (gradient_type)
        
        if(gradient_type == 'min_limit') then
            ! Standard min limiter, with theta [0,1] optionall pushing closer to first order (if < 1)
            edge_value = U_local + extrapolation_sign * theta * HALF_dp * &
                minmod(U_upper - U_local, U_local - U_lower)
        end if

        if(gradient_type == 'standard') then 

            !a = U_upper - U_local
            !b = U_local - U_lower
            d = minmod(U_upper - U_local, U_local - U_lower) * LIMITER_PAR_dp * theta

            c = merge(ZERO_dp, HALF_dp * (U_upper - U_lower), d == ZERO_dp) 

            ! NOTE: IF d /= 0, then clearly d, c have the same sign
            ! We exploit this to avoid a further minmod call (which seems
            ! expensive)

            edge_value = U_local + HALF_dp * extrapolation_sign * &
                merge(min(c, d), max(c, d), d > ZERO_dp)
        end if

        if(gradient_type == 'debug') then
            !edge_value = HALF_dp * merge((U_local + U_upper), (U_local+U_lower), extrapolation_sign > 0_dp)
            edge_value = U_local !+ HALF_dp * extrapolation_sign * HALF_dp * (U_upper - U_lower)
        end if

    END FUNCTION
    
    !
    ! Extrapolation to cell edges on a regular mesh without any limiting
    !
    ! edge = U_local + 0.25*(U_upper - U_lower) * theta * extrapolation_sign
    ! Typically theta = 1.0, but values closer to 0 can be used for first order in space extrapolation
    !
    ELEMENTAL FUNCTION extrapolate_edge_second_order_nolimit(U_local, U_lower, U_upper, theta, &
        extrapolation_sign) RESULT(edge_value)
        REAL(dp), INTENT(IN):: U_local, U_lower, U_upper, theta 
        INTEGER(ip), INTENT(IN):: extrapolation_sign
        REAL(dp) :: edge_value
        
        ! Standard min limiter
        edge_value = U_local + theta * extrapolation_sign * HALF_dp * HALF_dp * (U_upper - U_lower)

    END FUNCTION

    !
    ! Get bottom edge values for the j'th North-South row. Because this
    ! is for a finite volume method, values at edges are discontinuous. The edge has a
    ! 'positive' and 'negative' value, viewed from cell j and j-1 respectively
    ! 
    ! @param domain the domain
    ! @param j get the j'th row
    ! @param nx, ny domain number of cells along x and y directions
    ! @param theta_wd_neg_B 'theta' parameter controlling extrapolation as viewed from cell on negative side of edge
    ! @param theta_wd_pos_B 'theta' parameter controlling extrapolation as viewed from cell on positive side of edge
    ! @param stage_neg_B stage edge value as viewed from cell on negative side of edge
    ! @param stage_pos_B stage edge value as viewed from cell on positive side of edge
    ! @param depth_neg_B depth edge value as viewed from cell on negative side of edge
    ! @param depth_pos_B depth edge value as viewed from cell on positive side of edge
    ! @param u_neg_B x-velocity edge value as viewed from cell on negative side of edge
    ! @param u_pos_B x-velocity edge value as viewed from cell on positive side of edge
    ! @param v_neg_B y-velocity edge value as viewed from cell on negative side of edge
    ! @param v_pos_B y-velocity edge value as viewed from cell on positive side of edge
    ! 
    SUBROUTINE get_bottom_edge_values(domain, j, nx, ny, &
        theta_wd_neg_B, theta_wd_pos_B, &
        stage_neg_B, stage_pos_B, &
        depth_neg_B, depth_pos_B, &
        u_neg_B, u_pos_B, &
        v_neg_B, v_pos_B)

        CLASS(domain_type), INTENT(IN):: domain
        INTEGER(ip), INTENT(IN) :: j, ny, nx
        REAL(dp), INTENT(OUT):: theta_wd_neg_B(nx), theta_wd_pos_B(nx)
        REAL(dp), INTENT(OUT):: depth_neg_B(nx), depth_pos_B(nx)
        REAL(dp), INTENT(OUT):: stage_neg_B(nx), stage_pos_B(nx)
        REAL(dp), INTENT(OUT):: u_neg_B(nx), u_pos_B(nx)
        REAL(dp), INTENT(OUT):: v_neg_B(nx), v_pos_B(nx)

        ! Bottom edge, negative side
        IF(j > 2) THEN
            ! Extrapolate using points j-2, j-1, j
            ! theta = max(min( 5 * ( [min(depth) - eps]/[max(depth) + 1000 eps] - 0.1), 1), 0)
            !     
            theta_wd_neg_B = 5.0_dp*( &
                (min(domain%depth(:,j-1), domain%depth(:,j-2), domain%depth(:,j)) &
                    - minimum_allowed_depth) / &
                (max(domain%depth(:,j-1), domain%depth(:,j-2), domain%depth(:,j)) &
                    + 1000.0_dp*minimum_allowed_depth) &
                - 0.1_dp)
            theta_wd_neg_B = max(min(domain%theta, theta_wd_neg_B), ZERO_dp)

            ! FIXME: test only
            !theta_wd_neg_B = 1.0_dp

            stage_neg_B = extrapolate_edge_second_order(domain%U(:, j-1,STG), &
            !stage_neg_B = extrapolate_edge_second_order_nolimit(domain%U(:, j-1,STG), &
                domain%U(:, j-2,STG), domain%U(:, j, STG), theta_wd_neg_B, 1_ip)
            depth_neg_B = extrapolate_edge_second_order(domain%depth(:, j-1), &
            !depth_neg_B = extrapolate_edge_second_order_nolimit(domain%depth(:, j-1), &
                domain%depth(:, j-2), domain%depth(:, j), theta_wd_neg_B, 1_ip)
            u_neg_B = extrapolate_edge_second_order(domain%velocity(:, j-1, UH), &
            !u_neg_B = extrapolate_edge_second_order_nolimit(domain%velocity(:, j-1, UH), &
                domain%velocity(:, j-2, UH), domain%velocity(:, j, UH), theta_wd_neg_B, 1_ip)
            v_neg_B = extrapolate_edge_second_order(domain%velocity(:, j-1, VH), &
            !v_neg_B = extrapolate_edge_second_order_nolimit(domain%velocity(:, j-1, VH), &
                domain%velocity(:, j-2, VH), domain%velocity(:, j, VH), theta_wd_neg_B, 1_ip)

        ELSE
            ! j == 2, cannot extrapolate so just use the lower neighbour value (first order accurate)
            theta_wd_neg_B = ZERO_dp
            stage_neg_B = domain%U(:, j-1, STG)
            depth_neg_B = domain%depth(:, j-1)
            u_neg_B = domain%velocity(:,j-1, UH)
            v_neg_B = domain%velocity(:,j-1, VH)
        END IF

        IF (j < ny) THEN
            ! Extrapolate using points j-1, j, j+1
            theta_wd_pos_B = 5.0_dp*( &
                (min(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1)) &
                    - minimum_allowed_depth) / &
                (max(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1)) &
                    + 1000.0_dp*minimum_allowed_depth) &
                - 0.1_dp)
            theta_wd_pos_B = max(min(domain%theta, theta_wd_pos_B), ZERO_dp)
            
            ! FIXME: test only
            !theta_wd_pos_B = 1.0_dp

            stage_pos_B = extrapolate_edge_second_order(domain%U(:, j,STG), &
            !stage_pos_B = extrapolate_edge_second_order_nolimit(domain%U(:, j,STG), &
                domain%U(:, j-1,STG), domain%U(:, j+1, STG), theta_wd_pos_B, -1_ip)
            depth_pos_B = extrapolate_edge_second_order(domain%depth(:, j), &
            !depth_pos_B = extrapolate_edge_second_order_nolimit(domain%depth(:, j), &
                domain%depth(:, j-1), domain%depth(:, j+1), theta_wd_pos_B, -1_ip)
            u_pos_B = extrapolate_edge_second_order(domain%velocity(:, j, UH), &
            !u_pos_B = extrapolate_edge_second_order_nolimit(domain%velocity(:, j, UH), &
                domain%velocity(:, j-1,UH), domain%velocity(:, j+1, UH), theta_wd_pos_B, -1_ip)
            v_pos_B = extrapolate_edge_second_order(domain%velocity(:, j, VH), &
            !v_pos_B = extrapolate_edge_second_order_nolimit(domain%velocity(:, j, VH), &
                domain%velocity(:, j-1, VH), domain%velocity(:, j+1, VH), theta_wd_pos_B, -1_ip)

        ELSE
            ! j == ny, cannot extrapolate, so use the j value (first order accurate)
            theta_wd_pos_B = ZERO_dp
            stage_pos_B = domain%U(:, j, STG)
            depth_pos_B = domain%depth(:, j)
            u_pos_B = domain%velocity(:, j, UH)
            v_pos_B = domain%velocity(:, j, VH)
        END IF

    END SUBROUTINE

    !
    ! Get left edge values for the j'th North-South row. Because this
    ! is for a finite volume method, values at edges are discontinuous. 
    ! The edge has a 'positive' and 'negative' value, viewed from 
    ! cell i and i-1 respectively
    !
    ! @param domain the domain
    ! @param j get the j'th row
    ! @param nx domain number of cells along x direction
    ! @param theta_wd_neg_L 'theta' parameter controlling extrapolation as viewed from cell on negative side of edge
    ! @param theta_wd_pos_L 'theta' parameter controlling extrapolation as viewed from cell on positive side of edge
    ! @param stage_neg_L stage edge value as viewed from cell on negative side of edge
    ! @param stage_pos_L stage edge value as viewed from cell on positive side of edge
    ! @param depth_neg_L depth edge value as viewed from cell on negative side of edge
    ! @param depth_pos_L depth edge value as viewed from cell on positive side of edge
    ! @param u_neg_L x-velocity edge value as viewed from cell on negative side of edge
    ! @param u_pos_L x-velocity edge value as viewed from cell on positive side of edge
    ! @param v_neg_L y-velocity edge value as viewed from cell on negative side of edge
    ! @param v_pos_L y-velocity edge value as viewed from cell on positive side of edge
    !
    SUBROUTINE get_left_edge_values(domain, j, nx, &
        theta_wd_neg_L, theta_wd_pos_L, &
        stage_neg_L, stage_pos_L, &
        depth_neg_L, depth_pos_L, &
        u_neg_L, u_pos_L, &
        v_neg_L, v_pos_L)

        CLASS(domain_type), INTENT(IN):: domain
        INTEGER(ip), INTENT(IN) :: j, nx
        REAL(dp), INTENT(OUT):: theta_wd_neg_L(nx), theta_wd_pos_L(nx)
        REAL(dp), INTENT(OUT):: depth_neg_L(nx), depth_pos_L(nx)
        REAL(dp), INTENT(OUT):: stage_neg_L(nx), stage_pos_L(nx)
        REAL(dp), INTENT(OUT):: u_neg_L(nx), u_pos_L(nx)
        REAL(dp), INTENT(OUT):: v_neg_L(nx), v_pos_L(nx)

        ! Note: In practice the value in index (1) is not subsequently used for any variables

        ! View from negative side of edge

        theta_wd_neg_L(2:(nx-1)) = 5.0_dp * ( &
            (min(domain%depth(2:(nx-1), j), domain%depth(1:(nx-2), j), domain%depth(3:nx, j)) &
                - minimum_allowed_depth) / &
            (max(domain%depth(2:(nx-1), j), domain%depth(1:(nx-2), j), domain%depth(3:nx, j)) &
                + 1000.0_dp * minimum_allowed_depth) &
            - 0.1_dp) 
        theta_wd_neg_L(1) = ZERO_dp
        theta_wd_neg_L(nx) = ZERO_dp
        theta_wd_neg_L = max(min(domain%theta, theta_wd_neg_L), ZERO_dp)

        ! FIXME: test only
        !theta_wd_neg_L = 1.0_dp

        stage_neg_L(3:nx) = extrapolate_edge_second_order(domain%U(2:(nx-1), j, STG), &
        !stage_neg_L(3:nx) = extrapolate_edge_second_order_nolimit(domain%U(2:(nx-1), j, STG), &
            domain%U(1:(nx-2), j, STG), domain%U(3:nx, j, STG), theta_wd_neg_L(2:(nx-1)), 1_ip)
        stage_neg_L(2) = domain%U(1, j, STG)

        depth_neg_L(3:nx) = extrapolate_edge_second_order(domain%depth(2:(nx-1), j), &
        !depth_neg_L(3:nx) = extrapolate_edge_second_order_nolimit(domain%depth(2:(nx-1), j), &
            domain%depth(1:(nx-2), j), domain%depth(3:nx, j), theta_wd_neg_L(2:(nx-1)), 1_ip)
        depth_neg_L(2) = domain%depth(1, j)

        u_neg_L(3:nx) = extrapolate_edge_second_order(domain%velocity(2:(nx-1), j, UH), &
        !u_neg_L(3:nx) = extrapolate_edge_second_order_nolimit(domain%velocity(2:(nx-1), j, UH), &
            domain%velocity(1:(nx-2), j, UH), domain%velocity(3:nx, j, UH), theta_wd_neg_L(2:(nx-1)), 1_ip)
        u_neg_L(2) = domain%velocity(1, j, UH)

        v_neg_L(3:nx) = extrapolate_edge_second_order(domain%velocity(2:(nx-1), j, VH), &
        !v_neg_L(3:nx) = extrapolate_edge_second_order_nolimit(domain%velocity(2:(nx-1), j, VH), &
            domain%velocity(1:(nx-2), j, VH), domain%velocity(3:nx, j, VH), theta_wd_neg_L(2:(nx-1)), 1_ip)
        v_neg_L(2) = domain%velocity(1, j, VH)


        ! View from positive side of edge
 
        theta_wd_pos_L = theta_wd_neg_L

        stage_pos_L(2:(nx-1)) = extrapolate_edge_second_order(domain%U(2:(nx-1), j, STG), &
        !stage_pos_L(2:(nx-1)) = extrapolate_edge_second_order_nolimit(domain%U(2:(nx-1), j, STG), &
            domain%U(1:(nx-2), j, STG), domain%U(3:nx, j, STG), theta_wd_pos_L(2:(nx-1)), -1_ip)
        stage_pos_L(nx) = domain%U(nx, j, STG)

        depth_pos_L(2:(nx-1)) = extrapolate_edge_second_order(domain%depth(2:(nx-1), j), &
        !depth_pos_L(2:(nx-1)) = extrapolate_edge_second_order_nolimit(domain%depth(2:(nx-1), j), &
            domain%depth(1:(nx-2), j), domain%depth(3:nx, j), theta_wd_pos_L(2:(nx-1)), -1_ip)
        depth_pos_L(nx) = domain%depth(nx, j)
       
        u_pos_L(2:(nx-1)) = extrapolate_edge_second_order(domain%velocity(2:(nx-1), j, UH), &
        !u_pos_L(2:(nx-1)) = extrapolate_edge_second_order_nolimit(domain%velocity(2:(nx-1), j, UH), &
            domain%velocity(1:(nx-2), j, UH), domain%velocity(3:nx, j, UH), theta_wd_pos_L(2:(nx-1)), -1_ip)
        u_pos_L(nx) = domain%velocity(nx, j, UH)
         
        v_pos_L(2:(nx-1)) = extrapolate_edge_second_order(domain%velocity(2:(nx-1), j, VH), &
        !v_pos_L(2:(nx-1)) = extrapolate_edge_second_order_nolimit(domain%velocity(2:(nx-1), j, VH), &
            domain%velocity(1:(nx-2), j, VH), domain%velocity(3:nx, j, VH), theta_wd_pos_L(2:(nx-1)), -1_ip)
        v_pos_L(nx) = domain%velocity(nx, j, VH)
        
    END SUBROUTINE

    !
    ! Compute the fluxes, and other things, in preparation for an update of U 
    !
    ! Update values of:
    ! domain%flux_NS, domain%flux_EW, domain%max_dt, domain%explicit_source, 
    ! domain%explicit_source_VH_j_minus_1, domain%boundary_flux_store
    !
    ! @param domain the model domain type for which fluxes etc will be computed
    !
    SUBROUTINE compute_fluxes(domain)
        ! Compute fluxes for 2D shallow water equations on structured grid
        ! Use an Audusse type method, like ANUGA, but structured

        CLASS(domain_type), INTENT(INOUT):: domain

        INTEGER(ip):: i, j, nx, ny
        ! wavespeeds
        REAL(dp):: s_max, s_min, gs_pos, gs_neg, sminsmax
        ! stage/depth at + and - side of edge
        REAL(dp):: stage_pos, stage_neg, stage_pos_star, stage_neg_star, &
                   depth_pos, depth_pos_star, depth_pos_c, depth_pos_inv,&
                   depth_neg, depth_neg_star, depth_neg_c, depth_neg_inv,&
                   z_half, z_neg, z_pos
        ! velocities and momenta at + and - sides of edge
        REAL(dp):: u_pos, v_pos, u_neg, v_neg, ud_neg, vd_neg, ud_pos, vd_pos, vel_beta_neg, vel_beta_pos
        ! convenience variables
        REAL(dp):: denom, inv_denom, max_speed, max_dt, min_dt_inv, dx_cfl_inv(2)
        REAL(dp):: pressure_s, pressure_n, pressure_e, pressure_w
        REAL(dp):: depth_local(3), depth_inv
        REAL(dp):: half_cfl, max_dt_inv, source_tmp
        CHARACTER(len=charlen):: timer_name
        REAL(dp), PARAMETER :: EPS = 1.0e-06_dp, diffusion_scale = 1.0_dp
        REAL(dp), PARAMETER:: half_gravity = HALF_dp * gravity
        INTEGER(ip):: masscon_error, n_ext
        REAL(dp):: masscon_error_neg_depth
        
        ! Bottom edge values on 'positive' and 'negative' side (i.e. viewed from j and j-1 respectively)
        ! theta_wd controls the limiting
        REAL(dp) :: theta_wd_neg_B(domain%nx(1)), theta_wd_pos_B(domain%nx(1)), &
            stage_neg_B(domain%nx(1)), stage_pos_B(domain%nx(1)), &
            depth_neg_B(domain%nx(1)), depth_pos_B(domain%nx(1)), &
            u_neg_B(domain%nx(1)), u_pos_B(domain%nx(1)), &
            v_neg_B(domain%nx(1)), v_pos_B(domain%nx(1))
       
        ! Left edge values on 'positive' and 'negative' side (i.e. viewed from i and i-1 respectively) 
        ! theta_wd controls the limiting
        REAL(dp) :: theta_wd_neg_L(domain%nx(1)), theta_wd_pos_L(domain%nx(1)), &
            stage_neg_L(domain%nx(1)), stage_pos_L(domain%nx(1)), &
            depth_neg_L(domain%nx(1)), depth_pos_L(domain%nx(1)), &
            u_neg_L(domain%nx(1)), u_pos_L(domain%nx(1)), &
            v_neg_L(domain%nx(1)), v_pos_L(domain%nx(1))

        REAL(dp) :: bed_j_minus_1(domain%nx(1))

        half_cfl = HALF_dp * domain%cfl
        ny = domain%nx(2)
        nx = domain%nx(1)
        n_ext = domain%exterior_cells_width


        ! Set dt to a high value (it will drop)
        max_dt = domain%maximum_timestep
        ! By computing the inverse we avoid division in the loop
        max_dt_inv = ONE_dp/max_dt

        ! Zero the explicit source, recompute depth and velocity
        ! Must be updated before the main loop when done in parallel
        masscon_error = 0_ip
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, half_cfl, nx, ny) REDUCTION(MAX:max_dt_inv, masscon_error)
        !$OMP DO SCHEDULE(GUIDED)
        DO j = 1, domain%nx(2) !ny 
            domain%explicit_source(:,j,:) = ZERO_dp
            DO i = 1, domain%nx(1) !nx
                domain%depth(i,j) = domain%U(i,j,STG) - domain%U(i,j,ELV)
                IF(domain%depth(i,j) > minimum_allowed_depth) THEN
                    depth_inv = ONE_dp/domain%depth(i,j)
                    domain%velocity(i,j,UH) = domain%U(i,j,UH) * depth_inv
                    domain%velocity(i,j,VH) = domain%U(i,j,VH) * depth_inv
                ELSE
                    IF(domain%depth(i,j) < ZERO_dp) then
                        masscon_error = 1_ip
                    END IF
                    domain%velocity(i,j,UH) = ZERO_dp
                    domain%velocity(i,j,VH) = ZERO_dp
                END IF
            END DO
        END DO
        !$OMP END DO
    
        if(masscon_error > 0_ip) then
            print*, 'stage < bed --> mass conservation error'
            print*, minval(domain%depth), domain%nsteps_advanced
            call generic_stop()
        end if

        !Main loop

        !$OMP DO SCHEDULE(GUIDED)
        DO j = 2, domain%nx(2) !ny
            ! Assume distance_left_edge is constant but distance_bottom_edge
            ! might change with y (spherical coordinates). For cartesian coordinates
            ! we could move this outside the loop.
            dx_cfl_inv(1) = ONE_dp/(domain%distance_bottom_edge(j) * half_cfl)
            dx_cfl_inv(2) = ONE_dp/(domain%distance_left_edge(1) * half_cfl)

            ! Get bed at j-1 (might improve cache access) 
            bed_j_minus_1 = domain%U(:,j-1,ELV)

            ! Get bottom edge values
            CALL domain%get_bottom_edge_values(j, nx, ny, &
                theta_wd_neg_B, theta_wd_pos_B, &
                stage_neg_B, stage_pos_B, &
                depth_neg_B, depth_pos_B, &
                u_neg_B, u_pos_B, &
                v_neg_B, v_pos_B)

            ! Get left edge values
            CALL domain%get_left_edge_values(j, nx, &
                theta_wd_neg_L, theta_wd_pos_L, &
                stage_neg_L, stage_pos_L, &
                depth_neg_L, depth_pos_L, &
                u_neg_L, u_pos_L, &
                v_neg_L, v_pos_L)

            ! Could potentially vectorize this loop
            !!$OMP SIMD 
            !DO i = 2, domain%nx(1) !nx
            DO CONCURRENT (i = 2:nx)
                !
                ! North-South flux computation
                ! Compute flux_NS(i,j) = flux(i,j-1/2)
                !
                ! Cell i,j has North-South flux(i,j+1/2) at flux_NS index i,j+1, 
                !          and North-South flux(i,j-1/2) at flux_NS index i,j
                !

                ! Negative/positive bottom edge values
                stage_neg = stage_neg_B(i)
                stage_pos = stage_pos_B(i)
                depth_neg = depth_neg_B(i)
                depth_pos = depth_pos_B(i)
                u_neg = u_neg_B(i)
                u_pos = u_pos_B(i)
                v_neg = v_neg_B(i)
                v_pos = v_pos_B(i)
                depth_neg_c = domain%depth(i, j-1)
                depth_pos_c = domain%depth(i, j)

                ! Bed elevation
                z_neg = stage_neg - depth_neg
                z_pos = stage_pos - depth_pos
                z_half = max(z_neg, z_pos)

                stage_pos_star = max(stage_pos, z_half)
                stage_neg_star = max(stage_neg, z_half)

                ! Audusse Depths
                depth_neg_star = stage_neg_star - z_half
                depth_pos_star = stage_pos_star - z_half

                ! Velocity (in NS direction)
                IF(depth_neg_star == ZERO_dp) THEN
                    v_neg = ZERO_dp
                    ud_neg = ZERO_dp
                    vd_neg = ZERO_dp

                    gs_neg = ZERO_dp
                ELSE
                    ! Audusse type depth_integrated_velocity correction
                    vd_neg = v_neg * depth_neg_star
                    ud_neg = u_neg * depth_neg_star

                    ! Gravity wave celerity
                    gs_neg = sqrt(gravity * depth_neg_star)
                END IF
                
                ! Velocity (in NS direction)
                IF(depth_pos_star == ZERO_dp) THEN
                    v_pos = ZERO_dp
                    ud_pos = ZERO_dp
                    vd_pos = ZERO_dp

                    gs_pos = ZERO_dp
                ELSE
                    ! Correct depth-integrated_velocity (Audusse type approach)
                    vd_pos = v_pos * depth_pos_star
                    ud_pos = u_pos * depth_pos_star

                    ! Gravity wave celerity
                    gs_pos = sqrt(gravity * depth_pos_star)
                END IF

                ! Wave-celerities
                s_max = max(max(v_neg + gs_neg, v_pos + gs_pos), ZERO_dp)
                s_min = min(min(v_neg - gs_neg, v_pos - gs_pos), ZERO_dp)

                denom = s_max - s_min

                IF (denom > EPS) THEN
                    inv_denom = domain%distance_bottom_edge(j) / denom 
                    sminsmax = s_min * s_max * diffusion_scale
                    vel_beta_neg = v_neg * advection_beta
                    vel_beta_pos = v_pos * advection_beta

                    domain%flux_NS(i, j, STG) = (s_max * vd_neg - &
                                               s_min * vd_pos + &
                                               sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
                    domain%flux_NS(i, j, UH) = (s_max * ud_neg * vel_beta_neg - &
                                               s_min * ud_pos * vel_beta_pos + &
                                               sminsmax *(ud_pos - ud_neg)) * inv_denom
                    domain%flux_NS(i, j, VH) = (s_max * (vd_neg * vel_beta_neg + half_gravity * depth_neg_star*depth_neg_star) - & 
                                               s_min * (vd_pos * vel_beta_pos + half_gravity * depth_pos_star*depth_pos_star) + &
                                               sminsmax *(vd_pos - vd_neg)) * inv_denom
                ELSE
                    domain%flux_NS(i,j,STG) = ZERO_dp
                    domain%flux_NS(i,j,UH) = ZERO_dp
                    domain%flux_NS(i,j,VH) = ZERO_dp
                END IF

                !! Pressure gradient term at i,j-1, north side
                pressure_n = half_gravity * (depth_neg_star*depth_neg_star - depth_neg*depth_neg - &
                    (depth_neg + depth_neg_c)*(z_neg - bed_j_minus_1(i)))
               
                !! NOTE regarding OPENMP!!
                !! We cannot naively do:
                ! domain%explicit_source(i, j - 1 ,VH) = domain%explicit_source(i, j-1, VH) + pressure_n
                !! Since j-1 might also be updated on another OMP thread
                !! The solution is below
                !! Other pressure gradient cases refer to i,j, or i-1,j, so should be ok
                domain%explicit_source_VH_j_minus_1(i,j) = pressure_n * domain%distance_bottom_edge(j) ! Add to explicit_source later

                !! Pressure gradient term at i,j, south-side

                pressure_s = half_gravity * (depth_pos_star*depth_pos_star - depth_pos*depth_pos - &
                    (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)))
                
                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) - &
                    pressure_s * domain%distance_bottom_edge(j) 


                ! Timestep
                max_speed = max(s_max, -s_min)
                IF (max_speed > EPS) THEN
                    max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(2))
                END IF

        !
        !
        ! EW Flux computation
        !
        !
                ! left edge variables 
                stage_neg = stage_neg_L(i)
                stage_pos = stage_pos_L(i)
                depth_neg = depth_neg_L(i)
                depth_pos = depth_pos_L(i)
                u_neg = u_neg_L(i)
                u_pos = u_pos_L(i)
                v_neg = v_neg_L(i)
                v_pos = v_pos_L(i)
                depth_neg_c = domain%depth(i-1, j)
                depth_pos_c = domain%depth(i, j)

                ! Bed elevation
                z_neg = stage_neg - depth_neg
                z_pos = stage_pos - depth_pos
                z_half = max(z_neg, z_pos)

                stage_neg_star = max(stage_neg, z_half)
                stage_pos_star = max(stage_pos, z_half)

                ! Audusse Depths
                depth_neg_star = stage_neg_star - z_half
                depth_pos_star = stage_pos_star - z_half

                ! Velocity (in EW direction)
                IF(depth_neg_star == ZERO_dp) THEN
                    u_neg = ZERO_dp
                    ud_neg = ZERO_dp
                    vd_neg = ZERO_dp

                    gs_neg = ZERO_dp 
                ELSE
                    ! Audusse type depth-integrated-velocity correction
                    ud_neg = u_neg * depth_neg_star
                    vd_neg = v_neg * depth_neg_star

                    gs_neg = sqrt(gravity * depth_neg_star)
                END IF

                ! Velocity (in NS direction)
                IF(depth_pos_star == ZERO_dp) THEN
                    u_pos = ZERO_dp
                    ud_pos = ZERO_dp
                    vd_pos = ZERO_dp

                    gs_pos = ZERO_dp
                ELSE
                    ud_pos = u_pos * depth_pos_star
                    vd_pos = v_pos * depth_pos_star

                    gs_pos = sqrt(gravity * depth_pos_star)
                END IF

                ! Wave-celerities
                s_max = max(max(u_neg + gs_neg, u_pos + gs_pos), ZERO_dp)
                s_min = min(min(u_neg - gs_neg, u_pos - gs_pos), ZERO_dp)

                denom = s_max - s_min

                IF (denom > EPS) THEN
                    inv_denom = domain%distance_left_edge(i) / denom 
                    sminsmax = s_min * s_max * diffusion_scale
                    vel_beta_neg = u_neg * advection_beta
                    vel_beta_pos = u_pos * advection_beta

                    domain%flux_EW(i,j,STG) = (s_max * ud_neg - &
                                             s_min * ud_pos + &
                                             sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
                    domain%flux_EW(i,j,UH) = (s_max * (ud_neg * vel_beta_neg + half_gravity * depth_neg_star*depth_neg_star) - &
                                             s_min * (ud_pos * vel_beta_pos + half_gravity * depth_pos_star*depth_pos_star) + &
                                             sminsmax * (ud_pos - ud_neg)) * inv_denom
                    domain%flux_EW(i,j,VH) = (s_max * vd_neg * vel_beta_neg - &
                                             s_min * vd_pos * vel_beta_pos + &
                                             sminsmax * (vd_pos - vd_neg)) * inv_denom
                ELSE
                    domain%flux_EW(i,j,STG) = ZERO_dp
                    domain%flux_EW(i,j,UH) = ZERO_dp
                    domain%flux_EW(i,j,VH) = ZERO_dp
                END IF

                !! Pressure gradient term at i-1,j -- east side of cell i-1,j

                pressure_e = half_gravity * (depth_neg_star*depth_neg_star - depth_neg*depth_neg -&
                                                 (depth_neg + depth_neg_c)*(z_neg - domain%U(i-1 , j, ELV)))

                domain%explicit_source(i-1, j, UH) = domain%explicit_source(i-1, j ,UH) + &
                    pressure_e * domain%distance_left_edge(i)

                !! Pressure gradient term at i,j -- west_side

                pressure_w = half_gravity * (depth_pos_star*depth_pos_star - depth_pos*depth_pos -&
                                                 (depth_pos + depth_pos_c)*(z_pos - domain%U(i, j, ELV)))
                domain%explicit_source(i, j, UH) = domain%explicit_source(i, j, UH) - &
                    pressure_w * domain%distance_left_edge(i)

#ifdef SPHERICAL
                !! Source term associated with spherical coordinates for the uvd advection term
                !!
                !! Using the product rule for derivatives we can write
                ! R cos(lat) dlon * [(uvd)_{lat+} - (uvd)_{lat-}] = flux - source
                !! where
                ! (flux)  = [ (uvd R cos(lat) dlon)_{lat+} - (uvd R cos(lat) dlon)_{lat-} ] 
                ! (source) = uvd [(R cos(lat) dlon)_{lat+} - (R cos(lat) dlon)_{lat-}]
                !! This puts the equations in flux-conservative form as required for finite volumes.
                !! So here we add the source term

                domain%explicit_source(i, j, UH) = domain%explicit_source(i, j, UH) + &
                    advection_beta * domain%velocity(i, j, UH) * domain%U(i, j, VH) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))

                !
                ! Similarly for the v^2h + gh^2/2 term in the VH momentum equation
                !

                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) + &
                    advection_beta * domain%velocity(i, j, VH) * domain%U(i, j, VH) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))

                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) + &
                    half_gravity * domain%depth(i, j) * domain%depth(i, j) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))

#endif

                ! Timestep
                max_speed = max(s_max, -s_min)
                IF (max_speed > EPS) THEN
                    max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(1))
                END IF

            END DO
        END DO
        !$OMP END DO

        !$OMP END PARALLEL
        
        max_dt = ONE_dp/max_dt_inv
        domain%max_dt = max_dt

        !
        ! Spatially integrate boundary fluxes. Order is N, E, S, W. 
        ! Note we integrate in the INTERIOR (i.e. cut the edge rows/columns, then
        ! integrate the boundary of what remains). At the moment the 'edge' is
        ! only 1 cell wide, but one could imagine that changing (e.g. nested grid)
        !
        ! Outward boundary flux over the north
        domain%boundary_flux_store(1) = sum(domain%flux_NS((1+n_ext):(nx-n_ext),ny+1-n_ext,1))
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%flux_EW(nx+1-n_ext,(1+n_ext):(ny-n_ext),1))
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%flux_NS((1+n_ext):(nx-n_ext),1+n_ext,1))
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%flux_EW(1+n_ext,(1+n_ext):(ny-n_ext),1))

    END SUBROUTINE

    !
    ! Update the values of domain%U (i.e. the main flow variables), based on the fluxes and sources in domain
    !
    ! @param domain the domain to be updated
    ! @param dt timestep to advance
    !
    SUBROUTINE update_U(domain, dt)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), INTENT(IN):: dt

        INTEGER(ip) :: nx

        REAL(dp):: inv_cell_area_dt, depth, implicit_factor, dt_gravity, fs
        INTEGER(ip):: j, i, k, kk


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(GUIDED)
        DO j = 1, domain%nx(2)
        !DO concurrent (j = 1:domain%nx(2))
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !DO i = 1, nx 
            DO CONCURRENT (i = 1:domain%nx(1))
                depth = domain%depth(i,j)
        
                !! Fluxes
                DO kk = 1, 3
                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
                END DO

                ! Velocity clipping
                IF (depth <= minimum_allowed_depth) THEN
                    domain%U(i,j,UH) = ZERO_dp 
                    domain%U(i,j,VH) = ZERO_dp 
                ELSE
                    ! Implicit friction slope update
                    ! U_new = U_last + U_explicit_update - dt*depth*friction_slope_multiplier*U_new

                    ! If we multiply this by UH or VH, we get the associated friction slope term
                    fs = domain%manning_squared(i,j) * &
                        sqrt(domain%velocity(i,j,UH) * domain%velocity(i,j,UH) + &
                             domain%velocity(i,j,VH) * domain%velocity(i,j,VH)) * &
                        (max(depth, minimum_allowed_depth)**(NEG_SEVEN_ON_THREE_dp))

                    implicit_factor = ONE_dp/(ONE_dp + dt_gravity*depth*fs)

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
                END IF
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + dt
    
    END SUBROUTINE

    ! 
    ! Standard forward euler 1st order timestepping scheme 
    !
    ! This is also used as a component of more advanced timestepping schemes
    ! Argument 'timestep' is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced
    !
    SUBROUTINE one_euler_step(domain, timestep)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), OPTIONAL, INTENT(IN):: timestep

        CHARACTER(len=charlen):: timer_name
        REAL(dp):: ts

        CALL domain%update_boundary()

        TIMER_START('flux')
        CALL domain%compute_fluxes()
        TIMER_STOP('flux')

        TIMER_START('update')
        IF(PRESENT(timestep)) THEN
            ts = timestep
        ELSE
            ts = domain%max_dt
        END IF

        CALL domain%update_U(ts)

        ! Track flux through boundaries
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral + &
            ts * sum(domain%boundary_flux_store)

        TIMER_STOP('update')

        ! Coarray communication, if required
        !TIMER_START('comms')
        !CALL domain%comms%communicate(domain%U)
        !TIMER_STOP('comms')


    END SUBROUTINE

    !
    ! Standard 2-step second order timestepping runge-kutta scheme
    ! Argument timestep is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced
    !
    !
    SUBROUTINE one_rk2_step(domain, timestep)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), OPTIONAL, INTENT(IN):: timestep

        REAL(dp):: backup_time, dt_first_step
        INTEGER(ip):: j
        CHARACTER(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        CALL domain%backup_quantities()
        
        ! First euler step
        IF(PRESENT(timestep)) THEN
            CALL domain%one_euler_step(timestep)
            dt_first_step = timestep
        ELSE
            CALL domain%one_euler_step()
            dt_first_step = domain%max_dt
        END IF
       
        ! Second euler step with the same timestep 
        CALL domain%one_euler_step(dt_first_step)

        TIMER_START('average')

        ! Take average (but allow for openmp)
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        DO j = 1, domain%nx(2)
            domain%U(:, j, [STG, UH, VH]) = HALF_dp * (domain%U(:, j, [STG, UH, VH]) +&
                domain%backup_U(:, j, [STG, UH, VH]))
        END DO
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time (since we updated twice) and boundary flux integral
        domain%time = backup_time + dt_first_step
        domain%boundary_flux_evolve_integral = HALF_dp * domain%boundary_flux_evolve_integral
        
        TIMER_STOP('average')

    END SUBROUTINE
    

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
    ! @param domain the domain to be updated
    ! @param timestep the timestep, by which the solution is advanced (n-1) times
    !
    SUBROUTINE one_rk2n_step(domain, timestep)
        ! Advance (n-1) * timesteps in this routine, with 2nd order in time accuracy
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), OPTIONAL, INTENT(IN):: timestep

        INTEGER(ip), PARAMETER:: n = 5 ! Number of substeps to take, must be > 2
        REAL(dp), PARAMETER:: n_inverse = ONE_dp / (ONE_dp * n)
        REAL(dp):: backup_time, dt_first_step, reduced_dt, backup_flux_integral
        INTEGER(ip):: j,k
        CHARACTER(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        CALL domain%backup_quantities()
     
        IF(PRESENT(timestep)) THEN 
            ! First step 
            CALL domain%one_euler_step(timestep)
            ! Store timestep
            dt_first_step = timestep
        ELSE
            ! First step 
            CALL domain%one_euler_step()
            ! Store timestep
            dt_first_step = domain%max_dt
        END IF

        ! Steps 2, n-1
        DO j = 2, n-1        
            CALL domain%one_euler_step(dt_first_step)
        END DO

        ! Store 1/n * (original_u - u)

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        DO k = 1, 3
            DO j = 1, domain%nx(2)
                domain%backup_U(:,j,k) = (domain%backup_U(:,j, k) - domain%U(:,j,k)) * n_inverse
            END DO
        END DO
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Store this for boundary flux tracking
        backup_flux_integral = domain%boundary_flux_evolve_integral*n_inverse

        ! Now take one step of duration (n-1)/n * dt
        reduced_dt = (ONE_dp * n - ONE_dp) * n_inverse * dt_first_step
        CALL domain%one_euler_step(reduced_dt)
        
        TIMER_START('final_update')

        ! Final update
        ! domain%U = domain%U + domain%backup_U
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        DO k = 1, 3
            DO j = 1, domain%nx(2)
                domain%U(:, j, k) = domain%U(:, j, k) + domain%backup_U(:, j, k)
            END DO
        END DO
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Get final boundary flux integral
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral - &
            backup_flux_integral

        ! Fix time and timestep (since we updated (n-1)*dt regular timesteps)
        domain%time = backup_time + dt_first_step * (n-1) 

        TIMER_STOP('final_update')

    END SUBROUTINE


    !
    ! Another 2nd order timestepping scheme (like the trapezoidal rule) 
    ! Argument timestep is optional, but if provided should satisfy the CFL condition
    !
    ! @param domain the domain to be updated
    ! @param timestep the timestep by which the solution is advanced
    !
    SUBROUTINE one_midpoint_step(domain, timestep)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), OPTIONAL, INTENT(IN) :: timestep

        REAL(dp):: backup_time, dt_first_step
        INTEGER(ip):: j


        ! Backup quantities
        backup_time = domain%time
        CALL domain%backup_quantities()
        
        CALL domain%update_boundary()

        ! First euler sub-step
        CALL domain%compute_fluxes()

        IF(PRESENT(timestep)) THEN
            dt_first_step = timestep
            CALL domain%update_U(dt_first_step*HALF_dp)
        ELSE
            dt_first_step = domain%max_dt
            CALL domain%update_U(dt_first_step*HALF_dp)
        END IF

        !TIMER_START('comms')
        !CALL domain%comms%communicate(domain%U)
        !TIMER_STOP('comms')

        CALL domain%update_boundary()
        
        ! Compute fluxes 
        CALL domain%compute_fluxes()

        ! Set U back to backup_U
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        DO j = 1, domain%nx(2)
            domain%U(:, j, [STG,UH,VH]) = domain%backup_U(:, j, [STG,UH,VH])
        END DO
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time
        domain%time = backup_time

        ! Update U
        CALL domain%update_U(dt_first_step)

        !TIMER_START('comms')
        !CALL domain%comms%communicate(domain%U)
        !TIMER_STOP('comms')

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store)*dt_first_step

    END SUBROUTINE

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
    SUBROUTINE one_linear_leapfrog_step(domain, dt)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), INTENT(IN):: dt

        INTEGER(ip) :: nx, ny

        REAL(dp):: inv_cell_area_dt, inv_cell_area_dt_vh_g, inv_cell_area_dt_g
        REAL(dp):: dw_j(domain%nx(1)), h_jph_vec(domain%nx(1)), h_iph_vec(domain%nx(1))
        INTEGER(ip):: j, i, xL, xU, yL, yU, n_ext

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        LOGICAL, PARAMETER:: truely_linear = .TRUE.

        !
        ! idea: U(i, j, UH) = UH_{i+1/2, j}
        !     : U(i, j, VH) = UH_{i, j+1/2}
        ! We rely on boundary conditions for ALL boundary stage values, and
        ! apply the boundary condition after the stage update. We then update
        ! all uh/vh values which can be updated using those stage values [so avoid
        ! having to specify boundary conditions for them]

        TIMER_START('LF_update')

        nx = domain%nx(1)
        ny = domain%nx(2)
      
        xL = domain%xL
        xU = domain%xU
        yL = domain%yL
        yU = domain%yU
        n_ext = domain%exterior_cells_width

        !CALL domain%update_boundary()
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU)
        
        !
        ! Update stage
        !
        !$OMP DO SCHEDULE(STATIC)
        !DO j = 2, ny-1
        DO j = (yL+1),(yU-1)
        !DO concurrent (j = 1:domain%nx(2))
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !DO i = 1, nx 
            !DO CONCURRENT (i = 2:(nx-1))
            DO CONCURRENT (i = (xL+1):(xU-1))
                ! dstage/dt = - 1/(R cos (lat)) [ d uh / dlon + dvhcos(lat)/dlat ]
                domain%U(i, j, STG) = domain%U(i, j, STG) - inv_cell_area_dt * &
                    ((domain%U(i, j, UH) - domain%U(i-1, j, UH))*domain%distance_left_edge(i) + &
                    (domain%U(i, j, VH)*domain%distance_bottom_edge(j+1) - &
                        domain%U(i, j-1, VH)*domain%distance_bottom_edge(j)))
        
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        
        domain%time = domain%time + dt * HALF_dp
        CALL domain%update_boundary()

        !
        ! Update uh, vh
        !

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt, nx, ny, xL, xU, yL, yU)
        !$OMP DO SCHEDULE(STATIC)
        !DO j = 1, ny - 1
        DO j = yL, yU - 1
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            ! For leap-frog, the area associated with the NS momentum term is different
            inv_cell_area_dt_g = gravity * dt / domain%area_cell_y(j)
            inv_cell_area_dt_vh_g = gravity * dt / (HALF_dp * (domain%area_cell_y(j) + domain%area_cell_y(j+1)))
       
            ! 
            ! Try to keep control-flow and non-local memory jumps out of inner loop
            ! This improves speed on my machine with gfortran (11/08/2016)
            !
            dw_j = domain%U(:, j+1, STG) - domain%U(:, j, STG)

            IF(truely_linear) THEN
                !
                ! In the g * d * dStage/dx type term, let d be constant 
                !

                ! Depth at j-plus-half
                h_jph_vec(xL:xU) = merge(domain%MSL_linear - HALF_dp * (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:xU,j+1,ELV) < -minimum_allowed_depth + domain%MSL_linear).AND. &
                     ( domain%U(xL:xU,j,ELV)   < -minimum_allowed_depth + domain%MSL_linear)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    domain%MSL_linear - HALF_dp * (domain%U((xL+1):xU, j, ELV) + domain%U((xL):(xU-1), j, ELV)), &
                    ZERO_dp, &
                    (( domain%U(xL:(xU-1),j,ELV) < -minimum_allowed_depth + domain%MSL_linear).AND.&
                     ( domain%U((xL+1):xU,j,ELV) < -minimum_allowed_depth + domain%MSL_linear)))  
            ELSE
                !
                ! In the g * d * dStage/dx type term, let d vary. This means
                ! the equations are not actually linear!
                !

                ! Depth at j-plus-half
                h_jph_vec(xL:xU) = merge(&
                    HALF_dp * ((domain%U(xL:xU,j+1,STG) + domain%U(xL:xU,j,STG)) - &
                               (domain%U(xL:xU,j+1,ELV) + domain%U(xL:xU,j,ELV))), &
                    ZERO_dp, &
                    ((domain%U(xL:xU,j+1,STG) - domain%U(xL:xU,j+1,ELV) > minimum_allowed_depth).AND. &
                     (domain%U(xL:xU,j,STG) - domain%U(xL:xU,j,ELV) > minimum_allowed_depth)))

                ! Depth at i-plus-half
                h_iph_vec(xL:(xU-1)) = merge( &
                    HALF_dp * ((domain%U((xL+1):xU, j, STG)  + domain%U(xL:(xU-1), j, STG)) -&
                               (domain%U((xL+1):xU, j, ELV) - domain%U(xL:(xU-1), j, ELV))), &
                    ZERO_dp, &
                    ((domain%U(xL:(xU-1), j, STG) - domain%U(xL:(xU-1),j,ELV) > minimum_allowed_depth).AND.&
                        (domain%U((xL+1):xU, j, STG) - domain%U((xL+1):xU,j,ELV) > minimum_allowed_depth)))  


            END IF

            !DO CONCURRENT (i = 1:(nx-1))
            DO CONCURRENT (i = xL:(xU-1))

                ! duh/dt = - g * h0/(R cos (lat)) [ d stage / dlon ]
                domain%U(i, j, UH) = domain%U(i, j, UH) - inv_cell_area_dt_g * h_iph_vec(i) *&
                    (domain%U(i+1, j, STG) - domain%U(i, j, STG)) * domain%distance_left_edge(i+1)

                ! dvh/dt = - g * h0/(R) [ d stage / dlat ]
                domain%U(i, j, VH) = domain%U(i, j, VH) - inv_cell_area_dt_vh_g * h_jph_vec(i) *&
                    dw_j(i) * domain%distance_bottom_edge(j+1)
        
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        domain%time = domain%time + HALF_dp*dt

        !! Boundary flux integration
        ! note : U(i, j, UH) = UH_{i+1/2, j}
        !      : U(i, j, VH) = UH_{i, j+1/2}
        ! Outward boundary flux over the north -- integrate over the 'interior' cells
        domain%boundary_flux_store(1) = sum(domain%U((1+n_ext):(nx-n_ext),ny-n_ext,VH)) * domain%distance_bottom_edge(ny-n_ext+1)
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%U(nx-n_ext,(1+n_ext):(ny-n_ext),UH)) * domain%distance_left_edge(nx-n_ext+1)
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%U((1+n_ext):(nx-n_ext),n_ext,VH)) * domain%distance_bottom_edge(n_ext+1)
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%U(n_ext,(1+n_ext):(ny-n_ext),UH)) * domain%distance_left_edge(n_ext+1)

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store) * dt

        TIMER_STOP('LF_update')

        !TIMER_START('comms')
        !CALL domain%comms%communicate(domain%U)
        !TIMER_STOP('comms')

    END SUBROUTINE
    
    !
    ! Routine to run all boundary conditions
    !
    SUBROUTINE update_boundary(domain)
        CLASS(domain_type), INTENT(INOUT):: domain 

        TIMER_START('boundary_update')

        IF(associated(domain%boundary_subroutine)) THEN
            CALL domain%boundary_subroutine(domain)
        END IF

        TIMER_STOP('boundary_update')

    END SUBROUTINE

    !
    ! Stub routine to allow the user to not provide a boundary condition
    !
    SUBROUTINE default_boundary_subroutine(domain)
        TYPE(domain_type), INTENT(INOUT):: domain
        ! Do nothing (default case)
    END SUBROUTINE

    !
    ! Copy domain%U to domain%backup_U
    !
    SUBROUTINE backup_quantities(domain)
        
        CLASS(domain_type), INTENT(INOUT):: domain
        INTEGER(ip):: j, k

        TIMER_START('backup')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        DO k = 1, 3
            DO j = 1, domain%nx(2)
                domain%backup_U(:, j, k) = domain%U(:, j, k)
            END DO
        END DO
        !$OMP END DO 
        !$OMP END PARALLEL

        TIMER_STOP('backup')

    END SUBROUTINE

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
    SUBROUTINE evolve_one_step(domain, timestep)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), OPTIONAL, INTENT(IN):: timestep

        REAL(dp):: time0


        time0 = domain%time

        ! Reset the boundary fluxes integrated over the evolve step to zero
        domain%boundary_flux_evolve_integral = ZERO_dp
   
        SELECT CASE (domain%timestepping_method) 
        
        CASE ('linear')
            IF(PRESENT(timestep)) THEN
                CALL domain%one_linear_leapfrog_step(timestep)
            ELSE
                print*, 'ERROR: timestep must be provided for linear evolve_one_step'
                call generic_stop()
            END IF
        CASE DEFAULT
            print*, 'ERROR: we are only supporting the linear solver herein'
            call generic_stop()
        END SELECT

        ! For some problems updating max U can take a significant fraction of the time,
        ! so we allow it to only be done occasionally
        if(mod(domain%nsteps_advanced, domain%max_U_update_frequency) == 0) then
            CALL domain%update_max_quantities()
        end if

        ! Record the timestep here
        domain%evolve_step_dt = domain%time - time0

        ! Update the boundary flux time integral
        domain%boundary_flux_time_integral = domain%boundary_flux_time_integral + &
            domain%boundary_flux_evolve_integral

        domain%nsteps_advanced = domain%nsteps_advanced + 1

    END SUBROUTINE

    !
    ! Convenience function to compute the volume of water in the 'interior' of the domain
    ! This involves all parts of the domain that are more than domain%exterior_cells_width
    ! from the edge.
    !
    FUNCTION volume_interior(domain) RESULT(domain_volume)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp) :: domain_volume

        INTEGER(ip) :: j, n_ext
        
        ! Volume on the interior. At the moment the interior is all
        ! but the outer cells of the domain, but that could change.

        TIMER_START('volume_interior')
        n_ext = domain%exterior_cells_width
        domain_volume = ZERO_dp
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, n_ext) REDUCTION(+:domain_volume)
        !$OMP DO SCHEDULE(STATIC)
        DO j = (1+n_ext), (domain%nx(2) - n_ext) 
            domain_volume = domain_volume + domain%area_cell_y(j) * &
                (sum(domain%U((1+n_ext):(domain%nx(1)-n_ext),j,STG) - &
                     domain%U((1+n_ext):(domain%nx(1)-n_ext),j,ELV)))
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        TIMER_STOP('volume_interior')

    END FUNCTION

    !
    ! Compute the volume on interior cells and add to
    ! the integrated boundary flux. This should sum to a constant
    ! in the absence of mass sources.
    !
    ! Note this relies on the time-stepping routines correctly
    ! computing the boundary_flux_time_integral
    !
    FUNCTION mass_balance_interior(domain) RESULT(mass_balance)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp) :: mass_balance

        mass_balance = domain%volume_interior() + domain%boundary_flux_time_integral 

    END FUNCTION

    !
    ! Function to compute the max timestep allowed for the linear shallow water equations,
    ! using the provided CFL.
    !
    FUNCTION linear_timestep_max(domain) RESULT(timestep)
        CLASS(domain_type), INTENT(IN):: domain
        REAL(dp) :: timestep

        REAL(dp) :: ts_max
        INTEGER(ip) :: i, j

        ts_max = HUGE(1.0_dp)
    
        DO j = 1, domain%nx(2)
            DO i = 1, domain%nx(1)
                ! max timestep = Distance along latitude / wave-speed <= dt
                ! Beware -- might need to use a different CFL number?
                ts_max = min(ts_max, &
                    0.5_dp * min(&
                        (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)),& 
                        (domain%distance_left_edge(i+1)   + domain%distance_left_edge(i)  ) ) / &
                    sqrt(gravity * max(-domain%U(i,j,ELV), minimum_allowed_depth)) )
            END DO
        END DO
        timestep = ts_max * domain%cfl

    END FUNCTION

    !
    ! Initialise output files
    !
    SUBROUTINE create_output_files(domain)
        CLASS(domain_type), INTENT(INOUT):: domain

        CHARACTER(len=charlen):: mkdir_command, cp_command, t1, t2, t3, &
                                 output_folder_name
        INTEGER(ip):: i, metadata_unit

        ALLOCATE(domain%output_variable_unit_number(domain%nvar))

        ! Create output directory
        CALL DATE_AND_TIME(t1, t2, t3)
        ! Get domain id as a character
        WRITE(t3, '(I6)') 100000 + domain%myid

        output_folder_name = TRIM(domain%output_basedir) // '/RUN_ID' // TRIM(t3) // &
            '_' // TRIM(t1) // '_' // TRIM(t2)
        mkdir_command = 'mkdir -p ' // TRIM(output_folder_name)
        !CALL EXECUTE_COMMAND_LINE(TRIM(mkdir_command))
        CALL SYSTEM(TRIM(mkdir_command))


        ! Make a filename to hold domain metadata, and write the metadata
        t1 = TRIM(output_folder_name) // '/' // 'Domain_info_ID' // TRIM(t3) // '.txt'
        domain%metadata_ascii_filename = t1
        OPEN(newunit = metadata_unit, file=domain%metadata_ascii_filename)
        WRITE(metadata_unit, *) 'id :', domain%myid
        WRITE(metadata_unit, *) 'nx :', domain%nx
        WRITE(metadata_unit, *) 'dx :', domain%dx
        WRITE(metadata_unit, *) 'lower_left_corner: ', domain%lower_left
        WRITE(metadata_unit, *) 'dp_precision: ', dp
        WRITE(metadata_unit, *) 'ip_precision: ', ip
        WRITE(metadata_unit, *) 'output_precision: ', output_precision
        CLOSE(metadata_unit)

        DO i=1, domain%nvar
            ! Make output_file_name
            WRITE(t1,'(I1)') i
            WRITE(t2, *) 'Var_', TRIM(t1), '_ID', TRIM(t3)
            t1 = TRIM(output_folder_name) // '/' // ADJUSTL(TRIM(t2))

            ! Open the files

            ! Ascii
            !OPEN(newunit = domain%output_variable_unit_number(i), file = t1)

            ! Binary. Using 'stream' access makes it easy to read in R
            OPEN(newunit = domain%output_variable_unit_number(i), file = t1, &
                access='stream', form='unformatted')
        END DO

        ! Make a time file_name. Store as ascii
        t1 = TRIM(output_folder_name) // '/' // 'Time_ID' // TRIM(t3) // '.txt'
        OPEN(newunit = domain%output_time_unit_number, file = t1)

        ! Copy code to the output directory

        mkdir_command = 'mkdir -p ' // TRIM(output_folder_name) // '/Code'
        cp_command = 'cp *.f* *.R *.c make* ' // TRIM(output_folder_name) // '/Code'
        !CALL EXECUTE_COMMAND_LINE(TRIM(mkdir_command))
        CALL SYSTEM(TRIM(mkdir_command))
        !CALL EXECUTE_COMMAND_LINE(TRIM(cp_command))
        CALL SYSTEM(TRIM(cp_command))

        domain%output_folder_name = output_folder_name

    END SUBROUTINE create_output_files

    !
    ! Write output to files (must call create_output_files first).
    ! If time_only=.TRUE., only write the time. This is used to avoid
    ! writing the main model grids. Typically useful when values at point-gauges
    ! are being recorded, and we want to store the time too, but it would
    ! take too much disk to store the model grids
    !
    SUBROUTINE write_to_output_files(domain, time_only)
        CLASS(domain_type), INTENT(INOUT):: domain
        LOGICAL, OPTIONAL, INTENT(IN) :: time_only
        INTEGER(ip):: i, j
        LOGICAL:: to

        if(present(time_only)) then
            to = time_only
        else
            to = .FALSE.
        end if

        TIMER_START('fileIO')
       
        if(to .eqv. .FALSE.) then 
            DO i = 1, domain%nvar
                DO j = 1, domain%nx(2)
                    ! Binary
                    WRITE(domain%output_variable_unit_number(i)) REAL(domain%U(:,j,i), output_precision)
                END DO
            END DO
        end if

        ! Time too, as ascii
        WRITE(domain%output_time_unit_number, *) domain%time

        TIMER_STOP('fileIO')

    END SUBROUTINE write_to_output_files

    !
    ! Set the domain%logfile_unit to point to an actual file
    !
    SUBROUTINE divert_logfile_unit_to_file(domain)
        CLASS(domain_type), INTENT(INOUT):: domain
        CHARACTER(len=charlen):: logfile_name, domain_ID

        WRITE(domain_ID, '(I6)') 100000 + domain%myid

        if(domain%output_folder_name == '') then
            logfile_name = './logfile_' // TRIM(domain_ID) // '.log'
        else
            logfile_name = TRIM(domain%output_folder_name) // '/logfile_' // TRIM(domain_ID) // '.log'
        end if

        open(newunit=domain%logfile_unit, file=logfile_name)

    END SUBROUTINE

    !
    ! Keep track of the maxima of stage (only) 
    !
    SUBROUTINE update_max_quantities(domain)
        CLASS(domain_type), INTENT(INOUT):: domain
        INTEGER(ip):: j, k, i

        if(domain%record_max_U) then
            TIMER_START('update_max_quantities')

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)

            !! Previously we stored all quantities in max_U. However, generally
            !! max uh and max vh don't have much meaning by themselves. So now
            !! we just store max stage
            !!$OMP DO SCHEDULE(GUIDED), COLLAPSE(2)
            !DO k = 1, domain%nvar

            !! UPDATE: Only record max stage
            !$OMP DO SCHEDULE(GUIDED)
                DO j = domain%yL, domain%yU !1, domain%nx(2)
                    DO i = domain%xL, domain%xU !1, domain%nx(1)
                        domain%max_U(i,j,1) = max(domain%max_U(i,j,1), domain%U(i,j,1))
                    END DO
                END DO
            !END DO
            !$OMP END DO
            !$OMP END PARALLEL
            
            TIMER_STOP('update_max_quantities')
        end if

    END SUBROUTINE update_max_quantities

    !
    ! Write max quantities to a file (usually just called once at the end of a
    ! simulation). 
    !
    ! Currently we only write stage, followed by the elevation. Although the
    ! latter usually doesn't evolve, we generally want both the max stage and the
    ! elevation for plotting purposes, so it is saved here too.
    !
    SUBROUTINE write_max_quantities(domain)
        CLASS(domain_type), INTENT(IN):: domain
        CHARACTER(len=charlen):: max_quantities_filename
        INTEGER:: i, j, k

        if(domain%record_max_U) then

            max_quantities_filename = TRIM(domain%output_folder_name) // '/Max_quantities'

            OPEN(newunit=i, file=max_quantities_filename, access='stream', form='unformatted')
        
            !DO k = 1, domain%nvar
            !! Update: Only record max stage
            k = 1
                DO j = 1, domain%nx(2)
                    WRITE(i) REAL(domain%max_U(:,j,k), output_precision)
                END DO
            !END DO

            ! Also store elevation since it is typically useful, and we might not store it otherwise
            DO j = 1, domain%nx(2)
                WRITE(i) REAL(domain%U(:,j,ELV), output_precision)
            END DO
            
            CLOSE(i)

        end if
 
    END SUBROUTINE

    !
    ! Set up point gauges, which record values of variables at cells nearest
    ! the given xy points.
    ! 
    ! @param domain The domain within which we will record outputs (in domain%U)
    ! @param xy_coords numeric array of x,y coordinates with 2 rows and as many
    !  columns as points. These must all be inside the extents of the domain
    ! @param time_series_var (optional) array with the indices of U to store at
    !  gauges each timestep, Default is [1, 2, 3] to store stage, uh, vh
    ! @param static_var (optional) array with the indices of U to store at
    !  gauges only once. Default is [4] to store elevation
    ! @param gauge_ids (optional) an REAL ID for each gauge. Default gives
    !  (1:size(xy_coords(1,:))) * 1.0. Even though integers are natural we store
    !   as REAL to avoid precision loss issues
    !
    SUBROUTINE setup_point_gauges(domain, xy_coords, time_series_var, static_var, gauge_ids,&
        attribute_names, attribute_values)
        CLASS(domain_type), INTENT(INOUT):: domain
        REAL(dp), INTENT(IN) :: xy_coords(:,:)
        INTEGER(ip), OPTIONAL, INTENT(IN):: time_series_var(:), static_var(:)
        REAL(dp), OPTIONAL, INTENT(IN):: gauge_ids(:)
        CHARACTER(charlen), OPTIONAL, INTENT(IN):: attribute_names(:), attribute_values(:)

        CHARACTER(charlen) :: netcdf_gauge_output_file, t3
        INTEGER(ip), ALLOCATABLE:: tsv(:), sv(:)
        REAL(dp), ALLOCATABLE:: gauge_ids_local(:)
        INTEGER(ip):: i

        if(present(time_series_var)) then
            ALLOCATE(tsv(size(time_series_var)))
            tsv = time_series_var
        else
            ! Default case -- store stage/uh/vh every output step
            ALLOCATE(tsv(3))
            tsv = [STG, UH, VH]
        end if

        if(present(static_var)) then
            ALLOCATE(sv(size(static_var)))
            sv = static_var
        else 
            ! Default case -- store elevation once
            ALLOCATE(sv(1))
            sv = [ELV]
        end if

        if(present(gauge_ids)) then
            if(size(gauge_ids) /= size(xy_coords(1,:))) then
                print*, 'Number of gauge ids does not equal number of coordinates' 
            end if
            ALLOCATE(gauge_ids_local(size(gauge_ids)))
            gauge_ids_local = gauge_ids
        else
            ! Default case -- give sequential integer ids
            ALLOCATE(gauge_ids_local(size(xy_coords(1,:))))
            DO i = 1, size(gauge_ids_local)
                gauge_ids_local(i) = i*1.0_dp
            END DO
        end if
    
        ! Get domain id as a character and put it in the output file name
        WRITE(t3, '(I6)') 100000 + domain%myid
        netcdf_gauge_output_file = TRIM(domain%output_folder_name) // '/' // &
            'Gauges_data_ID' // TRIM(t3) // '.nc'

        call domain%point_gauges%allocate_gauges(xy_coords, tsv, sv, gauge_ids_local)

        if((present(attribute_names)).AND.(present(attribute_values))) then
            call domain%point_gauges%initialise_gauges(domain%lower_left, domain%dx, &
                domain%nx, domain%U, netcdf_gauge_output_file, &
                attribute_names=attribute_names, attribute_values=attribute_values)
        else
            call domain%point_gauges%initialise_gauges(domain%lower_left, domain%dx, &
                domain%nx, domain%U, netcdf_gauge_output_file)
        end if

    END SUBROUTINE

    !
    ! Write the values of point gauges to a file
    !
    SUBROUTINE write_gauge_time_series(domain)
        CLASS(domain_type), INTENT(INOUT):: domain

        if(allocated(domain%point_gauges%time_series_values)) then
            call domain%point_gauges%write_current_time_series(domain%U, domain%time)
        end if

    END SUBROUTINE

    !
    ! Routine to call once we no longer need the domain. One case where this is
    ! important is when using netcdf output -- since if the files are not closed,
    ! then they may not be completely written out.
    !
    SUBROUTINE finalise_domain(domain)
        CLASS(domain_type), INTENT(INOUT):: domain

        ! Close the gauges netcdf file -- since otherwise it might not finish
        ! writing.
        call domain%point_gauges%finalise()
    
        ! Flush all open file units
        call flush()

    END SUBROUTINE

END MODULE domain_mod
