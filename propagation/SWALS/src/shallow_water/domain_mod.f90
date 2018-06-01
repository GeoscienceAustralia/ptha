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
                          minimum_allowed_depth, &
                          default_timestepping_method, extrapolation_theta,&
                          wall_elevation, &
                          default_output_folder
    use timer_mod, only: timer_type
    use point_gauge_mod, only: point_gauge_type
    use coarray_utilities_mod, only: partitioned_domain_nesw_comms_type
    use nested_grid_comms_mod, only: two_way_nesting_comms_type
    use stop_mod, only: generic_stop
    use iso_fortran_env, only: output_unit
    use iso_c_binding, only: C_DOUBLE

#ifdef SPHERICAL    
    ! Compile with -DSPHERICAL to get the code to run in spherical coordinates
    !
    ! Can only have coriolis if we have spherical. However, the code
    ! will still work if we have spherical and not coriolis. But by default,
    ! we include coriolis with spherical. Add the flag 'NOCORIOLIS' to remove it
    !
    ! Note coriolis is only implemented for linear solver.
#ifndef NOCORIOLIS
#define CORIOLIS
#endif
    use spherical_mod, only: area_lonlat_rectangle, deg2rad, earth_angular_freq
    use global_mod, only: radius_earth
#endif

#ifndef NOOPENMP
    use omp_lib
#endif

    implicit none

    ! Make everything private, except domain_type, which has its own methods,
    ! and the domain_metadata_type, which is sometimes more convenient to
    ! work with than the domain [since it is 'lightweight']
    private
    public:: domain_type !, domain_metadata_type

    ! Indices for arrays: Stage, depth-integrated-x-velocity,
    ! depth-integrated-v-velocity, elevation. So e.g. stage
    ! is in domain%U(:,:,STG), and elevation is in domain%U(:,:ELV)
    integer(4), parameter, public:: STG=1, UH=2, VH=3, ELV=4

    real(dp), parameter :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp
    real(dp), parameter:: NEG_SEVEN_ON_THREE_dp = -7.0_dp/3.0_dp

    type :: domain_type
    !
    ! Type to only hold domain metadata
    !
    ! This can be useful when reading in domains, prior to allocation. For example,
    ! we can manipulate the domain bounding box to support nesting buffers, etc
    !
    ! However, using non-type-bound procedures with the CLASS attribute of
    ! dummy arguments (as required for polymorphism) seems to have a bug
    ! with older gfortran (e.g. 4.8, which is default on ubuntu 14.04). So 
    ! be careful about using such subroutines/functions in the code.
    !
    !type :: domain_metadata_type

        ! [Length,width] of domain in x/y units
        real(dp):: lw(2) 
        ! grid size [number of x-cells, number of y-cells]
        integer(ip):: nx(2) 
        ! cell size [dx,dy]
        real(dp):: dx(2)
        ! Absolute lower-left coordinate (bottom left corner of cell(1,1))
        real(dp):: lower_left(2)

        ! Parameter controlling extrapolation for finite volume methods
        real(dp) :: theta = extrapolation_theta
        ! CFL number
        real(dp):: cfl = cfl

        ! This information is useful in the context of nesting, where interior
        ! domains have cell sizes that are an integer divisor of the
        ! coarse-domain dx. Likewise one domain might take integer
        real(dp) :: dx_refinement_factor = -1.0_dp
        real(dp) :: timestepping_refinement_factor = 1.0_dp

        ! The width of the nesting layer for this domain -- which will depend
        ! on the dx value relative to the neighbouring domains, on the timestepping_method,
        ! and on the timestepping_refinement_factor. See 'get_domain_nesting_layer_thickness'
        integer(ip) :: nest_layer_width = -1_ip

        ! Useful to hold the interior bounding box [ = originally provided bounding box]
        ! since the actual bounding box might be changed to accommodate nesting 
        real(dp) :: interior_bounding_box(4,2) = 0.0_dp
        
        ! Domain ID, which is useful if multiple domains are running
        integer(4):: myid = 1

        ! Flag to denote boundaries at which nesting occurs: order is N, E, S, W.
        logical :: is_nesting_boundary(4) != .FALSE. ![.FALSE., .FALSE., .FALSE., .FALSE.]

        ! timestepping_method determines the choice of solver
        character(len=charlen):: timestepping_method = default_timestepping_method

        real(dp) :: max_parent_dx_ratio

    !end type

    !
    ! Main type used in the program. It holds the model arrays, information on the domain,
    ! etc
    !type, extends(domain_metadata_type) :: domain_type

        ! Number of quantities (stage, uh, vh, elevation)
        integer(ip):: nvar = 4 !global_nvar

        ! Name of ascii file where we output metadata
        character(charlen):: metadata_ascii_filename

        ! The domain 'interior' is surrounded by 'exterior' cells which are
        ! updated by boundary conditions, or copied from other domains. When
        ! tracking mass conservation, we only want to record inflows/outflows to
        ! interior cells. The interior cells are:
        !    [(1+exterior_cells_width):(domain%nx(1)-exterior_cells_width), &
        !     (1+exterior_cells_width):(domain%nx(2)-exterior_cells_width)]
        ! NOTE: exterior_cells_width should only be used to influence mass conservation tracking calculations.
        ! It is not strictly related to the halo width for parallel computations, although we may well want 
        ! that to be the same.
        integer(ip):: exterior_cells_width = 2
    
        ! Count number of time steps (useful if we want to do something only
        ! every n'th timestep)
        integer(ip):: nsteps_advanced = 0
        
        real(dp):: time = ZERO_dp
        ! dt computed from CFL condition (e.g. during flux computation)
        real(dp):: max_dt = ZERO_dp
        ! the time-step evolved by domain%evolve_one_step.
        real(dp):: evolve_step_dt = ZERO_dp
        ! the maximum allowed timestep [used to prevent high timesteps on dry domains]
        real(dp):: maximum_timestep = maximum_timestep
        
        ! Parameter which determines how the 'static depth' is computed in
        ! linear solver [corresponding to 'h0' in 'g h_0 d(free_surface)/dx']
        real(dp):: msl_linear = 0.0_dp

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
        character(len=charlen):: boundary_type = ''
        procedure(boundary_fun), pointer, nopass:: boundary_function => NULL()
        procedure(boundary_subroutine), pointer, nopass:: boundary_subroutine => default_boundary_subroutine
        !
        ! Flag whether the boundary is exterior (TRUE) or interior (FALSE). 
        ! Order is North (1), East (2), South (3), West (4) 
        ! The interior boundaries are 'overlaps between domains', and are dealt
        ! with by communication. 
        ! The numerical boundary conditions (e.g. reflective, etc) are applied
        ! to the exterior boundaries
        logical :: boundary_exterior(4) = [.TRUE., .TRUE., .TRUE., .TRUE.]

        ! 
        ! Mass conservation tracking 
        !
        ! Store the flux through the N, E, S, W boundaries 
        real(C_DOUBLE):: boundary_flux_store(4)  = ZERO_dp
        ! Time integrate the boundary fluxes
        real(C_DOUBLE):: boundary_flux_time_integral = ZERO_dp
        ! We need an intermediate variable to take care of time-stepping
        ! This integrates the boundary fluxes within the evolve step only
        real(C_DOUBLE):: boundary_flux_evolve_integral = ZERO_dp
        
        ! Lower/upper x and y indices over which the SWE computation takes place
        ! These might restrict the SWE update to a fraction of the domain (e.g.
        ! beginning of an earthquake-tsunami run where only a fraction of the domain is
        ! active.
        ! Currently implemented for linear SWE only
        integer(ip):: xL, xU, yL, yU

        ! Output folder units
        ! (FIXME: Consider making a 'output_data_type' which
        !    could hold the variables prefixed with 'output_'. Not essential but
        !    might reduce proliferation of variables in domain.)

        integer(ip), allocatable :: output_variable_unit_number(:)
        integer(ip):: output_time_unit_number
        character(len=charlen):: output_basedir = default_output_folder
        character(len=charlen):: output_folder_name = ''
        integer(ip):: logfile_unit = output_unit

        ! Type to record CPU timings
        type(timer_type):: timer

        ! Type to do single-grid coarray communication
        type(partitioned_domain_nesw_comms_type):: partitioned_comms

        ! Opportunity to do nesting 
        type(two_way_nesting_comms_type), allocatable :: nesting_boundaries(:)
 
        ! We don't have to store the max_U. For some problems (linear solver) that can
        ! take a significant fraction of the total time, or use too much memory.
        logical:: record_max_U = .true.
        integer(ip):: max_U_update_frequency = 1
        
        ! Spatial coordinates, dx/dy distances (useful for spherical coordinates)
        real(dp), allocatable :: x(:), y(:), distance_bottom_edge(:), distance_left_edge(:)
        real(dp), allocatable :: area_cell_y(:)
#ifdef SPHERICAL
        real(dp), allocatable :: coslat(:), coslat_bottom_edge(:)
#endif
#ifdef CORIOLIS
        real(dp), allocatable :: coriolis(:), coriolis_bottom_edge(:)
#endif

        ! Tide gauges
        type(point_gauge_type) :: point_gauges

        ! Big arrays
        !
        ! U holds quantities
        ! First 2 dimensions = space, 3rd = number of quantities
        real(dp), allocatable :: U(:,:,:) ! Needed always
        real(dp), allocatable :: max_U(:,:,:) ! Needed if max_U is to be output, but can reduce precision
        !
        ! Multi-dimensional arrays below are only required for nonlinear. Would be possible
        ! to further reduce memory usage, but not completely trivial
        real(dp), allocatable :: flux_NS(:,:,:) ! Could avoid
        real(dp), allocatable :: flux_EW(:,:,:) ! Could avoid
        real(dp), allocatable :: depth(:,:) ! Could avoid
        real(dp), allocatable :: explicit_source(:,:,:) ! Could avoid
        real(dp), allocatable :: explicit_source_VH_j_minus_1(:,:) ! Separate from explicit_source for OPENMP parallel logic
        real(dp), allocatable :: velocity(:,:, :) ! Could avoid
        real(dp), allocatable :: manning_squared(:,:) ! Needed for variable manning
        real(dp), allocatable :: backup_U(:,:,:) ! Needed

        ! nesting
        integer(ip), allocatable :: priority_domain_index(:,:)
        integer(ip), allocatable :: priority_domain_image(:,:)

        CONTAINS

        ! Initialisation
        procedure:: allocate_quantities => allocate_quantities

        ! Reporting
        procedure:: print => print_domain_statistics

        ! Core routines that occur within a timestep
        ! (consider making these not type bound -- since the user should not
        !  really call them)
        procedure:: get_bottom_edge_values => get_bottom_edge_values
        procedure:: get_left_edge_values => get_left_edge_values
        procedure:: compute_fluxes => compute_fluxes
        procedure:: update_U => update_U 
        procedure:: backup_quantities => backup_quantities

        ! Timestepping
        ! (consider making only 'evolve_one_step' type bound -- since the user should not
        !  really call others)
        procedure:: one_euler_step => one_euler_step
        procedure:: one_rk2_step => one_rk2_step 
        procedure:: one_rk2n_step => one_rk2n_step 
        procedure:: one_midpoint_step => one_midpoint_step
        procedure:: one_linear_leapfrog_step => one_linear_leapfrog_step
        procedure:: evolve_one_step => evolve_one_step

        ! boundary conditions. This just calls whatever domain%boundary_subroutine points to
        ! (consider making not type bound)
        procedure:: update_boundary => update_boundary

        ! io
        procedure:: create_output_files => create_output_files
        procedure:: write_to_output_files => write_to_output_files
        procedure:: update_max_quantities => update_max_quantities
        procedure:: write_max_quantities => write_max_quantities
        procedure:: log_outputs => divert_logfile_unit_to_file

        ! mass conservation tracking
        procedure:: mass_balance_interior => mass_balance_interior
        procedure:: volume_interior => volume_interior

        ! time-step for linear solver (a constant)
        procedure:: linear_timestep_max => linear_timestep_max

        ! gauges
        procedure:: setup_point_gauges => setup_point_gauges
        procedure:: write_gauge_time_series => write_gauge_time_series

        ! finalization (e.g. close netcdf files)
        procedure:: finalise => finalise_domain

    end type domain_type

    interface

        ! This gives the interface for a generic 'boundary function' which
        ! returns an array of 4 output values (stage/uh/vh/elevation)
        !
        ! It is used in conjunction with a number of different types of
        ! boundary conditions below
        !
        function boundary_fun(domain, t, x, y) RESULT(stage_uh_vh_elev)
            import dp, domain_type
            implicit none
            type(domain_type), intent(in):: domain 
            real(dp), intent(in):: t, x, y
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
    subroutine print_domain_statistics(domain)
        class(domain_type), intent(inout):: domain
        real(dp):: maxstage, minstage
        integer:: i,j, ecw
        real(dp):: dry_depth_threshold, energy_total, energy_potential, energy_kinetic
        real(dp):: depth, depth_iplus, depth_jplus
        logical, parameter:: report_energy_statistics=.FALSE.

        TIMER_START('printing_stats')

        dry_depth_threshold = minimum_allowed_depth
            
        ! Min/max stage in wet areas
        maxstage = -huge(1.0_dp)
        minstage = huge(1.0_dp)
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
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

        ! Optionally report integrated energy
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
            else
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
                            energy_kinetic = energy_kinetic + domain%area_cell_y(j) * &
                                (domain%U(i,j,UH)**2 + domain%U(i,j,VH)**2)/depth
                        end if 
                    end do
                end do

            end if

            ! Rescale energy statistics appropriately 
            energy_potential = energy_potential * gravity * HALF_dp
            energy_kinetic = energy_kinetic * HALF_dp
            energy_total = energy_potential + energy_kinetic

            write(domain%logfile_unit, *) 'Energy Total / rho: ', energy_total
            write(domain%logfile_unit, *) 'Energy Potential / rho: ', energy_potential
            write(domain%logfile_unit, *) 'Energy Kinetic / rho: ', energy_kinetic

        end if

        ! Flux of mass through model boundaries
        !write(domain%logfile_unit, *) 'Boundary flux time integral: ', domain%boundary_flux_time_integral
        ! Mass conservation check
        write(domain%logfile_unit, *) 'Mass Balance (domain interior): ', domain%mass_balance_interior()

        TIMER_STOP('printing_stats')
    end subroutine

    ! 
    ! Set up the full domain, allocate arrays, etc
    !
    ! @param global_lw real(dp) array size 2. length/width of domain in same
    !  units as x,y coordinates
    ! @param global_nx integer array size 2. number of x/y cells in the domain
    ! @param global_ll real(dp) array size 2. lower left x/y coordinate of the
    !  domain (at the corner of the lower left cell)
    ! @param create_output_files optional. If .TRUE. or not provided, then make
    ! @param co_size_xy Split up domain into sub-tiles of this dimension, using coarrays
    ! @param ew_periodic Use EW periodic boundaries [coarray only]
    ! @param ns_periodic Use NS periodic boundaries [coarray only]
    !  output files.
    !
    subroutine allocate_quantities(domain, global_lw, global_nx, global_ll, create_output_files,&
        co_size_xy, ew_periodic, ns_periodic)

        class(domain_type), target, intent(inout):: domain
        real(dp), intent(in):: global_lw(2), global_ll(2)
        integer(ip), intent(in):: global_nx(2)
        logical, optional, intent(in) :: create_output_files
        integer(ip), optional, intent(in):: co_size_xy(2)
        logical, optional, intent(in) :: ew_periodic, ns_periodic

        integer(ip), pointer:: nx, ny, nvar
        integer(ip) :: i
        logical :: create_output, use_partitioned_comms, ew_periodic_, ns_periodic_
        real(dp):: local_lw(2), local_ll(2)
        integer(ip):: local_nx(2)

        if(present(create_output_files)) then
            create_output = create_output_files
        else
            create_output = .TRUE.
        end if

        ! Only split the domain if we provided a co_size with > 1 image
        if (present(co_size_xy)) then
            if(maxval(co_size_xy) > 1) then
                use_partitioned_comms = .TRUE.

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

        ! Use coarrays if co_size_xy is provided
        if(use_partitioned_comms) then

            ! Compute the ll/lw/nx for this sub-domain
            call domain%partitioned_comms%initialise(&
                co_size_xy, global_ll, global_lw, global_nx, &
                local_ll, local_lw, local_nx, &
                ew_periodic=ew_periodic_, ns_periodic=ns_periodic_)
            domain%lower_left = local_ll
            domain%lw = local_lw 
            domain%nx = local_nx

            call domain%partitioned_comms%print()
        
            ! Make sure that communication boundaries are not
            ! numerical boundaries
            do i = 1,4
                if(domain%partitioned_comms%neighbour_images(i) > 0) domain%boundary_exterior(i) = .FALSE.
            end do

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

        nx => domain%nx(1)
        ny => domain%nx(2)
        nvar => domain%nvar 

        ! x/y coordinates (only stored along domain edge, assumed constant)
        allocate(domain%x(nx), domain%y(ny))
        do i = 1, nx
            domain%x(i) = domain%lower_left(1) + (i - HALF_dp)/(nx*ONE_dp)*domain%lw(1)
        end do
        do i = 1, ny
            domain%y(i) = domain%lower_left(2) + (i - HALF_dp)/(ny*ONE_dp)* domain%lw(2) 
        end do
     
#ifdef SPHERICAL
        ! For spherical coordinates it saves computation to have cos(latitude)
        ! at cells and edges
        allocate(domain%coslat(ny), domain%coslat_bottom_edge(ny+1))
        domain%coslat = cos(domain%y * deg2rad)
        domain%coslat_bottom_edge(1:ny) = cos((domain%y - HALF_dp * domain%dx(2))*deg2rad)
        domain%coslat_bottom_edge(ny+1) = cos((domain%y(ny) + HALF_dp*domain%dx(2))*deg2rad)

#endif

#ifdef CORIOLIS
        ! Coriolis parameter
        allocate(domain%coriolis_bottom_edge(ny+1), domain%coriolis(ny))
        domain%coriolis = 2.0_dp * sin(domain%y * deg2rad) * earth_angular_freq
        domain%coriolis_bottom_edge(1:ny) = 2.0_dp * earth_angular_freq * &
            sin((domain%y - HALF_dp * domain%dx(2))*deg2rad)
        domain%coriolis_bottom_edge(ny+1) = 2.0_dp * earth_angular_freq * &
            sin((domain%y(ny) + HALF_dp * domain%dx(2))*deg2rad)

#ifndef SPHERICAL
        stop('Cannot define preprocessing flag CORIOLIS without also defining SPHERICAL')
#endif    
#endif


        ! Distances along edges 
        ! For spherical coordinates, distance_xedge changes with y 
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

        allocate(domain%U(nx, ny, 1:nvar))
        domain%U = ZERO_dp
       
        if(domain%record_max_U) then 
            allocate(domain%max_U(nx, ny, 1))
            domain%max_U = -huge(1.0_dp)
        end if

        ! Many other variables are required for non-linear FV, but not for
        ! linear leap-frog
        if(domain%timestepping_method /= 'linear') then 

            allocate(domain%velocity(nx, ny, UH:VH))
            domain%velocity = ZERO_dp
            allocate(domain%depth(nx, ny))
            domain%depth = ZERO_dp
            
            ! NOTE: We don't need a flux for elevation 
            allocate(domain%flux_NS(nx, ny+1_ip, 1:3))
            domain%flux_NS = ZERO_dp
            allocate(domain%flux_EW(nx + 1_ip, ny, 1:3))
            domain%flux_EW = ZERO_dp


            ! Pressure gradient applies only to uh and vh (2 and 3 variables in U)
            allocate(domain%explicit_source(nx, ny, UH:VH))
            domain%explicit_source = ZERO_dp
            
            allocate(domain%explicit_source_VH_j_minus_1(nx, ny+1))
            domain%explicit_source_VH_j_minus_1 = ZERO_dp

            allocate(domain%manning_squared(nx, ny))
            domain%manning_squared = ZERO_dp

            ! If we use complex timestepping we need to back-up U for variables 1-3
            if (domain%timestepping_method /= 'euler') then
                allocate(domain%backup_U(nx, ny, 1:3))
                domain%backup_U = ZERO_dp
            end if

        endif

        
        if(create_output) then
            CALL domain%create_output_files()
        end if

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
    elemental function extrapolate_edge_second_order(U_local, U_lower, U_upper, theta, &
        extrapolation_sign) result(edge_value)
        real(dp), intent(in):: U_local, U_lower, U_upper, theta 
        integer(ip), intent(in):: extrapolation_sign

        real(dp) :: edge_value
        character(len=charlen), parameter :: gradient_type = 'standard'
        real(dp), parameter:: LIMITER_PAR_dp = 1.9_dp
        real(dp):: c, d

        !! This is much faster
        !FIXME: deliberate bug
        !edge_value = U_local
        !return

        !SELECT CASE (gradient_type)

        ! The following repeated 'if' statements seem to be optimized away when
        ! gradient_type is a parameter.
        
        if(gradient_type == 'min_limit') then
            ! Standard min limiter, with theta [0,1] optionall pushing closer to first order (if < 1)
            edge_value = U_local + extrapolation_sign * theta * HALF_dp * &
                minmod(U_upper - U_local, U_local - U_lower)
        end if

        if(gradient_type == 'standard') then 

            d = minmod(U_upper - U_local, U_local - U_lower) * LIMITER_PAR_dp * theta
            c = merge(ZERO_dp, HALF_dp * (U_upper - U_lower), d == ZERO_dp) 

            ! NOTE: IF d /= 0, then clearly d, c have the same sign
            ! We exploit this to avoid a further minmod call (which seems
            ! expensive)

            edge_value = U_local + HALF_dp * extrapolation_sign * &
                merge(min(c, d), max(c, d), d > ZERO_dp)
        end if

        if(gradient_type == 'debug') then
            edge_value = U_local 
        end if

    end function
    
    !
    ! Extrapolation to cell edges on a regular mesh without any limiting
    !
    ! edge = U_local + 0.25*(U_upper - U_lower) * theta * extrapolation_sign
    ! Typically theta = 1.0, but values closer to 0 can be used for first order in space extrapolation
    !
    elemental function extrapolate_edge_second_order_nolimit(U_local, U_lower, U_upper, theta, &
        extrapolation_sign) result(edge_value)
        real(dp), intent(in):: U_local, U_lower, U_upper, theta 
        integer(ip), intent(in):: extrapolation_sign
        real(dp) :: edge_value

        real(dp), parameter :: QUARTER_dp = HALF_dp * HALF_dp
        
        ! Standard min limiter
        edge_value = U_local + theta * extrapolation_sign * QUARTER_dp * (U_upper - U_lower)

    end function

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
    pure subroutine get_bottom_edge_values(domain, j, nx, ny, &
        theta_wd_neg_B, theta_wd_pos_B, &
        stage_neg_B, stage_pos_B, &
        depth_neg_B, depth_pos_B, &
        u_neg_B, u_pos_B, &
        v_neg_B, v_pos_B)

        class(domain_type), intent(in):: domain
        integer(ip), intent(in) :: j, ny, nx
        real(dp), intent(out):: theta_wd_neg_B(nx), theta_wd_pos_B(nx)
        real(dp), intent(out):: depth_neg_B(nx), depth_pos_B(nx)
        real(dp), intent(out):: stage_neg_B(nx), stage_pos_B(nx)
        real(dp), intent(out):: u_neg_B(nx), u_pos_B(nx)
        real(dp), intent(out):: v_neg_B(nx), v_pos_B(nx)

        ! Bottom edge, negative side
        if(j > 2) then
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

        else
            ! j == 2, cannot extrapolate so just use the lower neighbour value (first order accurate)
            theta_wd_neg_B = ZERO_dp
            stage_neg_B = domain%U(:, j-1, STG)
            depth_neg_B = domain%depth(:, j-1)
            u_neg_B = domain%velocity(:,j-1, UH)
            v_neg_B = domain%velocity(:,j-1, VH)
        end if

        if (j < ny) then
            ! Extrapolate using points j-1, j, j+1
            theta_wd_pos_B = 5.0_dp*( &
                (min(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1)) &
                    - minimum_allowed_depth) / &
                (max(domain%depth(:,j), domain%depth(:,j-1), domain%depth(:,j+1)) &
                    + 1000.0_dp*minimum_allowed_depth) &
                - 0.1_dp)
            theta_wd_pos_B = max(min(domain%theta, theta_wd_pos_B), ZERO_dp)
            

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

        else
            ! j == ny, cannot extrapolate, so use the j value (first order accurate)
            theta_wd_pos_B = ZERO_dp
            stage_pos_B = domain%U(:, j, STG)
            depth_pos_B = domain%depth(:, j)
            u_pos_B = domain%velocity(:, j, UH)
            v_pos_B = domain%velocity(:, j, VH)
        end if

    end subroutine

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
    pure subroutine get_left_edge_values(domain, j, nx, &
        theta_wd_neg_L, theta_wd_pos_L, &
        stage_neg_L, stage_pos_L, &
        depth_neg_L, depth_pos_L, &
        u_neg_L, u_pos_L, &
        v_neg_L, v_pos_L)

        class(domain_type), intent(in):: domain
        integer(ip), intent(in) :: j, nx
        real(dp), intent(out):: theta_wd_neg_L(nx), theta_wd_pos_L(nx)
        real(dp), intent(out):: depth_neg_L(nx), depth_pos_L(nx)
        real(dp), intent(out):: stage_neg_L(nx), stage_pos_L(nx)
        real(dp), intent(out):: u_neg_L(nx), u_pos_L(nx)
        real(dp), intent(out):: v_neg_L(nx), v_pos_L(nx)

        integer(ip) :: i

        ! Note: In practice the value in index (1) is not subsequently used for any variables

        ! View from negative side of edge
        theta_wd_neg_L(1) = ZERO_dp
        theta_wd_neg_L(2:(nx-1)) = 5.0_dp * ( &
            (min(domain%depth(2:(nx-1), j), domain%depth(1:(nx-2), j), domain%depth(3:nx, j)) &
                - minimum_allowed_depth) / &
            (max(domain%depth(2:(nx-1), j), domain%depth(1:(nx-2), j), domain%depth(3:nx, j)) &
                + 1000.0_dp * minimum_allowed_depth) &
            - 0.1_dp)
        theta_wd_neg_L(nx) = ZERO_dp
        theta_wd_neg_L = max(min(domain%theta, theta_wd_neg_L), ZERO_dp)

        !theta_wd_neg_L(1) = ZERO_dp
        !do i = 2, (nx-1)
        !    theta_wd_neg_L(i) = 5.0_dp * ( &
        !        (min(domain%depth(i,j), domain%depth(i-1,j), domain%depth(i+1,j)) &
        !            - minimum_allowed_depth)/&
        !        (max(domain%depth(i,j), domain%depth(i-1,j), domain%depth(i+1,j)) &
        !            + 1000.0_dp * minimum_allowed_depth) &
        !        -0.1_dp)
        !    theta_wd_neg_L(i) = max(min(domain%theta, theta_wd_neg_L(i)), ZERO_dp)
        !end do
        !theta_wd_neg_L(nx) = ZERO_dp


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
        
    end subroutine

    !
    ! Compute the fluxes, and other things, in preparation for an update of U 
    !
    ! Update values of:
    ! domain%flux_NS, domain%flux_EW, domain%max_dt, domain%explicit_source, 
    ! domain%explicit_source_VH_j_minus_1, domain%boundary_flux_store
    !
    ! @param domain the model domain type for which fluxes etc will be computed
    !
    subroutine compute_fluxes(domain)
        ! Compute fluxes for 2D shallow water equations on structured grid
        ! Use an Audusse type method, like ANUGA, but structured

        class(domain_type), intent(inout):: domain

        integer(ip):: i, j, nx, ny
        ! wavespeeds
        real(dp):: s_max, s_min, gs_pos, gs_neg, sminsmax
        ! stage/depth at + and - side of edge
        real(dp):: stage_pos, stage_neg, stage_pos_star, stage_neg_star, &
                   depth_pos, depth_pos_star, depth_pos_c, depth_pos_inv,&
                   depth_neg, depth_neg_star, depth_neg_c, depth_neg_inv,&
                   z_half, z_neg, z_pos
        ! velocities and momenta at + and - sides of edge
        real(dp):: u_pos, v_pos, u_neg, v_neg, ud_neg, vd_neg, ud_pos, vd_pos, vel_beta_neg, vel_beta_pos
        ! convenience variables
        real(dp):: denom, inv_denom, max_speed, max_dt, min_dt_inv, dx_cfl_inv(2)
        real(dp):: pressure_s, pressure_n, pressure_e, pressure_w
        real(dp):: depth_local(3), depth_inv
        real(dp):: half_cfl, max_dt_inv, source_tmp
        character(len=charlen):: timer_name
        real(dp), parameter :: EPS = 1.0e-06_dp, diffusion_scale = 1.0_dp
        real(dp), parameter:: half_gravity = HALF_dp * gravity
        integer(ip):: masscon_error, n_ext
        real(dp):: masscon_error_neg_depth
        
        ! Bottom edge values on 'positive' and 'negative' side (i.e. viewed from j and j-1 respectively)
        ! theta_wd controls the limiting
        real(dp) :: theta_wd_neg_B(domain%nx(1)), theta_wd_pos_B(domain%nx(1)), &
            stage_neg_B(domain%nx(1)), stage_pos_B(domain%nx(1)), &
            depth_neg_B(domain%nx(1)), depth_pos_B(domain%nx(1)), &
            u_neg_B(domain%nx(1)), u_pos_B(domain%nx(1)), &
            v_neg_B(domain%nx(1)), v_pos_B(domain%nx(1))
       
        ! Left edge values on 'positive' and 'negative' side (i.e. viewed from i and i-1 respectively) 
        ! theta_wd controls the limiting
        real(dp) :: theta_wd_neg_L(domain%nx(1)), theta_wd_pos_L(domain%nx(1)), &
            stage_neg_L(domain%nx(1)), stage_pos_L(domain%nx(1)), &
            depth_neg_L(domain%nx(1)), depth_pos_L(domain%nx(1)), &
            u_neg_L(domain%nx(1)), u_pos_L(domain%nx(1)), &
            v_neg_L(domain%nx(1)), v_pos_L(domain%nx(1))

        real(dp) :: bed_j_minus_1(domain%nx(1))

        half_cfl = HALF_dp * domain%cfl
        ny = domain%nx(2)
        nx = domain%nx(1)


        ! Set dt to a high value (it will drop)
        max_dt = domain%maximum_timestep
        ! By computing the inverse we avoid division in the loop
        max_dt_inv = ONE_dp/max_dt

        ! Zero the explicit source, recompute depth and velocity
        ! Must be updated before the main loop when done in parallel
        masscon_error = 0_ip
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, half_cfl, nx, ny) REDUCTION(MAX:max_dt_inv, masscon_error)
        !$OMP DO SCHEDULE(GUIDED)
        do j = 1, domain%nx(2) !ny 
            domain%explicit_source(:,j,:) = ZERO_dp
            do i = 1, domain%nx(1) !nx
                domain%depth(i,j) = domain%U(i,j,STG) - domain%U(i,j,ELV)
                if(domain%depth(i,j) > minimum_allowed_depth) then
                    depth_inv = ONE_dp/domain%depth(i,j)
                    domain%velocity(i,j,UH) = domain%U(i,j,UH) * depth_inv
                    domain%velocity(i,j,VH) = domain%U(i,j,VH) * depth_inv
                else
                    if(domain%depth(i,j) < ZERO_dp) then
                        masscon_error = 1_ip
                    end if
                    domain%velocity(i,j,UH) = ZERO_dp
                    domain%velocity(i,j,VH) = ZERO_dp
                end if
            end do
        end do
        !$OMP END DO
    
        if(masscon_error > 0_ip) then
            masscon_error_neg_depth = minval(domain%depth)
            print*, 'stage < bed --> mass conservation error'
            print*, masscon_error_neg_depth, domain%nsteps_advanced
            do j = 1, domain%nx(2)
                do i = 1, domain%nx(1)
                    if(domain%depth(i,j) == masscon_error_neg_depth) then
                        print*, i, j, domain%x(i) + domain%lower_left(1), &
                            domain%y(j) + domain%lower_left(2)
                    end if
                end do 
            end do
    
            call generic_stop()
        end if

        !Main loop

        !$OMP DO SCHEDULE(GUIDED)
        do j = 2, domain%nx(2) !ny
            ! Assume distance_left_edge is constant but distance_bottom_edge
            ! might change with y (spherical coordinates). For cartesian coordinates
            ! we could move this outside the loop.
            dx_cfl_inv(1) = ONE_dp/(domain%distance_bottom_edge(j) * half_cfl)
            dx_cfl_inv(2) = ONE_dp/(domain%distance_left_edge(1) * half_cfl)

            ! Get bed at j-1 (might improve cache access) 
            bed_j_minus_1 = domain%U(:,j-1,ELV)

            ! Get bottom edge values
            call domain%get_bottom_edge_values(j, nx, ny, &
                theta_wd_neg_B, theta_wd_pos_B, &
                stage_neg_B, stage_pos_B, &
                depth_neg_B, depth_pos_B, &
                u_neg_B, u_pos_B, &
                v_neg_B, v_pos_B)

            ! Get left edge values
            call domain%get_left_edge_values(j, nx, &
                theta_wd_neg_L, theta_wd_pos_L, &
                stage_neg_L, stage_pos_L, &
                depth_neg_L, depth_pos_L, &
                u_neg_L, u_pos_L, &
                v_neg_L, v_pos_L)

            ! Could potentially vectorize this loop
            !!$OMP SIMD 
            !DO i = 2, domain%nx(1) !nx
            do concurrent (i = 2:nx)
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
                if(depth_neg_star == ZERO_dp) then
                    v_neg = ZERO_dp
                    ud_neg = ZERO_dp
                    vd_neg = ZERO_dp

                    gs_neg = ZERO_dp
                else
                    ! Audusse type depth_integrated_velocity correction
                    vd_neg = v_neg * depth_neg_star
                    ud_neg = u_neg * depth_neg_star

                    ! Gravity wave celerity
                    gs_neg = sqrt(gravity * depth_neg_star)
                end if
                
                ! Velocity (in NS direction)
                if(depth_pos_star == ZERO_dp) then
                    v_pos = ZERO_dp
                    ud_pos = ZERO_dp
                    vd_pos = ZERO_dp

                    gs_pos = ZERO_dp
                else
                    ! Correct depth-integrated_velocity (Audusse type approach)
                    vd_pos = v_pos * depth_pos_star
                    ud_pos = u_pos * depth_pos_star

                    ! Gravity wave celerity
                    gs_pos = sqrt(gravity * depth_pos_star)
                end if

                ! Wave-celerities
                s_max = max(max(v_neg + gs_neg, v_pos + gs_pos), ZERO_dp)
                s_min = min(min(v_neg - gs_neg, v_pos - gs_pos), ZERO_dp)

                denom = s_max - s_min

                if (denom > EPS) then
                    inv_denom = domain%distance_bottom_edge(j) / denom 
                    sminsmax = s_min * s_max * diffusion_scale
                    vel_beta_neg = v_neg * advection_beta
                    vel_beta_pos = v_pos * advection_beta

                    domain%flux_NS(i, j, STG) = &
                        (s_max * vd_neg - &
                         s_min * vd_pos + &
                        sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
                    domain%flux_NS(i, j, UH) = &
                        (s_max * ud_neg * vel_beta_neg - &
                         s_min * ud_pos * vel_beta_pos + &
                        sminsmax *(ud_pos - ud_neg)) * inv_denom
                    domain%flux_NS(i, j, VH) = &
                        (s_max * (vd_neg * vel_beta_neg + half_gravity * depth_neg_star*depth_neg_star) - & 
                         s_min * (vd_pos * vel_beta_pos + half_gravity * depth_pos_star*depth_pos_star) + &
                        sminsmax *(vd_pos - vd_neg)) * inv_denom
                else
                    domain%flux_NS(i,j,STG) = ZERO_dp
                    domain%flux_NS(i,j,UH) = ZERO_dp
                    domain%flux_NS(i,j,VH) = ZERO_dp
                end if

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
                if (max_speed > EPS) then
                    max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(2))
                end if

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
                if(depth_neg_star == ZERO_dp) then
                    u_neg = ZERO_dp
                    ud_neg = ZERO_dp
                    vd_neg = ZERO_dp

                    gs_neg = ZERO_dp 
                else
                    ! Audusse type depth-integrated-velocity correction
                    ud_neg = u_neg * depth_neg_star
                    vd_neg = v_neg * depth_neg_star

                    gs_neg = sqrt(gravity * depth_neg_star)
                end if

                ! Velocity (in NS direction)
                if(depth_pos_star == ZERO_dp) then
                    u_pos = ZERO_dp
                    ud_pos = ZERO_dp
                    vd_pos = ZERO_dp

                    gs_pos = ZERO_dp
                else
                    ud_pos = u_pos * depth_pos_star
                    vd_pos = v_pos * depth_pos_star

                    gs_pos = sqrt(gravity * depth_pos_star)
                end if

                ! Wave-celerities
                s_max = max(max(u_neg + gs_neg, u_pos + gs_pos), ZERO_dp)
                s_min = min(min(u_neg - gs_neg, u_pos - gs_pos), ZERO_dp)

                denom = s_max - s_min

                if (denom > EPS) then
                    inv_denom = domain%distance_left_edge(i) / denom 
                    sminsmax = s_min * s_max * diffusion_scale
                    vel_beta_neg = u_neg * advection_beta
                    vel_beta_pos = u_pos * advection_beta

                    domain%flux_EW(i,j,STG) = &
                        (s_max * ud_neg - &
                         s_min * ud_pos + &
                        sminsmax * (stage_pos_star - stage_neg_star)) * inv_denom
                    domain%flux_EW(i,j,UH) = &
                        (s_max * (ud_neg * vel_beta_neg + half_gravity * depth_neg_star*depth_neg_star) - &
                         s_min * (ud_pos * vel_beta_pos + half_gravity * depth_pos_star*depth_pos_star) + &
                         sminsmax * (ud_pos - ud_neg)) * inv_denom
                    domain%flux_EW(i,j,VH) = &
                        (s_max * vd_neg * vel_beta_neg - &
                         s_min * vd_pos * vel_beta_pos + &
                        sminsmax * (vd_pos - vd_neg)) * inv_denom
                else
                    domain%flux_EW(i,j,STG) = ZERO_dp
                    domain%flux_EW(i,j,UH) = ZERO_dp
                    domain%flux_EW(i,j,VH) = ZERO_dp
                end if

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

                !! FIXME: Consider adding coriolis

                !! Source term associated with spherical coordinates for the pressure gradient term
                !!
                domain%explicit_source(i, j, VH) = domain%explicit_source(i, j, VH) + &
                    half_gravity * domain%depth(i, j) * domain%depth(i, j) * &
                    (domain%distance_bottom_edge(j+1) - domain%distance_bottom_edge(j))

#endif

                ! Timestep
                max_speed = max(s_max, -s_min)
                if (max_speed > EPS) then
                    max_dt_inv = max(max_dt_inv, max_speed * dx_cfl_inv(1))
                end if

            end do
        end do
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
        n_ext = domain%exterior_cells_width
        ! Outward boundary flux over the north
        domain%boundary_flux_store(1) = sum(domain%flux_NS( (1+n_ext):(nx-n_ext), ny+1-n_ext, STG))
        ! Outward boundary flux over the east
        domain%boundary_flux_store(2) = sum(domain%flux_EW( nx+1-n_ext, (1+n_ext):(ny-n_ext), STG))
        ! Outward boundary flux over the south
        domain%boundary_flux_store(3) = -sum(domain%flux_NS( (1+n_ext):(nx-n_ext), 1+n_ext, STG))
        ! Outward boundary flux over the west
        domain%boundary_flux_store(4) = -sum(domain%flux_EW( 1+n_ext, (1+n_ext):(ny-n_ext), STG))

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

        integer(ip) :: nx

        real(dp):: inv_cell_area_dt, depth, implicit_factor, dt_gravity, fs
        integer(ip):: j, i, k, kk


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, dt)
        dt_gravity = dt * gravity
        !$OMP DO SCHEDULE(GUIDED)
        do j = 1, domain%nx(2)
        !do concurrent (j = 1:domain%nx(2))
            ! For spherical coordiantes, cell area changes with y.
            ! For cartesian coordinates this could be moved out of the loop
            inv_cell_area_dt = dt / domain%area_cell_y(j)

            !do i = 1, nx 
            do concurrent (i = 1:domain%nx(1))
                depth = domain%depth(i,j)
        
                !! Fluxes
                do kk = 1, 3
                    domain%U(i,j,kk) = domain%U(i,j,kk) - inv_cell_area_dt * ( & 
                        (domain%flux_NS(i, j+1, kk) - domain%flux_NS(i, j, kk)) + &
                        (domain%flux_EW(i+1, j, kk) - domain%flux_EW(i, j, kk) ))
                end do

                ! Velocity clipping
                if (depth <= minimum_allowed_depth) then
                    domain%U(i,j,UH) = ZERO_dp 
                    domain%U(i,j,VH) = ZERO_dp 
                else
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
                end if
            end do
        end do
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
    subroutine one_euler_step(domain, timestep)
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep

        character(len=charlen):: timer_name
        real(dp):: ts

        call domain%update_boundary()

        TIMER_START('flux')
        call domain%compute_fluxes()
        TIMER_STOP('flux')

        TIMER_START('update')
        if(present(timestep)) then
            ts = timestep
        else
            ts = domain%max_dt
        end if

        call domain%update_U(ts)

        ! Track flux through boundaries
        domain%boundary_flux_evolve_integral = domain%boundary_flux_evolve_integral + &
            ts * sum(domain%boundary_flux_store)

        TIMER_STOP('update')

        ! Coarray communication, if required
        TIMER_START('partitioned_comms')
        call domain%partitioned_comms%communicate(domain%U)
        TIMER_STOP('partitioned_comms')


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

        real(dp):: backup_time, dt_first_step
        integer(ip):: j
        character(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
        
        ! First euler step
        if(present(timestep)) then
            call domain%one_euler_step(timestep)
            dt_first_step = timestep
        else
            call domain%one_euler_step()
            dt_first_step = domain%max_dt
        end if
       
        ! Second euler step with the same timestep 
        call domain%one_euler_step(dt_first_step)

        TIMER_START('average')

        ! Take average (but allow for openmp)
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            domain%U(:, j, [STG, UH, VH]) = HALF_dp * (domain%U(:, j, [STG, UH, VH]) +&
                domain%backup_U(:, j, [STG, UH, VH]))
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time (since we updated twice) and boundary flux integral
        domain%time = backup_time + dt_first_step
        domain%boundary_flux_evolve_integral = HALF_dp * domain%boundary_flux_evolve_integral
        
        TIMER_STOP('average')

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
    ! @param domain the domain to be updated
    ! @param timestep the timestep, by which the solution is advanced (n-1) times
    !
    subroutine one_rk2n_step(domain, timestep)
        ! Advance (n-1) * timesteps in this routine, with 2nd order in time accuracy
        class(domain_type), intent(inout):: domain
        real(dp), optional, intent(in):: timestep

        integer(ip), parameter:: n = 5 ! number of substeps to take, must be > 2
        real(dp), parameter:: n_inverse = ONE_dp / (ONE_dp * n)
        real(dp):: backup_time, dt_first_step, reduced_dt, backup_flux_integral
        integer(ip):: j,k
        character(len=charlen):: timer_name

        ! Backup quantities
        backup_time = domain%time
        call domain%backup_quantities()
     
        if(present(timestep)) then 
            ! first step 
            call domain%one_euler_step(timestep)
            ! store timestep
            dt_first_step = timestep
        else
            ! first step 
            call domain%one_euler_step()
            ! store timestep
            dt_first_step = domain%max_dt
        end if

        ! Steps 2, n-1
        do j = 2, n-1        
            call domain%one_euler_step(dt_first_step)
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

        ! Now take one step of duration (n-1)/n * dt
        reduced_dt = (ONE_dp * n - ONE_dp) * n_inverse * dt_first_step
        call domain%one_euler_step(reduced_dt)
        
        TIMER_START('final_update')

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

        ! Fix time and timestep (since we updated (n-1)*dt regular timesteps)
        domain%time = backup_time + dt_first_step * (n-1) 

        TIMER_STOP('final_update')

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
        
        call domain%update_boundary()

        ! First euler sub-step
        call domain%compute_fluxes()

        if(present(timestep)) then
            dt_first_step = timestep
            call domain%update_U(dt_first_step*HALF_dp)
        else
            dt_first_step = domain%max_dt
            call domain%update_U(dt_first_step*HALF_dp)
        end if

        TIMER_START('partitioned_comms')
        call domain%partitioned_comms%communicate(domain%U)
        TIMER_STOP('partitioned_comms')

        call domain%update_boundary()
        
        ! Compute fluxes 
        call domain%compute_fluxes()

        ! Set U back to backup_U
        !
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            domain%U(:, j, [STG,UH,VH]) = domain%backup_U(:, j, [STG,UH,VH])
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        ! Fix time
        domain%time = backup_time

        ! Update U
        call domain%update_U(dt_first_step)

        TIMER_START('partitioned_comms')
        call domain%partitioned_comms%communicate(domain%U)
        TIMER_STOP('partitioned_comms')

        domain%boundary_flux_evolve_integral = sum(domain%boundary_flux_store)*dt_first_step

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

        ! Do we represent pressure gradients with a 'truely' linear term g * depth0 * dStage/dx,
        ! or with a nonlinear term g * depth * dStage/dx (i.e. where the 'depth' varies)?
        logical, parameter:: truely_linear = .TRUE.

        ! The linear solver code has become complex [in order to reduce memory footprint, and
        ! include coriolis, while still having openmp work]. So it is moved here. 
#include "domain_linear_solver_include.inc"        

    end subroutine
    
    !
    ! Routine to run all boundary conditions
    !
    subroutine update_boundary(domain)
        class(domain_type), intent(inout):: domain 

        TIMER_START('boundary_update')

        if(associated(domain%boundary_subroutine)) then
            CALL domain%boundary_subroutine(domain)
        end if

        TIMER_STOP('boundary_update')

    end subroutine

    !
    ! Stub routine to allow the user to not provide a boundary condition
    !
    subroutine default_boundary_subroutine(domain)
        type(domain_type), intent(inout):: domain
        ! Do nothing (default case)
    end subroutine

    !
    ! Copy domain%U to domain%backup_U
    !
    subroutine backup_quantities(domain)
        
        class(domain_type), intent(inout):: domain
        integer(ip):: j, k

        TIMER_START('backup')

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)
        !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
        do k = 1, 3
            do j = 1, domain%nx(2)
                domain%backup_U(:, j, k) = domain%U(:, j, k)
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        TIMER_STOP('backup')

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


        time0 = domain%time

        ! Reset the boundary fluxes integrated over the evolve step to zero
        domain%boundary_flux_evolve_integral = ZERO_dp
   
        select case (domain%timestepping_method) 
        
        !case ('euler')
        !    if(present(timestep)) then
        !        call domain%one_euler_step(timestep)
        !    else
        !        call domain%one_euler_step()
        !    end if
        !case ('rk2')
        !    if(present(timestep)) then
        !        call domain%one_rk2_step(timestep)
        !    else
        !        call domain%one_rk2_step()
        !    end if
        !case('rk2n')
        !    if(present(timestep)) then
        !        call domain%one_rk2n_step(timestep)
        !    else
        !        call domain%one_rk2n_step()
        !    end if
        !case ('midpoint')
        !    if(present(timestep)) then
        !        call domain%one_midpoint_step(timestep)
        !    else
        !        call domain%one_midpoint_step()
        !    end if
        case ('linear')
            if(present(timestep)) then
                call domain%one_linear_leapfrog_step(timestep)
            else
                print*, 'ERROR: timestep must be provided for linear evolve_one_step'
                call generic_stop()
            end if
        case default
            ! Deliberately only support linear solver for now
            print*, 'ERROR: domain%timestepping_method is not == linear'
            call generic_stop()
        end select

        ! For some problems updating max U can take a significant fraction of the time,
        ! so we allow it to only be done occasionally
        if(mod(domain%nsteps_advanced, domain%max_U_update_frequency) == 0) then
            call domain%update_max_quantities()
        end if

        ! Record the timestep here
        domain%evolve_step_dt = domain%time - time0

        ! Update the boundary flux time integral
        domain%boundary_flux_time_integral = domain%boundary_flux_time_integral + &
            domain%boundary_flux_evolve_integral

        domain%nsteps_advanced = domain%nsteps_advanced + 1

    end subroutine

    !
    ! Convenience function to compute the volume of water in the 'interior' of the domain
    ! This involves all parts of the domain that are more than domain%exterior_cells_width
    ! from the edge.
    !
    function volume_interior(domain) result(domain_volume)
        class(domain_type), intent(in):: domain
        real(dp) :: domain_volume

        integer(ip) :: j, n_ext
        
        ! Volume on the interior. At the moment the interior is all
        ! but the outer cells of the domain, but that could change.

        !TIMER_START('volume_interior')
        n_ext = domain%exterior_cells_width
        domain_volume = ZERO_dp
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, n_ext) REDUCTION(+:domain_volume)
        !$OMP DO SCHEDULE(STATIC)
        do j = (1+n_ext), (domain%nx(2) - n_ext) 
            domain_volume = domain_volume + domain%area_cell_y(j) * &
                (sum(domain%U((1+n_ext):(domain%nx(1)-n_ext),j,STG) - &
                     domain%U((1+n_ext):(domain%nx(1)-n_ext),j,ELV)))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        !TIMER_STOP('volume_interior')

    end function

    !
    ! Compute the volume on interior cells and add to
    ! the integrated boundary flux. This should sum to a constant
    ! in the absence of mass sources.
    !
    ! Note this relies on the time-stepping routines correctly
    ! computing the boundary_flux_time_integral
    !
    function mass_balance_interior(domain) result(mass_balance)
        class(domain_type), intent(in):: domain
        real(dp) :: mass_balance

        mass_balance = domain%volume_interior() + domain%boundary_flux_time_integral 

    end function

    !
    ! Function to compute the max timestep allowed for the linear shallow water equations,
    ! using the provided CFL.
    !
    function linear_timestep_max(domain) result(timestep)
        class(domain_type), intent(in):: domain
        real(dp) :: timestep

        real(dp) :: ts_max
        integer(ip) :: i, j

        ts_max = huge(1.0_dp)
    
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                ! max timestep = Distance along latitude / wave-speed <= dt
                ! Beware -- might need to use a different CFL number?
                ts_max = min(ts_max, &
                    0.5_dp * min(&
                        (domain%distance_bottom_edge(j+1) + domain%distance_bottom_edge(j)),& 
                        (domain%distance_left_edge(i+1)   + domain%distance_left_edge(i)  ) ) / &
                    sqrt(gravity * max(-domain%U(i,j,ELV), minimum_allowed_depth)) )
            end do
        end do
        timestep = ts_max * domain%cfl

    end function

    !
    ! Initialise output files
    !
    subroutine create_output_files(domain)
        class(domain_type), intent(inout):: domain

        character(len=charlen):: mkdir_command, cp_command, t1, t2, t3, &
                                 output_folder_name
        integer(ip):: i, metadata_unit

        allocate(domain%output_variable_unit_number(domain%nvar))

        ! Create output directory
        call date_and_time(t1, t2, t3)
        ! Get domain id as a character
        write(t3, '(I6)') 100000 + domain%myid

        output_folder_name = trim(domain%output_basedir) // '/RUN_ID' // trim(t3) // &
            '_' // trim(t1) // '_' // trim(t2)
        mkdir_command = 'mkdir -p ' // trim(output_folder_name)
        !call execute_command_line(trim(mkdir_command))
        call system(trim(mkdir_command))


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

        do i=1, domain%nvar
            ! Make output_file_name
            write(t1,'(I1)') i
            write(t2, *) 'Var_', trim(t1), '_ID', trim(t3)
            t1 = trim(output_folder_name) // '/' // adjustl(trim(t2))

            ! Open the files

            ! Ascii
            !OPEN(newunit = domain%output_variable_unit_number(i), file = t1)

            ! Binary. Using 'stream' access makes it easy to read in R
            open(newunit = domain%output_variable_unit_number(i), file = t1, &
                access='stream', form='unformatted')
        end do

        ! Make a time file_name. Store as ascii
        t1 = trim(output_folder_name) // '/' // 'Time_ID' // trim(t3) // '.txt'
        open(newunit = domain%output_time_unit_number, file = t1)

        ! Copy code to the output directory

        mkdir_command = 'mkdir -p ' // trim(output_folder_name) // '/Code'
        cp_command = 'cp *.f* *.R *.c make* ' // trim(output_folder_name) // '/Code'
        !call execute_command_line(trim(mkdir_command))
        call system(trim(mkdir_command))
        !call execute_command_line(trim(cp_command))
        call system(trim(cp_command))

        domain%output_folder_name = output_folder_name

    end subroutine create_output_files

    !
    ! Write output to files (must call create_output_files first).
    ! If time_only=.TRUE., only write the time. This is used to avoid
    ! writing the main model grids. Typically useful when values at point-gauges
    ! are being recorded, and we want to store the time too, but it would
    ! take too much disk to store the model grids
    !
    subroutine write_to_output_files(domain, time_only)
        class(domain_type), intent(inout):: domain
        logical, optional, intent(in) :: time_only
        integer(ip):: i, j
        logical:: to

        if(present(time_only)) then
            to = time_only
        else
            to = .FALSE.
        end if

        TIMER_START('fileIO')
       
        if(to .eqv. .FALSE.) then 
            do i = 1, domain%nvar
                do j = 1, domain%nx(2)
                    ! Binary
                    write(domain%output_variable_unit_number(i)) real(domain%u(:,j,i), output_precision)
                end do
            end do
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

        write(domain_ID, '(I6)') 100000 + domain%myid

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
            TIMER_START('update_max_quantities')

            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain)

            !! Previously we stored all quantities in max_U. However, generally
            !! max uh and max vh don't have much meaning by themselves. So now
            !! we just store max stage
            !!$OMP DO SCHEDULE(GUIDED), COLLAPSE(2)
            !DO k = 1, domain%nvar

            !! UPDATE: Only record max stage
            !$OMP DO SCHEDULE(GUIDED)
                do j = domain%yL, domain%yU !1, domain%nx(2)
                    do i = domain%xL, domain%xU !1, domain%nx(1)
                        domain%max_U(i,j,1) = max(domain%max_U(i,j,1), domain%U(i,j,1))
                    end do
                end do
            !END DO
            !$OMP END DO
            !$OMP END PARALLEL
            
            TIMER_STOP('update_max_quantities')
        end if

    end subroutine update_max_quantities

    !
    ! Write max quantities to a file (usually just called once at the end of a
    ! simulation). 
    !
    ! Currently we only write stage, followed by the elevation. Although the
    ! latter usually doesn't evolve, we generally want both the max stage and the
    ! elevation for plotting purposes, so it is saved here too.
    !
    subroutine write_max_quantities(domain)
        class(domain_type), intent(in):: domain
        character(len=charlen):: max_quantities_filename
        integer:: i, j, k

        if(domain%record_max_U) then

            max_quantities_filename = trim(domain%output_folder_name) // '/Max_quantities'

            open(newunit=i, file=max_quantities_filename, access='stream', form='unformatted')
        
            !DO k = 1, domain%nvar
            !! Update: Only record max stage
            k = 1
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
 
    end subroutine

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
        integer(ip):: i

        if(present(time_series_var)) then
            allocate(tsv(size(time_series_var)))
            tsv = time_series_var
        else
            ! Default case -- store stage/uh/vh every output step
            allocate(tsv(3))
            tsv = [STG, UH, VH]
        end if

        if(present(static_var)) then
            allocate(sv(size(static_var)))
            sv = static_var
        else 
            ! Default case -- store elevation once
            allocate(sv(1))
            sv = [ELV]
        end if

        if(present(gauge_ids)) then
            if(size(gauge_ids) /= size(xy_coords(1,:))) then
                print*, 'Number of gauge ids does not equal number of coordinates' 
            end if
            allocate(gauge_ids_local(size(gauge_ids)))
            gauge_ids_local = gauge_ids
        else
            ! Default case -- give sequential integer ids
            allocate(gauge_ids_local(size(xy_coords(1,:))))
            do i = 1, size(gauge_ids_local)
                gauge_ids_local(i) = i*1.0_dp
            end do
        end if
    
        ! Get domain id as a character and put it in the output file name
        write(t3, '(I6)') 100000 + domain%myid
        netcdf_gauge_output_file = trim(domain%output_folder_name) // '/' // &
            'Gauges_data_ID' // trim(t3) // '.nc'

        ! Pass the bounding box
        bounding_box(1,1:2) = domain%lower_left
        bounding_box(2,1:2) = domain%lower_left + domain%lw

        ! Allocate the gauges
        call domain%point_gauges%allocate_gauges(xy_coords, tsv, sv, gauge_ids_local, &
            bounding_box=bounding_box)

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

        if(allocated(domain%point_gauges%time_series_values)) then
            call domain%point_gauges%write_current_time_series(domain%U, domain%time)
        end if

    end subroutine

    !
    ! Routine to call once we no longer need the domain. One case where this is
    ! important is when using netcdf output -- since if the files are not closed,
    ! then they may not be completely written out.
    !
    subroutine finalise_domain(domain)

        class(domain_type), intent(inout):: domain

        ! Close the gauges netcdf file -- since otherwise it might not finish
        ! writing.
        call domain%point_gauges%finalise()
    
        ! Flush all open file units
        call flush()

    end subroutine

end module domain_mod
