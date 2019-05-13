! Compile with -DTIMER to add timing to the code 
#ifdef TIMER
#   define TIMER_START(tname) call md%timer%timer_start(tname)
#   define TIMER_STOP(tname)  call md%timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

!
! Contains a multidomain_type, to hold multiple rectangular domains which
! communicate with each other. 
! - Finer resolution domains may overlap coarser domains. 
! - Domains with the same cell size cannot overlap each other, because where domains
!   overlap, we need to decide which one has 'priority', and this is currently
!   done by selecting the one with finest cell area. Overlaps would lead to an ambiguous choice,
!   and so are not currently supported.
! - Internally, the code will extend each domain with 'communication buffers',
!   where they take values copied from neighbouring 'priority' domains. Flux correction
!   is also applied (by default, see a flag in global_mod) so that the multidomain
!   is mass conservative. 
! - When determining the size of the communication buffers, each domain figures out
!   which priority domains are 'just outside' its own initially provided domain extent.
!   It then extends communication buffers, which completely cover enough cells in the neighbouring
!   domain [so that waves do not propagate through the communication buffer in-between communications].
!
! Note that most of the 'setup' type functions in this module work
! with domains that have not run domain%allocate_quantities. 
!
module multidomain_mod

    use global_mod, only: dp, ip, charlen, default_timestepping_method, &
        wall_elevation, minimum_allowed_depth, force_double, force_long_double, &
        default_output_folder, send_boundary_flux_data, &
        real_bytes, force_double_bytes, integer_bytes, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV 
    use stop_mod, only: generic_stop
    use points_in_poly_mod, only: point_in_poly
    use ragged_array_mod, only: ragged_array_2d_dp_type, ragged_array_2d_ip_type
    use which_mod, only: which, rle_ip, bind_arrays_ip, remove_rows_ip, cumsum_ip
    use qsort_mod, only: match
    use nested_grid_comms_mod, only: process_received_data ! Import this to avoid using type bound procedure in omp region
    use coarray_point2point_comms_mod, only: allocate_p2p_comms, deallocate_p2p_comms, &
        linked_p2p_images, communicate_p2p, size_of_send_recv_buffers !, send_buffer, recv_buffer
    use iso_fortran_env, only: int64
    use logging_mod, only: log_output_unit, send_log_output_to_file
    use file_io_mod, only: mkdir_p, read_ragged_array_2d_ip
    use timer_mod, only: timer_type
#ifdef COARRAY_PROVIDE_CO_ROUTINES
    use coarray_intrinsic_alternatives, only: co_broadcast, co_max, co_sum, co_min
#endif

    implicit none

    private
    public :: multidomain_type
    public :: setup_multidomain, test_multidomain_mod

    ! Key parameter affecting parallel communication:
    ! - If TRUE, send halos as soon as they've been computed. This might give more time to overlap computation and comms
    ! - If FALSE, then send all at once. This allows us to do only one send to each image, which may have other efficiencies
    logical, parameter :: send_halos_immediately = .false. ! Must be .false. with MPI

    ! If we are doing coarray communication, the puts are non-blocking, and we
    ! sync at the start/end of we md%recv_halos. But if we use MPI for the main
    ! communication, this is not required (because alltoallv is blocking)
#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    logical, parameter :: sync_before_recv = .false., sync_after_recv=.false.
#else
    logical, parameter :: sync_before_recv = .true., sync_after_recv=.true.
#endif


    ! Parameter affecting default halo approach.
    !
    ! Set this to an integer > 0 so that "parts of another domain that I receive from"
    ! do not directly neighbour "parts of my domain that I send to that domain". 
    integer(ip), parameter :: extra_halo_buffer_default = 0_ip !
    ! This can be overridden in the multidomain setup stage as e.g. " call md%setup(extra_halo_buffer=1_ip) "
    ! This could help with stability with the "old" evolve_multidomain_one_step approach.
    ! Other references suggest such an approach can help with stability (e.g. Debreu et al 2012 Ocean modelling, paper on ROMS
    ! nesting). But it seems not required with the 'revised' nesting technique in SWALS
    ! It is ignored if we only communicate with domains that have the same domain%dx, as the benefit is really for coarse-to-fine
    ! communication

    ! This is another "padding" factor for the halos. By adding an extra-pad, we can ensure that e.g. gradient calculations are
    ! valid, when otherwise they might involve 'out-of-date' cells.
    integer(ip), parameter :: extra_cells_in_halo_default = 1_ip

    ! Local timestepping of domains
    ! This affects the distributed-memory version of the model, where we partition the 
    ! larger domains in parallel. Doing that allows "in principle" for some domains to take
    ! shorter timesteps than others (if the depth/speed vary significantly). We can exploit
    ! this to reduce model run-times (generally load-balancing will be required)
#ifdef LOCAL_TIMESTEP_PARTITIONED_DOMAINS
    logical, parameter :: local_timestep_partitioned_domains = .true.
#else
    logical, parameter :: local_timestep_partitioned_domains = .false.
#endif

    ! (co)array to store the interior bounding-box of ALL domains (i.e. not 
    ! including their nesting layer), and their dx(1:2) values, and some 
    ! metadata describing inter-domain communication needs. 
    ! Also store info on whether each domain uses a staggered grid, because this
    ! should affect details of how we do nesting communication (interpolation)
    ! It would seem cleaner to have this in a derived type, but for now coarray
    ! in derived type does not sound well supported.
#ifdef COARRAY
    real(dp), allocatable :: all_bbox(:,:,:,:)[:]
    real(dp), allocatable :: all_dx(:,:,:)[:]
    integer(ip), allocatable :: all_timestepping_refinement_factor(:,:)[:]
    real(dp), allocatable :: all_recv_metadata(:,:,:,:)[:]
    real(dp), allocatable :: all_send_metadata(:,:,:,:)[:]
    integer(ip), allocatable :: all_is_staggered_grid(:,:)[:]
#else
    real(dp), allocatable :: all_bbox(:,:,:,:)
    real(dp), allocatable :: all_dx(:,:,:)
    integer(ip), allocatable :: all_timestepping_refinement_factor(:,:)
    real(dp), allocatable :: all_recv_metadata(:,:,:,:)
    real(dp), allocatable :: all_send_metadata(:,:,:,:)
    integer(ip), allocatable :: all_is_staggered_grid(:,:)
#endif

    ! Store "this_image()", "num_images()"
    integer :: ti = -1 ! Dummy values which will cause error if not set later
    integer :: ni = -1


    ! Number of variables used to describe send/recv_metadata. We often need 
    ! this constant, so put it here to reduce 'magic numbers'
    integer(ip), parameter :: srm = 6

    ! Useful to have a very large integer, which is still smaller than huge(1_int64)
    integer(int64), parameter :: large_64_int = 10000000000_int64

    ! Main 'array of domains'. Idea is for this to have methods with names
    ! matching domain as much as possible, for easy translation of existing
    ! scripts.
    type :: multidomain_type

        ! If coarse and fine domains overlap, their dx's are either: identical,
        ! or; differ by some multiple of this value
        !integer(ip) :: nesting_ratio

        ! Array of domains 
        type(domain_type), allocatable :: domains(:)

        ! If we split domains in parallel (COARRAYS), then
        ! this will hold the 'initial, unallocated' domains
        type(domain_type), allocatable :: domain_metadata(:)

        ! Domain volume tracking
        real(force_double), allocatable :: volume_initial(:), volume(:)

        ! Output directory
        character(len=charlen) :: output_basedir = default_output_folder

        type(timer_type) :: timer

        ! The lower/upper extents of periodic domains. For example, on the earth,
        ! we could have EW periodic spherical coordinates with 
        ! periodic_xs = [-180.0_dp, 180.0_dp]
        real(dp) :: periodic_xs(2) = [-HUGE(1.0_dp), HUGE(1.0_dp)]
        real(dp) :: periodic_ys(2) = [-HUGE(1.0_dp), HUGE(1.0_dp)]

        ! File with information on load balancing
        character(len=charlen) :: load_balance_file = ''
        ! If load_balancing_file is provided, we can use this
        type(ragged_array_2d_ip_type) :: load_balance_part

        ! Control additional halo width. This can enable separation of send/recv halos,
        ! which can reduce nesting artefacts in some cases
        integer(ip) :: extra_halo_buffer = extra_halo_buffer_default
        integer(ip) :: extra_cells_in_halo = extra_cells_in_halo_default

        ! Store all_dx in the multidomain
        real(dp), allocatable :: all_dx_md(:,:,:)

        ! Store 'all_is_staggered_grid' in the multidomain
        integer(ip), allocatable :: all_is_staggered_grid_md(:,:)

        contains

        ! Main initialisation routine
        procedure :: setup => setup_multidomain
        procedure :: partition_domains => partition_domains
        ! Main time-stepper
        procedure :: evolve_one_step => evolve_multidomain_one_step
        ! Utilities to ensure consistency of nesting areas prior to main computation
        procedure :: set_null_regions_to_dry => set_null_regions_to_dry
        procedure :: use_constant_wetdry_send_elevation => use_constant_wetdry_send_elevation
        procedure :: make_initial_conditions_consistent => make_initial_conditions_consistent
        ! Mass tracking
        procedure :: get_flow_volume => get_flow_volume
        procedure :: record_initial_volume => record_initial_volume
        procedure :: report_mass_conservation_statistics => report_mass_conservation_statistics
        ! Inter-domain communication
        procedure :: send_halos => send_multidomain_halos
        procedure :: recv_halos => receive_multidomain_halos
        procedure :: separate_halos => separate_fine_halos_from_their_domain
        ! Memory tracking
        procedure :: memory_summary => memory_summary
        ! Convenience
        procedure :: print => print_multidomain
        procedure :: finalise_and_print_timers => finalise_and_print_timers

    end type

    contains 

    !
    ! Given an x, y point, and size=2 arrays periodic_xs, periodic_ys defining
    ! the x and y extent of a periodic domain, compute x/y values inside the periodic domain.
    ! If this leads to a change in the x/y coordinate (i.e. if the initial x or y was
    ! outside the ranges provided by periodic_xs and periodic_ys), then set periodic_point
    ! to .true., otherwise set it to .false.
    !
    ! @param x real x coordinate value
    ! @param y real y coordinate value
    ! @param periodic_x array with (xmin, xmax) giving the range of x values in the periodic domain
    ! @param periodic_y array with (ymin, ymax) giving the range of y values in the periodic domain
    ! @param periodic_point logical output variable
    ! @param adjust_coordinates if true then change x/y coordinate values to be inside the main domain, assuming x and y are "just
    ! outside" the main domain (e.g. if periodic_xs = [0, 360], then x can be in the range (-360, 720), but we can't have more
    ! coordinate wrapping.).
    subroutine check_periodic(x, y, periodic_xs, periodic_ys, periodic_point, adjust_coordinates)
        real(dp), intent(inout) :: x, y
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)
        logical, intent(out) :: periodic_point
        logical, intent(in) :: adjust_coordinates

        if(periodic_xs(2) <= periodic_xs(1)) then
            print*, 'periodic_xs(2) should be > periodic_xs(1)'
            call generic_stop
        end if
        if(periodic_ys(2) <= periodic_ys(1)) then
            print*, 'periodic_ys(2) should be > periodic_ys(1)'
            call generic_stop
        end if

        periodic_point = .false.

        ! Assess 'x'
        if(x < periodic_xs(1)) then
            ! Be careful about precision loss
            if(adjust_coordinates) x = real(periodic_xs(2) - periodic_xs(1), force_long_double) + real(x, force_long_double)
            periodic_point = .true.
        end if
        if(x > periodic_xs(2)) then
            ! Be careful about precision loss
            if(adjust_coordinates) x = real(periodic_xs(1) - periodic_xs(2), force_long_double) + real(x, force_long_double)
            periodic_point = .true.
        end if
        if(adjust_coordinates .and. (x < periodic_xs(1) .or. x > periodic_xs(2))) then
            print*, 'Error in setting up periodic boundaries'
            print*, 'Make sure periodic_xs covers the input x-range of the multidomain'
        end if

        ! Assess 'y'
        if(y < periodic_ys(1)) then
            ! Be careful about precision loss
            if(adjust_coordinates) y = real(periodic_ys(2) - periodic_ys(1), force_long_double) + real(y, force_long_double)
            periodic_point = .true.
        end if
        if(y > periodic_ys(2)) then
            ! Be careful about precision loss
            if(adjust_coordinates) y = real(periodic_ys(1) - periodic_ys(2), force_long_double) + real(y, force_long_double)
            periodic_point = .true.
        end if
        if(adjust_coordinates .and. (y < periodic_ys(1) .or. y > periodic_ys(2))) then
            print*, 'Error in setting up periodic boundaries'
            print*, 'Make sure periodic_ys covers the input y-range of the multidomain'
        end if

    end subroutine

    !
    ! For each domain in domains, figure out how thick the nesting layer
    ! has to be. Also determine which boundaries have a nesting buffer, 
    ! which are physical boundaries, etc.
    !
    ! This is called before domains have been buffered to include nesting
    ! regions -- so that we know how thick to make the buffer.
    ! 
    ! @param domains array of domains 
    ! @param verbose logical
    !
    subroutine compute_multidomain_nesting_layer_width(domains, verbose,&
        periodic_xs, periodic_ys, extra_halo_buffer, extra_cells_in_halo)
        type(domain_type), intent(inout) :: domains(:)
        logical, optional, intent(in) :: verbose
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)
        integer(ip), intent(in) :: extra_halo_buffer
        integer(ip), intent(in) :: extra_cells_in_halo

        integer(ip):: i, j, ii, nd_local, nest_layer_width, boundary_flag, n1

        type(ragged_array_2d_ip_type), allocatable :: nbr_domain_ind(:), nbr_image_ind(:)

        real(dp), allocatable:: nbr_domain_dx_local(:,:)
        logical, allocatable :: is_nesting_boundary(:,:)

        !logical :: inb
        real(dp) :: max_parent_dx_ratio
        !real(dp) :: max_parent_timestepping_refinement_ratio, rts
        logical :: verbose1

        if(present(verbose)) then
            verbose1 = verbose
        else
            verbose1 = .true.
        end if

        nd_local = size(domains)
        allocate(is_nesting_boundary(4,nd_local))

        ! Allocate arrays of ragged arrays (one entry for each domain on this_image), which
        ! will hold metadata on the other domains just outside its boundary. We use that 
        ! to help determine nesting relationships.
        allocate(nbr_domain_ind(nd_local), &
            nbr_image_ind(nd_local))

        do j = 1, nd_local
            ! Allocate ragged arrays to store data for each boundary (N, E, S, W)
            allocate(nbr_domain_ind(j)%i2(4), &
                nbr_image_ind(j)%i2(4))
        end do

        ! Set the coarray-related module variables ni, ti 
#ifdef COARRAY
        ni = num_images()
        ti = this_image()
#else
        ni = 1
        ti = 1
#endif        

        ! Store the dx and 'interior' bounding box of all domains
        call create_all_bbox_and_dx(domains, periodic_xs, periodic_ys)

        !
        ! Loop over north(1), east(2), south(3), west(4) domain interior 
        ! boundaries and find which other domains (if any) cells just outside 
        ! our interior boundary would lie inside.
        !
        ! We want to receive data from the 'best' domain available. If we need 
        ! data from a region containing multiple overlapping domains, then we 
        ! take data from the one having finest dx.
        !
        ! If cells just outside our boundary are not inside any domain, then
        ! that boundary is a 'physical' boundary [i.e. better be treated with 
        ! a boundary condition]. 
        !
        ! 
        do j = 1, nd_local

            if(verbose1) write(log_output_unit, *) 'Domain ', j

            ! Initial value to be refined in loop            
            max_parent_dx_ratio = 1.0_dp

            ! Loop over all boundaries
            do i = 1, 4 ! N, E, S, W

                boundary_flag = i

                ! Preliminary only
                nest_layer_width = 1_ip

                ! For points one cell width outside our domain, find which
                ! other domain they are inside [if any]. If they are inside
                ! multiple domains, then select the one with finest dx
                call find_priority_nesting_domains_along_a_boundary(&
                    domains(j), boundary_flag, & 
                    nest_layer_width, &
                    nbr_domain_ind(j)%i2(i)%i1, &
                    nbr_image_ind(j)%i2(i)%i1, &
                    nbr_domain_dx_local,&
                    periodic_xs = periodic_xs,&
                    periodic_ys = periodic_ys)

                ! Find the max ratio between parent cell sizes and the current
                ! domain's cell size. Note that nbr_domain_dx_local will hold a
                ! large negative number where there is no neighbour domain
                !
                ! Ultimately this will be (Maximum value of 'parent domain dx' / 'my dx')
                !
                ! x-direction
                max_parent_dx_ratio = max(max_parent_dx_ratio, real(nint(&
                    maxval(nbr_domain_dx_local(:,1)/domains(j)%dx(1))), dp) )
                ! y-direction
                max_parent_dx_ratio = max(max_parent_dx_ratio, real(nint(&
                    maxval(nbr_domain_dx_local(:,2)/domains(j)%dx(2))), dp) )

            end do

            if (verbose1) write(log_output_unit, *) 'max_parent_dx_ratio: ', max_parent_dx_ratio
            domains(j)%max_parent_dx_ratio = max_parent_dx_ratio

            ! Record whether this boundary does nesting (or not)
            do i = 1, 4
                if(verbose1) write(log_output_unit,*) '.. ', i, any(nbr_image_ind(j)%i2(i)%i1 > 0)
                is_nesting_boundary(i,j) = any(nbr_image_ind(j)%i2(i)%i1 > 0)
            end do
            if(verbose1) write(log_output_unit,*) 'Nesting boundaries: ', is_nesting_boundary(:,j)
            domains(j)%is_nesting_boundary = is_nesting_boundary(:,j)
    
            ! Now we know the dx of all parent domains, we can compute  
            ! the required nest_layer_width for domains(j)
            nest_layer_width = get_domain_nesting_layer_thickness(&
                domains(j), max_parent_dx_ratio, extra_halo_buffer, extra_cells_in_halo)
            domains(j)%nest_layer_width = nest_layer_width
            !domains(j)%exterior_cells_width = nest_layer_width
            if(verbose1) write(log_output_unit,*) 'Nesting layer width: ', domains(j)%nest_layer_width

            ! Logical checks
            do i = 1, 4
                ! Ensure that there are no 'half-nested' boundaries
                ! Only check cells right next to the original domain boundary,
                ! since it is possible that other cells are not nested [e.g. if this domain
                ! lies on the physical boundary of the "full multidomain extent" ]
                do ii = 2, size(nbr_image_ind(j)%i2(i)%i1) - 1 
                    if((nbr_image_ind(j)%i2(i)%i1(ii) > 0)) then
                        if((nbr_image_ind(j)%i2(i)%i1(ii) < 0)) then
                            write(log_output_unit,*) 'Boundary cannot be half nested and half un-nested'
                            write(log_output_unit,*) ' but this occurs on image ', ti, ' on domains ', j , &
                                'on boundary ', i, '. The domain coordinates are: '
                            write(log_output_unit,*) all_bbox(:,:, j, ti)
                            call generic_stop()
                        end if
                    end if
                end do
            end do
        end do

    end subroutine

    !
    ! Convert the arrays 'priority_domain_index' and 'priority_domain_image' to a set
    ! of boxes which are well suited to passing to communication routines
    !
    !
    subroutine convert_priority_domain_info_to_boxes(priority_domain_index, &
        priority_domain_image, halo_width, myimage, myindex, box_metadata, &
        xs, ys, periodic_xs, periodic_ys, max_parent_dx_ratio)

        integer(ip), intent(in) :: priority_domain_index(:,:), &
            priority_domain_image(:,:)
        integer(ip), intent(in) :: myindex, myimage, halo_width
        integer(ip), allocatable, intent(inout) :: box_metadata(:,:)
        real(dp), intent(in):: xs(:), ys(:)
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)
        integer(ip), intent(in) :: max_parent_dx_ratio

        integer(ip) :: dims(2), i, j, ii, i1, i2, j1, j2, i1n, j1n, mpdx, ind1, im1
        integer(ip), allocatable :: nest_ind(:), nest_image(:), match_ind(:)
        integer(ip), allocatable :: run_indices(:), run_values(:), &
            run_lengths(:), ends(:), starts(:)
        integer(ip), allocatable :: candidate_boxes(:,:), final_boxes(:,:), &
            candidate_boxes_new(:,:)
        integer(ip), allocatable :: new_boxes(:), end_boxes(:)
        logical, allocatable :: in_active(:), in_periodic_x(:), in_periodic_y(:)
        character(len=charlen), allocatable :: box_i_string(:), cbox_i_string(:)

        integer(ip), allocatable :: buffer_priority_domain_index(:,:)
        integer(ip), allocatable :: buffer_priority_domain_image(:,:)

        i1n = size(priority_domain_image, 1)
        j1n = size(priority_domain_image, 2)
        allocate(buffer_priority_domain_index(i1n, j1n))
        allocate(buffer_priority_domain_image(i1n, j1n))

        ! ID cells that are in the periodic boundary regions. This is necessary because
        ! we must avoid the 'boxes' crossing the periodic boundary
        allocate(in_periodic_x(size(xs)), in_periodic_y(size(ys)))
        in_periodic_x = (xs < periodic_xs(1) .or. xs > periodic_xs(2))
        in_periodic_y = (ys < periodic_ys(1) .or. ys > periodic_ys(2))

        !
        ! Make a 'buffered' priority domain index/image, to move nesting regions
        ! away from the priority domain. This prevents send and recv regions
        ! sharing boundaries. We only move recv regions for 'fine' domains.
        ! For coarse domains, for mass conservation purposes it's important
        ! to keep its recv regions right on its boundary
        !
        buffer_priority_domain_index = priority_domain_index
        buffer_priority_domain_image = priority_domain_image

        !! Define buffer size
        !mpdx = max_parent_dx_ratio * halo_separation_factor 
        !! Update the buffer_priority_domain variables
        !do j = 1, size(priority_domain_index,2)
        !    do i = 1, size(priority_domain_index,1)
        !        ! Buffer by (1 * max_parent_dx_ratio)   
        !        !
        !        ! If any cell within "max_parent_dx_ratio" cells is in {my_index, my_image}
        !        ! then we do not want to have a send region here
        !        do j1 = -mpdx, mpdx
        !            do i1 = -mpdx, mpdx
        !                ! Avoid out-of-bounds errors
        !                if( ((i1 + i) < 1) .or. ((j1 + j) < 1) .or. & 
        !                    ((i1 + i) > size(priority_domain_index, 1)) .or. &
        !                    ((j1 + j) > size(priority_domain_index, 2)) ) cycle 

        !                ! Avoid changes in periodic regions
        !                if(in_periodic_x(i) .or. in_periodic_y(j)) cycle

        !                ind1 = priority_domain_index(i, j)
        !                im1  = priority_domain_image(i, j)

        !                ! If the associated domain cell is finer or equal to the domain of interest, do not expand
        !                ! Use a factor 1.1 to avoid round-off, considering the areas should be integer
        !                ! multiples of each other
        !                if(product(all_dx(1:2, ind1, im1)) <= 1.1_dp*product(all_dx(1:2, myindex, myimage))) cycle

        !                !! Only expand if at least one of the domains is staggered
        !                !if(.not.(all_is_staggered_grid(ind1, im1) == 1 .or. &
        !                !         all_is_staggered_grid(myindex, myimage) == 1)) cycle

        !                ! "Grow" the region with {my_index, my_image} to prevent receiving there
        !                if( (priority_domain_index(i1 + i, j1 + j) == myindex) .and. &
        !                    (priority_domain_image(i1 + i, j1 + j) == myimage) ) then
        !                    buffer_priority_domain_index(i,j) = priority_domain_index(i1+i, j1+j)
        !                    buffer_priority_domain_image(i,j) = priority_domain_image(i1+i, j1+j)
        !                end if
        !            end do
        !        end do 
        !    end do
        !end do 

        ! box_metadata = main output
        if(allocated(box_metadata)) deallocate(box_metadata)

        dims = shape(buffer_priority_domain_index)

        ! Work arrays
        allocate(in_active(dims(1)), nest_ind(dims(1)), nest_image(dims(1)), &
            run_indices(dims(1)))

        run_indices = (/ (i, i=1, dims(1)) /)

    
    
        ! Loop over domain y-coordinate
        do j = 1, dims(2)

            ! For each x coordinate, find whether we are in the active part of
            ! the domain. This let's us exclude inactive 'NULL' areas, which are
            ! internal parts of the domain covered by finer domains, which are not
            ! even required for communication
            do i = 1, dims(1)
                ! FIXME: This would have to be changed to allow a range of
                ! halo_width's, if we do 'multi-time-level-comms'
                call is_in_active_region(buffer_priority_domain_index, &
                    buffer_priority_domain_image, myindex, myimage, &
                    halo_width, i, j, in_active(i))
            end do

            ! Record the domain-index and image of each point in the active
            ! region -- non-active points are set to zero
            nest_ind = merge(buffer_priority_domain_index(:,j), &
                0 * buffer_priority_domain_index(:,j), in_active)
            nest_image = merge(buffer_priority_domain_image(:,j), &
                0 * buffer_priority_domain_image(:,j), in_active)
            ! Find 'start' and 'end' indices corresponding to x values with a
            ! contiguous priority domain/image, AND never crossing periodic_xs or periodic_ys. 
            ! Note that 'rle_ip' is similar
            ! to the function 'rle' in the R language, which gives 'run-lengths'
            ! of repeated variables. 
            call rle_ip(run_indices, run_values, run_lengths, &
                equality_function=equal_domain_image_and_index)
            
            ! 'ends' gives the end index of the chunk
            if(allocated(ends)) deallocate(ends)
            ends = run_lengths
            call cumsum_ip(ends) 
            ! 'starts' gives the start index of the chunk
            if(allocated(starts)) deallocate(starts)
            allocate(starts(size(ends)))
            starts = ends - (run_lengths - 1) 

            ! Sanity check
            if(maxval(ends) > dims(1)) then
                write(log_output_unit,*) 'ERROR: j = ', j
                write(log_output_unit,*) 'ends: ', ends
                write(log_output_unit,*) 'run_values: ', run_values
                write(log_output_unit,*) 'run_lengths: ', run_lengths
                stop
            end if

            if(j == 1) then
                ! Setup 'candidate_boxes' and 'final_boxes' (the latter will
                ! eventually be box_metadata)
                !
                ! Candidate_boxes is a table (integer matrix) with 6 variables:
                ! domain_index, domain_image, box_left_i, box_right_i, 
                ! box_bottom_j, box_top_j
                allocate(candidate_boxes(srm, size(run_values)))
                candidate_boxes(1,:) = nest_ind(run_values)
                candidate_boxes(2,:) = nest_image(run_values)
                candidate_boxes(3,:) = starts
                candidate_boxes(4,:) = ends
                candidate_boxes(5,:) = j + 0*starts
                ! The upper-right y index is unknown at present, use -1 to
                ! represent NA. We will find the upper-right index when we scan
                ! the j row just past the box. At that time, the box is can be
                ! finalised
                candidate_boxes(6,:) = -1 + 0*starts 

                ! Final_boxes has the same structure, but the 6th column is
                ! known. Once boxes are finalised, we move them to 'final_boxes'
                allocate(final_boxes(srm,0))

            else
                
                if(allocated(box_i_string)) deallocate(box_i_string)
                allocate( box_i_string( size(starts) ) )
                if(allocated(cbox_i_string)) deallocate(cbox_i_string)
                allocate( cbox_i_string( size(candidate_boxes(1,:)) ) )

                ! Perform a 'match' of starts/ends on the current row with 
                ! starts/ends of existing candidate boxes.
                ! When there is no match of starts/ends/nest_index/nest_image, 
                ! we know the candiate box has finished, and so we can populate 
                ! its upper_right y value, and move it to the final boxes.
                ! The match is done using strings. For current domain y-value, 
                ! and the candidate boxes, we make a string containing 
                ! domain_index, domain_image, box_left_i, box_right_i. We put 'in_periodic_y(j)'
                ! on the end (to prevent boxes from crossing 'y' periodic' boundaries). Note
                ! we already prevented it from crossing x periodic boundaries by the equality function
                ! passed to rle_ip
                do i = 1, size(run_values)
                    ii = run_values(i)
                    write(box_i_string(i), *) nest_ind(ii), ',', &
                        nest_image(ii), ',', starts(i), ',', ends(i), ',', in_periodic_y(j)
                end do
                do i = 1, size(cbox_i_string)
                    write(cbox_i_string(i), *) candidate_boxes(1,i), ',', &
                        candidate_boxes(2,i), ',', candidate_boxes(3,i), ',', &
                        candidate_boxes(4,i), ',', in_periodic_y(j-1)
                end do

                ! Boxes without a match are 'new'
                if(allocated(match_ind)) deallocate(match_ind)
                allocate(match_ind(size(box_i_string)))
                call match(box_i_string, cbox_i_string, match_ind)
                call which(match_ind == -1_ip, new_boxes)
                
                ! Candidate boxes without a match are 'finished' 
                deallocate(match_ind)
                allocate(match_ind(size(cbox_i_string)))
                call match(cbox_i_string, box_i_string, match_ind)
                call which(match_ind == -1_ip, end_boxes)

                ! Transfer completed candidate boxes to the final_boxes array,
                ! and remove them from the candidate boxes array
                if( size(end_boxes) > 0) then
                    ! We now know the upper-right y index of the box
                    candidate_boxes(srm,end_boxes) = j - 1 
                    call bind_arrays_ip(final_boxes, &
                        candidate_boxes(:,end_boxes), rowbind=.false.)
                    call remove_rows_ip(candidate_boxes, end_boxes, &
                        apply_to_rows = .false.)
                end if

                !
                ! Add new boxes to the candidate_boxes array
                !
                if( size(new_boxes) > 0) then
                    if(allocated(candidate_boxes_new)) deallocate(candidate_boxes_new)
                    allocate( candidate_boxes_new(srm, size(new_boxes)) )
                    candidate_boxes_new(1,:) = nest_ind(run_values(new_boxes))
                    candidate_boxes_new(2,:) = nest_image(run_values(new_boxes))
                    candidate_boxes_new(3,:) = starts(new_boxes)
                    candidate_boxes_new(4,:) = ends(new_boxes)
                    candidate_boxes_new(5,:) = j + 0*new_boxes
                    ! For unknown upper-right y-index, use -1 
                    candidate_boxes_new(6,:) = -1 + 0*new_boxes 
                    call bind_arrays_ip(candidate_boxes, candidate_boxes_new, &
                        rowbind=.false.)
                end if

                ! On the last row, any candidate boxes should be finalised
                if( j == dims(2) .and. size(candidate_boxes) > 0) then
                    candidate_boxes(6,:) = j
                    call bind_arrays_ip(final_boxes, candidate_boxes, &
                        rowbind=.false.)
                end if

            end if

        end do

        dims = shape(final_boxes) 
        allocate(box_metadata(dims(1), dims(2))) 
        box_metadata = final_boxes

        ! Final check
        do j = 1, size(box_metadata(1,:))
            ! The logic below won't work for boxes matching 0
            ! because we don't distinguish that for priority_domain_index
            if(box_metadata(1,j) == 0) cycle
            i1 = box_metadata(3,j)
            i2 = box_metadata(4,j)
            j1 = box_metadata(5,j)
            j2 = box_metadata(6,j)
            if( any( buffer_priority_domain_index(i1:i2,j1:j2) /= box_metadata(1,j)) .or. &
                any( buffer_priority_domain_image(i1:i2,j1:j2) /= box_metadata(2,j)) .or. &
                any( in_periodic_x(i1:i2) .neqv. in_periodic_x(i1) ).or. &
                any( in_periodic_y(j1:j2) .neqv. in_periodic_y(j1) ) ) then
                write(log_output_unit,*) 'image :', myimage, ' index: ', myindex
                write(log_output_unit,*) 'invalid: ', box_metadata(:,j)
                write(log_output_unit,*) '.......: ', final_boxes(:,j)
                stop 'Error in creating box_metadata for priority domain information'
            end if
        end do

        deallocate(buffer_priority_domain_index)
        deallocate(buffer_priority_domain_image)

        contains
            ! Subroutine used in 'rle_ip' call to detect changes in
            ! priority_domain, or crossing of periodic boundaries
            function equal_domain_image_and_index(i1, i2) result(istrue)
                integer(ip), intent(in) :: i1, i2
                logical :: istrue

                istrue = ( (nest_ind(i1) == nest_ind(i2)) .and. &
                           (nest_image(i1) == nest_image(i2)) .and. &
                            ! This condition ensures no boxes cross the
                            ! periodic boundary along the x-axis. We use
                            ! a different constraint (in string matching of candidate boxes)
                            ! to ensure boxes do not cross the periodic boundary along
                            ! the y axis
                           (in_periodic_x(i1) .eqv. in_periodic_x(i2)) ) 

            end function

    end subroutine

    !
    ! Update the values of all_bbox and all_dx, so they contain
    ! info on the bounding box's and dx values of all domains on all images.
    !
    ! For example, all_bbox(1:4, 1:2, i, j) will contain the bounding box
    ! of domains(i) on image j
    !
    ! @param domains the array of domains
    ! 
    subroutine create_all_bbox_and_dx(domains, periodic_xs, periodic_ys)

        type(domain_type), intent(inout) :: domains(:)
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)

        integer(ip):: nd, i, j, k, nd_local

        nd_local = size(domains)

#ifdef COARRAY
        nd = nd_local
        call co_max(nd)
#else
        nd = nd_local
#endif        

        if(allocated(all_bbox)) deallocate(all_bbox)
        if(allocated(all_dx)) deallocate(all_dx)
        if(allocated(all_timestepping_refinement_factor)) deallocate(all_timestepping_refinement_factor)
        if(allocated(all_is_staggered_grid)) deallocate(all_is_staggered_grid)

        ! Get interior box of other boundaries
#ifdef COARRAY
        allocate(all_bbox(4, 2, nd, ni)[*])
        allocate(all_dx(2, nd, ni)[*])
        allocate(all_timestepping_refinement_factor(nd, ni)[*])
        allocate(all_is_staggered_grid(nd, ni)[*])
#else
        allocate(all_bbox(4, 2, nd, ni))
        allocate(all_dx(2, nd, ni))
        allocate(all_timestepping_refinement_factor(nd, ni))
        allocate(all_is_staggered_grid(nd, ni))
#endif
        all_bbox = 0.0_dp
        all_dx = 0.0_dp
        all_timestepping_refinement_factor = 0_ip
        all_is_staggered_grid = -1_ip

        ! Save the interior bounding box. If the latter has not yet been defined,
        ! assume it is the same as the bounding box defined by the lower left and length/width.
        do i = 1, nd_local
            if(all(domains(i)%interior_bounding_box == 0.0_dp)) then

                domains(i)%interior_bounding_box(1, 1:2) = domains(i)%lower_left(1:2)
                domains(i)%interior_bounding_box(2, 1:2) = domains(i)%lower_left(1:2) + &
                    domains(i)%lw * [1.0_dp, 0.0_dp]
                domains(i)%interior_bounding_box(3, 1:2) = domains(i)%lower_left(1:2) + &
                    domains(i)%lw
                domains(i)%interior_bounding_box(4, 1:2) = domains(i)%lower_left(1:2) + &
                    domains(i)%lw * [0.0_dp, 1.0_dp]
            end if
            ! Store interior bounding box
            all_bbox(:,:,i,ti) = domains(i)%interior_bounding_box(1:4, 1:2)
            ! Store dx
            all_dx(1:2, i, ti) = domains(i)%dx
            all_timestepping_refinement_factor(i,ti) = domains(i)%timestepping_refinement_factor
            ! 1 if our scheme uses a staggered grid, 0 otherwise. Currently only the linear solver does
            all_is_staggered_grid(i,ti) = merge(1_ip, 0_ip, domains(i)%timestepping_method == 'linear')
        end do

#ifdef COARRAY
        ! Need a synt to ensure that every domain has already set its all_bbox
        ! to zero a few lines above. Otherwise, below we could communicate our data, but then
        ! have it set to zero elsewhere!
        sync all 
        ! Broadcast the current 'domains' metadata to all other images.
        ! ni-to-ni one-sided communication.
        do k = 1, ni
            all_bbox(1:4, 1:2, 1:nd, ti)[k] = all_bbox(1:4, 1:2, 1:nd, ti)
            all_dx(1:2, 1:nd, ti)[k] = all_dx(1:2, 1:nd, ti)
            all_timestepping_refinement_factor(1:nd,ti)[k] = all_timestepping_refinement_factor(1:nd, ti) 
            all_is_staggered_grid(1:nd, ti)[k] = all_is_staggered_grid(1:nd, ti)
        end do
        sync all
        ! At this point we know the interior bounding box for every other domain
#endif
        !write(log_output_unit, *) 'All_bbox: '
        do k = 1, ni
            do i = 1, nd
                !
                ! Double check that periodic_xs/periodic_ys do not cross the input multidomain.
                ! If this happens, there is an input error
                if( (any(all_bbox(:,1,i,k) < periodic_xs(1)) .or. any(all_bbox(:,1,i,k) > periodic_xs(2)) .or. &
                     any(all_bbox(:,2,i,k) < periodic_ys(1)) .or. any(all_bbox(:,2,i,k) > periodic_ys(2))) .and. &
                    (.not. all(all_bbox(:,:,i,k) == 0.0_dp)) ) then
                    write(log_output_unit,*) 'Error: one of the input bounding boxes extends outside the periodic domain'
                    write(log_output_unit,*) 'This is not permitted. If periodic boundary conditions are desired, then'
                    write(log_output_unit,*) 'the extremes of the bounding boxes in the multidomain should touch the '
                    write(log_output_unit,*) 'box defined by periodic_xs and periodic_ys, but not go outside it'
                    write(log_output_unit,*) ''
                    write(log_output_unit,*) '.... the offending bounding box is ... '
                    write(log_output_unit,*) all_bbox(:,:,i,k)
                    write(log_output_unit,*) 'with index ', i, k
                    write(log_output_unit,*) '.....and the periodic extents are ...'
                    write(log_output_unit,*) 'periodic_xs: ', periodic_xs, '; periodic_ys: ', periodic_ys
                    call generic_stop
                end if
            end do
        end do
        write(log_output_unit, *) '    '

    end subroutine


    ! Timestep all grids (by the global time-step).
    ! 
    ! This method only includes communication at each global time-step. While that's
    ! good in terms of reducing comms, it does mean comms buffers might need
    ! to be very large. To reduce the comms buffer size, we could have more
    ! communication in sub-steps of the global timestep 
    !
    subroutine evolve_multidomain_one_step(md, dt)
        class(multidomain_type), intent(inout) :: md
        real(dp), intent(in) :: dt

        integer(ip) :: j, i, nt
        real(dp) :: dt_local, max_dt_store, tend

        ! All the domains will evolve to this time 
        tend = md%domains(1)%time + dt
        ! However they may take different numbers of steps to get 
        ! there, leading to tiny numerical differences in the time.
        ! To avoid any issues, we explicitly set them to the same number

        ! Evolve every domain, sending the data to communicate straight after
        ! the evolve
        do j = 1, size(md%domains)

            TIMER_START('domain_evolve')

            ! Re-set the nesting boundary flux integrals
            ! These will be used to determine how much flux-correction to apply
            ! following these evolve steps
            call md%domains(j)%nesting_boundary_flux_integral_multiply(c=0.0_dp)

            ! Evolve each domain one or more steps, for a total time of dt
            nt = md%domains(j)%timestepping_refinement_factor

            if(local_timestep_partitioned_domains .and. md%domains(j)%timestepping_method /= 'linear') then
                ! For nonlinear domains, allow fewer timesteps, if it won't cause blowup. 
                !
                ! This is most important in distributed-memory applications where the partitioned
                ! domain could support substantially different time-steps in different parts of the
                ! "big domain". In combination with load balancing, we can get large speedups.
                !
                ! NOTE: This will lead to different timestepping with different numbers of cores,
                ! so the results should depend on the number of cores.
                nt = max(1, min(nt, ceiling(dt/md%domains(j)%max_dt) ))
            end if

            dt_local = dt/(1.0_dp * nt)

            ! Step once
            call md%domains(j)%evolve_one_step(dt_local)

            ! Report max_dt as the peak dt during the first timestep
            ! In practice max_dt may cycle between time-steps
            max_dt_store = md%domains(j)%max_dt

            ! Do the remaining sub-steps
            do i = 2, nt ! Loop never runs if nt < 2
                call md%domains(j)%evolve_one_step(dt_local)
            end do
            md%domains(j)%max_dt = max_dt_store

            TIMER_STOP('domain_evolve')

            ! Send the halos only for domain j
            call md%send_halos(domain_index=j, send_to_recv_buffer = send_halos_immediately)

            ! Ensure numerically identical time for all domains
            md%domains(j)%time = tend
        end do

        if(.not. send_halos_immediately) then
            ! Do all the coarray communication in one hit
            ! This lets us 'collapse' multiple sends to a single image,
            ! and is more efficient in a bunch of cases.
            TIMER_START('comms1')
            call communicate_p2p
            TIMER_STOP('comms1')
        end if
        
        ! Get the halo information from neighbours
        ! For coarray comms, we need to sync before to ensure the information is sent, and
        ! also sync after to ensure it is not overwritten before it is used
        call md%recv_halos(sync_before=sync_before_recv, sync_after=sync_after_recv)

        if(send_boundary_flux_data) then
            ! Do flux correction
            do j = 1, size(md%domains)
                call md%domains(j)%nesting_flux_correction_everywhere(md%all_dx_md, &
                    md%all_is_staggered_grid_md, fraction_of = 1.0_dp)
            end do
        end if

    end subroutine

    !! DELETE THIS AT SOME POINT
    !!
    !! THIS OLD VERSION OF evolve_multidomain_one_step REQUIRES MULTIPLE COMMUNICATIONS,
    !! AND DOES NOT HAVE SUCH GOOD REPRODUCIBILITY OF OMP/COARRAY VERSIONS, COMPARED WITH
    !! THE NEW VERSION.
    !! FOR STABILITY IT REQUIRES EXTRA_HALO_BUFFER=1 ON COMPLICATED NESTING PROBLEMS.
    !! IT MAY ALSO REQUIRE THE WET/DRY "iterative correction" TREATMENT IN domain_mod.f90 (compute_depth_and_velocity).
    !! THUS THE NEW VERSION SEEMS MUCH BETTER -- BUT KEEP THIS FOR NOW IN CASE WE DISCOVER ANY ISSUES.
    !! 
    !!
    !! Timestep all grids (by the global time-step).
    !! 
    !! This method only includes communication at each global time-step. While that's
    !! good in terms of reducing comms, it does mean comms buffers might need
    !! to be very large. To reduce the comms buffer size, we could have more
    !! communication in sub-steps of the global timestep 
    !!
    !subroutine evolve_multidomain_one_step(md, dt)
    !    class(multidomain_type), intent(inout) :: md
    !    real(dp), intent(in) :: dt

    !    integer(ip) :: j, i, nt
    !    real(dp) :: dt_local, max_dt_store
    !    logical :: lagged_flux_correction_local

    !    ! If send_boundary_flux_data is .false., then we should never do lagged_flux_correction
    !    ! (because no flux corrections are sent in any case!)
    !    lagged_flux_correction_local = (md%lagged_flux_correction .and. send_boundary_flux_data)

    !    ! Evolve every domain, sending the data to communicate straight after
    !    ! the evolve
    !    do j = 1, size(md%domains)

    !        TIMER_START('domain_evolve')

    !        ! Re-set the nesting boundary flux integrals
    !        ! These will be used to determine how much flux-correction to apply
    !        ! following these evolve steps
    !        call md%domains(j)%nesting_boundary_flux_integral_multiply(c=0.0_dp)

    !        ! Evolve each domain one or more steps, for a total time of dt
    !        nt = md%domains(j)%timestepping_refinement_factor
    !        dt_local = dt/(1.0_dp * nt)

    !        ! Step once
    !        call md%domains(j)%evolve_one_step(dt_local)
    !        if(lagged_flux_correction_local) call md%domains(j)%nesting_flux_correction(1.0_dp/nt)
    !        ! Report max_dt as the peak dt during the first timestep
    !        ! In practice max_dt may cycle between time-steps
    !        max_dt_store = md%domains(j)%max_dt

    !        ! Do the remaining sub-steps
    !        do i = 2, nt ! Loop never runs if nt < 2
    !            call md%domains(j)%evolve_one_step(dt_local)
    !            if(lagged_flux_correction_local) call md%domains(j)%nesting_flux_correction(1.0_dp/nt)
    !        end do
    !        md%domains(j)%max_dt = max_dt_store

    !        TIMER_STOP('domain_evolve')

    !        ! Send the halos only for domain j
    !        call md%send_halos(domain_index=j, send_to_recv_buffer = send_halos_immediately)

    !    end do

    !    if(.not. send_halos_immediately) then
    !        ! Do all the coarray communication in one hit
    !        ! This lets us 'collapse' multiple sends to a single image,
    !        ! and is more efficient in a bunch of cases.
    !        TIMER_START('comms1')
    !        call communicate_p2p
    !        TIMER_STOP('comms1')
    !    end if

    !    ! Get the halo information from neighbours
    !    ! For coarray comms, we need to sync before to ensure the information is sent, and
    !    ! also sync after to ensure it is not overwritten before it is used
    !    call md%recv_halos(sync_before=sync_before_recv, sync_after=sync_after_recv)

    !    if(send_boundary_flux_data .and. (.not. lagged_flux_correction_local)) then
    !        ! Do flux correction, and send the 'possibly corrected' data. 
    !        !
    !        ! Note that in many cases, this 'second send' might not be of practical importance.
    !        !    Indeed, even flux correction itself might often not practically matter, often.
    !        !    A good, cheap compromise might be to use md%lagged_flux_correction = .true., which
    !        !    should stop the mass errors from accumulating too much
    !        !
    !        ! FIXME: We might be able to reduce the amount of data sent (e.g. by only sending from domains
    !        !       that have had some correction applied, or by not sending boundary fluxes)
    !        do j = 1, size(md%domains)
    !            call md%domains(j)%nesting_flux_correction()
    !            call md%send_halos(domain_index=j, send_to_recv_buffer = send_halos_immediately)
    !        end do

    !        if(.not. send_halos_immediately) then
    !            ! Do all the coarray communication in one hit
    !            ! This lets us 'collapse' multiple sends to a single image,
    !            ! and is more efficient in a bunch of cases.
    !            TIMER_START('comms2')
    !            call communicate_p2p
    !            TIMER_STOP('comms2')
    !        end if

    !        ! Get the (possibly corrected) halo data.
    !        ! We need to sync before to ensure the information is sent, and
    !        ! also sync after to ensure it is not overwritten before it is used
    !        call md%recv_halos(sync_before=sync_before_recv, sync_after=sync_after_recv)

    !    end if

    !end subroutine
   
   

    !
    ! Given some x, y points, find a domain bounding box which contains them
    ! (if any). If multiple domains contain them, then pick the one with smallest
    ! cell area.
    !
    ! This is a useful workhorse routine for nesting setup
    !
    ! @param xs x-coordinates of test points
    ! @param ys y-coordinates of test points
    ! @param nbr_domain_ind integer array with the same length as xs. Will
    !   be filled with the index of the domain that contains the points, 
    !   or (-1) for points not inside any domain
    ! @param nbr_image_ind integer array with the same length as xs. Will
    !   be filled with the image_index of the domain that contains the points
    !   or (-1) for points not inside any domain
    ! @param nbr_domain_dx real rank2 array with first dimension the same
    !   length as xs. Will be filled with the dx values of the domain that
    !   contains the points, or a large negative number for points not in any domain
    ! @param error_domain_overlap_same_dx
    !
    subroutine find_priority_domain_containing_xy(xs, ys, nbr_domain_ind, &
        nbr_image_ind, nbr_domain_dx, error_domain_overlap_same_dx, &
        periodic_xs, periodic_ys)

        real(dp), intent(in) :: xs(:), ys(:)
        integer(ip), intent(inout) :: nbr_domain_ind(:), nbr_image_ind(:)
        real(dp), intent(inout) :: nbr_domain_dx(:,:)
        logical, intent(in) :: error_domain_overlap_same_dx
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)

        integer(ip) :: di, ii, xi, nd
        real(dp) :: xs1, ys1
        logical :: is_inside, periodic_point
        real(dp) :: nbr_dx(2), bbox(4,2), cell_area1, cell_area2, err_tol

        ! Tolerence when checking for equality of neighbouring domain dx values.
        ! Since all dx should attain a discrete set of values that
        ! differ by at least a factor of 2 between them, a relatively crude eps
        ! is fine.
        real(dp), parameter :: eps = 1.0e-02_dp

        ! Default values for 'not contained by any other bbox' points
        nbr_domain_ind = -1_ip
        nbr_image_ind = -1_ip
        ! Default neighbour cell area = very large negative number
        nbr_domain_dx = -sqrt(HUGE(1.0_dp) - 1.0e+06_dp)

        !
        ! Loop over all bounding boxes and find which one xs, ys are inside.
        ! If they are inside multiple, then give preference to the ones
        ! which have the smallest cell area
        !
        !
        !ni = size(all_bbox(1,1,1,:))
        nd = size(all_bbox, 3)
        do ii = 1, ni
            do di = 1, nd

                ! Properties of neighbour domain
                bbox = all_bbox(:, :, di, ii)
                nbr_dx = all_dx(1:2, di, ii)

                ! Quick exit -- by default all_bbox = 0.0_dp
                if(all(bbox(:,1) == 0.0_dp) .and. &
                    all(bbox(:,2) == 0.0_dp)) cycle

                do xi = 1, size(xs) 
                    ! Find whether xs(xi), ys(x) is inside the di'th domain on
                    ! image ii

                    xs1 = xs(xi)
                    ys1 = ys(xi)

                    ! If xs1 or ys1 are outside the periodic boundaries, deal with it
                    call check_periodic(xs1, ys1, periodic_xs, periodic_ys, periodic_point, adjust_coordinates=.TRUE.)

                    call point_in_poly(4_ip, bbox(:, 1), bbox(:,2), xs1, ys1, is_inside)

                    if(is_inside) then
                        !
                        ! We want to take the point from the highest-resolution
                        ! domain -- i.e. smallest cell 'area'. 
                        !
                        ! Note that the computation of cell_area1 and 
                        ! cell_area (product(nbr_dx)) may not lead to physical
                        ! areas [e.g. spherical coordinates] but because all
                        ! domains have dx/dy as some integer multiple of the
                        ! parent dx/dy, the relations '>' '=' '<' will be
                        ! preserved.
                        !
                        cell_area1 = product(nbr_dx)
                        cell_area2 = product(nbr_domain_dx(xi,1:2))

                        if(cell_area1 < cell_area2) then
                            ! Update priority domain
                            nbr_image_ind(xi) = ii
                            nbr_domain_ind(xi) = di
                            nbr_domain_dx(xi, 1:2) = nbr_dx
                        else
                            ! In this case we expect cell_area2 < cell_area1
                            !
                            ! If cell2 == cell_area1, then we have some error,
                            ! since the multi-domain is not allowed overlapping
                            ! areas with the same resolution.
                            !

                            ! Check they are not equal to within err_tol
                            err_tol = eps*min(cell_area1, cell_area2)

                            if(error_domain_overlap_same_dx .and. &
                                (abs(cell_area1 - cell_area2) < err_tol)) then

                                !! Update priority domain
                                !nbr_image_ind(xi) = ii
                                !nbr_domain_ind(xi) = di
                                !nbr_domain_dx(xi, 1:2) = nbr_dx

                                write(log_output_unit,*) 'Overlapping domains with the same', &
                                    ' cell size exist around: ', &
                                    xs(xi), ys(xi), &
                                    '. This is not permitted.', &
                                    cell_area1, cell_area2,&
                                    abs(cell_area1 - cell_area2), &
                                    sqrt(cell_area1)

                                write(log_output_unit,*) 'ii: ', ii, ' di: ', di
                                write(log_output_unit,*) 'xy: ', xs(xi), ys(xi)
                                write(log_output_unit,*) 'bbox: ' 
                                write(log_output_unit,*) bbox(1,:)
                                write(log_output_unit,*) bbox(2,:)
                                write(log_output_unit,*) bbox(3,:)
                                write(log_output_unit,*) bbox(4,:)
                                stop

                            end if
                        end if
                    end if
                end do

            end do
        end do

    end subroutine

    ! Given a domain, and a given boundary (N/E/S/W), find which other domains
    ! each point along its boundary should nest with.
    !
    ! @param domain the domain
    ! @param boundary_flag integer flag [1,2,3,4] corresponding to [N,E,S,W]
    !     boundary
    ! @param nest_layer_width Number of cells in the nesting layer
    ! @param nbr_domain_ind rank 1 array with indices of the neighbour domain
    !     we nest with
    ! @param nbr_image_ind rank 1 array with image_indices holding the
    !     neighbour domain we nest with
    ! @param nbr_domain_dx rank 2 array with dx(1:2) values for the domains we
    !     nest with
    !
    subroutine find_priority_nesting_domains_along_a_boundary(&
        domain, &
        boundary_flag, &
        nest_layer_width, &
        nbr_domain_ind, &
        nbr_image_ind, &
        nbr_domain_dx,&
        periodic_xs,&
        periodic_ys)

        type(domain_type), intent(in) :: domain
        integer(ip), intent(in) :: boundary_flag, nest_layer_width
        integer(ip), allocatable, intent(out):: nbr_domain_ind(:), &
            nbr_image_ind(:)
        real(dp), allocatable, intent(out):: nbr_domain_dx(:,:)
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)

        ! Local vars
        real(dp), allocatable :: domain_xs(:), domain_ys(:), xs(:), ys(:)

        integer(ip) :: i, n1, n2
        real(dp) :: dx, dy

        call get_domain_xs_ys(domain, domain_xs, domain_ys) 

        ! Make points just outside boundaries
        select case(boundary_flag)
            case(1) 
                ! Make row of points just outside the north boundary, extending
                ! 'nest_layer_width' over both edges
                n1 = domain%nx(1)
                n2 = n1 + 2 * nest_layer_width
                dx = domain%dx(1)
                dy = domain%dx(2)
                allocate(xs(n2), ys(n2))
                xs((nest_layer_width+1):(n2-nest_layer_width)) = domain_xs
                do i = 1, nest_layer_width
                    xs(i) = domain_xs(1) - (nest_layer_width - (i-1))*dx
                    xs(n2 - (nest_layer_width - i)) = domain_xs(n1) + i * dx
                end do
                ys = domain_ys(domain%nx(2)) + dy * nest_layer_width

            case(2) 
                ! Make row of points just outside the east boundary, extending
                ! 'nest_layer_width' over both edges
                n1 = domain%nx(2)
                n2 = n1 + 2 * nest_layer_width
                dx = domain%dx(1)
                dy = domain%dx(2)
                allocate(xs(n2), ys(n2))
                ys((nest_layer_width+1):(n2-nest_layer_width)) = domain_ys
                do i = 1, nest_layer_width
                    ys(i) = domain_ys(1) - (nest_layer_width - (i-1))*dy
                    ys(n2 - (nest_layer_width - i)) = domain_ys(n1) + i * dy
                end do
                xs = domain_xs(domain%nx(1)) + dx * nest_layer_width

            case(3) 
                ! Make row of points just outside the south boundary, extending
                ! 'nest_layer_width' over both edges
                n1 = domain%nx(1)
                n2 = n1 + 2 * nest_layer_width
                dx = domain%dx(1)
                dy = domain%dx(2)
                allocate(xs(n2), ys(n2))
                xs((nest_layer_width+1):(n2-nest_layer_width)) = domain_xs
                do i = 1, nest_layer_width
                    xs(i) = domain_xs(1) - (nest_layer_width - (i-1))*dx
                    xs(n2 - (nest_layer_width - i)) = domain_xs(n1) + i * dx
                end do
                ys = domain_ys(1) - dy * nest_layer_width

            case(4) 
                ! Make row of points just outside the west boundary, extending
                ! 'nest_layer_width' over both edges
                n1 = domain%nx(2)
                n2 = n1 + 2 * nest_layer_width
                dx = domain%dx(1)
                dy = domain%dx(2)
                allocate(xs(n2), ys(n2))
                ys((nest_layer_width+1):(n2-nest_layer_width)) = domain_ys
                do i = 1, nest_layer_width
                    ys(i) = domain_ys(1) - (nest_layer_width - (i-1))*dy
                    ys(n2 - (nest_layer_width - i)) = domain_ys(n1) + i * dy
                end do
                xs = domain_xs(1) - dx * nest_layer_width

            case default
                stop "Bad boundary_flag value"
        end select

        ! Store neighbour info
        if(allocated(nbr_domain_ind)) deallocate(nbr_domain_ind)
        if(allocated(nbr_image_ind)) deallocate(nbr_image_ind)
        if(allocated(nbr_domain_dx)) deallocate(nbr_domain_dx)

        ! Make space to record the domain/image index we nest with
        allocate(nbr_domain_ind(size(xs)), &
            nbr_image_ind(size(xs)), &
            nbr_domain_dx(size(xs), 2))

        call find_priority_domain_containing_xy(xs, ys, nbr_domain_ind, &
            nbr_image_ind, nbr_domain_dx, error_domain_overlap_same_dx=.true.,&
            periodic_xs=periodic_xs, periodic_ys=periodic_ys)

    end subroutine


    !
    ! Determine the required thickness of cells in the nesting region,
    !  in which the domain receives data. 
    ! It has to be thick enough so we can take one_evolve_step() without
    !  the full nesting_region information propagating into the interior. 
    ! It also has to be an integer multiple of the ratio between the coarsest
    !  parent domain's dx size and the current domain's dx -- so that the 
    !  nesting region represents 'full' coarse cells
    !
    function get_domain_nesting_layer_thickness(domain, &
        max_parent_dx_ratio, extra_halo_buffer, extra_cells_in_halo) result(thickness)

        type(domain_type), intent(in) :: domain
        real(dp), intent(in) :: max_parent_dx_ratio
        integer(ip), intent(in) :: extra_halo_buffer, extra_cells_in_halo
        integer(ip) :: thickness

        integer(ip) :: required_cells_ts_method

        ! This is useful for checking the effect of 'extra' padding (e.g. for nesting-debug)
        !! e.g. if we underestimate cells required e.g. due to use of gradients in nesting,
        !! increasing this can catch it.
        !integer(ip), parameter :: extra_cells = extra_cells_in_halo_default

        ! Find the minimum nesting layer thickness required by the
        ! time-stepping algorithm, assuming we only take one evolve_one_step
        ! in-between nesting updates
        select case(domain%timestepping_method)
            case('euler')
                required_cells_ts_method = 2
            case('rk2')
                required_cells_ts_method = 4
            case('rk2n')
                required_cells_ts_method = 10

                ! Deliberately fail for now.
                write(log_output_unit, *) " ERROR: Currently rk2n doesn't have the mass-conservation-tracking ", &
                    "infrastucture for multidomain. Until implemented, we deliberately STOP."
                call generic_stop()

            case('midpoint')
                required_cells_ts_method = 4
            case('linear')
                required_cells_ts_method = 2
            case default
                write(log_output_unit,*) 'timestepping_method not recognized'
                error stop
        end select

        ! Now multiply by the number of evolve_one_step calls actually taken
        thickness = (required_cells_ts_method ) * &
            domain%timestepping_refinement_factor + extra_cells_in_halo

        if(nint(max_parent_dx_ratio) > 1) then
            ! Now round up to be a multiple of the parent cell-sizes
            thickness = (ceiling(thickness*1.0_dp/max_parent_dx_ratio) + extra_halo_buffer) * nint(max_parent_dx_ratio)
        else
            ! Now round up to be a multiple of the parent cell-sizes
            thickness = (ceiling(thickness*1.0_dp/max_parent_dx_ratio)                    ) * nint(max_parent_dx_ratio)
        end if
     
    end function

    ! 
    ! Compute the coordinate x/y arrays associated with a given domain
    !
    ! The domain is in this subroutine is typically used prior to a call to domain%allocate_quantities()
    ! Thus it is missing the x/y coordinates. Sometimes we need the x and y coordinates implied by the
    ! domain, so we get them here.
    !
    subroutine get_domain_xs_ys(domain, xs, ys)

        type(domain_type), intent(in) :: domain
        real(dp), allocatable, intent(inout) :: xs(:), ys(:)

        integer(ip) :: i

            if(allocated(xs)) deallocate(xs)
            if(allocated(ys)) deallocate(ys)

            ! Compute xs, ys for the domain. These are currently not stored, so
            ! we directly calculate them from the metadata
            allocate(xs(domain%nx(1)), ys(domain%nx(2)))

            do i = 1, domain%nx(1)
                xs(i) = domain%lower_left(1) + &
                    (i - 0.5_dp)*domain%lw(1) / domain%nx(1)
            end do
            do i = 1, domain%nx(2)
                ys(i) = domain%lower_left(2) + &
                    (i - 0.5_dp)*domain%lw(2) / domain%nx(2)
            end do

    end subroutine


    ! Compute volume in 'active areas' for all domains (or optionally, in just
    ! one domain). This also updates the volume data stored in md
    !
    ! @param md multidomain
    ! @param vol output goes here
    ! @param domain_index optional index of only one domain for which we want the volume. If
    !     not present, sum over all domains
    !
    subroutine get_flow_volume(md, vol, domain_index)
        class(multidomain_type), intent(inout) :: md
        real(force_double), intent(inout) :: vol
        integer(ip), intent(in), optional :: domain_index

        integer(ip) :: i, j, k, n_ext, ny, nx, kstart, kend
        real(force_double) :: vol_k, local_sum

        ! If domain_index is provided, then only get volume in that domain
        if(present(domain_index)) then
            kstart = domain_index
            kend = domain_index
        else
            kstart = 1
            kend = size(md%domains)
        end if

        vol = 0.0_force_double

        do k = kstart, kend 
   
            vol_k = 0.0_force_long_double

            n_ext = md%domains(k)%exterior_cells_width
            ny = md%domains(k)%nx(2)
            nx = md%domains(k)%nx(1)
            ! Integrate over the area with this domain = priority domain, which is not in periodic regions.
            ! Be careful about precision ( see calls to real(...., force_double) )
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(md, n_ext, nx, ny, ti, k) REDUCTION(+:vol_k)
            vol_k = 0.0_force_long_double
            !$OMP DO SCHEDULE(STATIC)
            do j = (n_ext+1), ny - n_ext
                local_sum = sum(&
                        real(md%domains(k)%U( (n_ext+1):(nx-n_ext) ,j ,STG), force_long_double) - &
                        real(md%domains(k)%U( (n_ext+1):(nx-n_ext) ,j ,ELV), force_long_double ), &
                    mask=(&
                        (md%domains(k)%nesting%priority_domain_index( (n_ext+1):(nx-n_ext) ,j) == k) .and. &
                        (md%domains(k)%nesting%priority_domain_image( (n_ext+1):(nx-n_ext) ,j) == ti) .and. &
                        ! Avoid periodic regions (because they might not be excluded by the previous criteria)
                        (md%domains(k)%x( (n_ext+1):(nx-n_ext)) > md%periodic_xs(1)) .and.&
                        (md%domains(k)%x( (n_ext+1):(nx-n_ext)) < md%periodic_xs(2))))
                ! Avoid periodic regions 
                if(md%domains(k)%y(j) < md%periodic_ys(1) .or. md%domains(k)%y(j) > md%periodic_ys(2)) local_sum = 0.0_dp
                vol_k = vol_k + real(md%domains(k)%area_cell_y(j) * local_sum, force_long_double)
            end do 
            !$OMP END DO
            !$OMP END PARALLEL

            md%volume(k) = vol_k
            vol = vol + vol_k

        end do

    end subroutine


    ! Make sure initial flow/elev values are consistent between parallel domains, by doing one halo exchange
    !
    ! If separate_halos = .true., then call md%separate_halos() to make recv_weights zero in parts of fine-domain nesting regions
    ! that are close to other nesting regions
    !
    subroutine make_initial_conditions_consistent(md, separate_halos)
        class(multidomain_type), intent(inout) :: md
        logical, optional, intent(in) :: separate_halos

        logical :: separate_halos_local

        if(present(separate_halos)) then
            separate_halos_local = separate_halos
        else
            separate_halos_local = .false.
        end if

        call md%send_halos()
#ifdef COARRAY
        sync all
#endif
        call md%recv_halos()

        if(separate_halos_local) then
            call md%separate_halos()
        end if

        ! For mass conservation tracking it is good to do this
        call md%record_initial_volume()

    end subroutine

    subroutine memory_summary(md) 
        class(multidomain_type), intent(in) :: md

        integer(int64) :: domain_U_size, nesting_buffer_size, big_tmp
        integer(ip) ::  i, p2p_buffer_size, tmp
        real(dp) :: buffer_size_on_domain_U_size
        

        domain_U_size = 0
        nesting_buffer_size = 0
        tmp = 0

        do i = 1, size(md%domains)
            ! This can overflow if we do not force it to be a 'long' integer!
            big_tmp = 1_int64 * md%domains(i)%nx(1) * md%domains(i)%nx(2) * &
                md%domains(i)%nvar * real_bytes
            domain_U_size = domain_U_size + big_tmp
            call md%domains(i)%nesting%memory_size(tmp)
            nesting_buffer_size = nesting_buffer_size + tmp
        end do 

        p2p_buffer_size = size_of_send_recv_buffers()

        buffer_size_on_domain_U_size = ((p2p_buffer_size+nesting_buffer_size)*1.0_dp)/(domain_U_size*1.0_dp)

        write(log_output_unit, *) 'Multidomain has ', domain_U_size, ' U bytes and ', &
            p2p_buffer_size , ' p2p bytes and ', nesting_buffer_size , &
            ' nesting buffer bytes , ratio= ', buffer_size_on_domain_U_size
        write(log_output_unit, *) 'Nonlinear domains will include other large arrays ', &
            '(e.g. for fluxes) with total size larger than U, which are not accounted for above'

    end subroutine

    subroutine report_mass_conservation_statistics(md)
        class(multidomain_type), intent(inout) :: md

        real(force_double) :: vol, vol0, bfi, dvol, vol_and_bfi(4)
        integer(ip) :: i

            ! Compute the volume in the domain
            vol0 = sum(md%volume_initial)
            !vol = md%get_flow_volume()
            call md%get_flow_volume(vol)

            dvol = sum(md%volume - md%volume_initial)

            ! Compute the 'exterior' boundary flux integral
            bfi = 0.0_dp
            do i = 1, size(md%domains)
                bfi = bfi + md%domains(i)%boundary_flux_time_integral_exterior
            end do

#ifdef COARRAY
            ! Sum volume and flux-integral over all cores. Pack into an array, so
            ! we only need one collective call
            vol_and_bfi(1) = vol
            vol_and_bfi(2) = bfi
            vol_and_bfi(3) = vol0
            vol_and_bfi(4) = dvol
#ifndef COARRAY_PROVIDE_CO_ROUTINES
            call co_sum(vol_and_bfi) 
#else
            ! The MPI-replacement version of co_sum does not operate on a
            ! vector
            do i = 1, 4
                call co_sum(vol_and_bfi(i)) 
            end do
#endif
            vol = vol_and_bfi(1)
            bfi = vol_and_bfi(2)
            vol0 = vol_and_bfi(3)
            dvol = vol_and_bfi(4) 
#endif
            write(log_output_unit,*) 'Volume statistics (m^3) integrated over all domains and images:'
            write(log_output_unit,*) '  Multidomain volume       : ', vol
            write(log_output_unit,*) '              volume change: ', dvol
            write(log_output_unit,*) '     boundary flux integral: ', bfi
            write(log_output_unit,*) '         unexplained change: ', dvol + bfi

    end subroutine

    !
    ! Find if cell 'i,j' on domains(myindex) on 'this_image() == myimage'
    ! is within the active region
    !
    ! The 'active region' contains cells which we need to get from other images, AND
    ! cells which have priority domain = [ domains(myindex) on 'this_image() == myimage'],
    ! (although in the latter case, no communication is required, except if using periodic boundaries). 
    ! An example of cells which are not in the 'active region' is coarse domain cells 
    ! which are completely covered by finer domains (and not within the communication buffers).
    ! Such 'inactive' cells do not affect the solution in priority domains
    !
    ! @param priority_domain_index gives index of priority domain for all cells
    !     on [ domains(myindex) on 'this_image() == myimage']
    ! @param priority_domain_image gives image index of priority domain for all
    !     cells on [ domains(myindex) on 'this_image() == myimage']
    ! @param myindex see above
    ! @param myimage see above
    ! @param halo_width the width of the communication buffer
    ! @param i x-index to be queried
    ! @param j y-index to be queried
    ! @param is_in_active output (logical)
    !
    pure subroutine is_in_active_region(priority_domain_index, priority_domain_image, &
        myindex, myimage, halo_width, i, j, is_in_active)

        integer(ip), intent(in) :: priority_domain_index(:,:), priority_domain_image(:,:)        
        integer(ip), intent(in) :: myindex, myimage, halo_width, i, j
        logical, intent(inout) :: is_in_active

        integer(ip) :: ii, jj, dims(2)

        dims = shape(priority_domain_index)

        is_in_active = .FALSE.

        do jj = max(1, j-halo_width), min(dims(2), j + halo_width)
            do ii = max(1, i-halo_width), min(dims(1), i + halo_width)

                ! If ii,jj is within the original domain, then quick exit
                if(priority_domain_index(ii,jj) == myindex .and. &
                   priority_domain_image(ii,jj) == myimage) then
                    is_in_active = .TRUE.
                    return
                end if

            end do
        end do

    end subroutine

    !
    ! This routine can be used to partition domains among coarray images
    ! Basically we move the provided md%domains to md%domain_metadata, and
    ! make a new md%domains, with the right piece of the domain on each image
    !
    subroutine partition_domains(md)
        class(multidomain_type), intent(inout) :: md

        integer(ip) :: nd, next_d, i, j, ni, ti, ii, i0, i1
        integer(ip) :: local_ti, local_ni, local_co_size_xy(2), local_co_index(2)
        integer(ip) :: domain_nx(2), nx(2), lower_left_nx(2), upper_right_nx(2) 
        integer(ip) :: domain_dx_refinement_factor(2), dx_refine_X_co_size_xy(2)
        integer(ip) :: parent_domain_dx_refinement_factor(2), dx_refine_ratio(2)
        integer(ip) :: best_co_size_xy(2), fileID, extra_halo_padding(2)
        real(dp) :: nboundary, best_nboundary

#ifndef COARRAY
        ti = 1
        ni = 1
#else
        ! Coarray case, we split the images
        ti = this_image()
        ni = num_images()
#endif

        ! Name domains as " image_index * 1e+10 + original_domain_index "
        ! The naming will be implemented below -- for now, check that the naming will not overflow.
        if(ni * large_64_int + size(md%domains) > HUGE(md%domains(1)%myid)) then
            write(log_output_unit, *) 'The number of images and local domains is too high for the naming convention'
            call generic_stop()
        end if

        ! Clear md%domains
        call move_alloc(md%domains, md%domain_metadata)

        ! Read load balancing info if provided
        if(md%load_balance_file /= '') then

            call read_ragged_array_2d_ip(md%load_balance_part, md%load_balance_file)

            if(size(md%load_balance_part%i2) /= size(md%domains)) &
                stop 'Rows in load_balance_file not equal to number of domains '

            ! No values < 1
            do i = 1, size(md%load_balance_part%i2)
                if( minval(md%load_balance_part%i2(i)%i1) < 1 ) then
                    stop 'load_balance_file contains numbers < 1'
                end if
            end do

            ! If there are values > ni, convert to the range 1-ni, with a warning
            do i = 1, size(md%load_balance_part%i2)
                if( maxval(md%load_balance_part%i2(i)%i1) > ni ) then
                    write(log_output_unit, *) &
                        ' WARNING: load_balance_file contains values > num_images: converting to smaller values.'
                    md%load_balance_part%i2(i)%i1 = mod((md%load_balance_part%i2(i)%i1 - 1), ni) + 1_ip
                end if
            end do

            ! Final check
            do i = 1, size(md%load_balance_part%i2)
                if( minval(md%load_balance_part%i2(i)%i1) < 1 .or.  &
                    maxval(md%load_balance_part%i2(i)%i1) > ni) then
                    stop 'hit uncorrectable error in md%load_balance_part'
                end if
            end do
            !! DEBUG -- print out the file
            !if(this_image() == 1) then
            !    do i = 1, size(md%domain_metadata)
            !        print*, 'Row ', i, ':', md%load_balance_partition(i,:)
            !    end do
            !end if

        else

            ! By default, spread all images over all domains
            allocate(md%load_balance_part%i2(size(md%domain_metadata)))
            do i = 1, size(md%domain_metadata)
                allocate(md%load_balance_part%i2(i)%i1(ni))
                do j = 1, ni
                    md%load_balance_part%i2(i)%i1(j) = j
                end do
            end do

        end if

        ! Count how many domains we need on the current image
        nd = 0
        do i = 1, size(md%load_balance_part%i2)
            nd = nd + count( md%load_balance_part%i2(i)%i1 == ti )
        end do
        ! Make space for nd domains
        allocate(md%domains(nd))


        ! Split the originally provided domains 
        next_d = 0
        do i = 1, size(md%domain_metadata)
            ! Number of pieces that we split the i'th domain into
            local_ni = size(md%load_balance_part%i2(i)%i1)

            ! Loop over each piece
            do j = 1, local_ni

                ! If this piece is not assigned to the current image, do nothing
                if(ti /= md%load_balance_part%i2(i)%i1(j)) cycle 

                ! Index of the 'piece' of the i'th domain we are working on
                local_ti = j

                ! Index in md%domains where we put this domain
                next_d = next_d + 1

                ! Key info about the domain that will be split
                domain_nx = md%domain_metadata(i)%nx
                domain_dx_refinement_factor = int(nint(md%domain_metadata(i)%dx_refinement_factor), ip)
                parent_domain_dx_refinement_factor = int(nint(md%domain_metadata(i)%dx_refinement_factor_of_parent_domain), ip)

                !write(log_output_unit, *) 'domain(', i, '), dx_refinement_factor = ', domain_dx_refinement_factor, &
                !    ' parent_domain_dx_refinement_factor = ', parent_domain_dx_refinement_factor


                ! Split domain into (nx , ny) tiles where " nx*ny == local_ni "
                ! Approach: Loop over all possible combinations of
                ! local_co_size_xy(1:2) that have product = local_ni.
                ! Find the 'best' one (i.e. smallest total boundary to minimise comms)
                best_nboundary = HUGE(1.0_dp)
                best_co_size_xy = -1
                do ii = local_ni, 1, -1

                    local_co_size_xy(1) = ii
                    local_co_size_xy(2) = local_ni/local_co_size_xy(1) ! Deliberate integer division

                    if(product(local_co_size_xy) /= local_ni) cycle

                    ! Approximate total number of 'boundary cells' in all tiles
                    nboundary = 2 * local_ni * sum(domain_nx/(1.0_dp*local_co_size_xy)) 
                    if(nboundary < best_nboundary) then
                        ! We have a new 'best' local_co_size_xy
                        best_co_size_xy = local_co_size_xy 
                        best_nboundary = nboundary
                    end if
                end do
                local_co_size_xy = best_co_size_xy
                ! Check we got some sensible result
                if(any(local_co_size_xy < 0)) then
                    write(log_output_unit,*) 'Error in splitting up images into 2d', local_co_size_xy, local_ni
                    call generic_stop
                end if

                ! Get the 'co-index' equivalent for the current image
                local_co_index(1) = mod(local_ti, local_co_size_xy(1))
                if(local_co_index(1) == 0) local_co_index(1) = local_co_size_xy(1)
                local_co_index(2) = (local_ti-local_co_index(1))/local_co_size_xy(1) + 1 ! Deliberate integer division

                !write(log_output_unit,*) 'ti :', ti, ' local_ti :', local_ti, ' local_ni: ', local_ni, &
                !    ' local_co_size_xy: ', local_co_size_xy, ' local_co_index:', local_co_index    


                ! Check that the domain is large enough to be split
                !dx_refine_ratio = domain_dx_refinement_factor ! Assuming the domain's parent is the global domain
                dx_refine_ratio = ((domain_dx_refinement_factor/parent_domain_dx_refinement_factor))
                if(.not. all(dx_refine_ratio * parent_domain_dx_refinement_factor == domain_dx_refinement_factor)) then
                    write(log_output_unit,*) 'dx_refine_ratio should be an exact integer upon creation'
                    call generic_stop()
                end if
                dx_refine_X_co_size_xy = dx_refine_ratio * local_co_size_xy
                if(.not. all( (domain_nx/(dx_refine_X_co_size_xy)) > 1)) then
                    write(log_output_unit, *) 'Cannot split domain ', i, ' with nx=', domain_nx, &
                        ' into local_co_size_xy=', local_co_size_xy, ' when domain_dx_refinement_factor=', &
                        domain_dx_refinement_factor, ' and dx_refine_ratio = ', dx_refine_ratio
                    write(log_output_unit,*) 'Consider enlarging the domain, or specifying the decomposition more thoroughly'
                    call generic_stop()
                end if

                ! Find the boundaries of the tile
                ! Compute the lower-left/upper-right cells in the domain that will be assigned to this image
                ! The splitted domains must all have edges that align perfectly with the parent domains
                ! We exploit integer division here to ensure that
                lower_left_nx  = (((local_co_index - 1) * domain_nx)/(dx_refine_X_co_size_xy)) * &
                    dx_refine_ratio + 1 ! Deliberate integer division
                upper_right_nx = (((local_co_index    ) * domain_nx)/(dx_refine_X_co_size_xy)) * &
                    dx_refine_ratio ! Deliberate integer division

                ! Now we need to set variables, like
                ! md%domains(1)%lw = global_lw
                md%domains(next_d)%lw = ((upper_right_nx - lower_left_nx + 1)*&
                    real(md%domain_metadata(i)%lw, force_long_double))/ & 
                    md%domain_metadata(i)%nx
                ! md%domains(1)%lower_left =global_ll
                md%domains(next_d)%lower_left = real(md%domain_metadata(i)%lower_left, force_long_double) + &
                    real( ((lower_left_nx - 1)*real(md%domain_metadata(i)%lw, force_long_double))&
                    /md%domain_metadata(i)%nx, force_long_double)
                ! md%domains(1)%nx = global_nx
                md%domains(next_d)%nx = upper_right_nx - lower_left_nx + 1
                ! md%domains(1)%dx = md%domains(1)%lw/md%domains(1)%nx
                md%domains(next_d)%dx = md%domain_metadata(i)%dx
                ! md%domains(1)%timestepping_method = 'linear'
                md%domains(next_d)%timestepping_method = md%domain_metadata(i)%timestepping_method
                md%domains(next_d)%linear_solver_is_truely_linear = md%domain_metadata(i)%linear_solver_is_truely_linear
                ! md%domains(1)%timestepping_refinement_factor = 1_ip
                md%domains(next_d)%timestepping_refinement_factor = md%domain_metadata(i)%timestepping_refinement_factor
                ! md%domains(1)%dx_refinement_factor = 1.0_dp
                md%domains(next_d)%dx_refinement_factor = md%domain_metadata(i)%dx_refinement_factor
                md%domains(next_d)%dx_refinement_factor_of_parent_domain = &
                    md%domain_metadata(i)%dx_refinement_factor_of_parent_domain

                ! Copy over the netcdf down-sampling information
                md%domains(next_d)%nc_grid_output%spatial_stride = md%domain_metadata(i)%nc_grid_output%spatial_stride
                ! Let the nc_grid_output know it is part of a larger domain with the given lower-left
                md%domains(next_d)%nc_grid_output%spatial_ll_full_domain = md%domain_metadata(i)%lower_left

                ! Time before which we do trivial evolves
                md%domains(next_d)%static_before_time = md%domain_metadata(i)%static_before_time

                ! CFL and theta
                md%domains(next_d)%cfl = md%domain_metadata(i)%cfl
                md%domains(next_d)%theta = md%domain_metadata(i)%theta

                ! NAME THE DOMAIN, using the "original" domain index, rather than the one
                ! after partitioning
                md%domains(next_d)%myid = ti * large_64_int + i 
                md%domains(next_d)%local_index = local_ti

                write(log_output_unit,*) 'i: ', i, &
                    ' local_ti: ', local_ti, ' local_ni: ', local_ni, ' ti: ', ti, ' ni: ', ni, &
                    ' lower_left_nx: ', lower_left_nx, &
                    ' upper_right_nx: ', upper_right_nx, ' lower-left: ', md%domains(next_d)%lower_left, &
                    ' dx: ', md%domains(next_d)%dx, ' lw: ', md%domains(next_d)%lw
            end do
        end do 
#ifdef COARRAY
        sync all
#endif
        flush(log_output_unit)

    end subroutine


    !
    ! Convenience print routine
    !
    subroutine print_multidomain(md, global_stats_only)
        class(multidomain_type), intent(inout) :: md ! inout allows for timer to run
        logical, optional, intent(in) :: global_stats_only

        integer(ip) :: i, j, k
        real(dp) :: minstage, maxstage, minspeed, maxspeed, stg1, speed1, depth1
        logical :: is_nesting, only_global_stats
        real(dp) :: global_max_stage, global_min_stage, global_max_speed, global_min_speed
        real(dp) :: wrk(2)

        if(present(global_stats_only)) then
            only_global_stats = global_stats_only
        else
            only_global_stats = .false.
        end if


        TIMER_START('printing_stats')

        ! Mark out the new time-step
        write(log_output_unit, *) ''
        write(log_output_unit, *) '###################################'

        ! Variables to track extremes over all domains
        global_max_stage = -HUGE(1.0_dp)
        global_min_stage = HUGE(1.0_dp)
        global_max_speed = 0.0_dp 
        global_min_speed = 0.0_dp

        do k = 1, size(md%domains)

            if(.not. only_global_stats) then
                write(log_output_unit,*) ''
                write(log_output_unit,*) '-----------'
                write(log_output_unit,*) 'domain ', k
            end if

            ! Reset stats of interest
            maxstage = -HUGE(1.0_dp)
            maxspeed = 0.0_dp ! speed is always > 0
            minstage = HUGE(1.0_dp)
            minspeed = 0.0_dp ! speed is always > 0

            ! If nesting is occurring, the above stats are only computed
            ! in the priority domain area
            if(md%domains(k)%nesting%my_index > 0) then
                is_nesting = .true.
            else
                is_nesting = .false.
            end if

            ! Compute stats of interest
            do j = 2, md%domains(k)%nx(2) - 1
                do i = 2, md%domains(k)%nx(1) - 1
    
                    if(.not. is_nesting) cycle

                    ! Only act inside the priority domain
                    if( (md%domains(k)%nesting%priority_domain_index(i,j) /= md%domains(k)%nesting%my_index) .or. &
                        (md%domains(k)%nesting%priority_domain_image(i,j) /= md%domains(k)%nesting%my_image)) cycle

                    ! Actual stats    
                    if( (md%domains(k)%U(i,j,STG) > md%domains(k)%U(i,j,ELV) + minimum_allowed_depth)) then
                        stg1 = md%domains(k)%U(i,j,STG)
                        depth1 = stg1 - md%domains(k)%U(i,j,ELV)
        
                        maxstage = max(maxstage, stg1)
                        minstage = min(minstage, stg1)

                        ! Compute 'velocity**2', and update later 
                        speed1 = (md%domains(k)%U(i,j,UH) * md%domains(k)%U(i,j,UH) + &
                                  md%domains(k)%U(i,j,VH) * md%domains(k)%U(i,j,VH)   )/ &
                                 (depth1*depth1)

                        maxspeed = max(maxspeed, speed1)
                        minspeed = min(minspeed, speed1)
                    end if

                end do
            end do

            maxspeed = sqrt(maxspeed)
            minspeed = sqrt(minspeed)
            global_max_stage = max(global_max_stage, maxstage)
            global_min_stage = min(global_min_stage, minstage)
            global_max_speed = max(global_max_speed, maxspeed)
            global_min_speed = min(global_min_speed, minspeed)

            if(only_global_stats) cycle

            ! Print main statistics
            write(log_output_unit, *) ''
            write(log_output_unit, *) 'Domain ID: '
            write(log_output_unit, *) '        ', md%domains(k)%myid
            write(log_output_unit, *) 'Time: '
            write(log_output_unit, *) '        ', md%domains(k)%time
            write(log_output_unit, *) 'nsteps_advanced:'
            write(log_output_unit, *) '        ', md%domains(k)%nsteps_advanced
            write(log_output_unit, *) 'dt: '
            write(log_output_unit, *) '        ', md%domains(k)%max_dt
            write(log_output_unit, *) 'evolve_step_dt: '
            write(log_output_unit, *) '        ', md%domains(k)%evolve_step_dt
            write(log_output_unit, *) 'Stage: '
            write(log_output_unit, *) '        ', maxstage
            write(log_output_unit, *) '        ', minstage
            write(log_output_unit, *) 'Speed: '
            write(log_output_unit, *) '        ', maxspeed
            write(log_output_unit, *) '        ', minspeed
            write(log_output_unit, *) 'Negative_depth_clip_counter: '
            write(log_output_unit, *) '        ', md%domains(k)%negative_depth_fix_counter

        end do

        ! Even if reporting global stats only, we'll still want to know the time
        if(only_global_stats) then
            write(log_output_unit, *) 'Time: '
            write(log_output_unit, *) '        ', md%domains(1)%time
        end if

#ifdef COARRAY
        call co_max(global_max_stage)
        call co_max(global_max_speed)
        call co_min(global_min_stage)
        call co_min(global_min_speed)
#endif

        write(log_output_unit, *) ''
        write(log_output_unit,*) '-----------'
        write(log_output_unit, *) 'Global stage range (over all domains and images): '
        write(log_output_unit, *) '        ', global_max_stage
        write(log_output_unit, *) '        ', global_min_stage
        write(log_output_unit, *) 'Global speed range (over all domains and images) -- ', &
            'note linear domains can have very high velocities, even when stable: '
        write(log_output_unit, *) '        ', global_max_speed
        write(log_output_unit, *) '        ', global_min_speed
        write(log_output_unit,*) '-----------'
        call md%report_mass_conservation_statistics()


        TIMER_STOP('printing_stats')

    end subroutine

    !
    ! Convenience routine to get nested grid halo data from communication buffer in a multidomain
    !
    subroutine receive_multidomain_halos(md, sync_before, sync_after)
        class(multidomain_type), intent(inout) :: md
        logical, optional, intent(in) :: sync_before, sync_after

        integer(ip) :: i, j, ii, jj, iL, iU, jL, jU, ip1, jp1
        real(dp) :: elev_lim
        logical :: sync_before_local, sync_after_local

        ! By default do not sync before or after
        if(present(sync_before)) then
            sync_before_local = sync_before
        else
            sync_before_local = .false.
        end if

        if(present(sync_after)) then
            sync_after_local = sync_after
        else
            sync_after_local = .false.
        end if

#ifdef COARRAY
        if(sync_before_local) then
            TIMER_START('sync_before_recv')
            sync images(linked_p2p_images)
            TIMER_STOP('sync_before_recv')
        end if
#endif


        TIMER_START('receive_multidomain_halos')

        ! Note: This loop can go in OMP, but apparently they do not support 
        ! type bound procedures, so we avoid calling like that. 
        ! {Fixed in openmp 5 standard?}

        !!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(md)
        !!!$OMP DO SCHEDULE(DYNAMIC)
        do j = 1, size(md%domains)
#ifdef TIMER
            call md%domains(j)%timer%timer_start('receive_halos')
#endif
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(md, j)
            !$OMP DO SCHEDULE(DYNAMIC)
            do i = 1, size(md%domains(j)%nesting%recv_comms)

                !! Avoid use of type-bound-procedure in openmp region
                !call md%domains(j)%nesting%recv_comms(i)%process_received_data(md%domains(j)%U)
                call process_received_data(md%domains(j)%nesting%recv_comms(i), md%domains(j)%U)

                ! If a linear domain is receiving, need to avoid getting non-zero UH/VH at wet/dry
                ! boundaries
                if(md%domains(j)%nesting%recv_comms(i)%recv_active .and. &
                   md%domains(j)%nesting%recv_comms(i)%use_wetdry_limiting .and. &
                   md%domains(j)%timestepping_method == 'linear') then

                    ! Indices of the region that received data
                    iL = md%domains(j)%nesting%recv_comms(i)%recv_inds(1,1)
                    iU = md%domains(j)%nesting%recv_comms(i)%recv_inds(2,1)
                    jL = md%domains(j)%nesting%recv_comms(i)%recv_inds(1,2)
                    jU = md%domains(j)%nesting%recv_comms(i)%recv_inds(2,2)

                    ! Elevation above which flux = 0 for linear domain
                    elev_lim = (md%domains(j)%msl_linear - minimum_allowed_depth)

                    ! Loop over the received region, and ensure UH/VH on
                    ! wet-dry boundaries are 0
                    do jj = jL, jU

                        ! jj+1, avoiding out-of-bounds
                        jp1 = min(jj+1, size(md%domains(j)%U, 2))

                        do ii = iL, iU

                            ! ii+1, avoiding out-of-bounds
                            ip1 = min(ii+1, size(md%domains(j)%U, 1))

                            ! Condition for zero EW flux on staggered grid
                            if( (md%domains(j)%U(ii,jj,ELV) >= elev_lim) .or. &
                                (md%domains(j)%U(ip1,jj,ELV) >= elev_lim) ) then
                                md%domains(j)%U(ii,jj,UH) = 0.0_dp
                            end if

                            ! Condition for zero NS flux on staggered grid
                            if( (md%domains(j)%U(ii,jj,ELV) >= elev_lim) .or. &
                                (md%domains(j)%U(ii,jp1,ELV) >= elev_lim) ) then
                                md%domains(j)%U(ii,jj,VH) = 0.0_dp
                            end if
                            
                        end do
                    end do
                    
                end if

            end do
            !$OMP END DO
            !$OMP END PARALLEL

#ifdef TIMER            
            call md%domains(j)%timer%timer_end('receive_halos')
#endif
        end do
        !!$OMP END DO
        !!$OMP END PARALLEL

        TIMER_STOP('receive_multidomain_halos')

#ifdef COARRAY
        if(sync_after_local) then
            TIMER_START('sync_after_recv')
            sync images(linked_p2p_images)
            TIMER_STOP('sync_after_recv')
        end if
#endif

    end subroutine

    !
    ! Record the flow volume in all domains at the 'start' of a simulation
    ! Useful for mass conservation tracking
    !
    subroutine record_initial_volume(md)
        class(multidomain_type), intent(inout) :: md

        integer(ip):: i

        do i = 1, size(md%volume_initial)        
            call md%get_flow_volume(md%volume_initial(i), i)
        end do

    end subroutine

    !
    ! Convenience routine to send nested grid halo data to neighbour images in a multidomain
    !
    ! @param md multidomain object
    ! @param domain_index optional integer. If provided, send halos of
    !   md%domains(domain_index). Otherwise, send halos for all domains.
    ! @param send_to_recv_buffer optional logical. Default TRUE. If FALSE, we
    !   do not do a coarray communication, but only copy data to the send buffer.
    !
    subroutine send_multidomain_halos(md, domain_index, send_to_recv_buffer) !, coarser_domains_only)
        class(multidomain_type), intent(inout) :: md
        integer(ip), intent(in), optional:: domain_index
        logical, intent(in), optional:: send_to_recv_buffer

        integer(ip) :: i, j, jmn, jmx
        logical::send_to_recv_buffer_local

        TIMER_START('send_multidomain_halos')

        if(present(domain_index)) then
            ! Only send domain_index
            jmn = domain_index
            jmx = domain_index
        else
            jmn = 1
            jmx = size(md%domains)
        end if

        if(present(send_to_recv_buffer)) then
            send_to_recv_buffer_local = send_to_recv_buffer
        else
            send_to_recv_buffer_local = .TRUE.
        end if
#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
        if(send_halos_immediately) then
            write(log_output_unit, *) 'Error: Cannot do send_halos_immediately with COARRAY_USE_MPI_FOR_INTENSIVE_COMMS', &
                __LINE__, __FILE__
            call generic_stop
        end if
#endif

        do j = jmn, jmx
#ifdef TIMER
            call md%domains(j)%timer%timer_start('send_halos')
#endif
            do i = 1, size(md%domains(j)%nesting%send_comms)

                call md%domains(j)%nesting%send_comms(i)%process_data_to_send(md%domains(j)%U)
                call md%domains(j)%nesting%send_comms(i)%send_data(send_to_recv_buffer=send_to_recv_buffer_local)
            end do
#ifdef TIMER
            call md%domains(j)%timer%timer_end('send_halos')
#endif
        end do

        TIMER_STOP('send_multidomain_halos')
    end subroutine

    ! In nesting communication, when a fine domain receives from a coarse domain, 
    ! we may want to 'not update' fine halo cells that are very close to their domain.
    !
    ! This can be beneficial for numerical stability, although it can also cause problems. 
    !
    ! This routine sets the nesting recv_weights to 0.0_dp, for fine-domain cells that 
    ! are close to a coarse domain.
    ! 
    ! Beware that other parameters in nested_grid_comms_mod.f90 control whether these weights
    ! are even used -- which depends also on the depth and depth-gradients.

    ! In any case, halos must be 'fatter than required'.
    !
    subroutine separate_fine_halos_from_their_domain(md)
        class(multidomain_type), intent(inout) :: md
       
        integer(ip) :: i, j, nx, ny, irecv_L, irecv_U, ir, jrecv_L, jrecv_U, jr
        integer(ip) :: ia, ib, ja, jb, i1, j1
        integer(ip) :: sep


        ! Search all domains
        do j = 1, size(md%domains)

            ! Quick exit
            if(.not. allocated(md%domains(j)%nesting%recv_comms)) cycle

            nx = size(md%domains(j)%nesting%priority_domain_index, 1)
            ny = size(md%domains(j)%nesting%priority_domain_index, 2)

            if(md%domains(j)%max_parent_dx_ratio > 1_ip) then
                ! Separate this much from halos -- must be a multiple of md%domains(j)%max_parent_dx_ratio
                sep = md%extra_halo_buffer*md%domains(j)%max_parent_dx_ratio 
            else
                ! No extra separation if we only receive data from domains of the same or finer size
                cycle
            end if

            ! Search the recv comms
            do i = 1, size(md%domains(j)%nesting%recv_comms)

                ! Look for recv_comms where a fine domain receives
                if(.not. (md%domains(j)%nesting%recv_comms(i)%recv_active .and. &
                          md%domains(j)%nesting%recv_comms(i)%my_domain_is_finer)) cycle
       
                ! Indices of the recv domain where the recv happens 
                irecv_L = md%domains(j)%nesting%recv_comms(i)%recv_inds(1,1)
                irecv_U = md%domains(j)%nesting%recv_comms(i)%recv_inds(2,1)
                jrecv_L = md%domains(j)%nesting%recv_comms(i)%recv_inds(1,2)
                jrecv_U = md%domains(j)%nesting%recv_comms(i)%recv_inds(2,2)

                ! Loop over every cell in the recv region, and see if it is within 'sep' cells
                ! of a cell having priority_domain_index/image equal to the recv domain.
                ! If such cells are found, then the halo is 'near', so we set the recv_weight to 0
                do j1 = jrecv_L, jrecv_U
                    do i1 = irecv_L, irecv_U
                        ! Avoid out-of-bounds
                        ia = max(1, i1-sep)
                        ib = min(nx, i1+sep)
                        ja = max(1, j1-sep)
                        jb = min(ny, j1+sep)

                        ! Check if the halo cell is near. If so, ignore the received data
                        if(any(&
                            (md%domains(j)%nesting%priority_domain_image(ia:ib, ja:jb) == &
                                md%domains(j)%nesting%recv_comms(i)%my_domain_image_index) .and. &
                            (md%domains(j)%nesting%priority_domain_index(ia:ib, ja:jb) == &
                                md%domains(j)%nesting%recv_comms(i)%my_domain_index) ) ) then

                                md%domains(j)%nesting%recv_comms(i)%recv_weights(i1, j1) = 0.0_dp
                            
                        end if
                    end do
                end do
                !do j1 = jrecv_L, jrecv_U
                !    print*, ' col : ', j1
                !    print*, md%domains(j)%nesting%recv_comms(i)%recv_weights(:,j1)
                !end do
            end do
        end do

    end subroutine


    ! Set the flow variables in 'null regions' of a domain to 'high and dry' values
    !
    ! The 'null regions' are areas where nesting%recv_metadata(:,1) == 0 
    ! Such regions have no impact on the computed flow solution
    ! in 'priority_domain' regions. 
    !
    ! @param md
    ! @param ignore_linear logical. Default FALSE. If TRUE, do not make any changes to null
    !     regions of domains with timestepping_method == 'linear'. Previously this was
    !     used to avoid minor wet-dry artefacts (now fixed), caused by neighbouring domains communicating
    !     non-zero U(:,:,UH) or U(:,:,VH) along wet-dry boundaries. The linear solver
    !     should always have zero flux terms along such boundaries, but that can be
    !     violated by nesting communication. 
    subroutine set_null_regions_to_dry(md, ignore_linear)

        class(multidomain_type), intent(inout) :: md
        logical, intent(in), optional :: ignore_linear

        ! Null regions have particular values of the recv_metadata index and image
        integer(ip), parameter :: null_index = 0, null_image = 0

        ! Values used to denote 'null regions'
        real(dp), parameter :: stage_null = wall_elevation
        real(dp), parameter :: elev_null = wall_elevation
        real(dp), parameter :: uh_null = 0.0_dp , vh_null = 0.0_dp

        integer(ip) :: i, j, x1, x2, y1, y2
        logical :: ignore_linear_local

        if(present(ignore_linear)) then
            ignore_linear_local = ignore_linear
        else
            ignore_linear_local = .false.
        end if 

        do j = 1, size(md%domains)

            ! Quick exit
            if(.not. allocated(md%domains(j)%nesting%recv_metadata)) cycle

            ! For wetting and drying reasons, we may wish to not apply this to domains
            ! solving the linear equations
            if(md%domains(j)%timestepping_method == 'linear' .and. ignore_linear_local) cycle

            ! Loop over recv_metadata boxes
            do i = 1, size(md%domains(j)%nesting%recv_metadata, 2)

                ! Check for null region
                if( (md%domains(j)%nesting%recv_metadata(1, i) == null_index) .and. &
                    (md%domains(j)%nesting%recv_metadata(2, i) == null_image) ) then

                    ! Inside null region, set domains%U to the NULL values
                    x1 = md%domains(j)%nesting%recv_metadata(3,i)
                    x2 = md%domains(j)%nesting%recv_metadata(4,i)
                    y1 = md%domains(j)%nesting%recv_metadata(5,i)
                    y2 = md%domains(j)%nesting%recv_metadata(6,i)

                    md%domains(j)%U(x1:x2, y1:y2, STG) = stage_null
                    md%domains(j)%U(x1:x2, y1:y2, UH) = uh_null
                    md%domains(j)%U(x1:x2, y1:y2, VH) = vh_null
                    md%domains(j)%U(x1:x2, y1:y2, ELV) = elev_null

                end if
            end do
        end do

    end subroutine

    ! Routine to put all domains outputs in a single folder
    subroutine setup_multidomain_output_folder(md)

        class(multidomain_type), intent(inout) :: md
        integer:: i, date_time_values(8) ! For date_and_time
        character(len=charlen) :: datetime_char, mkdir_command

        call date_and_time(values=date_time_values)
#ifdef COARRAY
        call co_broadcast(date_time_values, source_image = 1)
#endif
        write(datetime_char, '(I4,I2.2,I2.2,A,I2.2,I2.2,I2.2,I3.3)') &
            date_time_values(1), date_time_values(2), date_time_values(3), '_', date_time_values(5),&
            date_time_values(6), date_time_values(7), date_time_values(8)

        md%output_basedir = trim(md%output_basedir) // '/RUN_' // trim(datetime_char)

        ! Make the directory
#ifndef COARRAY
        call mkdir_p(md%output_basedir)
#else
        if(this_image() == 1) call mkdir_p(md%output_basedir)
        sync all ! Ensure the directory creation is finished
#endif
    end subroutine


    ! Initial entry-point to setup multidomain
    !
    ! Setup defaults, partition the domains if required, and call setup_multidomain_domains
    !
    ! @param md multidomain
    ! @param verbose control print-out verbosity
    ! @param use_wetdry_limiting_nesting If .true., then in wet-dry regions the nesting communication differs. This is important to
    ! avoid wet-dry artefacts related to nesting.
    ! @param capture_log log outputs in a multidomain specific file
    ! @param extra_halo_buffer integer controlling the 'extra' width of halos beyond that strictly required for the solver.
    ! Increased buffer size can be used to separate neighbouring send halo regions, which can be important for numerical stability,
    ! although this also requires the use of nesting receive weights (controlled by other parts of the code)
    ! @param allocate_comms Call allocate_p2p_comms to setup the parallel communication structures
    !
    subroutine setup_multidomain(md, verbose, use_wetdry_limiting_nesting, capture_log, extra_halo_buffer, allocate_comms)

        class(multidomain_type), intent(inout) :: md
        logical, optional, intent(in) :: verbose
        logical, optional, intent(in) :: use_wetdry_limiting_nesting
        logical, optional, intent(in) :: capture_log
        integer(ip), optional, intent(in) :: extra_halo_buffer
        logical, optional, intent(in) :: allocate_comms

        logical:: verbose1, use_wetdry_limiting, capture_log_local, allocate_comms_local
        integer(ip) :: i, extra_halo_buffer_local
        character(len=charlen) :: log_filename
#ifdef TIMER
        call md%timer%timer_start('setup')
#endif

        if(present(verbose)) then
            verbose1 = verbose
        else
            verbose1 = .true.
        end if

        ! By default, the nesting method does use wet-dry limiting
        if(present(use_wetdry_limiting_nesting)) then
            use_wetdry_limiting = use_wetdry_limiting_nesting
        else
            use_wetdry_limiting = .true.
        end if

        if(present(capture_log)) then
            capture_log_local = capture_log
        else
            capture_log_local = .true.
        end if

        if(present(allocate_comms)) then
            allocate_comms_local = allocate_comms
        else
            allocate_comms_local = .true.
        end if

        ! Allow additional padding of halos -- some experiments to improve numerical stability require this
        if(present(extra_halo_buffer)) then
            extra_halo_buffer_local = extra_halo_buffer
        else
            ! This should be 0, unless we find extra halos are important
            extra_halo_buffer_local = extra_halo_buffer_default
        end if
        md%extra_halo_buffer = extra_halo_buffer_local

        ! Make separate output folder
        call setup_multidomain_output_folder(md)
   
        if(capture_log_local) then     
            ! Log in that folder
            log_filename = trim(md%output_basedir) // '/multidomain_log'
            call send_log_output_to_file(log_filename)
        end if

        ! Split up domains among images, and create all md%domains(i)%myid
        call md%partition_domains()

        ! Make sure all domains use the new output_basedir 
        do i = 1, size(md%domains)
            md%domains(i)%output_basedir = md%output_basedir
        end do


        ! Main setup routine
        call setup_multidomain_domains(domains=md%domains, &
            verbose=verbose1, &
            use_wetdry_limiting_nesting=use_wetdry_limiting,&
            periodic_xs = md%periodic_xs, &
            periodic_ys = md%periodic_ys, &
            extra_halo_buffer = md%extra_halo_buffer, &
            extra_cells_in_halo = md%extra_cells_in_halo, &
            all_dx_md = md%all_dx_md, &
            all_is_staggered_grid_md = md%all_is_staggered_grid_md)

        ! Storage space for mass conservation
        allocate(md%volume_initial(size(md%domains)), md%volume(size(md%domains)))
    
        md%volume_initial = 0.0_dp
        md%volume = 0.0_dp

        if(allocate_comms_local) call allocate_p2p_comms
#ifdef TIMER
        call md%timer%timer_end('setup')
#endif
    end subroutine

    ! Main routine for setting up the multidomain 
    !
    ! Determine which domains we need to have a two-way-nesting-comms relation
    ! with, adjust domain bboxes, setup the nesting_boundaries, etc
    ! 
    ! What this routine does:
    !   Step 1: For each domain, compute its required nesting layer width. This is done
    !           by looking 'just outside' the outer boundary of each domain, to see if another
    !           domain is there. If so, we determine the nesting layer width so that
    !           a) cells in overlapping areas are completely overlapping [i.e. no half-covered coarse cells]
    !           b) We can sync only once per global time-step, without having
    !              the flow propagate over the full nesting layer width in that time
    !           c) There may be other constraints too (e.g. separate buffers ?)
    !   Step 2: For each domain, make a 'priority_domain_metadata', which tells us what
    !           the priority domain for each cell is. We blank out cells which are more than
    !           one nesting layer width away from that domain's priority region.
    !   Step 3: Collapse the priority_domain_metadata to a set of rectangles, represented in a 
    !           tabular format. These DEFINE the regions that the domain needs to receive data from,
    !           and so we call them the 'receive metadata'
    !   Step 4: Broadcast the 'receive metadata', and create the 'send metadata', based on the
    !           regions that need to be received
    !   Step 5: Broadcast the 'send metadata'. This is done because the receives need to know
    !           part of the send_metadata (i.e. the comms index of domains(j)%nesting%send_comms that they will be
    !           getting data from) to build the communication data structure. Further, having the full send metadata
    !           is useful for debugging.
    !   Step 6: Set-up the domains(j)%nesting%send_comms, and domains(j)%nesting%recv_comms. This is 
    !           straightforward given the information above.
    !
    ! @param domains array of domains 
    ! @param verbose control print-out verbosity
    ! @param use_wetdry_limiting_nesting If .true., then in wet-dry regions the nesting communication differs. This is important to
    ! avoid wet-dry artefacts related to nesting.
    ! @param periodic_xs cells with x-coordinate < periodic_xs(1) or > periodic_xs(2) are halos that will receive from the other
    ! side of the domain.
    ! @param periodic_ys Analogous to periodic_xs
    ! @param extra_halo_buffer integer controlling the 'extra' width of halos beyond that strictly required for the solver.
    ! @param extra_cells_in_halo another integer controlling the 'extra' width of halos.
    ! @param all_dx_md stores the dx values for all domains in multidomain (i.e. copy of all_dx)
    ! @param all_is_staggered_grid_md stores the all_is_staggered_grid variable in the multidomain.
    ! Increased buffer size can be used to separate neighbouring send halo regions, which can be important for numerical stability,
    ! although this also requires the use of nesting receive weights (controlled by other parts of the code)
    !
    subroutine setup_multidomain_domains(domains, verbose, use_wetdry_limiting_nesting,&
        periodic_xs, periodic_ys, extra_halo_buffer, extra_cells_in_halo, &
        all_dx_md, all_is_staggered_grid_md)

        type(domain_type), intent(inout) :: domains(:)
        logical, optional, intent(in) :: verbose
        logical, optional, intent(in) :: use_wetdry_limiting_nesting
        real(dp), intent(in) :: periodic_xs(2), periodic_ys(2)
        integer(ip), intent(in) :: extra_halo_buffer, extra_cells_in_halo
        real(dp), allocatable, intent(inout) :: all_dx_md(:,:,:)
        integer(ip), allocatable, intent(inout) :: all_is_staggered_grid_md(:,:)

        integer(ip):: i, j, k, jj, ii, nd_local, nest_layer_width, tmp2(2), nbox_max, counter, nd_global

        real(dp), allocatable:: nbr_domain_dx_local(:,:), yc_tmp(:), xc_tmp(:), y_tmp(:)
        real(dp), allocatable:: local_metadata(:,:)

        integer(ip) :: ijk_to_send(2,3), ijk_to_recv(2,3)

        real(dp) :: d_lw(2), d_dx(2), d_ll(2), a_row(srm+2), tmp_xs(2), tmp_ys(2), tmp_dxs(2)
        real(dp) :: box_diff(4)
        integer(ip) :: d_nx(2), send_count, recv_count, n_ind, n_img, n_comms, n_row, tmp
        logical :: verbose1, use_wetdry_limiting, in_periodic_region1, in_periodic_region2, in_periodic_region
        character(len=charlen) :: msg

        !! Index labels for send/recv metadata 
        !! The following parameters give the 'row-index' at which we store various types of metadata.
        ! IND --> index in domains(:)
        ! IMG --> image holding the domain 
        ! XLO, XHI --> indices for x coordinates in communication region. These must be consecutive
        ! YLO, YHI --> indices for y coordinates in communication region. These must be consecutive
        !
        integer(ip), parameter :: IND = 1,  IMG = 2, XLO = 3, XHI = 4, YLO = 5, YHI = 6

        if(present(verbose)) then
            verbose1 = verbose
        else
            verbose1 = .true.
        end if

        ! By default, the nesting method does use wet-dry limiting (i.e. when
        ! doing coarse-to-fine interpolation, gradients are limited near wet-dry areas).
        if(present(use_wetdry_limiting_nesting)) then
            use_wetdry_limiting = use_wetdry_limiting_nesting
        else
            use_wetdry_limiting = .true.
        end if

        nd_local = size(domains)
        nd_global = nd_local
#ifdef COARRAY
        call co_max(nd_global)
#endif

        !
        ! Figure out how thick the nesting layer buffer has to be for each domain
        ! in the multi domain. 
        ! Note this is ONLY based on checking around the 'exterior boundary', one cell
        ! width outside of each domain. There might be 'internal' nesting
        ! regions, but they are computed later
        ! FIXME: This will need editing if we are to do 'multi-time-level-communication'.
        !        Or, we might be able to back-calculate the required thickness
        call compute_multidomain_nesting_layer_width(domains=domains, verbose=verbose1, &
            periodic_xs=periodic_xs, periodic_ys=periodic_ys, extra_halo_buffer=extra_halo_buffer, &
            extra_cells_in_halo=extra_cells_in_halo)

        !
        ! Extend domain bounding boxes to include nesting layers
        !
        do j = 1, nd_local
            do i = 1, 4 
                !
                ! i = {1, 2, 3, 4} <==> boundary = {N, E, S, W}
                !
                if(domains(j)%is_nesting_boundary(i)) then

                    d_lw = domains(j)%lw
                    d_dx = domains(j)%dx
                    d_ll = domains(j)%lower_left
                    d_nx = domains(j)%nx
                    nest_layer_width = domains(j)%nest_layer_width

                    ! Case where we expand in x direction
                    if( (i == 2) .or. (i == 4) ) then
                        domains(j)%lw(1) = d_lw(1) + nest_layer_width * d_dx(1)
                        domains(j)%nx(1) = d_nx(1) + nest_layer_width
                    end if
                    ! Case where we expand in y direction
                    if( (i == 1) .or. (i == 3) ) then
                        domains(j)%lw(2) = d_lw(2) + nest_layer_width * d_dx(2)
                        domains(j)%nx(2) = d_nx(2) + nest_layer_width
                    end if

                    ! Southern and western nesting boundaries also imply the
                    ! lower_left is shifted
                    if(i == 3) then
                        domains(j)%lower_left(2) = d_ll(2) - nest_layer_width * d_dx(2)
                    end if

                    if(i == 4) then
                        domains(j)%lower_left(1) = d_ll(1) - nest_layer_width * d_dx(1)
                    end if

                end if
            end do
        end do

        !
        ! Convert nesting information to 'box' metadata
        !
        nbox_max = 0 ! Store the maximum number of rows in the box metadata
        do j = 1, nd_local

            d_nx = domains(j)%nx
            ! For every point in the domain, store the index/image of the priority domain
            ! which contains it. Often, this will just be domains(j). However, if the
            ! domain is overlapped by finer domains, then they will be recorded here
            allocate( domains(j)%nesting%priority_domain_index(d_nx(1), d_nx(2)) )
            allocate( domains(j)%nesting%priority_domain_image(d_nx(1), d_nx(2)) )

            ! Below we call 'find_priority_domain_containing_xy', and need to
            ! pass an array which can store the neighbour dx.
            if(allocated(nbr_domain_dx_local)) deallocate(nbr_domain_dx_local)
            allocate(nbr_domain_dx_local(d_nx(1), 2))

            ! Also convenient to have a x/y row-wise coordinate vector
            call get_domain_xs_ys(domains(j), xc_tmp, yc_tmp)
            if(allocated(y_tmp)) deallocate(y_tmp)
            allocate(y_tmp(size(xc_tmp)))

            ! Loop over each column (i.e. with fixed y coordinate) and for each
            ! cell, find the priority domain index and corresponding image in
            ! which it is contained
            do jj = 1, d_nx(2) 
                y_tmp = xc_tmp * 0.0_dp + yc_tmp(jj)
                call find_priority_domain_containing_xy(&
                    xc_tmp, y_tmp, &
                    domains(j)%nesting%priority_domain_index(:,jj), &
                    domains(j)%nesting%priority_domain_image(:,jj), &
                    nbr_domain_dx_local, error_domain_overlap_same_dx=.true.,&
                    periodic_xs=periodic_xs, periodic_ys=periodic_ys)

                ! Useful to store the current domain index and image inside the nesting structure
                domains(j)%nesting%my_index = j
                domains(j)%nesting%my_image = ti
          
                if(any(domains(j)%nesting%priority_domain_index(:,jj) < 0)) then
                    write(log_output_unit, *) &
                        'Error: priority_domain_index contains areas that are not inside any domain. ', &
                        'This can happen if domains are small and have nesting buffers large enough ', &
                        'to spill outside their neighbours (e.g. using too much parallel refinement)'
                    call generic_stop
                end if
 
                ! Logical check 
                do ii = 1, 2
                    if(any(nbr_domain_dx_local(:,ii) > (1.0001_dp * domains(j)%max_parent_dx_ratio*domains(j)%dx(ii)))) then
                        write(log_output_unit,*) &
                            'ERROR: nesting boundary overlaps with a domain coarser than domains on its immediate boundary.'
                        write(log_output_unit,*) &
                            'This may occur if the nesting buffers are very fat, and there are close domains that do not touch'
                        write(log_output_unit,*) 'domain_index: ', j
                        write(log_output_unit,*) domains(j)%max_parent_dx_ratio
                        write(log_output_unit,*) domains(j)%dx(ii)
                        write(log_output_unit,*) domains(j)%dx(ii) * domains(j)%max_parent_dx_ratio
                        write(log_output_unit,*) maxval(nbr_domain_dx_local(:,ii))
                        stop
                    end if
                end do
            end do

            ! Convert the priority domain information to a set of boxes, which is
            ! well designed for nested comms
            !
            ! The box_metadata in domains(j)%nesting%recv_metadata contains one column
            ! for each box, with the variables:
            !    domain_index, domain_image, box_left_i, box_right_i, box_bottom_j, box_top_j
            !
            ! FIXME: This would have to be changed to allow 'multi-time-level-communication',
            !        to deal with variable halow_width
            call get_domain_xs_ys(domains(j), xc_tmp, yc_tmp)
            call convert_priority_domain_info_to_boxes(&
                priority_domain_index = domains(j)%nesting%priority_domain_index, &
                priority_domain_image = domains(j)%nesting%priority_domain_image, &
                halo_width = domains(j)%nest_layer_width, &
                myimage = ti, &
                myindex = j, &
                box_metadata = domains(j)%nesting%recv_metadata, &
                xs=xc_tmp, &
                ys=yc_tmp, &
                periodic_xs=periodic_xs,&
                periodic_ys=periodic_ys,&
                max_parent_dx_ratio=nint(domains(j)%max_parent_dx_ratio))

            ! Keep track of the maximum number of 'boxes' in the box metadata
            tmp2 = shape(domains(j)%nesting%recv_metadata)
            nbox_max = max(nbox_max, tmp2(2))
  
        end do
#ifdef COARRAY
        call co_max(nbox_max)
#endif

        !
        ! Allocate arrays in each domain
        !
        do j = 1, nd_local

            if(verbose1) write(log_output_unit,*) ' Domain ID: ', domains(j)%myid

            domains(j)%boundary_exterior = (.NOT. domains(j)%is_nesting_boundary)

            d_nx = domains(j)%nx
            d_lw = domains(j)%lw
            d_ll = domains(j)%lower_left
            call domains(j)%allocate_quantities(global_lw = d_lw, global_nx = d_nx, global_ll = d_ll, verbose=verbose1)

        end do

        !
        ! Broadcast all domains(j)%nesting%recv_metadata, but replace x/y coordinate
        ! indices with actual x-y values. For simplicity make the whole thing contain reals,
        ! even though both reals and integers are stored
        !

        ! Dimensions of all_recv_metadata correspond to: 
        ! recv_from_domain_index, recv_from_image_index, xlo, xhi, ylo, yhi, send_metadata_row_index, recv_comms_index
        if(allocated(all_recv_metadata)) deallocate(all_recv_metadata)
#ifdef COARRAY
        allocate( all_recv_metadata(srm+2, nbox_max, nd_global, ni )[*])
#else     
        allocate( all_recv_metadata(srm+2, nbox_max, nd_global, ni ) )
#endif

        ! Fill the local_metadata (which will be communicated to all_recv_metadata)
        if(allocated(local_metadata)) deallocate(local_metadata)
        allocate(local_metadata(srm+2, nbox_max))

        do j = 1, nd_local

            local_metadata = -1.0_dp ! Default 'empty box' value
            counter = 0 ! Keep track of 'real' receive indices -- i.e. those that communicate with another domain/image
    
            do i = 1, size(domains(j)%nesting%recv_metadata, 2)

                ! Domain index and image that is needed for box i [as reals, even though the values are integer]
                local_metadata([IND,IMG],i) = 1.0_dp * domains(j)%nesting%recv_metadata([IND,IMG], i)

                ! Range of x/y values in communication zone
                ! Firstly, get centroid  coordinates, and determine if we are in the periodic region
                ! These were constructed to never cross a periodic boundary
                tmp_xs = domains(j)%x( domains(j)%nesting%recv_metadata(3:4,i) )
                tmp_ys = domains(j)%y( domains(j)%nesting%recv_metadata(5:6,i) ) 
                call check_periodic(tmp_xs(1), tmp_ys(1), periodic_xs, periodic_ys, in_periodic_region, &
                    adjust_coordinates=.FALSE.)
                call check_periodic(tmp_xs(2), tmp_ys(2), periodic_xs, periodic_ys, in_periodic_region, &
                    adjust_coordinates=.FALSE.)
                ! Range of the communication regions, extended to cell edges
                local_metadata(3:4,i) = tmp_xs + [-0.5_dp, 0.5_dp]*domains(j)%dx(1)
                local_metadata(5:6,i) = tmp_ys +  [-0.5_dp, 0.5_dp]*domains(j)%dx(2)
    
                ! Sequential index for 'true' communication rows
                if( ( (.not. in_periodic_region).and.&
                      (local_metadata(IND,i) == j .and. local_metadata(IMG,i) == ti)) .or. &
                      (local_metadata(IND,i) < 1) ) then
                    ! Not 'true' communication, because we are "not in the periodic region, and we are
                    ! communicating with our own image and index", or we're not in an active region 
                    ! The code design ensures this corresponds to regions which do not need nesting information.
                else
                    ! True communication.

                    ! For now we skip index 7 in local_metadata
 
                    counter = counter + 1
                    local_metadata(8,i) = counter
                end if

            end do
#ifdef COARRAY
            ! One-sided parallel communication to all images 
            do k = 1, ni
                all_recv_metadata(:,:,j,ti)[k] = local_metadata
            end do
#else
            all_recv_metadata(:,:,j,ti) = local_metadata
#endif
        end do
#ifdef COARRAY        
        sync all
#endif

        !
        ! Set up the 'send' metadata 
        !

        nbox_max = 0
        do j = 1, nd_local

            ! Make space for send_metadata, based on the number of domains
            ! which receive data from domain-index j on image ti
            send_count = count( (all_recv_metadata(IND,:,:,:) == j) .and. (all_recv_metadata(IMG,:,:,:) == ti))
            nbox_max = max(nbox_max, send_count)
            allocate(domains(j)%nesting%send_metadata(srm+2, send_count))
            domains(j)%nesting%send_metadata = 0_ip ! Pre-set
            if(send_count == 0) cycle

            ! Fill the send metadata
            counter = 0
            do k = 1, size(all_recv_metadata, 4) ! Number of images
                do jj = 1, size(all_recv_metadata,3) ! Number of domains
                    do i = 1, size(all_recv_metadata,2) ! Every row in the receive metadata

                        ! Find 'recvs' which need to get data from this index/image
                        if(all_recv_metadata(IND,i,jj,k) == j .and. all_recv_metadata(IMG,i,jj,k) == ti) then
                            counter = counter + 1
                            a_row = all_recv_metadata(:,i,jj,k)

                            ! Adjust the coordinates if we are in the periodic area
                            ! First shrink the boxes to avoid having any coordinate fall on a periodic boundary
                            ! when we check if the point is in a periodic region. This will be un-done shortly
                            tmp_dxs = all_dx(1:2, jj, k)
                            tmp_xs = [a_row(XLO), a_row(XHI)] + 0.5_dp * [1.0_dp, -1.0_dp] * tmp_dxs(1)
                            tmp_ys = [a_row(YLO), a_row(YHI)] + 0.5_dp * [1.0_dp, -1.0_dp] * tmp_dxs(2)
                            ! Then do the periodic correction
                            call check_periodic(tmp_xs(1), tmp_ys(1), periodic_xs, periodic_ys, &
                                in_periodic_region, adjust_coordinates=.TRUE.)
                            call check_periodic(tmp_xs(2), tmp_ys(2), periodic_xs, periodic_ys, &
                                in_periodic_region, adjust_coordinates=.TRUE.)
                            ! Then un-shrink the boxes
                            tmp_xs = tmp_xs - 0.5_dp * [1.0_dp, -1.0_dp] * tmp_dxs(1)
                            tmp_ys = tmp_ys - 0.5_dp * [1.0_dp, -1.0_dp] * tmp_dxs(2)

                            ! domain_index and image that we send to
                            domains(j)%nesting%send_metadata([IND, IMG] ,counter) = [jj, k]
                            ! x_indices that we will send
                            domains(j)%nesting%send_metadata(XLO,counter) = count(domains(j)%x < tmp_xs(1)) + 1
                            domains(j)%nesting%send_metadata(XHI,counter) = count(domains(j)%x < tmp_xs(2))
                            ! y_indices that we will send
                            domains(j)%nesting%send_metadata(YLO,counter) = count(domains(j)%y < tmp_ys(1)) + 1
                            domains(j)%nesting%send_metadata(YHI,counter) = count(domains(j)%y < tmp_ys(2))

                            ! Store the row indices of corresponding recv. Very handy for establishing the comms
                            domains(j)%nesting%send_metadata(7,counter) = i

                            ! In column 8, store whether or not we are doing a periodic nest
                            domains(j)%nesting%send_metadata(8, counter) = count([in_periodic_region])

                            ! Now, tell all_recv_metadata this row index in send_metadata
#ifdef COARRAY
                            do ii = 1, ni
                                ! One-sided all-to-all communication of the send_metadata index to the recv_metadata
                                ! This index will be the same as the one into all_send_metadata
                                all_recv_metadata(7,i,jj,k)[ii] = (1.0_dp * counter)
                            end do
#else
                            all_recv_metadata(7,i,jj,k) = (1.0_dp * counter)
#endif
                        end if
                    end do
                end do
            end do

        end do
#ifdef COARRAY
        sync all
#endif        

#ifdef MULTIDOMAIN_DEBUG
        ! DEBUG: Write out all_recv_metadata. Compile with -DMULTIDOMAIN_DEBUG to do this
        write(log_output_unit, *) ''
        write(log_output_unit, *) '## all_recv_metadata:'
        do k = 1, ni
            do j = 1, size(all_recv_metadata, 3) 
                write(log_output_unit, *) '    ', k, j
                do i = 1, size(all_recv_metadata, 2)
                    write(log_output_unit, *) '      ', all_recv_metadata(:,i,j,k)
                end do
            end do
        end do
        write(log_output_unit, *) ''
        flush(log_output_unit)
#endif
        !
        ! Broadcast the send metadata
        !

        if(allocated(all_send_metadata)) deallocate(all_send_metadata)
#ifdef COARRAY
        call co_max(nbox_max)
        allocate(all_send_metadata(srm+2, nbox_max, nd_global, ni)[*])
#else
        allocate(all_send_metadata(srm+2, nbox_max, nd_global, ni))
#endif
        if(allocated(local_metadata)) deallocate(local_metadata)
        allocate(local_metadata(srm+2, nbox_max))

        !
        ! For every domain, pack the metadata into a local array, and send to all_send_metadata
        !
        do j = 1, nd_local

            local_metadata = -1.0_dp ! Default 'empty box' value
            counter = 0 ! Keep track of 'true sends', as opposed to rows that don't need sends

            do i = 1, size(domains(j)%nesting%send_metadata, 2)

                ! Domain index and image that is needed for box i [as reals, even though the values are integer]
                local_metadata([IND, IMG],i) = 1.0_dp * domains(j)%nesting%send_metadata([IND, IMG], i)
                ! Range of x values in communication zone
                local_metadata([XLO, XHI],i) = domains(j)%x( domains(j)%nesting%send_metadata([XLO, XHI],i) ) + &
                    [-0.5_dp, 0.5_dp]*domains(j)%dx(1)
                ! Range of y values in communication zone
                local_metadata([YLO, YHI],i) = domains(j)%y( domains(j)%nesting%send_metadata([YLO, YHI],i) ) + &
                    [-0.5_dp, 0.5_dp]*domains(j)%dx(2)

                ! Is the receive in a periodic region? (Of course the send is not, by construction, since periodic regions
                ! cannot be in the priority domain)
                in_periodic_region = (domains(j)%nesting%send_metadata(8,i) == 1_ip)

                ! Sequential index for 'true' communication rows
                if((.not. in_periodic_region) .and. &
                   ((local_metadata(IND,i) == j).and.(local_metadata(IMG,i) == ti)).or.(local_metadata(IND,i) < 1)) then
                    ! Not real communication. Either the domain is sending to
                    ! itself, at a site that is not in a periodic region, OR the metadata is empty
                else

                    local_metadata(7,i) = 1.0_dp * domains(j)%nesting%send_metadata(7,i)

                    counter = counter+1
                    local_metadata(8,i) = 1.0_dp * counter
                end if

            end do
#ifdef COARRAY
            ! One-sided parallel communication to all images
            do k = 1, ni
                all_send_metadata(:,:,j,ti)[k] = local_metadata
            end do
#else
            all_send_metadata(:,:,j,ti) = local_metadata
#endif
        end do
#ifdef COARRAY        
        sync all
#endif

#ifdef MULTIDOMAIN_DEBUG
        ! DEBUG: Write out all_send_metadata
        write(log_output_unit, *) ''
        write(log_output_unit, *) '## all_send_metadata:'
        do k = 1, ni
            do j = 1, size(all_send_metadata, 3) 
                write(log_output_unit, *) '    ', k, j
                do i = 1, size(all_send_metadata, 2)
                    write(log_output_unit, *) '      ', all_send_metadata(:,i,j,k)
                end do
            end do
        end do
        write(log_output_unit, *) ''
        flush(log_output_unit)
#endif
        !
        ! Set up sends and recvs
        !

        do j = 1, nd_local
    
            !
            ! Setup sends
            !
            send_count = int( nint(max(0.0_dp , maxval(all_send_metadata(8,:,j,ti)))), ip)
            allocate(domains(j)%nesting%send_comms(send_count))

            ! Loop over rows of all_send_metadata(:,:,j,ti), and set up the
            ! send comms if required
            counter = 0
            do i = 1, size(all_send_metadata, 2)

                if(all_send_metadata(8,i,j,ti) < 1) cycle

                ! Neighbour domain index/image
                n_ind = nint(all_send_metadata(IND,i,j,ti))
                n_img = nint(all_send_metadata(IMG,i,j,ti))
                ! Row in all_recv_metadata corresponding to this send
                n_row = nint(all_send_metadata(7,i,j,ti)) 
                ! Send to domains(n_ind)%nesting%recv_comms(n_comms) on image n_img
                n_comms = nint(all_recv_metadata(8, n_row, n_ind, n_img)) 

                ! DEBUG: Check that the send/recv coordinate boxes are
                ! identical, to within tolerable rounding
                box_diff = all_send_metadata(3:6, i    ,     j,    ti) - &
                           all_recv_metadata(3:6, n_row, n_ind, n_img)
                if( any( abs(box_diff) > 100.0_dp*spacing(maxval(abs(all_send_metadata(3:6, i, j, ti)))) )) then
                    ! Looks like the boxes don't match up so well
                    ! But this could be due to periodic boundaries
                    if(any(abs(box_diff(1:2)) > 0.9999_dp*(periodic_xs(2) - periodic_xs(1))) .or. &
                       any(abs(box_diff(3:4)) > 0.9999_dp*(periodic_ys(2) - periodic_ys(1)))) then
                        ! Periodic boundaries. Skip them! 
                    else
                        msg = 'Error: Send/recv metadata do not appear to have the same boxes to within round-off'
                        write(log_output_unit,*) trim(msg)
                        write(log_output_unit,*) i, j, ti
                        write(log_output_unit,*) n_ind, n_img, n_row, n_comms
                        write(log_output_unit,*) all_send_metadata(3:6,i,j,ti)
                        write(log_output_unit,*) all_recv_metadata(3:6,n_row, n_ind, n_img)
                        write(log_output_unit,*) 100.0_dp*spacing(maxval(abs(all_send_metadata(3:6, i, j, ti))))
                    end if
                end if

                ! Array bounds to send
                ijk_to_send(1:2,1) = domains(j)%nesting%send_metadata([XLO, XHI],i) ! Xrange
                ijk_to_send(1:2,2) = domains(j)%nesting%send_metadata([YLO, YHI],i) ! Yrange
                ijk_to_send(1:2,3) = [1,4] ! Send all variables in domain%U

                ! Array bounds to receive -- NULL
                ijk_to_recv = -1
    
                counter = counter + 1
                call domains(j)%nesting%send_comms(counter)%initialise(&
                    my_dx = all_dx(1:2, j, ti), &
                    neighbour_dx = all_dx(1:2, n_ind, n_img), &
                    ijk_to_send = ijk_to_send, &
                    ijk_to_recv = ijk_to_recv, &
                    neighbour_domain_index = n_ind, &
                    my_domain_index = j, &
                    neighbour_domain_comms_index = n_comms , &
                    my_domain_comms_index = counter, &
                    neighbour_domain_image_index = n_img, &
                    my_domain_image_index = ti, &
                    send_only = .true., &
                    use_wetdry_limiting=use_wetdry_limiting, &
                    neighbour_domain_staggered_grid = all_is_staggered_grid(n_ind, n_img),&
                    my_domain_staggered_grid = all_is_staggered_grid(j, ti))

            end do

            !
            ! Set up recvs
            !
            recv_count = nint(max(0.0_dp, maxval(all_recv_metadata(8,:,j,ti))))
            allocate(domains(j)%nesting%recv_comms(recv_count))

            ! Loop over rows of all_recv_metadata(:,:,j,ti), and set up the
            ! recv comms if required
            counter = 0
            do i = 1, size(all_recv_metadata, 2)

                if(all_recv_metadata(8,i,j,ti) < 1) cycle

                ! Neighbour domain index/image
                n_ind = nint( all_recv_metadata(IND,i,j,ti))
                n_img = nint( all_recv_metadata(IMG,i,j,ti))
                n_row = nint( all_recv_metadata(7, i, j, ti)) 
                n_comms = nint(all_send_metadata(8, n_row, n_ind, n_img) )

                ! Check that boxes have the same x/y coordinates to within tolerable round-off
                box_diff = all_recv_metadata(3:6, i    ,     j,    ti) - &
                           all_send_metadata(3:6, n_row, n_ind, n_img)
                if( any( (abs(box_diff) > 100.0_dp*spacing(maxval(abs(all_recv_metadata(3:6, i, j, ti))))) ) ) then
                    ! Looks like the boxes don't match up so well

                    ! But this could be due to periodic boundaries
                    if(any(abs(box_diff(1:2)) > 0.9999_dp*(periodic_xs(2) - periodic_xs(1))) .or. &
                       any(abs(box_diff(3:4)) > 0.9999_dp*(periodic_ys(2) - periodic_ys(1)))) then
                        ! Periodic boundaries. Skip them! 
                    else
                        ! Looks like a genuine error
                        msg = 'Error: recv/send metadata do not appear to have the same boxes to within round-off'
                        write(log_output_unit,*) trim(msg)
                        write(log_output_unit,*) i, j, ti
                        write(log_output_unit,*) n_ind, n_img, n_row, n_comms
                        write(log_output_unit,*) all_recv_metadata(3:6,i,j,ti)
                        write(log_output_unit,*) all_send_metadata(3:6,n_row, n_ind, n_img)
                        write(log_output_unit,*) 100.0_dp*spacing(maxval(abs(all_recv_metadata(3:6, i, j, ti))))
                    end if
                end if

                ! Array bounds to send -- NULL
                ijk_to_send = -1

                ! Array bounds to receive
                ijk_to_recv(1:2,1) = domains(j)%nesting%recv_metadata([XLO, XHI],i)
                ijk_to_recv(1:2,2) = domains(j)%nesting%recv_metadata([YLO, YHI],i)
                ijk_to_recv(1:2,3) = [1,4] ! Recv all variables
    
                counter = counter + 1
                call domains(j)%nesting%recv_comms(counter)%initialise(&
                    my_dx = all_dx(1:2, j, ti), &
                    neighbour_dx = all_dx(1:2, n_ind, n_img), &
                    ijk_to_send = ijk_to_send, &
                    ijk_to_recv = ijk_to_recv, &
                    neighbour_domain_index = n_ind, &
                    my_domain_index = j, &
                    neighbour_domain_comms_index = n_comms , &
                    my_domain_comms_index = counter, &
                    neighbour_domain_image_index = n_img, &
                    my_domain_image_index = ti, &
                    recv_only = .true., &
                    use_wetdry_limiting=use_wetdry_limiting,&
                    neighbour_domain_staggered_grid = all_is_staggered_grid(n_ind, n_img),&
                    my_domain_staggered_grid = all_is_staggered_grid(j, ti))

            end do

        end do
        flush(log_output_unit)

        ! Free memory we don't need. Note some of these variables
        ! will grow in size rapidly with num_images, so good to be careful.
        ! If memory becomes a problem, could consider removing some of the big arrays
        ! in the domains until after the following are deallocated (?)
        if(allocated(all_recv_metadata)) deallocate(all_recv_metadata)
        if(allocated(all_send_metadata)) deallocate(all_send_metadata)
        if(allocated(all_timestepping_refinement_factor)) deallocate(all_timestepping_refinement_factor)

        ! Store the 'all_is_staggered_grid' variable -- in practice it will go inside the multidomain
        ! Cannot just have it in the multidomain from the start, because it might be a coarray
        if(allocated(all_is_staggered_grid)) then
            if(allocated(all_is_staggered_grid_md)) deallocate(all_is_staggered_grid_md)
            allocate(all_is_staggered_grid_md(size(all_is_staggered_grid, 1), size(all_is_staggered_grid,2)))
            all_is_staggered_grid_md = all_is_staggered_grid
            deallocate(all_is_staggered_grid)
        end if

        ! Store all dx in a variable all_dx_md, which in practice will be inside the multidomain
        ! Cannot just have it in the multidomain from the start, because it might be a coarray
        if(allocated(all_dx)) then
            if(allocated(all_dx_md)) deallocate(all_dx_md)  
            allocate( all_dx_md(2, size(all_dx, 2), size(all_dx,3)) )
            all_dx_md = all_dx
            deallocate(all_dx)
        end if

        ! Currently we use all_bbox in a test, but should get rid of that too
        !!if(allocated(all_bbox)) deallocate(all_bbox)

    end subroutine


    ! Make elevation constant in send_regions that go to a single coarser cell,
    ! if the maximum elevation is above elevation_threshold
    ! 
    ! This is done to avoid wet-dry instabilities, caused by aggregating over
    ! wet-and-dry cells on a finer domain, which is then sent to a coarser domain.
    ! Such an operation will break the hydrostatic balance, unless the elevation
    ! in the fine cells is constant
    !
    ! NOTE: CURRENTLY DEFUNCT -- we just send a central fine cell instead of aggregating,
    ! simpler solution to the issues
    !
    ! @param md multidomain
    ! @param elevation_threshold constant -- only do the aggregation if the
    !     elevation is above this threshold. The idea is to make this small enough 
    !     to encompass all potentially wet-dry regions, while not affecting deep
    !     cells. [e.g, a typical value might be -10.0 m or similar]
    !
    subroutine use_constant_wetdry_send_elevation(&
        md, elevation_threshold)

        class(multidomain_type), intent(inout) :: md
        ! Only enforce consistency in regions with (elevation > elevation_threshold)
        ! This should help to better focus on possibly wet-dry regions
        real(dp), intent(in) :: elevation_threshold

        integer(ip) :: nd, j

        nd = size(md%domains)

        ! Loop over domains
        do j = 1, nd
            call md%domains(j)%use_constant_wetdry_send_elevation(elevation_threshold)
        end do

    end subroutine

    ! Print all domain timers, as well as the multidomain timer itself, finalise the domains,
    ! and write max quantities. This is a typical step at the end of a program
    ! 
    subroutine finalise_and_print_timers(md)
        class(multidomain_type), intent(inout) :: md

        integer(ip) :: i

        ! Print out timing info for each
        do i = 1, size(md%domains)
            write(log_output_unit,*) ''
            write(log_output_unit,*) 'Timer of md%domains(', i, ')'
            write(log_output_unit, *) trim(md%domains(i)%output_folder_name)
            write(log_output_unit, *) mod(md%domains(i)%myid, large_64_int), md%domains(i)%local_index
            write(log_output_unit,*) ''
            call md%domains(i)%timer%print(log_output_unit)
            call md%domains(i)%write_max_quantities()
            call md%domains(i)%finalise()
        end do

        write(log_output_unit, *) ''
        write(log_output_unit, *) 'Multidomain timer'
        write(log_output_unit, *) ''
        call md%timer%print(log_output_unit)

    end subroutine


    !
    ! Test codes
    !
    subroutine test_multidomain_mod1()

        integer(ip):: j, i, k, ri(srm), jj
        real(dp):: xmx, xmn, ymx, ymn, val, err, c1, c2

        ! type holding all domains 
        type(multidomain_type) :: md

        ! length/width
        real(dp), parameter, dimension(2):: global_lw = [641000.0_dp, 367000.0_dp]
        ! lower-left corner coordinate
        real(dp), parameter, dimension(2):: global_ll = [185985.0_dp, 8795845.0_dp]
        ! grid size (number of x/y cells)
        integer(ip), parameter, dimension(2):: global_nx = [360_ip, 210_ip] 
        !integer(ip):: global_nx(2)

        ! inner domain in the 'middle 3rd' of the outer domain inside
        integer(ip), parameter :: nest_ratio = 3_ip
        real(dp), parameter, dimension(2) :: inner_lw = global_lw * (1.0_dp /nest_ratio)
        real(dp), parameter, dimension(2) :: inner_ll = global_ll + inner_lw
        integer(ip), parameter, dimension(2) :: inner_nx = global_nx

        integer(ip), allocatable :: reference_results(:,:)
        logical :: in_active

        real(dp) :: bboxes(2, 2, 4)

        ! Clear out any preliminary comms setup
        call deallocate_p2p_comms 

        ! 4 domains in this model
        allocate(md%domains(4))

        ! This is required to match the reference results below (because they were developed when
        ! it was the default approach). 
        md%extra_cells_in_halo = 0_ip

        !
        ! Parent domain
        ! 
        md%domains(1)%lw = global_lw
        md%domains(1)%lower_left =global_ll
        md%domains(1)%nx = global_nx
        md%domains(1)%dx = md%domains(1)%lw / md%domains(1)%nx
        md%domains(1)%timestepping_method = 'euler' !'linear'
        md%domains(1)%timestepping_refinement_factor = 1_ip
        md%domains(1)%dx_refinement_factor = 1.0_dp

        bboxes(1:2, 1, 1) = global_lw
        bboxes(1:2, 2, 1) = global_lw + global_ll

        !
        ! A detailed domain which shares a physical boundary with the global
        ! domain. NOTE: The flow solver doesn't work with this configuration.
        !
        md%domains(2)%lower_left = global_ll 
        md%domains(2)%lw = nint(inner_lw/md%domains(1)%dx)*md%domains(1)%dx ! Is a multiple of parent domain dx
        md%domains(2)%dx_refinement_factor = real(nest_ratio, dp)
        md%domains(2)%dx = md%domains(1)%dx / md%domains(2)%dx_refinement_factor
        md%domains(2)%nx = nint(md%domains(2)%lw / md%domains(2)%dx) ! Is a multiple of nest_ratio
        md%domains(2)%timestepping_method = 'rk2'
        md%domains(2)%timestepping_refinement_factor = 1_ip

        bboxes(1:2, 1, 2) = md%domains(2)%lower_left
        bboxes(1:2, 2, 2) = md%domains(2)%lower_left + md%domains(2)%lw

        !
        ! A detailed domain, fully inside the global domain, with double time-stepping
        !
        md%domains(3)%lower_left = global_ll + nint((inner_ll - global_ll)/md%domains(1)%dx)*md%domains(1)%dx ! Is a multiple of parent domain dx
        md%domains(3)%lw = nint(inner_lw/md%domains(1)%dx)*md%domains(1)%dx ! Is a multiple of parent domain dx
        md%domains(3)%dx_refinement_factor = real(nest_ratio, dp)
        md%domains(3)%dx = md%domains(1)%dx / md%domains(3)%dx_refinement_factor 
        md%domains(3)%nx = nint(md%domains(3)%lw /  md%domains(3)%dx)
        md%domains(3)%timestepping_method = 'rk2'
        md%domains(3)%timestepping_refinement_factor = 2_ip

        bboxes(1:2, 1, 3) = md%domains(3)%lower_left
        bboxes(1:2, 2, 3) = md%domains(3)%lower_left + md%domains(3)%lw

        !
        ! Yet another finer domain
        !
        md%domains(4)%lw = nint(md%domains(3)%lw * [0.3_dp, 0.5_dp] / md%domains(1)%dx)*md%domains(1)%dx ! Is a multiple of parent domain dx
        md%domains(4)%lower_left = global_ll + nint(1.1*md%domains(2)%lw/md%domains(1)%dx)*md%domains(1)%dx  ! Is a multiple of parent domain dx
        md%domains(4)%dx_refinement_factor = real(nest_ratio*nest_ratio, dp)
        md%domains(4)%dx = md%domains(1)%dx / md%domains(4)%dx_refinement_factor
        md%domains(4)%nx = nint(md%domains(4)%lw / md%domains(4)%dx)
        md%domains(4)%timestepping_method = 'rk2'
        md%domains(4)%timestepping_refinement_factor = 5_ip

        bboxes(1:2, 1, 4) = md%domains(4)%lower_left
        bboxes(1:2, 2, 4) = md%domains(4)%lower_left + md%domains(4)%lw


        ! For the tests, we turn off wet-dry limiting in nesting extrapolation, because
        ! this makes it easier to test things analytically
        call md%setup(verbose=.false., use_wetdry_limiting_nesting=.false., capture_log=.false., extra_halo_buffer=0_ip)

        ! Here we 'regression test' against results without coarray features
        if(ni == 1) then
            !
            ! Compare nesting metadata with reference results, domains(1)
            !

            reference_results = reshape( &
                [0, 0,   1, 118,  1,  68, &
                 2, 1, 119, 120,  1,  68, &
                 1, 1, 121, 360,  1,  70, &
                 2, 1,   1, 120, 69,  70, &
                 3, 1, 121, 240, 71,  72, &
                 3, 1, 121, 122, 73, 138, &
                 0, 0, 123, 238, 73, 138, &
                 3, 1, 239, 240, 73, 138, &
                 1, 1,   1, 120, 71, 140, &
                 1, 1, 241, 360, 71, 140, &
                 3, 1, 121, 240,139, 140, &
                 1, 1,   1, 360,141, 210],&
                [6, 12])

            if(all(md%domains(1)%nesting%recv_metadata(1:srm,:) == reference_results)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,& 
__FILE__
            end if
            deallocate(reference_results)

            !
            ! Compare nesting metadata with reference results, domains(2)
            !

            reference_results = reshape( &
                [2, 1,   1, 360,   1, 210, &
                 1, 1, 361, 366,   1, 210, &
                 1, 1,   1, 360, 211, 216, &
                 3, 1, 361, 366, 211, 216], &
                [6, 4])


            if(all(md%domains(2)%nesting%recv_metadata(1:srm,:) == reference_results)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
            end if
            deallocate(reference_results)

            reference_results = reshape(&    
                [2, 1,    1,   9,   1,   9, &
                 1, 1,   10, 378,   1,   9, &
                 3, 1,   10, 369,  10,  30, &
                 4, 1,   46, 153,  31,  39, &
                 4, 1,   46,  54,  40, 126, &
                 0, 0,   55, 144,  40, 126, &
                 4, 1,  145, 153,  40, 126, &
                 3, 1,   10,  45,  31, 135, &
                 3, 1,  154, 369,  31, 135, &
                 4, 1,   46, 153, 127, 135, &
                 1, 1,    1,   9,  10, 219, &
                 1, 1,  370, 378,  10, 219, &
                 3, 1,   10, 369, 136, 219, &
                 1, 1,    1, 378, 220, 228], &
                [6, 14])

            if(all(md%domains(3)%nesting%recv_metadata == reference_results)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
            end if
            deallocate(reference_results)

            reference_results = reshape(&
                [3, 1,   1, 366,   1,  21, &
                 3, 1,   1,  21,  22, 336, &
                 4, 1,  22, 345,  22, 336, &
                 3, 1, 346, 366,  22, 336, &
                 3, 1,   1, 366, 337, 357], &
                [6, 5])

            if(all(md%domains(4)%nesting%recv_metadata(1:srm,:) == reference_results)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
            end if
            deallocate(reference_results)

            ! Loop over all domains
            do i = 1, 4
                ! Loop over all rows in nesting%recv_metadata
                do j = 1, size(md%domains(i)%nesting%recv_metadata(1,:))

                    ! Row i
                    ri = md%domains(i)%nesting%recv_metadata(:,j)                
                    if( ri(1) == 0) cycle
       
                    ! Bounding box of the priority domain associated with row i 
                    xmn = minval(all_bbox(:,1,ri(1),ri(2)))
                    xmx = maxval(all_bbox(:,1,ri(1),ri(2)))
                    ymn = minval(all_bbox(:,2,ri(1),ri(2)))
                    ymx = maxval(all_bbox(:,2,ri(1),ri(2)))

                    ! For each box in the priority domain metadata, check that
                    ! the x/y values on domain i really are inside it.
                    if( all((md%domains(i)%x(ri(3):ri(4)) >= xmn) .and. & 
                            (md%domains(i)%x(ri(3):ri(4)) <= xmx)) .and. &
                        all((md%domains(i)%y(ri(5):ri(6)) >= ymn) .and. &
                            (md%domains(i)%y(ri(5):ri(6)) <= ymx ))) then
                        write(log_output_unit,*) 'PASS'
                    else
                        write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
                    end if

                    ! Now test that priority domain points are NOT inside higher resolution domains
                    ! Loop over all images
                    do k = 1, ni


                        if( ri(1) /= 4) then
                            ! ri is a box with priority domain index = 1, 2, or 3
                            !
                            ! Check that priority_domain points are NEVER inside
                            ! domains(4), which has higher resolution

                            ! Bounding box of the priority domain associated with domains(jj)
                            xmn = minval(all_bbox(:,1,4,k))
                            xmx = maxval(all_bbox(:,1,4,k))
                            ymn = minval(all_bbox(:,2,4,k))
                            ymx = maxval(all_bbox(:,2,4,k))

                            if(any((md%domains(i)%x(ri(3):ri(4)) >= xmn) .and. & 
                                   (md%domains(i)%x(ri(3):ri(4)) <= xmx)) .and. &
                               any((md%domains(i)%y(ri(5):ri(6)) >= ymn) .and. &
                                   (md%domains(i)%y(ri(5):ri(6)) <= ymx )) ) then
                                write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
                            else
                                write(log_output_unit,*) 'PASS'
                            end if

                        end if
                        
                        if( ri(1) == 1 ) then
                            ! Check that priority_domain points are NEVER in
                            ! domains(2:4), which have higher resolution

                            do jj = 2, 4
                                ! Bounding box of the priority domain associated with domains(jj) 
                                xmn = minval(all_bbox(:,1,jj,k))
                                xmx = maxval(all_bbox(:,1,jj,k))
                                ymn = minval(all_bbox(:,2,jj,k))
                                ymx = maxval(all_bbox(:,2,jj,k))

                                if(any((md%domains(i)%x(ri(3):ri(4)) >= xmn) .and. & 
                                       (md%domains(i)%x(ri(3):ri(4)) <= xmx)) .and. &
                                   any((md%domains(i)%y(ri(5):ri(6)) >= ymn) .and. &
                                       (md%domains(i)%y(ri(5):ri(6)) <= ymx )) ) then
                                    write(log_output_unit,*) 'FAIL', __LINE__,&
__FILE__
                                else
                                    write(log_output_unit,*) 'PASS'
                                end if

                            end do
                        end if

                    end do !k
                end do ! j
            end do ! i
        end if

        !
        ! Set up comms
        !
        !call allocate_p2p_comms

        ! Make a planar solution, offset by the domain_index. If we
        ! communicate, then the interopolation / aggregation should not lead
        ! to any error (aside from round-off), so this is a good test case.
        c1 = 0.0001_dp
        c2 = 0.00003_dp
        do j = 1, 4
            do jj = 1, md%domains(j)%nx(2)
                do i = 1, md%domains(j)%nx(1)
                    val = (md%domains(j)%x( i) - global_ll(1))*c1 + &
                          (md%domains(j)%y(jj) - global_ll(2))*c2 + &
                          j
                    md%domains(j)%U(i,jj,1:4) = val + [0, 3, 6, 9]
                end do
            end do

        end do

        call md%send_halos()

#ifdef COARRAY
        sync all
#endif

        call md%recv_halos()

        !
        ! Test the planar solution
        !
        xmx = 0.0_dp ! Store indices with largest error
        ymx = 0.0_dp
        do j = 1, 4
            err = 0.0_dp
            do jj = 1, md%domains(j)%nx(2)
                do i = 1, md%domains(j)%nx(1)

                    !
                    ! Points outside the 'active' region will not have been updated, so
                    ! do not test them
                    !
                    call is_in_active_region(md%domains(j)%nesting%priority_domain_index, &
                        md%domains(j)%nesting%priority_domain_image, j, ti, &
                        md%domains(j)%nest_layer_width, i, jj, in_active)

                    if(.not. in_active) cycle

                    do k = 1, 4
                        val = (md%domains(j)%x( i) - global_ll(1))*c1 + &
                              (md%domains(j)%y(jj) - global_ll(2))*c2 + &
                               md%domains(j)%nesting%priority_domain_index(i,jj) + 3*(k-1) ! 
                        val = abs(md%domains(j)%U(i,jj,k) - val )

                        if(val > err) then
                            xmx = i
                            ymx = jj 
                            err = val
                        end if

                        if(val > 1.0e-04_dp) then
                             write(log_output_unit,*) 'domain: ', j, ti, ', Max err: ', val, i, jj, k,&
                                'Priority-domain: ', md%domains(j)%nesting%priority_domain_index(i,jj),&
                                md%domains(j)%nesting%priority_domain_image(i,jj)
                        end if
                    end do

                end do
            end do
#ifdef COARRAY
            call flush(log_output_unit)
            sync all
#endif
            ! Here if it's working, the 'err' will depend on the precision.
            ! A pragmatic guide to the allowed error is the square root of the distance between
            ! numbers at 1.0.
            if(err*err > spacing(1.0_dp)) then
                write(log_output_unit,*) 'Error summary data:'
                write(log_output_unit,*) 'domain: ', j, ti, ', Max err: ', err, xmx, ymx, global_ll(1), minval(md%domains(j)%x)
                write(log_output_unit,*) 'U:', md%domains(j)%U(nint(xmx), nint(ymx), 1)
                write(log_output_unit,*) 'Var: ', &
                    (md%domains(j)%x(nint(xmx))-global_ll(1))*c1, &
                    (md%domains(j)%y(nint(ymx))-global_ll(2))*c2
                write(log_output_unit, *) 'FAIL', __LINE__, &
__FILE__
                call generic_stop
            else
                write(log_output_unit,*) 'PASS'
            end if
        end do

        ! Some info on storage
        !write(log_output_unit,*) 'Length send buffer: ', size(send_buffer)
        !write(log_output_unit,*) 'Length recv buffer: ', size(recv_buffer)
        !j = 0
        !do i = 1, size(md%domains)
        !    j = j + size(md%domains(i)%U)
        !end do
        !write(log_output_unit,*) 'Sizeof domains: ', j

        call deallocate_p2p_comms

    end subroutine

    !
    ! Check if we can make the multidomain be periodic
    !
    subroutine test_multidomain_mod_periodic()

        ! type holding all domains 
        type(multidomain_type) :: md

        integer(ip) :: nd, i, j, k, ii, region, sgn

        real(dp) :: lw(2) = [360.0_dp, 170.0_dp]
        real(dp) :: ll(2) = [-10.0_dp, -80.0_dp]
        integer(ip) :: nx(2) = [200_ip, 200_ip]
        integer(ip), allocatable :: inds(:)
        real(dp) :: trueval, x, y, err, err_max, xchange, ychange, err_tol
        logical :: passed, periodic_point

        call deallocate_p2p_comms
        
        ! Try to make a periodic domain
        md%periodic_xs = [-10.0_dp, 350.0_dp]
        md%periodic_ys = [-80.0_dp, 90.0_dp]

        nd = 1
        allocate(md%domains(nd))

        md%domains(1)%lw = lw
        md%domains(1)%nx = nx
        md%domains(1)%lower_left = ll
        md%domains(1)%dx = lw/(1.0_dp * nx)
        md%domains(1)%timestepping_method = 'euler'
        md%domains(1)%dx_refinement_factor = 1.0_dp
        md%domains(1)%timestepping_refinement_factor = 1_ip

        call md%setup()
        !call allocate_p2p_comms()

        ! Set initial conditions
        do k = 1, size(md%domains)
            do j = 1, md%domains(k)%nx(2)
                do i = 1, md%domains(k)%nx(1)
                    ! Make up some data
                    x = md%domains(k)%x(i)
                    y = md%domains(k)%y(j)
                    md%domains(k)%U(i,j,1) = x + y*y + 1.0_dp
                    md%domains(k)%U(i,j,2) = x + y*y + 2.0_dp
                    md%domains(k)%U(i,j,3) = x + y*y + 3.0_dp
                    md%domains(k)%U(i,j,4) = x + y*y + 4.0_dp
                end do
            end do
        end do

        !
        ! Make sure flow/elev values are consistent between domains, by doing one
        ! halo exchange
        !
        call md%send_halos()
#ifdef COARRAY
        sync all
#endif
        call md%recv_halos()

        ! After the halo exchange, we should have updated domain%U in the halo region.
        ! Check it worked!
        passed = .true.
        err_max = -1.0_dp
        k = 1
        do j = 1, md%domains(k)%nx(2)
            do i = 1, md%domains(k)%nx(1)
                x = md%domains(k)%x(i)
                y = md%domains(k)%y(j)
                call check_periodic(x, y, md%periodic_xs, md%periodic_ys, periodic_point, adjust_coordinates=.true.)
                trueval = x + y*y

                ! The error
                err = maxval(abs(abs(md%domains(k)%U(i,j,1:4) - (trueval + [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]) )))

                ! Record the difference in the xs, so we note periodic points
                xchange = abs(x-md%domains(k)%x(i))
                ychange = abs(y-md%domains(k)%y(j))

                ! Ideally the error should be zero. In practice, in periodic halo regions the
                ! 'periodic' x and y will slightly differ (due to floating
                ! point precision issues) from the x and y at point which sends data
                !
                ! This is a heuristic tolerance, allowing higher errors in association with changes in x/y due to
                ! periodicity. Instead of (x+y^2) we have ( (x+dx) + (y+dy)^2), suggesting an error of the form 
                ! dx + 2*y*dy (and we use 'spacing' to get an idea of dx, dy). There can also be error due to the
                ! size of trueval. 
                err_tol = 1.0e-10_dp + &
                    20.0_dp*(spacing(abs(x) + abs(y*y)) + (spacing(xchange) + 2.0_dp * ychange * spacing(ychange))) !10.0_dp*spacing(trueval + xchange + ychange) 

    
                ! Track for debugging
                err_max = max(err, err_max)
                !print*, x, y, md%domains(k)%x(i), md%domains(k)%y(j), trueval, err, spacing(trueval + xchange + ychange)
                if(err <= err_tol) then
                    ! This should happen
                else
                    print*, 'err exceeds tol: ', x, y, xchange, ychange, err, err_tol
                    passed = .false.
                end if
            end do
        end do
        !print*, 'err_max: ', err_max, spacing(trueval)

        if(passed) then
            print*, 'PASS'
        else
            print*, 'FAIL' , __LINE__, &
__FILE__
        end if
    
        call deallocate_p2p_comms

    end subroutine

    !
    ! Test for a global multidomain. This has complicated nesting, the refinement is
    ! not an integer divisor of some global dx.
    !
    !
    subroutine test_multidomain_global


        ! Type holding all domains 
        type(multidomain_type) :: md

        ! Change this to decrease the cell size by mesh_refine (i.e. for convergence testing)
        real(dp), parameter :: mesh_refine = 2.0_dp ! 1 --> 4 arc minute

        ! Length/width
        real(dp), parameter, dimension(2):: global_lw = [360.0_dp, 180.0_dp]
        ! Lower-left corner coordinate
        real(dp), parameter, dimension(2):: global_ll = [-180.0_dp, -90.0_dp]

        ! Inner domain 
        integer(ip), parameter :: nest_ratio = 3_ip

        ! Useful misc variables
        integer(ip):: j, i, k, i0, j0, centoff, nd, test_set
        real(dp):: last_write_time, gx(4), gy(4), stage_err, max_residual, roundoff_tol, ci, cj
        character(len=charlen) :: md_file, ti_char
        logical :: has_passed
        real(dp), allocatable :: residual(:)

        ! Clean up
        call deallocate_p2p_comms

        ! Set periodic EW boundary condition
        md%periodic_xs = [global_ll(1), global_ll(1) + global_lw(1)]
        
        
        ! nd domains in this model
        nd = 6
        allocate(md%domains(nd))

        !
        ! Setup basic metadata
        !

        ! Main domain
        md%domains(1)%lw = [global_lw(1), 65.0_dp * 2]
        md%domains(1)%lower_left = [global_ll(1), -65.0_dp]
        md%domains(1)%dx = 1/60.0_dp * 4.0_dp * [1.0_dp, 1.0_dp] / mesh_refine 
        md%domains(1)%nx = nint(md%domains(1)%lw/md%domains(1)%dx)
        md%domains(1)%dx_refinement_factor = 3.0_dp
        md%domains(1)%timestepping_refinement_factor = 1_ip
        md%domains(1)%timestepping_method = 'euler'

        ! North of main domain
        md%domains(2)%lw = [global_lw(1), 82.0_dp - 65.0_dp]
        md%domains(2)%lower_left = [global_ll(1), 65.0_dp]
        md%domains(2)%dx = md%domains(1)%dx*[3,1]
        md%domains(2)%nx = int(nint(md%domains(2)%lw/md%domains(2)%dx))
        md%domains(2)%dx_refinement_factor = 3.0_dp
        md%domains(2)%timestepping_refinement_factor = 1_ip
        md%domains(2)%timestepping_method = 'euler'

        ! North again
        md%domains(3)%lw = [global_lw(1), 88.0_dp - 82.0_dp]
        md%domains(3)%lower_left = [global_ll(1), 82.0_dp]
        md%domains(3)%dx = md%domains(1)%dx*[9,1]
        md%domains(3)%nx = int(nint(md%domains(3)%lw/md%domains(3)%dx))
        md%domains(3)%dx_refinement_factor = 3.0_dp
        md%domains(3)%timestepping_refinement_factor = 2_ip
        md%domains(3)%timestepping_method = 'euler'

        ! Right near north pole
        md%domains(4)%lw = [global_lw(1), 90.0_dp - 88.0_dp]
        md%domains(4)%lower_left = [global_ll(1), 88.0_dp]
        md%domains(4)%dx = md%domains(1)%dx*[27,1]
        md%domains(4)%nx = int(nint(md%domains(4)%lw/md%domains(4)%dx))
        md%domains(4)%dx_refinement_factor = 1.0_dp
        md%domains(4)%timestepping_refinement_factor = 8_ip
        md%domains(4)%timestepping_method = 'euler'

        ! South of main domain. No need to go further south due to Antarctica
        md%domains(5)%lw = [global_lw(1), -65.0_dp - (-82.0_dp)]
        md%domains(5)%lower_left = [global_ll(1), -82.0_dp]
        md%domains(5)%dx = md%domains(1)%dx*[3,1]
        md%domains(5)%nx = int(nint(md%domains(5)%lw/md%domains(5)%dx))
        md%domains(5)%dx_refinement_factor = 1.0_dp
        md%domains(5)%timestepping_refinement_factor = 1_ip
        md%domains(5)%timestepping_method = 'euler'

        
        ! High-res domain fully inside domain 1
        call md%domains(6)%match_geometry_to_parent(&
            parent_domain=md%domains(1), &
            lower_left  = [150._dp, -36.0_dp], &
            upper_right = [152._dp, -34.0_dp], &
            dx_refinement_factor = 3_ip, &
            timestepping_refinement_factor = 2_ip,&
            rounding_method='nearest', &
            recursive_nesting=.false.)
        md%domains(6)%timestepping_method = 'rk2'

        ! Allocate domains and prepare comms
        call md%setup()

        ! Run 3 different tests
        do test_set = 1, 3

            select case(test_set)
                case(1)
                    ! Test 'x' alignment only
                    ci = 1.0_dp
                    cj = 0.0_dp
                case(2)
                    ! Test 'y' alignment only
                    ci = 0.0_dp
                    cj = 1.0_dp
                case(3)
                    ! Test 'xy' alignment
                    ci = 1.0_dp
                    cj = pi
                case default
                    print*, 'Invalid test_set value'
                    call generic_stop
            end select

            ! Set the domain variables based on the coordinates. We can later use 
            ! this to check that comms have been done correctly
            do j0 = 1, nd
                do k = 1, 4
                    do j = 1, md%domains(j0)%nx(2)

                        ! Make some 'random' equation that is a function of x/y coordinates and the variable
                        ! But ensure stage > elev, and also ensure (elev < 0) to prevent the multidomain receive
                        ! from clipping UH with the linear solver
                        md%domains(j0)%U(:,j,k) = ci*md%domains(j0)%x + cj*md%domains(j0)%y(j) - k*10 - 1000.0_dp

                    end do
                end do
            end do

            ! Make sure flow/elev values are consistent between domains, by doing one
            ! halo exchange
            call md%send_halos()
#ifdef COARRAY
            sync all
#endif
            call md%recv_halos()

            !print*, 'DEBUG: ', md%domains(1)%U(1,1,1:4)

            ! At this stage the domains should have consistent values, accounting for the periodicity
            ! The answer should be correct within some roundoff tolerance
            roundoff_tol = 10*spacing(1000.0_dp)

            ! Test domains one-by-one
            do j0 = 1, nd

                has_passed = .true.
                max_residual = -huge(1.0_dp)

                ! Make storage for the test result
                if(allocated(residual)) deallocate(residual)
                allocate(residual(md%domains(j0)%nx(1)))
                do k = 1, 4
                    do j = 1, md%domains(j0)%nx(2)

                        ! Subtract the same 'random' function defined above
                        residual = md%domains(j0)%U(:,j,k) - (ci*md%domains(j0)%x + cj*md%domains(j0)%y(j) - k*10 - 1000.0_dp) 
                        ! Account for periodic
                        where(md%domains(j0)%x < global_ll(1)) residual = residual - 360.0_dp*ci
                        where(md%domains(j0)%x > (global_ll(1) + 360.0_dp)) residual = residual + 360.0_dp*ci

                        !if(j0 == 1 .and. k == 1 .and. j == 1) then
                            !print*, 'AFTER'
                            do i = 1, md%domains(j0)%nx(1)
                                if(abs(residual(i)) < roundoff_tol) cycle
                                print*, test_set, j0, k, j, i, md%domains(j0)%x(i), md%domains(j0)%y(j), &
                                    md%domains(j0)%U(i,j,k), residual(i), roundoff_tol
                            end do
                            !print*, 'END AFTER'
                        !end if

                        if(any(abs(residual) > roundoff_tol )) then
                            has_passed = .false.
                            max_residual = max(max_residual, maxval(abs(residual)))
                        end if

                    end do
                end do
            
                if(has_passed) then
                    print*, 'PASS'
                else
                    print*, 'domain ', j0, max_residual, '; test_set ', test_set
                    print*, 'FAIL ', __LINE__, &
__FILE__
                    call generic_stop
                end if    
            end do

        end do
        
        call deallocate_p2p_comms 

    end subroutine
   
    !
    ! Main test subroutine
    ! 
    subroutine test_multidomain_mod
        call test_multidomain_mod1()
        call test_multidomain_mod_periodic()
        call test_multidomain_global()
    end subroutine

end module
