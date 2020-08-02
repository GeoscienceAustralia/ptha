
module nested_grid_comms_mod
    !!
    !! Module for nesting communication. It defines a type 'two_way_nesting_comms_type'
    !! which can send and/or receive data between corresponding rectangular regions in 2 domains. 
    !! In practice each domain in a multidomain contains two arrays of this type - one to manage
    !! the data it sends to other domains, and one to manage the data it receives from other domains.
    !!
    ! As of 09/2017, I have been using separate instances for 'sending' and
    ! 'receiving'. While not essential, this seems logically cleaner at a higher
    ! level in the code.
    !

    use global_mod, only: dp, ip, charlen, minimum_allowed_depth, send_boundary_flux_data, force_double, &
        real_bytes, integer_bytes, force_double_bytes
    use stop_mod, only: generic_stop
    ! Note: coarray_point2point_comms_mod works in serial or parallel 
    use coarray_point2point_comms_mod, only: include_in_p2p_send_buffer, &
        send_to_p2p_comms, recv_from_p2p_comms, linked_p2p_images, &
        allocate_p2p_comms, deallocate_p2p_comms, communicate_p2p
    use reshape_array_mod, only: repack_rank1_array
    use logging_mod, only: log_output_unit
    use coarray_intrinsic_alternatives, only: this_image2, num_images2
    !use mpi

    !use iso_c_binding, only: C_DOUBLE

    implicit none

    private
    public two_way_nesting_comms_type, domain_nesting_type, test_nested_grid_comms_mod
    public process_received_data

    logical, parameter :: use_averaging_for_fine_to_coarse_data_sends = .false.
    !!
    !! When sending from a fine grid to a coarse grid cell, should we average over all cells
    !! in the fine grid that are inside the coarse grid cell (.TRUE.), or just sent a 
    !! single 'central value' from the fine grid (.FALSE.). 
    !! While the former could be more accurarate, it also may have wet/dry challenges which
    !! are not an issue for the latter.
    !!

    !
    ! Useful constants
    !
    real(dp), parameter :: ONE_dp = 1.0_dp, ZERO_dp = 0.0_dp, HALF_dp = 0.5_dp
    real(dp), parameter :: EPS = 1.0e-06_dp

    integer(ip), parameter, public :: SPATIAL_DIM = 2
    !! Spatial dimensions of the problem -- we solve the 2D shallow water equations

    integer(ip), parameter :: number_of_gradients_coarse_to_fine = 4_ip
    !! In a coarse-to-fine send, we send the coarse data, and a number of gradient terms
    !! which the fine domain can use to interpolate the result. For instance we might send
    !! 2 x-gradients (forward and backward) and 2 y-gradients (forward and backward)

    ! Received data can partially weighted using "receive_weights". This might help
    ! with nesting stability, BUT can cause problems if applied in wet/dry areas. So 
    ! we need some parameters to detect "near-wet-dry" type regions 
    real(dp), parameter :: ignore_receive_weights_depth_shallower_than = 20.0_dp
    real(dp), parameter :: ignore_receive_weights_depth_ratio = 2.0_dp

    ! Indices into boundary flux tracking arrays
    integer(ip), parameter :: NORTH=1_ip, SOUTH=2_ip, EAST=3_ip, WEST=4_ip

    ! Use this to hold flux integrals around the edges of the two_way_nesting_comms
    type :: array_rank2_dp_type
        real(dp), allocatable :: x(:,:)
    end type
    ! As above, for the case we always need high precision
    type :: array_rank2_force_double_type
        real(force_double), allocatable :: x(:,:)
    end type

    type :: array_rank1_logical_type
        logical, allocatable :: x(:)
    end type

    !
    ! The main type of the module
    !
    type :: two_way_nesting_comms_type
        !!
        !! Send and/or receive data between rectangular sub-regions of two domains.
        !! In SWALS, multidomain nesting is achieved by splitting the nesting regions
        !! into rectangular send/receive regions, and managing each of those with this type.
        !!

        ! cell_ratios = my_dx/neighbour_dx
        real(dp) :: cell_ratios(SPATIAL_DIM)
        logical :: my_domain_is_finer 
        logical :: equal_cell_ratios


        !
        ! Variables used to send data
        !

        !
        ! send_inds =  min_i, min_j, min_k,
        !              max_i, max_j, max_k
        integer(ip) :: send_inds(2, SPATIAL_DIM + 1) 
        real(dp), allocatable :: send_buffer(:)
        integer(ip) :: nsend ! We only need to send send_buffer(1:nsend)
        integer(ip) :: nsend_interior ! Number of send values that are not edge fluxes
        logical :: send_active = .true. ! Can use to switch off 'send' (so only receives occur)
        ! Store the integrated flux integral around the send_inds bbox boundary
        ! length(4) = NORTH, SOUTH, EAST, WEST
        type(array_rank2_dp_type) :: send_box_flux_integral(4)
        ! This array will tell us if we **can use** the exterior cells 
        ! that are around the edges of the box. We might want to use
        ! them in derivative calculation.
        !type(array_rank1_logical_type) :: can_use_exterior_cells_send(4)

        !
        ! Variables used to receive data
        !

        !
        ! recv_inds = min_i, min_j, min_k,
        !             max_i, max_j, max_k
        integer(ip) :: recv_inds(2, SPATIAL_DIM + 1) 
        real(dp), allocatable :: recv_buffer(:)
        integer(ip) :: nrecv ! We only need to use recv_buffer(1:nrecv)
        integer(ip) :: nrecv_interior 
        logical :: recv_active = .true. ! Can use to switch off 'recv' (so only sends occur)
        ! Store the flux integral around the recv_inds bbox boundary
        ! length(4) = NORTH, SOUTH, EAST, WEST
        type(array_rank2_dp_type) :: recv_box_flux_integral(4)
        type(array_rank2_force_double_type) :: recv_box_flux_error(4)

        ! We might want to only partially weight the receive data, and partially weight the existing data
        ! This allows that
        real(dp), allocatable :: recv_weights(:,:)
        ! Work array that is useful
        real(dp), allocatable :: recv_work(:,:)
        integer(ip) :: recv_counter = 0_ip

        !
        ! Communication info -- the index of communicating 'domain' in the
        ! vector of 'domains', and the index of the two_way_nesting_comms_type
        ! inside the latter.
        ! We are communicating with domains(neighbour_domain_index) on image
        ! 'neighbour_domain_image_index', specifically with the
        ! 'neighbour_domain_comms_index' index of its two_way_nesting_comms array
        ! i.e., EITHER (if we recv)
        !     domains(neighbour_domain_index)%nesting%send_comms(neighbour_domain_comms_index) 
        !  OR (if we send)
        !     domains(neighbour_domain_index)%nesting%recv_comms(neighbour_domain_comms_index) 
        integer(ip) :: neighbour_domain_index = -1_ip
        integer(ip) :: neighbour_domain_comms_index = -1_ip
        integer(ip) :: my_domain_index = -1_ip
        integer(ip) :: my_domain_comms_index = -1_ip

#if defined(COARRAY) 
        integer(ip) :: neighbour_domain_image_index = -1_ip
        integer(ip) :: my_domain_image_index = -1_ip
#else
        integer(ip) :: neighbour_domain_image_index = 1_ip
        integer(ip) :: my_domain_image_index = 1_ip
#endif

        ! We might interpolate differently if one domain has a staggered grid
        ! and the other doesn't
        integer(ip) :: my_domain_staggered_grid = -1_ip
        integer(ip) :: neighbour_domain_staggered_grid = -1_ip

        ! ID's which are passed to send/recv communication routines
        character(len=charlen) :: send_ID, recv_ID

        logical :: use_wetdry_limiting = .true.

        contains

        procedure :: initialise => initialise_two_way_nesting_comms
        procedure :: process_data_to_send => process_data_to_send
        procedure :: send_data => send_data
        procedure :: process_received_data => process_received_data
        ! Routines to allow time-stepping of boundary flux integral terms
        procedure :: boundary_flux_integral_multiply => boundary_flux_integral_multiply
        procedure :: boundary_flux_integral_tstep => boundary_flux_integral_tstep
        ! Memory summary
        procedure :: memory_summary => two_way_nesting_comms_memory_summary 

    end type

    type :: domain_nesting_type
        !!
        !! Type which holds everything a domain needs to do nesting
        !!

        integer(ip), allocatable :: recv_metadata(:,:), send_metadata(:,:)
            !! 'Tables' with metadata describing the send/recv regions, and regions
            !! they communicate with.

        type(two_way_nesting_comms_type), allocatable :: send_comms(:), recv_comms(:)
            !! Arrays of communication types. Seems easiest to split into separate
            !! send and receive, although the two_way_nesting_comms_type
            !! can do both within a single instance. Each element of the array
            !! takes care of one rectangular send or recv region

        integer(ip), allocatable :: priority_domain_index(:,:), priority_domain_image(:,:)
            !! Arrays describing the 'priority domain' index and its image. The
            !! 'priority domain' is the domain holding the 'correct' representation of the flow at
            !! a point. Domains overlap [e.g. for nesting], and the flow
            !! in some of the overlapping regions might not be realistic (e.g.
            !! nesting buffer regions). However, on the priority domain, it will be
            !! realistic.
        integer(ip), allocatable :: is_priority_domain_not_periodic(:,:)
            !! This will be 1.0 where the domain is priority and not in a periodic region, and 0.0 otherwise.
        integer(ip) :: my_index = 0_ip
        integer(ip) :: my_image = 0_ip

        contains

        !procedure:: print => print_nesting_fluxes
        procedure:: memory_size => domain_nesting_type_memory_size

    end type

    contains

    !
    ! Report on memory use by domain nesting type
    !
    pure subroutine domain_nesting_type_memory_size(domain_nesting, buffer_size)
        class(domain_nesting_type), intent(in) :: domain_nesting
        integer(ip), intent(out) :: buffer_size

        integer(ip) :: i, local_buffer_size

        buffer_size = 0_ip

        if(allocated(domain_nesting%send_comms)) then
            do i = 1, size(domain_nesting%send_comms, kind=ip)
                call domain_nesting%send_comms(i)%memory_summary(local_buffer_size)
                buffer_size = buffer_size + local_buffer_size
            end do        
        end if
        if(allocated(domain_nesting%recv_comms)) then
            do i = 1, size(domain_nesting%recv_comms, kind=ip)
                call domain_nesting%recv_comms(i)%memory_summary(local_buffer_size)
                buffer_size = buffer_size + local_buffer_size
            end do        
        end if

        ! This information is useful, but it's not really 'buffers' -- it's metadata.
        ! That could be important, but consider reporting separately
        !if(allocated(domain_nesting%priority_domain_index)) then       
        !    buffer_size = buffer_size + integer_bytes * (size(domain_nesting%send_metadata) + &
        !        size(domain_nesting%recv_metadata) + size(domain_nesting%priority_domain_index) + &
        !        size(domain_nesting%priority_domain_image))
        !end if

    end subroutine


    ! Set up to way nesting communicator type
    !
    ! @param two_way_nesting_comms The two-way-nesting-comms-type to set up
    ! @param my_dx [dx,dy] of my domain
    ! @param neighbour_dx [dx, dy] of neighbour domain
    ! @param ijk_to_send array with send_inds = [[min_i, min_j, min_k],
    !                                            [max_i, max_j, max_k]]
    !    giving the ranges of indices from which we extract data to send
    !    Ignored if recv_only==.TRUE.
    ! @param ijk_to_recv array with recv_inds = [[min_i, min_j, min_k],
    !                                            [max_i, max_j, max_k]]
    !    giving the ranges of indices where which we will put received data.
    !    Ignored if send_only==.TRUE.
    ! @param neighbour_domain_index integer giving the index of the 'neighbour
    !    domain' [i.e. assuming all domains are in a vector 'domains', and
    !    'neighbour domain' = domains(neighbour_domain_index) ]
    ! @param my_domain_index integer giving the index of 'my domain'
    !    [assuming all domains are in a vector 'domains(:)', and 
    !     'my domain' = domains(my_domain_index)  ]
    ! @param neighbour_domain_image_index optional image_index holding
    !    the 'neighbour domain' (for coarrays)
    ! @param my_domain_image_index optional image_index holding 'my domain' (for
    !    coarrays)
    ! @param send_only logical. Only send data, do not receive any data [so ijk_to_recv is not used]
    ! @param recv_only logical. Only receive data, do not send any data [so ijk_to_send is not used]
    ! @param use_wetdry_limiting logical Limit extrapolation gradients near wet-dry boundaries if doing coarse-to-fine sends
    ! @param neighbour_domain_staggered_grid -- integer flag, 1 means neighbour domain has a staggered grid, 0 otherwise
    ! @param my_domain_staggered_grid -- integer flag, 1 means neighbour domain has a staggered grid, 0 otherwise
    !
    subroutine initialise_two_way_nesting_comms(&
        two_way_nesting_comms, &
        my_dx, neighbour_dx, &
        ijk_to_send, ijk_to_recv, &
        neighbour_domain_index, my_domain_index, &
        neighbour_domain_comms_index, my_domain_comms_index,&
        neighbour_domain_image_index, my_domain_image_index,&
        send_only, recv_only,&
        use_wetdry_limiting, &
        neighbour_domain_staggered_grid, my_domain_staggered_grid)

        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(in) :: my_dx(SPATIAL_DIM), &
                         neighbour_dx(SPATIAL_DIM)
        integer(ip), intent(in) :: ijk_to_send(2, SPATIAL_DIM+1), &
                                   ijk_to_recv(2, SPATIAL_DIM+1)
        integer(ip), intent(in) ::           neighbour_domain_index      , my_domain_index
        integer(ip), intent(in) ::           neighbour_domain_comms_index, my_domain_comms_index
        integer(ip), optional, intent(in) :: neighbour_domain_image_index, my_domain_image_index
        logical, optional, intent(in) :: send_only, recv_only, use_wetdry_limiting
        integer(ip), optional, intent(in) :: neighbour_domain_staggered_grid, my_domain_staggered_grid

        ! Local variables
        real(dp) :: cell_ratios(SPATIAL_DIM), prod_cell_ratios
        integer(ip):: iL, iU, jL, jU , kL, kU, sbs, rbs, bl, bw, bd, i, four_ip(4), dir_ip(4), &
            flux_integral_size
        character(charlen) :: n_char1, n_char2, n_char3, char1, char2, char3
        logical :: equal_cell_ratios

        !
        ! Allow 'send_only' or 'recv_only' instances of this type
        ! 
        if(present(send_only) .and. present(recv_only)) then
            write(log_output_unit,*) 'Error: Both send_only and recv_only were passed to initialise_two_way_nesting_comms' 
            write(log_output_unit,*) 'At most one of these arguments may be passed'
            call generic_stop()
        end if

        if(present(send_only)) then
            if(send_only) then
                two_way_nesting_comms%recv_active = .false.
                two_way_nesting_comms%send_active = .true.
            else
                write(log_output_unit,*) 'Warning: send_only = .false. has no effect.'
                write(log_output_unit,*) 'Use either send_only = .true., or do not pass send_only (it is an optional argument)'
            end if
        end if

        if(present(recv_only)) then
            if(recv_only) then
                two_way_nesting_comms%recv_active = .true.
                two_way_nesting_comms%send_active = .false.
            else
                write(log_output_unit,*) 'Warning: recv_only = .false. has no effect.'
                write(log_output_unit,*) 'Use either recv_only = .true., or do not pass recv_only (it is an optional argument)'
            end if
        end if

        ! By default when doing coarse-to-fine interpolation, limit gradients
        ! used in extrapolation near wet-dry boundaries
        if(present(use_wetdry_limiting)) then
            two_way_nesting_comms%use_wetdry_limiting = use_wetdry_limiting
        else
            two_way_nesting_comms%use_wetdry_limiting = .true.
        end if

        ! If one domain uses a staggered grid and the other doesn't, we might want
        ! to adjust the interpolation procedures
        if(present(my_domain_staggered_grid)) then
            two_way_nesting_comms%my_domain_staggered_grid = my_domain_staggered_grid
        end if
        if(present(neighbour_domain_staggered_grid)) then
            two_way_nesting_comms%neighbour_domain_staggered_grid = neighbour_domain_staggered_grid
        end if


      
        ! Find the ratios of the cell sizes 
        cell_ratios = my_dx/neighbour_dx
        equal_cell_ratios = all(abs(cell_ratios - 1.0_dp) < EPS)

        prod_cell_ratios = product(cell_ratios)

        if( prod_cell_ratios >= ONE_dp .or. equal_cell_ratios) then
            ! The current domain is 'coarser' than the other, or the same size
            two_way_nesting_comms%my_domain_is_finer = .FALSE.

            ! Perform various sanity checks

            ! Both dimensions should be coarser or equal
            if(any(cell_ratios < ONE_dp - EPS)) then
                write(log_output_unit,*) 'Nesting cell dimensions not consistently larger', cell_ratios
                call generic_stop()
            end if
           
            ! cell_ratios should be an integer (up to floating point) 
            if(any(abs(cell_ratios - nint(cell_ratios)) > EPS)) then
                write(log_output_unit,*) 'Apparent non-integer cell ratios (not finer) ', cell_ratios, &
                    my_dx, neighbour_dx
                call generic_stop()
            end if

        else
            ! The current domain is 'finer' than the other
            two_way_nesting_comms%my_domain_is_finer = .TRUE.

            ! Perform various sanity checks

            ! Both dimensions should be finer or equal
            if(cell_ratios(2) > ONE_dp + EPS) then
                write(log_output_unit,*) 'Nesting cell dimensions not consistently smaller'
                call generic_stop()
            end if

            ! 1/cell_ratios  should be an integer (up to floating point)
            if(any(abs(ONE_dp/cell_ratios - nint(ONE_dp/cell_ratios)) > EPS)) then
                write(log_output_unit,*) 'Apparent non-integer cell ratios (finer) ', cell_ratios, ONE_dp/cell_ratios, &
                    my_dx, neighbour_dx
                call generic_stop()
            end if

        end if

        two_way_nesting_comms%cell_ratios = cell_ratios
        two_way_nesting_comms%equal_cell_ratios = equal_cell_ratios

        ! Pass on the comms info
        two_way_nesting_comms%neighbour_domain_index = neighbour_domain_index
        two_way_nesting_comms%my_domain_index = my_domain_index
        two_way_nesting_comms%neighbour_domain_comms_index = neighbour_domain_comms_index
        two_way_nesting_comms%my_domain_comms_index = my_domain_comms_index

        if(present(neighbour_domain_image_index)) then
            two_way_nesting_comms%neighbour_domain_image_index = &
                neighbour_domain_image_index
        end if

        if(present(my_domain_image_index)) then
            two_way_nesting_comms%my_domain_image_index = &
                my_domain_image_index
        end if

        ! Construct character ID's to associated with the send/recv communication
        ! This will uniquely identify the communication
        write(n_char1, '(I0)') two_way_nesting_comms%neighbour_domain_index
        write(n_char2, '(I0)') two_way_nesting_comms%neighbour_domain_comms_index
        write(n_char3, '(I0)') two_way_nesting_comms%neighbour_domain_image_index
        write(char1, '(I0)') two_way_nesting_comms%my_domain_index
        write(char2, '(I0)') two_way_nesting_comms%my_domain_comms_index
        write(char3, '(I0)') two_way_nesting_comms%my_domain_image_index

        ! Set up ijk to send/recv
        two_way_nesting_comms%send_inds = ijk_to_send
        two_way_nesting_comms%recv_inds = ijk_to_recv

        if(two_way_nesting_comms%send_active) then
            ! Shorthand
            iL = ijk_to_send(1,1)
            iU = ijk_to_send(2,1)
            jL = ijk_to_send(1,2)
            jU = ijk_to_send(2,2)
            kL = ijk_to_send(1,3)
            kU = ijk_to_send(2,3)
 
            ! Allocate space to store fluxes around the boundary of the send region
            ! (iL:iU, jL:jU), for all variables from kL:kU

            bl = iU - iL + 1 ! Length of 'i' side of flux box
            bw = jU - jL + 1 ! Length of 'j' side of flux box
            bd = kU - kL + 1 ! Number of variables sent

            ! Allocate NORTH, SOUTH, EAST, WEST flux integral storage
            four_ip = [bl, bl, bw, bw]
            dir_ip = [NORTH, SOUTH, EAST, WEST]
            do i = 1, 4
                allocate(two_way_nesting_comms%send_box_flux_integral(dir_ip(i))%x(four_ip(i), bd))
                two_way_nesting_comms%send_box_flux_integral(dir_ip(i))%x = 0.0_dp
            end do

            ! Number of values to add to the send array, considering that we
            ! will transform them to match recv grid. 
            flux_integral_size =  count([send_boundary_flux_data]) * &
                (sum(int(nint(four_ip(1:2) * bd * cell_ratios(1))))  + &
                 sum(int(nint(four_ip(3:4) * bd * cell_ratios(2)))))

            ! Record whether we can extrapolate around the outer edge
            !do i = 1, 4
            !    allocate(two_way_nesting_comms%can_use_exterior_cells_send(dir_ip(i))%x(four_ip(i)))
            !end do
            
            ! Sanity check the send/recv indices
            if(two_way_nesting_comms%my_domain_is_finer) then 
                ! Finer grid should be perfectly contained in coarser grid
                ! (i.e. iU - iL + 1 should be a multiple of the relative cell sizes)
                if(mod(iU - iL + 1, int(nint(ONE_dp/cell_ratios(1)))) /= 0) then
                    write(log_output_unit,*) 'Imperfectly nested i'
                    write(log_output_unit,*) iU, iL, iU - iL + 1, int(nint(ONE_dp/cell_ratios(1))) 
                    write(log_output_unit,*) neighbour_domain_comms_index, my_domain_comms_index
                    call generic_stop()
                end if

                if(mod(jU - jL + 1, int(nint(ONE_dp/cell_ratios(2)))) /= 0) then
                    write(log_output_unit,*) 'Imperfectly nested j'
                    call generic_stop()
                end if
            end if

            if( (iL > iU).OR.(jL > jU).OR.(kL > kU) ) then
                write(log_output_unit,*) 'indices to send has min > max'
                call generic_stop()
            end if
            
            ! Figure out how much data we need to send 
            ! Make this the same length as the region we will send to
            if(two_way_nesting_comms%my_domain_is_finer .or. equal_cell_ratios) then
                sbs = (iU - iL + 1) * (jU - jL + 1) * (kU - kL + 1) / &
                    int(nint(product(ONE_dp/cell_ratios)))
            else
                ! Sends from a coarse domain to a fine domain should send
                ! the coarse data, and a number of 'gradient' terms that can
                ! be used by the recipient fine domain to reconstruct the data
                ! at all fine cells. 
                sbs = (iU - iL + 1) * (jU - jL + 1) * (kU - kL + 1) * &
                    (1_ip + number_of_gradients_coarse_to_fine)
            end if
            sbs = sbs + flux_integral_size 

            two_way_nesting_comms%nsend = sbs
            two_way_nesting_comms%nsend_interior = sbs - flux_integral_size

        else
            ! Don't send anything
            sbs = 0
            two_way_nesting_comms%nsend = sbs
            two_way_nesting_comms%nsend_interior = sbs
        end if

        if(two_way_nesting_comms%recv_active) then
            ! Make space to store boundary flux integral around recv region
            bl = ijk_to_recv(2,1) - ijk_to_recv(1,1) + 1
            bw = ijk_to_recv(2,2) - ijk_to_recv(1,2) + 1
            bd = ijk_to_recv(2,3) - ijk_to_recv(1,3) + 1

            ! Allocate NORTH, SOUTH, EAST, WEST flux integral storage
            four_ip = [bl, bl, bw, bw]
            dir_ip = [NORTH, SOUTH, EAST, WEST]
            do i = 1, 4
                allocate(two_way_nesting_comms%recv_box_flux_integral(dir_ip(i))%x(four_ip(i), bd))
                two_way_nesting_comms%recv_box_flux_integral(dir_ip(i))%x = 0.0_dp
                ! Also make space for the flux 'error'
                allocate(two_way_nesting_comms%recv_box_flux_error(dir_ip(i))%x(four_ip(i), bd))
                two_way_nesting_comms%recv_box_flux_error(dir_ip(i))%x = 0.0_force_double
            end do
    
            ! Number of values to add to recv array
            flux_integral_size = sum(four_ip*bd)*count([send_boundary_flux_data])

            if(two_way_nesting_comms%my_domain_is_finer .and. (.not. equal_cell_ratios)) then
                ! Receive one value and a number of gradients for every coarse cell
                rbs = int(nint(&
                    (product(ijk_to_recv(2,:) - ijk_to_recv(1,:) + 1) * product(cell_ratios)) * &
                    (1_ip + number_of_gradients_coarse_to_fine) ))
            else
                rbs = product(ijk_to_recv(2,:) - ijk_to_recv(1,:) + 1)
            endif

            ! The full recv buffer also includes the flux integral
            rbs = rbs + flux_integral_size

            ! Setup recv-weights. This will usually be 1.0, but other values can allow for
            ! complex updating in the recv region
            allocate(two_way_nesting_comms%recv_weights(&
                ijk_to_recv(1,1):ijk_to_recv(2,1), ijk_to_recv(1,2):ijk_to_recv(2,2)))
            two_way_nesting_comms%recv_weights = 1.0_dp

            allocate(two_way_nesting_comms%recv_work(&
                ijk_to_recv(1,1):ijk_to_recv(2,1), ijk_to_recv(1,2):ijk_to_recv(2,2)))
            two_way_nesting_comms%recv_work = -HUGE(1.0_dp)

        else
            ! Don't recv anything
            rbs = 0
            flux_integral_size = 0

            ! Setup recv-weights
            allocate(two_way_nesting_comms%recv_weights(0, 0))
            allocate(two_way_nesting_comms%recv_work(0, 0))
        end if

        two_way_nesting_comms%nrecv = rbs
        two_way_nesting_comms%nrecv_interior = rbs - flux_integral_size

        ! Allocate buffers    
        allocate(two_way_nesting_comms%send_buffer(sbs))
        allocate(two_way_nesting_comms%recv_buffer(rbs))

        ! Set boundary flux integral terms to zero
        call two_way_nesting_comms%boundary_flux_integral_multiply(0.0_dp)

        !
        ! Prepare communication
        !
        
        if(two_way_nesting_comms%send_active) then

            two_way_nesting_comms%send_ID = 'nesting_from_' // (trim(char1)) // '_' // &
                (trim(char2)) // '_' // (trim(char3)) // '_to_' // &
                (trim(n_char1)) // '_' // (trim(n_char2)) // '_' // &
                (trim(n_char3))
       
            ! Make sure the point2point coarray communicator has space for this variable
            call include_in_p2p_send_buffer(two_way_nesting_comms%send_buffer, &
                buffer_label = two_way_nesting_comms%send_ID, &
                receiver_image = two_way_nesting_comms%neighbour_domain_image_index)

        else
            two_way_nesting_comms%send_ID = 'inactive'
        end if

        if(two_way_nesting_comms%recv_active) then
            ! Note the ordering of this is reversed compared with above, so it will match the
            ! 'send_ID' of another two_way_nesting_comms instance which is sending the data 
            two_way_nesting_comms%recv_ID = 'nesting_from_' // (trim(n_char1)) // '_' // &
                (trim(n_char2)) // '_' // (trim(n_char3)) // '_to_' // &
                (trim(char1)) // '_' // (trim(char2)) // '_' // &
                (trim(char3))
        else
            two_way_nesting_comms%recv_ID = 'inactive'
        end if


    end subroutine



    !
    ! Report on amount of memory used in nesting
    !
    pure subroutine two_way_nesting_comms_memory_summary(two_way_nesting_comms, buffer_size)
        class(two_way_nesting_comms_type), intent(in) :: two_way_nesting_comms
        integer(ip), intent(out) :: buffer_size
        integer(ip) :: i

        buffer_size = 0

        if(allocated(two_way_nesting_comms%send_buffer)) then
            buffer_size = buffer_size + size(two_way_nesting_comms%send_buffer, kind=ip)*real_bytes
        end if

        if(allocated(two_way_nesting_comms%recv_buffer)) then
            buffer_size = buffer_size + size(two_way_nesting_comms%recv_buffer, kind=ip)*real_bytes
        end if

        do i = 1, size(two_way_nesting_comms%send_box_flux_integral)        
            if(.not. allocated(two_way_nesting_comms%send_box_flux_integral(i)%x)) continue
            buffer_size = buffer_size + size(two_way_nesting_comms%send_box_flux_integral(i)%x, kind=ip)*real_bytes
        end do

        do i = 1, size(two_way_nesting_comms%recv_box_flux_integral)        
            if(.not. allocated(two_way_nesting_comms%recv_box_flux_integral(i)%x)) continue
            buffer_size = buffer_size + size(two_way_nesting_comms%recv_box_flux_integral(i)%x, kind=ip)*real_bytes
        end do

        do i = 1, size(two_way_nesting_comms%recv_box_flux_error)        
            if(.not. allocated(two_way_nesting_comms%recv_box_flux_error(i)%x)) continue
            buffer_size = buffer_size + size(two_way_nesting_comms%recv_box_flux_error(i)%x, kind=ip)*force_double_bytes
        end do

        if(allocated(two_way_nesting_comms%recv_weights)) then
           buffer_size = buffer_size + size(two_way_nesting_comms%recv_weights, kind=ip)*real_bytes 
        end if

        if(allocated(two_way_nesting_comms%recv_work)) then
           buffer_size = buffer_size + size(two_way_nesting_comms%recv_work, kind=ip)*real_bytes 
        end if

    end subroutine

    
    ! Copy data to the send buffer, with appropriate averaging or interpolation
    !
    ! @param two_way_nesting_comms the communicator
    ! @param U real rank 3 array that we send some subset of
    !
    subroutine process_data_to_send(two_way_nesting_comms, U)
        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(in) :: U(:,:,:)

        ! Local variables
        integer(ip) :: i,j,k, nb, nvar, ii, jj
        integer(ip) :: iL, iU, jL, jU, kL, kU, n0, n1
        integer(ip) :: iCounter, jCounter, kCounter, ijkCounter
        integer(ip) :: iR, jR, kR
        integer(ip) :: inv_cell_ratios_ip(2), cell_ratios_ip(2)
        integer(ip) :: already_sent, send_ind_offset
        real(dp) :: product_cell_ratios, inv_product_cell_ratios
        real(dp) :: inv_cell_ratios(2), cell_ratios(2)
        real(dp) :: dU_di_p, dU_di_m, dU_dj_p, dU_dj_m, dep, uc, inv_stride
        integer(ip) :: imn, imx, jmn, jmx, iim1, iip1, jjm1, jjp1, dir_ip(4), stride(4), central_i, central_j
        integer(ip) :: ip1, jp1, im1, jm1, i1, j1, im2, jm2
        integer(ip) :: my_domain_staggered_grid, neighbour_domain_staggered_grid, jouter, jinner
        real(dp) :: gradient_scale, depth_min, depth_max, del
        real(dp) :: gradient_scale_x, gradient_scale_y, depth_min_x, depth_max_x, depth_min_y, depth_max_y
        logical :: equal_cell_ratios
        ! parameters controlling how gradients are limited when doing coarse-to-fine interpolation
        real(dp), parameter :: depth_limit_upper_threshold = 0.25_dp 
        real(dp), parameter :: depth_limit_lower_threshold = 0.05_dp
        real(dp), parameter :: smooth_scale_coarse2fine = 0.0_dp
        ! Indices for stage and elevation in U. This is used to do 'wet-dry-limiting'
        ! of gradients in coarse-to-fine sends (which involve interpolation)
        integer(ip), parameter :: STG=1, UH=2, VH=3, ELV=4

        ! Quick exit if the type is purely used to recv data
        if(.not. two_way_nesting_comms%send_active) return


        imn = 1
        jmn = 1
        imx = size(U, 1, kind=ip)
        jmx = size(U, 2, kind=ip)

        ! Lower/upper spatial indices involved in nesting
        iL = two_way_nesting_comms%send_inds(1, 1)
        iU = two_way_nesting_comms%send_inds(2, 1)
        jL = two_way_nesting_comms%send_inds(1, 2)
        jU = two_way_nesting_comms%send_inds(2, 2)

        ! Lower/upper variable indices involved in nesting.
        ! (NOTE, the code logic currently assumes a sequence of variables
        !  is sent, so e.g. if variable 2,4 are sent, then variable 3 should be
        !  as well). kL and kU should be in the range of STG-ELV (1 to 4)
        kL = two_way_nesting_comms%send_inds(1, 3)
        kU = two_way_nesting_comms%send_inds(2, 3)
        ! Number of variables that are sent
        kR = (kU - kL + 1)

        cell_ratios = two_way_nesting_comms%cell_ratios

        ! Optimized special case when cell ratios are equal
        equal_cell_ratios = two_way_nesting_comms%equal_cell_ratios

        ! Record whether 'my domain' and the 'neighbour domain' are staggered
        ! grids or not. Complexities arise if one is staggered, and the other not
        my_domain_staggered_grid = two_way_nesting_comms%my_domain_staggered_grid
        neighbour_domain_staggered_grid = two_way_nesting_comms%neighbour_domain_staggered_grid

        if(equal_cell_ratios) then
            ! Attempt to optimize the important 'equal cell size' case
            ! Define useful constants
            inv_cell_ratios = 1.0_dp
            inv_cell_ratios_ip = 1_ip
            cell_ratios_ip = 1_ip
            product_cell_ratios = 1.0_dp 
            inv_product_cell_ratios = 1.0_dp
            ! Number of i/j cells in receive zone
            iR = iU - iL + 1
            jR = jU - jL + 1
        else
            ! Define useful constants
            inv_cell_ratios = ONE_dp / cell_ratios 
            inv_cell_ratios_ip = int(nint(inv_cell_ratios))
            cell_ratios_ip = int(nint(cell_ratios))
            product_cell_ratios = cell_ratios(1) * cell_ratios(2) !product(cell_ratios)
            inv_product_cell_ratios = ONE_dp/product_cell_ratios
            ! Number of i/j cells in receive zone
            iR = int( nint((iU - iL + 1) * cell_ratios(1)) ) 
            jR = int( nint((jU - jL + 1) * cell_ratios(2)) ) 
        end if


        ! **Clear** the send buffer
        two_way_nesting_comms%send_buffer = ZERO_dp

        ! We use iCounter, jCounter, kCounter to track the index of the send
        ! buffer that needs modification. 
        iCounter = 0
        jCounter = 0
        kCounter = 0



        ! Loop over spatial dimensions, with operations depending on
        ! whether we do a 'fine-to-coarse' or a 'coarse-to-fine' data send

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(two_way_nesting_comms, U), &
        !$OMP SHARED(kL, jL, iL, kU, jU, iU, kR, jR, iR, imn, imx, jmn, jmx), &
        !$OMP SHARED(cell_ratios, inv_cell_ratios, inv_cell_ratios_ip, cell_ratios_ip), &
        !$OMP SHARED(inv_product_cell_ratios, product_cell_ratios, equal_cell_ratios), &
        !$OMP SHARED(my_domain_staggered_grid, neighbour_domain_staggered_grid)
        if(two_way_nesting_comms%my_domain_is_finer .or. equal_cell_ratios ) then
            ! Fine-to-coarse


            ! Loop over variables (third index of U)
            !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
            do k = kL, kU
                !
                ! Send the volume-integrated U values
                !
                !do j = jL, jU
                !! Note -- here we aim to achieve the same as "do j = jL, jU", but the loop is split
                !! so that a single openmp threat only ever updates an entry of two_way_nesting_comms%send_buffer.
                !! This helps avoid openmp non-reproducibility due to ordering of additions
                do jouter = jL, jU, inv_cell_ratios_ip(2)
                do jinner = 0, inv_cell_ratios_ip(2)-1
                    j = jouter + jinner
               
                    ! NOTE: Put this inside the loop, so we can collapse the j/k loops with OMP
                    kCounter = (k - kL) + 1

                    if(equal_cell_ratios) then
                        ! Attempt quick approach for important case with equal cell sizes
                        !
                        ! This assumes we do not have a staggered to non-staggered send, or vice versa.
                        !
                        jCounter = (j-jL) + 1
                        send_ind_offset = ( (jCounter - 1) + (kCounter - 1) * jR ) * iR
                        two_way_nesting_comms%send_buffer((send_ind_offset + 1):(send_ind_offset+1+(iU-iL))) = &
                            U(iL:iU, j, k)
                    else
                        !
                        ! Unequal cell sizes
                        !

                        ! Since the send_buffer is a rank 1 array, we have to
                        ! pack the data into it
                        ! Every time j-jL grows by 1/cell_ratios(2), we are
                        ! modifying another 'jCounter' index in the send buffer
                        ! Deliberate integer division
                        jCounter = (j-jL)/inv_cell_ratios_ip(2) + 1

                        ! Remove the following computation from the inner loop
                        send_ind_offset = ( (jCounter - 1) + (kCounter - 1) * jR ) * iR

                        do i = iL, iU

                            ! Since the send_buffer is a rank 1 array, we have to
                            ! pack the data into it
                            ! Every time i-iL grows by 1/cell_ratios(1), we are
                            ! modifying another 'iCounter' index in the send buffer
                            ! Deliberate integer division
                            iCounter = (i - iL)/inv_cell_ratios_ip(1) + 1
                            ijkCounter = iCounter + send_ind_offset

                            if(use_averaging_for_fine_to_coarse_data_sends) then
                                ! Average all 'fine' cells inside the coarse cell
                                ! This has not been adapted to deal with staggered-to-non-staggered transfers
                                two_way_nesting_comms%send_buffer(ijkCounter) = &
                                    two_way_nesting_comms%send_buffer(ijkCounter) + &
                                    product_cell_ratios * U(i,j,k)
                            else

                                ! Get indices of the fine cell 'inside the coarse cell'.
                                !   central_i = 0, 1, 2, ... inv_cell_ratios_ip(1) - 1, 0, 1, 2, ..
                                !   central_j = 0, 1, 2, ... inv_cell_ratios_ip(2) - 1, 0, 1, 2, ..
                                central_i = mod(i - iL , inv_cell_ratios_ip(1)) 
                                central_j = mod(j - jL , inv_cell_ratios_ip(2))

                                if(my_domain_staggered_grid == 1 .and. neighbour_domain_staggered_grid == 1) then
                                    !
                                    ! Staggered sends to staggered
                                    !

                                    if(k == STG .or. k == ELV) then
                                        if( (central_i == inv_cell_ratios_ip(1)/2) .and. &
                                            (central_j == inv_cell_ratios_ip(2)/2) ) then
                                            two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                        end if
                                    else if(k == UH) then
                                        ! East point. 
                                        if(central_i == inv_cell_ratios_ip(1)-1) then
                                            ! Average value of UH along 'y' direction inside cell
                                            ! Avoid using points on the edge of nest box (for stability reasons)
                                            !if(i < iU) then
                                                two_way_nesting_comms%send_buffer(ijkCounter) = &
                                                    two_way_nesting_comms%send_buffer(ijkCounter) + &
                                                    U(i,j,k) * cell_ratios(2)
                                            !else
                                            !    ! Extrapolation here is not always stable
                                            !    two_way_nesting_comms%send_buffer(ijkCounter) = &
                                            !        two_way_nesting_comms%send_buffer(ijkCounter) + &
                                            !        !(2.0_dp*U(i-1,j,k) - 1.0_dp*U(i-2,j,k)) * cell_ratios(2)
                                            !        U(i-1,j,k) * cell_ratios(2)
                                            !end if
                                        end if
                                    else if(k == VH) then
                                        ! North point.
                                        if(central_j == inv_cell_ratios_ip(2) - 1) then
                                            ! Average value of VH along 'x' direction inside cell
                                            ! Avoid using points on the edge of nest box (for stability reasons)
                                            !if(j < jU) then
                                                two_way_nesting_comms%send_buffer(ijkCounter) = &
                                                    two_way_nesting_comms%send_buffer(ijkCounter) + &
                                                    U(i,j,k) * cell_ratios(1)
                                            !else
                                            !    ! Extrapolation here is not always stable
                                            !    two_way_nesting_comms%send_buffer(ijkCounter) = &
                                            !        two_way_nesting_comms%send_buffer(ijkCounter) + &
                                            !        !(2.0_dp * U(i,j-1,k) - 1.0_dp*U(i,j-2,k)) * cell_ratios(1)
                                            !        U(i,j-1,k) * cell_ratios(1)
                                            !end if
                                        end if
                                    end if

                                else if(my_domain_staggered_grid == neighbour_domain_staggered_grid) then
                                    !
                                    ! Probably non-staggered to non-staggered (or else unspecified, default)
                                    !
                                    ! Just send central value
                                    if( (central_i == inv_cell_ratios_ip(1)/2) .and. &
                                        (central_j == inv_cell_ratios_ip(2)/2) ) then
                                        two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                    end if

                                else if(my_domain_staggered_grid == 0 .and. neighbour_domain_staggered_grid == 1) then
                                    !
                                    ! Non-staggered to staggered
                                    !

                                    if(k == STG .or. k == ELV) then
                                        if( (central_i == inv_cell_ratios_ip(1)/2) .and. &
                                            (central_j == inv_cell_ratios_ip(2)/2) ) then
                                            two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                        end if
                                    else if(k == UH) then
                                        ! East point. 
                                        if( (central_i == inv_cell_ratios_ip(1)-1) ) then
                                            im1 = max(i-1, 1) ! Safe (i-1) index
                                            ! Avoid using points on the edge of nest box (for stability regions)
                                            !if(i < iU) then
                                                two_way_nesting_comms%send_buffer(ijkCounter) = &
                                                    two_way_nesting_comms%send_buffer(ijkCounter) + &
                                                    ! Imprecise but try for stability
                                                    (1.0_dp * U(i,j,k) - 0.0_dp*U(im1,j,k))*cell_ratios(2)
                                            !else
                                            !    im2 = max(i-2, 1) ! Safe (i-2) index
                                            !    two_way_nesting_comms%send_buffer(ijkCounter) = &
                                            !        two_way_nesting_comms%send_buffer(ijkCounter) + &
                                            !        ! Imprecise but try for stability
                                            !        (1.0_dp * U(im1,j,k) - 0.0_dp*U(im2,j,k))*cell_ratios(2)
                                            !end if
                                        end if
                                        !if( (central_i == inv_cell_ratios_ip(1)-1) .and. & 
                                        !    (central_j == inv_cell_ratios_ip(2)/2) ) then
                                        !    two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                        !end if
                                    else
                                        ! North point.
                                        if( (central_j == inv_cell_ratios_ip(2)-1) ) then
                                            ! Extrapolation from the '-2' point. These choices are for stability
                                            jm1 = max(j-1, 1) ! Safe (j-1) index
                                            ! Avoid using points on the edge of nest box (for stability regions)
                                            !if(j < jU) then
                                                two_way_nesting_comms%send_buffer(ijkCounter) = &
                                                    two_way_nesting_comms%send_buffer(ijkCounter) + &
                                                    ! Imprecise but try for stability
                                                    (1.0_dp * U(i,j,k) - 0.0_dp*U(i,jm1,k))*cell_ratios(1)
                                            !else
                                            !    jm2 = max(j-2, 1) ! Safe (j-2) index
                                            !    two_way_nesting_comms%send_buffer(ijkCounter) = &
                                            !        two_way_nesting_comms%send_buffer(ijkCounter) + &
                                            !        ! Imprecise but try for stability
                                            !        (1.0_dp * U(i,jm1,k) - 0.0_dp*U(i,jm2,k))*cell_ratios(1)
                                            !end if
                                        end if

                                        !if( (central_i == inv_cell_ratios_ip(1)/2) .and. & 
                                        !    (central_j == inv_cell_ratios_ip(2)-1) ) then
                                        !    two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                        !end if
                                    end if

                                else
                                    ! Generic fallback case
                                    ! Pointwise elevation, averaging for UH/VH
                                    if( (k == STG) .or. (k == ELV)) then
                                        if( (central_i == inv_cell_ratios_ip(1)/2) .and. &
                                            (central_j == inv_cell_ratios_ip(2)/2) ) then
                                            two_way_nesting_comms%send_buffer(ijkCounter) = U(i,j,k)
                                        end if
                                    else
                                        ! Averaging -- more stable
                                        two_way_nesting_comms%send_buffer(ijkCounter) = &
                                            two_way_nesting_comms%send_buffer(ijkCounter) + &
                                            product_cell_ratios * U(i,j,k)
                                    end if

                                end if
                            
                            end if
                        end do !i
                    end if
                !end do !j
                end do !jouter
                end do !jinner
            end do !k
            !$OMP END DO
        else
            !
            ! Coarse-to-fine
            !

            ! Loop with 'k' as an inner variable, so that we only need to
            ! compute gradient_scale once for eack i,j pair 
            !$OMP DO SCHEDULE(STATIC), COLLAPSE(2)
            do j = jL, jU
                do i = iL, iU

                    iim1 = max(imn, i-1)
                    iip1 = min(imx, i+1)
                    jjm1 = max(jmn, j-1)
                    jjp1 = min(jmx, j+1)
                    gradient_scale = 1.0_dp 
                    gradient_scale_x = 1.0_dp 
                    gradient_scale_y = 1.0_dp 

                    if(two_way_nesting_comms%use_wetdry_limiting) then
                        ! Reduce "gradient_scale" if the depth is rapidly varying

                        depth_max = -HUGE(1.0_dp)
                        depth_min = HUGE(1.0_dp)
                        do jj = jjm1, jjp1
                            do ii = iim1, iip1
                                dep = U(ii,jj,STG) - U(ii,jj,ELV)  ! STAGE - ELEVATION  
                                depth_max = max(depth_max, dep)
                                depth_min = min(depth_min, dep)
                            end do
                        end do
                        !
                        ! If large depth gradients, then linearly
                        ! transition to 0 as depth_min/depth_max passes
                        ! between 2 pre-specified thresholds
                        if (depth_min < depth_limit_upper_threshold*depth_max) then
                            gradient_scale = &
                                max(depth_min/depth_max - depth_limit_lower_threshold, 0.0_dp) / &
                                (depth_limit_upper_threshold - depth_limit_lower_threshold)
                                
                        end if
                        ! negative depth should not occur, but this extrapolates gracefully in any case. 
                        if(depth_max < 0.0_dp .or. depth_min < 0.0_dp) gradient_scale = 0.0_dp

                        !
                        ! As above, but only in x-direction
                        !

                        depth_max_x = maxval(U(iim1:iip1, j, STG) - U(iim1:iip1, j, ELV))
                        depth_min_x = minval(U(iim1:iip1, j, STG) - U(iim1:iip1, j, ELV))
                        if (depth_min_x < depth_limit_upper_threshold*depth_max_x) then
                            gradient_scale_x = &
                                max(depth_min_x/depth_max_x - depth_limit_lower_threshold, 0.0_dp) / &
                                (depth_limit_upper_threshold - depth_limit_lower_threshold)
                                
                        end if
                        ! negative depth should not occur, but this extrapolates gracefully in any case. 
                        if(depth_max_x < 0.0_dp .or. depth_min_x < 0.0_dp) gradient_scale_x = 0.0_dp


                        !
                        ! As above, but only in y-direction
                        !

                        depth_max_y = maxval(U(i, jjm1:jjp1, STG) - U(i, jjm1:jjp1, ELV))
                        depth_min_y = minval(U(i, jjm1:jjp1, STG) - U(i, jjm1:jjp1, ELV))
                        if (depth_min_y < depth_limit_upper_threshold*depth_max_y) then
                            gradient_scale_y = &
                                max(depth_min_y/depth_max_y - depth_limit_lower_threshold, 0.0_dp) / &
                                (depth_limit_upper_threshold - depth_limit_lower_threshold)
                                
                        end if
                        ! negative depth should not occur, but this extrapolates gracefully in any case. 
                        if(depth_max_y < 0.0_dp .or. depth_min_y < 0.0_dp) gradient_scale_y = 0.0_dp


                    end if

                    !! Test for round-off related effects.
                    !gradient_scale_x = 0.0_dp
                    !gradient_scale_y = 0.0_dp

                    ! Put k in the inner loop, so we don't have to keep recomputing 'gradient_scale'
                    do k = kL, kU
                        !
                        ! There would seem to be logical issues with using points if we might get them from other domains.
                        if(i > iim1 .and. i < iip1) then
                        !if(i > iL .and. i < iU) then
                            ! Define positive and negative derivatives in x direction
                            dU_di_p = U(i+1, j, k) - U(i,j,k)
                            dU_di_m = U(i, j, k) - U(i-1,j,k)
                        else
                            ! Edge casees
                            if(i == iim1) then
                            !if(i == iL) then
                                dU_di_p = U(i+1, j, k) - U(i, j, k)
                            else
                                dU_di_p = U(i, j, k) - U(i-1, j, k)
                            end if
                            dU_di_m = dU_di_p
                        end if

                        if(j > jjm1 .and. j < jjp1) then
                        !if(j > jL .and. j < jU) then
                            ! Define positive and negative derivatives in y direction
                            dU_dj_p = U(i, j+1, k) - U(i,j,k)
                            dU_dj_m = U(i, j, k) - U(i,j-1,k)
                        else
                            ! Edge cases
                            if(j == jjm1) then
                            !if(j == jL) then
                                dU_dj_p = U(i, j+1, k) - U(i, j, k)
                            else
                                dU_dj_p = U(i, j, k) - U(i, j-1, k)
                            end if
                            dU_dj_m = dU_dj_p
                        end if

                        uc = U(i,j,k)

                        if(my_domain_staggered_grid == 1) then
                            ! Sending from a coarse staggered grid to a finer grid.
                            ! The finer grid will expect the UH/VH components to be nearer the middle of the cell

                            ! Make sure we only use points that really are in the 'send' zone. 
                            !! Potential instabilities / logic issues if we use points received from other domains (?).
                            ip1 = iip1 !min(iip1, iU) !min(i+1, iip1)
                            im1 = iim1 !max(iim1, iL) !max(i-1, iim1) 
                            jp1 = jjp1 !min(jjp1, jU) !min(j+1, jjp1)
                            jm1 = jjm1 !max(jjm1, jL) !max(j-1, jjm1) 

                            if(k == UH) then
                                ! UH is offset east of where the receive grid wants it 

                                if(neighbour_domain_staggered_grid == 1) then
                                    ! Slightly east of middle of the cell
                                    ! But not stable to do this
                                    del = 0.5_dp !- inv_cell_ratios(1)*0.5_dp
                                else
                                    ! Middle of the cell
                                    del = 0.5_dp
                                end if

                                ! Move uh back 
                                uc = uc - dU_di_m*del
                                ! There is only one relevent EW gradient -- the 'minus' one
                                dU_di_p = dU_di_m 

                                !if(neighbour_domain_staggered_grid /= 1) then
                                !    ! Add slight smoothing
                                !    uc = (1.0_dp - smooth_scale_coarse2fine * gradient_scale) * uc + &
                                !        smooth_scale_coarse2fine * gradient_scale * (0.25_dp * &
                                !        (U(i,jp1,k) + 0*U(im1, jp1, k) + 0*U(ip1, jp1, k) + &
                                !         U(i,jm1,k) + 0*U(im1, jm1, k) + 0*U(ip1, jm1, k) + &
                                !         U(ip1,j,k) + U(im1, j, k)))
                                !end if

                                ! NS gradients are also not quite right. Might not matter?
                                !if(im1 >= iL) then
                                    dU_dj_p = 0.5_dp * (U(i,jp1,k) - U(i,j  ,k)) + 0.5_dp * (U(im1,jp1,k) - U(im1,j  ,k))
                                    dU_dj_m = 0.5_dp * (U(i,j  ,k) - U(i,jm1,k)) + 0.5_dp * (U(im1,j  ,k) - U(im1,jm1,k))
                                !else
                                !    ! Avoid using a point outside of this region
                                !    ! Non-extrapolated derivative seems more stable
                                !    dU_dj_p = 1.0_dp * (U(i,jp1,k) - U(i,j  ,k)) !- 0.5_dp*(U(ip1,jp1,k) - U(ip1,j  ,k))
                                !    dU_dj_m = 1.0_dp * (U(i,j  ,k) - U(i,jm1,k)) !- 0.5_dp*(U(ip1,j  ,k) - U(ip1,jm1,k))
                                !end if

                            else if(k == VH) then
                                !
                                ! VH is offset north of where the receive grid wants it
                                !
                                if(neighbour_domain_staggered_grid == 1) then
                                    ! Slightly north of middle of the cell
                                    ! but not stable to do this
                                    del = 0.5_dp !- inv_cell_ratios(2)*0.5_dp
                                else
                                    ! Middle of the cell
                                    del = 0.5_dp
                                end if

                                ! Move VH back to where we want it
                                uc = uc - dU_dj_m*del
                                ! There is only one relevant gradient NS
                                dU_dj_p = dU_dj_m

                                !if(neighbour_domain_staggered_grid /= 1) then
                                !    ! Add slight smoothing
                                !    uc = (1.0_dp - smooth_scale_coarse2fine * gradient_scale) * uc + &
                                !        smooth_scale_coarse2fine * gradient_scale * (0.25_dp * &
                                !        (U(i,jp1,k) + 0*U(im1, jp1, k) + 0*U(ip1, jp1, k) + &
                                !         U(i,jm1,k) + 0*U(im1, jm1, k) + 0*U(ip1, jm1, k) + &
                                !         U(ip1,j,k) + U(im1, j, k)))
                                !end if

                                ! The EW gradients are not quite right. Might not matter?
                                ! Improve the approximation without using an edge jm1
                                !if(jm1 >= jL) then
                                    dU_di_p = 0.5_dp * (U(ip1,j,k) - U(i  ,j,k)) + 0.5_dp * (U(ip1,jm1,k) - U(i  ,jm1,k))
                                    dU_di_m = 0.5_dp * (U(i  ,j,k) - U(im1,j,k)) + 0.5_dp * (U(i  ,jm1,k) - U(im1,jm1,k))
                                !else
                                !    ! Avoid using a point outside of this region
                                !    ! Non-extrapolated derivative seems more stable
                                !    dU_di_p = 1.0_dp * (U(ip1,j,k) - U(i  ,j,k)) !- 0.5_dp*(U(ip1,jp1,k) - U(i  ,jp1,k))
                                !    dU_di_m = 1.0_dp * (U(i  ,j,k) - U(im1,j,k)) !- 0.5_dp*(U(i  ,jp1,k) - U(im1,jp1,k))
                                !end if
                            else
                                !
                                ! Stage and elevation
                                !
                                !if(neighbour_domain_staggered_grid /= 1) then
                                !    ! Add slight smoothing.
                                !    uc = (1.0_dp - smooth_scale_coarse2fine * gradient_scale) * uc + &
                                !        smooth_scale_coarse2fine * gradient_scale * 0.25_dp * &
                                !        (U(i,jp1,k) + 0*U(im1, jp1, k) + 0*U(ip1, jp1, k) + &
                                !         U(i,jm1,k) + 0*U(im1, jm1, k) + 0*U(ip1, jm1, k) + &
                                !         U(ip1,j,k) + U(im1,   j, k))
                                !end if

                            end if
                        end if

                        ! The definition of gradients above could cause issues at
                        ! wet-dry fronts, or shocks.
                        ! SUPPRESS GRADIENTS IF REQUIRED
                        dU_di_p = dU_di_p * gradient_scale_x
                        dU_dj_p = dU_dj_p * gradient_scale_y
                        dU_di_m = dU_di_m * gradient_scale_x
                        dU_dj_m = dU_dj_m * gradient_scale_y

                        !!if(my_domain_staggered_grid == 0 .and. neighbour_domain_staggered_grid == 1) then
                        !!    ! Sending from non-staggered to staggered
                        !!    if(k == UH) then
                        !!        !UH is offset slightly west of where the staggered grid wants it
                        !!        uc = uc + dU_di_p * 0.5_dp/cell_ratios(1)
                        !!    else if(k == VH) then
                        !!        ! VH is offset slightly south of where the staggered grid wants it
                        !!        uc = uc + dU_dj_p*0.5_dp/cell_ratios(2)
                        !!    end if
                        !!
                        !!end if

                        ! Index offset where we pack the send buffer
                        send_ind_offset = &
                            ! Number of cells for previous variables
                            (iU - iL + 1) * (jU - jL + 1) * (k-kL) * (1 + number_of_gradients_coarse_to_fine) + &
                            ! Number of cells for previous j, with this k
                            (iU - iL + 1) * (j - jL) * (1 + number_of_gradients_coarse_to_fine) + &
                            ! Cells for previous i, with this j and k
                            (i - iL) * (1 + number_of_gradients_coarse_to_fine) 

                        ! Send the variable value and a number of gradients
                        !
                        two_way_nesting_comms%send_buffer(send_ind_offset + 1) = uc
                        ! send (+-) x gradients
                        two_way_nesting_comms%send_buffer(send_ind_offset + 2) = dU_di_p
                        two_way_nesting_comms%send_buffer(send_ind_offset + 3) = dU_di_m
                        ! send (+-) y gradients
                        two_way_nesting_comms%send_buffer(send_ind_offset + 4) = dU_dj_p
                        two_way_nesting_comms%send_buffer(send_ind_offset + 5) = dU_dj_m

                    end do !k
                end do !i
            end do !j
            !$OMP END DO
        end if !my_domain_is_finer

        !$OMP END PARALLEL

        !
        ! Send boundary flux data as well
        !
        if(send_boundary_flux_data) then
   

            if(equal_cell_ratios) then
                !
                ! Optimized treatment of this special case
                !
                already_sent = iR * kR * jR

                ! Loop over send_box_flux_integral sides
                dir_ip = [NORTH, SOUTH, EAST, WEST]

                do nb = 1, size(dir_ip)
                    j = dir_ip(nb)
                    imx = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 1, kind=ip)
                    nvar = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 2, kind=ip)

                    do k = 1, nvar
                        two_way_nesting_comms%send_buffer((already_sent + 1):(already_sent+imx)) = &
                            two_way_nesting_comms%send_box_flux_integral(j)%x(1:imx, k)
                        already_sent = already_sent + imx
                    end do

                end do

            else
                ! Unequal cell sizes in send and receive regions

                if(two_way_nesting_comms%my_domain_is_finer) then
                    ! Fine-to-coarse
                    ! Pack each boundary flux into the send_buffer, by columns,
                    ! aggregating to the coarse recv grid 
                    already_sent = iR * kR * jR

                    ! Loop over send_box_flux_integral sides
                    dir_ip = [NORTH, SOUTH, EAST, WEST]
                    stride = [inv_cell_ratios_ip(1), inv_cell_ratios_ip(1), &
                              inv_cell_ratios_ip(2), inv_cell_ratios_ip(2)]

                    do nb = 1, size(dir_ip)

                        ! Boundary info
                        j = dir_ip(nb)
                        imx = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 1, kind=ip)
                        nvar = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 2, kind=ip)

                        do k = 1, nvar
                            ! Genuine fine to coarse case
                            do i = 1, imx
                                ! Aggregate cell values to match the coarser recv grid
                                two_way_nesting_comms%send_buffer(already_sent + 1) = &                
                                    two_way_nesting_comms%send_buffer(already_sent + 1) + &
                                    two_way_nesting_comms%send_box_flux_integral(j)%x(i, k)
                                if(mod(i, stride(nb)) == 0) already_sent = already_sent + 1
                            end do
                        end do
                    end do

                else
                    ! Coarse-to-fine
                    already_sent = (iU - iL + 1) * (jU - jL + 1) * (kU - kL + 1) * &
                        (1_ip + number_of_gradients_coarse_to_fine)

                    ! Pack each boundary flux into the send_buffer, by columns,
                    ! with redundancy to match the finer recv grid
                    ! Loop over send_box_flux_integral sides
                    dir_ip = [NORTH, SOUTH, EAST, WEST]
                    stride = [cell_ratios_ip(1), cell_ratios_ip(1), &
                              cell_ratios_ip(2), cell_ratios_ip(2)]
                    do nb = 1, size(dir_ip)

                        ! Boundary info
                        j = dir_ip(nb)
                        imx = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 1, kind=ip)
                        nvar = size(two_way_nesting_comms%send_box_flux_integral(j)%x, 2, kind=ip)
                        ! We need to 'spread' the flux integral over stride(nb) cells, hence this factor
                        inv_stride = 1.0_dp/stride(nb)
                        do k = 1, nvar
                            do i = 1, imx
                                ! Repeat values to match the coarser grid.
                                n0 = already_sent + 1
                                n1 = already_sent + stride(nb)
                                two_way_nesting_comms%send_buffer(n0:n1) = &
                                    two_way_nesting_comms%send_buffer(n0:n1) + &
                                    two_way_nesting_comms%send_box_flux_integral(j)%x(i, k)*inv_stride
                                already_sent = already_sent + stride(nb)
                            end do
                        end do ! k

                    end do ! nb

                end if ! my_domain_is_finer
            end if ! equal_cel_sizes
        end if !send_boundary_flux_data

    end subroutine

    
    ! Routine to copy the recv buffer into the main computational array, once
    ! it has been filled by the neighbour.
    ! 
    ! We keep this separate from the routine which copies the data, so that 
    ! in parallel we can do a sync between the two steps. This also gives us the flexibility
    ! to use the sync for other parallel comms processes.
    !
    ! @param two_way_nesting_comms 
    ! @param U the main computational array that we copy the received data to.
    ! 
    subroutine process_received_data(two_way_nesting_comms, U)
        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(inout) :: U(:,:,:)
  
        integer(ip) :: iL, iU, jL, jU, kL, kU, n0, n1, dir_ip(4), i, j, nvar, s1
        !real(dp) :: stage_mean_before, stage_mean_after
        real(dp) :: uc, dU_di_p, dU_di_m, dU_dj_p, dU_dj_m
        integer(ip), parameter :: coarse_factor = (1_ip + number_of_gradients_coarse_to_fine)
        integer(ip) :: iL_1, iU_1, jL_1, jU_1, i1, j1, k1
        integer(ip) :: inv_cell_ratios_ip(2), niC, njC, nC
        integer(ip) :: coarse_cell_i, coarse_cell_j, coarse_cell_k, coarse_cell_count, use_rw
        real(dp) :: middle_index(2), dU_di, dU_dj, dUj, dUi, rw, depth_low, depth_high
        integer(ip), parameter :: STG=1, UH=2, VH=3, ELV=4

        ! Quick exit if the type is purely used to send data
        if(.not. two_way_nesting_comms%recv_active) return

        ! Shorthand
        iL = two_way_nesting_comms%recv_inds(1,1)
        iU = two_way_nesting_comms%recv_inds(2,1)
        jL = two_way_nesting_comms%recv_inds(1,2)
        jU = two_way_nesting_comms%recv_inds(2,2)
        kL = two_way_nesting_comms%recv_inds(1,3)
        kU = two_way_nesting_comms%recv_inds(2,3)

        !stage_mean_before = sum(U(iL:iU, jL:jU, 1))/((iU-iL+1)*(jU-jL+1))

        call recv_from_p2p_comms(two_way_nesting_comms%recv_buffer(:), &
            buffer_label = two_way_nesting_comms%recv_ID)


        if(two_way_nesting_comms%my_domain_is_finer) then

            ! number of values that have been sent
            n1 = int(nint( (iU-iL+1)*(jU-jL+1)*(kU-kL+1) * &
                product(two_way_nesting_comms%cell_ratios) * coarse_factor))

            !write(log_output_unit, *) '@@DEBUG@@ ', n1, coarse_factor, iL, iU, jL, jU, kL, kU, &
            !    product(two_way_nesting_comms%cell_ratios)
            !flush(log_output_unit)

            ! Local storage
            two_way_nesting_comms%recv_work = U(iL:iU, jL:jU, STG) - U(iL:iU, jL:jU, ELV)

            !$OMP PARALLEL DEFAULT(PRIVATE), &
            !$OMP SHARED(two_way_nesting_comms, U, iL, iU, jL, jU, kL, kU, n1)

            !
            ! Unpack the coarse values and gradients, and interpolate onto the
            ! fine grid
            !

            ! Number of fine i/j cells per coarse cell
            inv_cell_ratios_ip = int(nint(1/two_way_nesting_comms%cell_ratios))

            ! 'index' of the fine cell in the middle of a coarse cell. 
            ! This is useful for interpolation, e.g.
            ! if inv_cell_ratios_ip = [3,3], then middle_index = [2,2],
            ! if inv_cell_ratios_ip = [4,4], then middle_index = [2.5,2.5]
            middle_index = 0.5_dp * (inv_cell_ratios_ip + 1_ip)

            ! Number of i/j cells in the coarser domain
            niC = (iU - iL + 1)/inv_cell_ratios_ip(1)
            njC = (jU - jL + 1)/inv_cell_ratios_ip(2) 
            nC = niC * njC

            !$OMP DO SCHEDULE(STATIC)
            do i = 1, (n1 - coarse_factor + 1), coarse_factor

                ! Unpack the coarse cell value and its gradients. The order
                ! should correspond to how the data was packed in 'process_data_to_send'

                uc = two_way_nesting_comms%recv_buffer(i)
                dU_di_p = two_way_nesting_comms%recv_buffer(i+1)
                dU_di_m = two_way_nesting_comms%recv_buffer(i+2)
                dU_dj_p = two_way_nesting_comms%recv_buffer(i+3)
                dU_dj_m = two_way_nesting_comms%recv_buffer(i+4)

                ! Compute C-style multidimensional indices corresponding to the
                ! coarser domain we are receiving from -- index starting from 0
                coarse_cell_count = (i-1)/coarse_factor
                coarse_cell_k = coarse_cell_count/nC 
                coarse_cell_j = (coarse_cell_count - (coarse_cell_k*nC))/niC
                coarse_cell_i = (coarse_cell_count - (coarse_cell_k*nC + coarse_cell_j*niC))

                ! Pack into the appropriate part of U. For this coarse cell, the
                ! fine-U indices range from (iL_1 ... iU_1), and (jL_1 ... jU_1)
                iL_1 = iL + coarse_cell_i * inv_cell_ratios_ip(1)
                iU_1 = iL_1 + inv_cell_ratios_ip(1) - 1
                jL_1 = jL + coarse_cell_j * inv_cell_ratios_ip(2)
                jU_1 = jL_1 + inv_cell_ratios_ip(2) - 1
                k1 = kL + coarse_cell_k

                ! Make sure we do not use the 'receive-weights' near wet-dry areas
                ! For mass conservation, we need to force exact matches.
                ! FIXME: Make these thresholds less ad-hoc
                depth_low  = minval(two_way_nesting_comms%recv_work(iL_1:iU_1, jL_1:jU_1))
                depth_high = maxval(two_way_nesting_comms%recv_work(iL_1:iU_1, jL_1:jU_1))
                if(depth_low <= ignore_receive_weights_depth_shallower_than .or. &
                    (depth_high > ignore_receive_weights_depth_ratio * depth_low)) then
                    ! We might be near a wet-dry boundary -- do not use recv-weights
                    use_rw = 0
                else
                    use_rw = 1
                end if

                ! Loop over all cells in the 'recv' grid that can be reconstructed with
                ! the data we just unpacked
                do j1 = jL_1, jU_1

                    ! Get the j gradient
                    if(j1 - jL_1 + 1 < middle_index(2)) then
                        dU_dj = dU_dj_m
                    else
                        dU_dj = dU_dj_p
                    end if

                    ! The 'j interpolation adjustment' of the variable
                    dUj = ((j1-jL_1+1) - middle_index(2)) * two_way_nesting_comms%cell_ratios(2) *  dU_dj

                    do i1 = iL_1, iU_1
                        ! Get the i gradient
                        if(i1 - iL_1 + 1 < middle_index(1)) then
                            dU_di = dU_di_m
                        else
                            dU_di = dU_di_p
                        end if

                        ! The i interpolation adjustment of the variable
                        dUi = ((i1-iL_1+1) - middle_index(1)) * two_way_nesting_comms%cell_ratios(1) *  dU_di

                        ! Recv weights (typically this is 1.0, but values in [0-1] give the option of
                        ! more complex updating
                        rw = 1.0_dp - (1.0_dp - two_way_nesting_comms%recv_weights(i1, j1))*use_rw
                        ! Interpolate
                        U(i1, j1, k1) = (uc + dUj + dUi)*rw + (1.0_dp - rw)*U(i1, j1, k1)
                    end do
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        else
            ! Simple case -- the processing was done on the recv side
            n1 = (iU-iL+1)*(jU-jL+1)*(kU-kL+1)

            call repack_rank1_array(two_way_nesting_comms%recv_buffer(1:n1), &
                U(iL:iU, jL:jU, kL:kU))
        end if

        if(send_boundary_flux_data) then

            ! Loop over boundaries, and compute the recv_box_flux_error
            dir_ip = [NORTH, SOUTH, EAST, WEST]

            do j = 1, size(dir_ip)
                nvar = size(two_way_nesting_comms%recv_box_flux_integral(dir_ip(j))%x, 2, kind=ip)
                s1 = size(two_way_nesting_comms%recv_box_flux_integral(dir_ip(j))%x, 1, kind=ip)
                do i = 1, nvar
                    ! Indices into recv_buffer
                    n0 = n1 + 1
                    n1 = n1 + s1
                    ! Get the difference with the recv'd values -- since this
                    ! gives the required flux correction
                    two_way_nesting_comms%recv_box_flux_error(dir_ip(j))%x(:,i) = &
                        real(two_way_nesting_comms%recv_box_flux_integral(dir_ip(j))%x(:,i), force_double) - &
                        real(two_way_nesting_comms%recv_buffer(n0:n1), force_double)
                end do
            end do
        end if

        !stage_mean_after = sum(U(iL:iU, jL:jU, 1))/((iU-iL+1)*(jU-jL+1))
        !
        !write(log_output_unit,*) ' .. ', stage_mean_before - stage_mean_after

    end subroutine

    ! Copy send_buffer to the right recv_buffer. This should be called AFTER process_data_to_send.
    !
    ! @param send_to_recv_buffer optional logical. If FALSE, then do not do a
    !   coarray put, but only pack the data in the send buffer [we can communicate later with communicate_p2p].
    subroutine send_data(two_way_nesting_comms, send_to_recv_buffer)
        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        logical, optional, intent(in) :: send_to_recv_buffer

        logical:: put_in_recv_buffer

        ! Quick exit if the type is purely used to recv data
        if(.not. two_way_nesting_comms%send_active) return

        if(present(send_to_recv_buffer)) then
            put_in_recv_buffer = send_to_recv_buffer
        else
            put_in_recv_buffer = .TRUE.
        end if
  
        ! Send to point2point coarray communication buffer 
        call send_to_p2p_comms(two_way_nesting_comms%send_buffer, &
            buffer_label=two_way_nesting_comms%send_ID, &
            put_in_recv_buffer=put_in_recv_buffer)
     
    end subroutine

    ! Multiply the boundary flux integral terms by some constant
    !
    pure subroutine boundary_flux_integral_multiply(two_way_nesting_comms, c)

        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(in) :: c

        integer(ip) :: i

        if(two_way_nesting_comms%send_active) then
            ! Do NORTH, SOUTH, EAST, WEST
            do i = 1, 4
                two_way_nesting_comms%send_box_flux_integral(i)%x = c * &
                    two_way_nesting_comms%send_box_flux_integral(i)%x
            end do
        end if

        if(two_way_nesting_comms%recv_active) then
            ! Do NORTH, SOUTH, EAST, WEST
            do i = 1, 4
                two_way_nesting_comms%recv_box_flux_integral(i)%x = c * &
                    two_way_nesting_comms%recv_box_flux_integral(i)%x
            end do

        end if

    end subroutine

    ! Replace each boundary_flux_integral_term with 
    !  (dt * current_fluxes * dx + boundary_flux_integral_term)
    !
    ! @param two_way_nesting_comms Communication object, which is also tracking
    !     fluxes through its boundaries.
    ! @param dt real constant in the equation above
    ! @param flux_NS rank 3 array with north-south fluxes
    ! @param flux_NS_lower_index integer. Assume flux_NS(:,1,:) contains the
    !   bottom edge flux for cells with j index = flux_NS_lower_index. For example,
    !   "flux_NS_lower_index=1" implies that flux_NS includes fluxes on the south boundary.
    !   For some of our solvers this is not true [e.g. linear leapfrog, because the 'mass flux'
    !   terms are effectively stored in the domain%U variable].
    ! @param distance_bottom_edge real rank 1 array with length = (1+number of 'j' cells in domain).
    !    Gives the length of the bottom edge of a cell with 'j' y index. This
    !    is an array because for spherical coordinates, the bottom edge cell 
    !    distance changes with latitude.
    ! @param flux_NS rank 3 array with east-west fluxes
    ! @param flux_NS_lower_index integer. Assume flux_EW(1,:,:) contains the left-edge
    !    flux for cells with i index = flux_EW_lower_index. For example,
    !   "flux_EW_lower_index=1" implies that flux_EW includes fluxes right to the boundary.
    !   For some of our solvers this is not true [e.g. linear leapfrog, because the 'mass flux'
    !   terms are effectively stored in the domain%U variable].
    ! @param distance_left_edge real rank 1 array with length = (1+number of 'i' cells in domain).
    !    Gives the length of the left edge of a cell with 'i' y index. This is usually constant.
    ! @param k_indices Only apply the update to a few variables with 'k_indices'.
    ! @param flux_already_multiplied_by_dx logical If TRUE, assume flux_NS and flux_EW have
    !   already been multiplied by their 'distance_bottom_edge' or 'distance_left_edge' term,
    !   so do not multiply by this again. Useful because some of our solvers store the fluxes
    !   as (flux . dx), while others do not.
    ! 
    pure subroutine boundary_flux_integral_tstep(two_way_nesting_comms, dt, &
        flux_NS, flux_NS_lower_index, distance_bottom_edge, &
        flux_EW, flux_EW_lower_index, distance_left_edge, &
        k_indices, flux_already_multiplied_by_dx)

        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(in) :: flux_NS(:,:,:), flux_EW(:,:,:), dt
        real(dp), intent(in) :: distance_bottom_edge(:), distance_left_edge(:)
        integer(ip), intent(in) :: flux_NS_lower_index, flux_EW_lower_index, k_indices(2)
        logical, intent(in) :: flux_already_multiplied_by_dx

        integer(ip) :: ks(2), ms(2), ns, dns, i, edge_ip(2), ind1_ip(2), ind2_ip(2), edge

        ! Range of variable indices we update is based on the 3rd dimension of flux_NS.
        ! This allows a few cases:
        !     - If flux_NS contains fluxes for [stg, uh, vh], then all those
        !       will be updated in the boundary flux integral. 
        !     - If flux_NS has only fluxes for [stg], then only that is updated.
        
        ks = k_indices

        if (two_way_nesting_comms%send_active) then

            !
            ! NORTH - SOUTH SIDE
            !

            edge_ip = [NORTH, SOUTH] ! Index into send_box_flux_integral(.)
            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%send_inds(.,2)
            ind2_ip = [1, 0] ! Offset required to get 'exterior edge flux' -- +1 for NORTH, +0 for SOUTH
            do i = 1, 2
                ! i indices in flux_NS
                ms(1:2) = two_way_nesting_comms%send_inds(:,1)
                ! j indices in flux_NS = ns-dns
                ns = two_way_nesting_comms%send_inds(ind1_ip(i),2) + ind2_ip(i)
                dns = (flux_NS_lower_index - 1_ip)

                edge = edge_ip(i)
                if(ns - dns >= 1_ip .and. ns - dns <= size(flux_NS,2, kind=ip)) then
                    two_way_nesting_comms%send_box_flux_integral(edge)%x(:,ks(1):ks(2)) = &
                        two_way_nesting_comms%send_box_flux_integral(edge)%x(:,ks(1):ks(2)) + &
                        dt * flux_NS(ms(1):ms(2), ns - dns, ks(1):ks(2)) * &
                        merge(1.0_dp, distance_bottom_edge(ns), flux_already_multiplied_by_dx)
                end if
            end do


            !
            ! EAST - WEST side
            !
            edge_ip = [EAST, WEST] ! Index into send_box_flux_integral(.)
            ind1_ip = [2, 1] ! Rank1 index into two_way_nesting_comms%send_inds(., 1)
            ind2_ip = [1, 0] ! Offset required to get 'exterior edge flux' -- +1 for EAST, +0 for WEST
            do i = 1, 2 
                ! i indices in flux_EW = ns-dns
                ns = two_way_nesting_comms%send_inds(ind1_ip(i),1) + ind2_ip(i)
                dns = (flux_EW_lower_index - 1_ip)
                ! j indices in flux_EW
                ms(1:2) = two_way_nesting_comms%send_inds(:,2)

                edge = edge_ip(i)
                if(ns - dns >= 1_ip .and. ns - dns <= size(flux_EW,1, kind=ip)) then
                    two_way_nesting_comms%send_box_flux_integral(edge)%x(:,ks(1):ks(2)) = &
                        two_way_nesting_comms%send_box_flux_integral(edge)%x(:,ks(1):ks(2)) + &
                        dt * flux_EW(ns-dns, ms(1):ms(2), ks(1):ks(2)) * &
                        merge(1.0_dp, distance_left_edge(ns), flux_already_multiplied_by_dx)
                end if
            end do

        end if

        !
        ! The  entire if statement below is identical to the one above, but
        ! with 'send' replaced with 'recv'
        !
        if (two_way_nesting_comms%recv_active) then

            !
            ! NORTH - SOUTH SIDE
            !
            edge_ip = [NORTH, SOUTH] ! Index into recv_box_flux_integral(.)
            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%recv_inds(.,2)
            ind2_ip = [1, 0] ! Offset required to get 'exterior edge flux' -- +1 for NORTH, +0 for SOUTH
            do i = 1, 2
                ! i indices in flux_NS
                ms(1:2) = two_way_nesting_comms%recv_inds(:,1)
                ! j indices in flux_NS = ns-dns
                ns = two_way_nesting_comms%recv_inds(ind1_ip(i),2) + ind2_ip(i)
                dns = (flux_NS_lower_index - 1_ip)

                edge = edge_ip(i)
                if(ns - dns >= 1_ip .and. ns - dns <= size(flux_NS,2, kind=ip)) then
                    two_way_nesting_comms%recv_box_flux_integral(edge)%x(:,ks(1):ks(2)) = &
                        two_way_nesting_comms%recv_box_flux_integral(edge)%x(:,ks(1):ks(2)) + &
                        dt * flux_NS(ms(1):ms(2), ns - dns, ks(1):ks(2)) * &
                        merge(1.0_dp, distance_bottom_edge(ns), flux_already_multiplied_by_dx)
                end if
            end do


            !
            ! EAST - WEST side
            !
            edge_ip = [EAST, WEST] ! Index into recv_box_flux_integral(.)
            ind1_ip = [2, 1] ! Rank1 index into two_way_nesting_comms%recv_inds(., 1)
            ind2_ip = [1, 0] ! Offset required to get 'exterior edge flux' -- +1 for EAST, +0 for WEST
            do i = 1, 2 
                ! i indices in flux_EW = ns-dns
                ns = two_way_nesting_comms%recv_inds(ind1_ip(i),1) + ind2_ip(i)
                dns = (flux_EW_lower_index - 1_ip)
                ! j indices in flux_EW
                ms(1:2) = two_way_nesting_comms%recv_inds(:,2)

                edge = edge_ip(i)
                if(ns - dns >= 1_ip .and. ns - dns <= size(flux_EW,1, kind=ip)) then
                    two_way_nesting_comms%recv_box_flux_integral(edge)%x(:,ks(1):ks(2)) = &
                        two_way_nesting_comms%recv_box_flux_integral(edge)%x(:,ks(1):ks(2)) + &
                        dt * flux_EW(ns-dns, ms(1):ms(2), ks(1):ks(2)) * &
                        merge(1.0_dp, distance_left_edge(ns), flux_already_multiplied_by_dx)
                end if
            end do

        end if


    end subroutine

    !! Report on fluxes between nesting comms type and target domain.
    !! FIXME: This was used in development but is likely out of date now.
    !!
    !! @param nesting nesting type
    !! @param xs optional vector of x coordinates, in domain%x
    !! @param ys optional vector of y coordinates, in domain%y
    !!
    !subroutine print_nesting_fluxes(nesting, xs, ys)

    !    class(domain_nesting_type), intent(in) :: nesting
    !    real(dp), intent(in), optional :: xs(:), ys(:)

    !    integer(ip) :: i, j, nbr_index, nbr_image, nbr_comms, n, my_index, my_image, my_comms, i1(2)
    !    real(dp) :: xr(2), yr(2)
    !    real(dp) :: flux_send(4), flux_recv(4), flux_err_coarse_recv(4)
    !    integer(ip) :: edge_ip(2), ind1_ip(2), ind2_ip(2), edge

    !    if(allocated(nesting%send_comms)) then

    !        ! Compute the mass fluxes that are going to the 'real' part of
    !        ! another domain (i.e. where that domain is the priority domain)

    !        do i = 1, size(nesting%send_comms, kind=ip)

    !            nbr_index = nesting%send_comms(i)%neighbour_domain_index
    !            nbr_image = nesting%send_comms(i)%neighbour_domain_image_index
    !            nbr_comms = nesting%send_comms(i)%neighbour_domain_comms_index
    !            my_index = nesting%send_comms(i)%my_domain_index
    !            my_image = nesting%send_comms(i)%my_domain_image_index
    !            my_comms = nesting%send_comms(i)%my_domain_comms_index
    !       
    !            !
    !            ! NORTH-SOUTH SIDE
    !            ! 
    !            edge_ip = [NORTH, SOUTH] ! Index into send_box_flux_integral(.)
    !            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%send_inds(.,2)
    !            ind2_ip = [1, -1] ! Offset required to get 'j index of cell just outside the bbox ' -- +1 for NORTH, -1 for SOUTH
    !            do j = 1, 2 

    !                n = nesting%send_comms(i)%send_inds(ind1_ip(j),2) + ind2_ip(j) ! j-index just outside the send bbox
    !                i1 = nesting%send_comms(i)%send_inds(:,1)

    !                edge = edge_ip(j)
    !                if(n <= size(nesting%priority_domain_index,2, kind=ip) .and. n >= 1) then
    !                    ! Find the mass flux northward to/from the neighbour domain (if any) 
    !                    flux_send(edge) = sum(&
    !                        nesting%send_comms(i)%send_box_flux_integral(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(i1(1):i1(2),n) == nbr_index) .and. &
    !                               (nesting%priority_domain_image(i1(1):i1(2),n) == nbr_image) ) )
    !                else
    !                    flux_send(edge) = 0.0_dp
    !                end if

    !            end do

    !            !
    !            ! Store x range of this send bounding box
    !            !
    !            if(present(xs)) then
    !                xr = xs(i1) + 0.5_dp * [xs(i1(1)) - xs(i1(1)+1), xs(i1(2))  - xs(i1(2)-1)]
    !            else
    !                xr = 0.0_dp
    !            end if

    !            !
    !            ! EAST-WEST SIDE
    !            !
    !            edge_ip = [EAST, WEST] ! Index into send_box_flux_integral(.)
    !            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%send_inds(.,2)
    !            ind2_ip = [1, -1] ! Offset required to get 'i index of cell just outside the bbox ' -- +1 for EAST, -1 for WEST 

    !            do j = 1, 2

    !                ! Find the mass flux going eastward to/from the neighbour domain (if any)
    !                n = nesting%send_comms(i)%send_inds(ind1_ip(j),1) + ind2_ip(j) ! i-index just outside the send bbox
    !                i1 = nesting%send_comms(i)%send_inds(:,2)

    !                edge = edge_ip(j)

    !                if(n <= size(nesting%priority_domain_index, 1, kind=ip) .and. n >= 1) then
    !                    flux_send(edge) = sum(&
    !                        nesting%send_comms(i)%send_box_flux_integral(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(n,i1(1):i1(2)) == nbr_index) .and. &
    !                               (nesting%priority_domain_image(n,i1(1):i1(2)) == nbr_image) ) )
    !                else
    !                    flux_send(edge) = 0.0_dp
    !                end if
    !            end do

    !            !
    !            ! Store y range of this send bounding box
    !            !
    !            if(present(ys)) then
    !                yr = ys(i1) + 0.5_dp * [ys(i1(1)) - ys(i1(1)+1), ys(i1(2)) - ys(i1(2) - 1)]
    !            else
    !                yr = 0.0_dp
    !            end if

    !            write(log_output_unit,*) 'send flux from ', &
    !                my_index, my_image, my_comms, ' to ', &
    !                nbr_index, nbr_image, nbr_comms, &
    !                ' N: ', flux_send(NORTH), ' S: ', flux_send(SOUTH), &
    !                ' E: ', flux_send(EAST), ' W: ', flux_send(WEST), &
    !                'bbox: ', xr, yr

    !        end do

    !    end if

    !    if(allocated(nesting%recv_comms)) then

    !        ! Compute the mass fluxes that are going to/from the 'real' part of
    !        ! my domain (i.e. where my domain is the priority domain)

    !        do i = 1, size(nesting%recv_comms, kind=ip)

    !            nbr_index = nesting%recv_comms(i)%neighbour_domain_index
    !            nbr_image = nesting%recv_comms(i)%neighbour_domain_image_index
    !            nbr_comms = nesting%recv_comms(i)%neighbour_domain_comms_index
    !            my_index = nesting%recv_comms(i)%my_domain_index
    !            my_image = nesting%recv_comms(i)%my_domain_image_index
    !            my_comms = nesting%recv_comms(i)%my_domain_comms_index


    !            !
    !            ! NORTH-SOUTH SIDE
    !            !
    !            edge_ip = [NORTH, SOUTH] ! Index into send_box_flux_integral(.)
    !            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%send_inds(.,2)
    !            ind2_ip = [1, -1] ! Offset required to get 'j index of cell just outside the bbox ' -- +1 for NORTH, -1 for SOUTH
    !            do j = 1, 2
    !                ! Find fluxes flowing north (or south) to/from the 'real' part of my domain
    !                n = nesting%recv_comms(i)%recv_inds(ind1_ip(j),2) + ind2_ip(j) ! j-index just outside the bbox
    !                i1 = nesting%recv_comms(i)%recv_inds(:,1)

    !                edge = edge_ip(j)
    !                if(n <= size(nesting%priority_domain_index,2, kind=ip) .and. n>=1) then
    !                    flux_recv(edge) = sum(&
    !                        nesting%recv_comms(i)%recv_box_flux_integral(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(i1(1):i1(2),n) == my_index) .and. &
    !                               (nesting%priority_domain_image(i1(1):i1(2),n) == my_image) ) )
    !                    ! On coarse domains which recv, record the flux error
    !                    flux_err_coarse_recv(edge) = &
    !                        count([send_boundary_flux_data])*&
    !                        count([(.not. nesting%recv_comms(i)%my_domain_is_finer)]) * &
    !                        sum(nesting%recv_comms(i)%recv_box_flux_error(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(i1(1):i1(2),n) == my_index) .and. &
    !                               (nesting%priority_domain_image(i1(1):i1(2),n) == my_image) ) )

    !                else
    !                    flux_recv(edge) = 0.0_dp
    !                    flux_err_coarse_recv(edge) = 0.0_dp
    !                end if

    !            end do

    !            !
    !            ! Store x range of this recv bounding box
    !            !
    !            if(present(xs)) then
    !                xr = xs(i1) + 0.5_dp * [xs(i1(1)) - xs(i1(1)+1), xs(i1(2))  - xs(i1(2)-1)]
    !            else
    !                xr = 0.0_dp
    !            end if

    !            !
    !            ! EAST-WEST SIDE
    !            !
    !            edge_ip = [EAST, WEST] ! Index into recv_box_flux_integral(.)
    !            ind1_ip = [2, 1] ! Rank1 index in two_way_nesting_comms%recv_inds(.,2)
    !            ind2_ip = [1, -1] ! Offset required to get 'i index of cell just outside the bbox ' -- +1 for EAST, -1 for WEST 

    !            do j = 1, 2
    !                ! Find the mass flux eastward/westward to/from the priority part of my domain (if any) 
    !                n = nesting%recv_comms(i)%recv_inds(ind1_ip(j),1) + ind2_ip(j) ! i-index just outside the recv bbox
    !                i1 = nesting%recv_comms(i)%recv_inds(:,2)

    !                edge = edge_ip(j)
    !                if(n <= size(nesting%priority_domain_index, 1, kind=ip) .and. n>=1) then
    !                    flux_recv(edge) = sum(&
    !                        nesting%recv_comms(i)%recv_box_flux_integral(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(n,i1(1):i1(2)) == my_index) .and. &
    !                               (nesting%priority_domain_image(n,i1(1):i1(2)) == my_image) ) )

    !                    ! On coarse domains which recv, record the flux error
    !                    flux_err_coarse_recv(edge) = &
    !                        count([send_boundary_flux_data])*&
    !                        count([(.not. nesting%recv_comms(i)%my_domain_is_finer)]) * &
    !                        sum(nesting%recv_comms(i)%recv_box_flux_error(edge)%x(:,1), &
    !                        mask=( (nesting%priority_domain_index(n,i1(1):i1(2)) == my_index) .and. &
    !                               (nesting%priority_domain_image(n,i1(1):i1(2)) == my_image) ) )
    !                else
    !                    flux_recv(edge) = 0.0_dp 
    !                    flux_err_coarse_recv(edge) = 0.0_dp
    !                end if

    !            end do

    !            ! Store the y-range of the recv bbox
    !            if(present(ys)) then
    !                yr = ys(i1) + 0.5_dp * [ys(i1(1)) - ys(i1(1)+1), ys(i1(2)) - ys(i1(2) - 1)]
    !            else
    !                yr = 0.0_dp
    !            end if

    !            write(log_output_unit,*) 'recv flux from ', &
    !                nbr_index, nbr_image, nbr_comms, ' to ', &
    !                my_index, my_image, my_comms, &
    !                ' N: ', flux_recv(NORTH), '(', flux_err_coarse_recv(NORTH), ')', &
    !                ' S: ', flux_recv(SOUTH), '(', flux_err_coarse_recv(SOUTH), ')', &
    !                ' E: ', flux_recv(EAST),  '(', flux_err_coarse_recv(EAST), ')', &
    !                ' W: ', flux_recv(WEST),  '(', flux_err_coarse_recv(WEST), ')', &
    !                'bbox: ', xr, yr
    !        end do

    !    end if

    !end subroutine

    !
    ! Workhorse unit-test routine.
    ! Allows different levels of nesting
    ! @param nest_ratio integer. The child grid dx is '1/nest_ratio' times the parent grid dx
    subroutine test_nested_grid_comms_mod_workhorse(nest_ratio)

        ! Relative grid sizes
        integer(ip), INTENT(IN) :: nest_ratio != 3

        ! Where communicating, each domain needs to receive a layer of cells at
        ! least this thick along its boundary. 
        integer(ip) :: min_exchange_halo = 2

        ! Number of variables in the problem
        integer(ip), parameter :: nvar = 4

        ! Domain 1 -- outer domain
        integer(ip), parameter :: nx1(2) = [20_ip, 15_ip]
        real(dp), parameter :: ll1(2) = [ZERO_dp, ZERO_dp]
        real(dp), parameter :: dx1(2) = [ONE_dp, ONE_dp]

        ! Domain 2 -- inner nested domain -- cells of size dx1/nest_ratio
        ! Starts at [3, 4], finishes at [12, 12]
        real(dp), parameter :: ll2(2) = [3.0_dp, 4.0_dp]
        integer(ip), parameter :: nx2_on_nest_ratio(2) = [9_ip, 8_ip]
        integer(ip) :: nx2(2), b4_flux_data_p, b4_flux_data_c
        real(dp) :: dx2(2)
        ! For testing
        real(dp), parameter :: max_parent_send_val = 32.0_dp, min_parent_send_val = 8.0_dp
        real(dp) :: max_parent_send_val_offset

        integer(ip), parameter :: parent_comms_index = 3, child_comms_index = 2

        ! Type which serves the role of the array holding the domains, and something
        ! to copy the above into, to ensure no changes due to communication
        type array_of_arrays_r3_dp_type
            real(dp), allocatable :: U(:,:,:)
            type(two_way_nesting_comms_type), allocatable :: two_way_nesting_comms(:)
        end type

        type(array_of_arrays_r3_dp_type) :: Us(2), Us_store(2)

        ! Local variables
        integer(ip) :: i, j, k, xi, yi, this_image_local, num_images_local !, ierr
        integer(ip) :: parent_image_index, child_image_index
        real(dp) :: x, y, z, err_tol, test_max, test_min
        integer(ip), parameter :: strd = 1_ip + number_of_gradients_coarse_to_fine

#ifdef COARRAY
        this_image_local = this_image2()
        num_images_local = num_images2()
        !! Approach using MPI
        !call mpi_comm_rank(MPI_COMM_WORLD, this_image_local, ierr)
        !this_image_local = this_image_local + 1
        !call mpi_comm_size(MPI_COMM_WORLD, num_images_local, ierr)
#else
        this_image_local = 1
        num_images_local = 1
#endif    

        ! Do communication in parallel case
        child_image_index = this_image_local - 1
        parent_image_index = this_image_local + 1
        if(parent_image_index > num_images_local) parent_image_index = 1
        if(child_image_index < 1) child_image_index = num_images_local

        ! Number of cells in child domain 
        nx2 = nx2_on_nest_ratio * nest_ratio
        ! dx in child domain
        dx2 = dx1/nest_ratio

        allocate(Us(1)%U(nx1(1), nx1(2), nvar))
        allocate(Us_store(1)%U(nx1(1), nx1(2), nvar))

        allocate(Us(2)%U(nx2(1), nx2(2), nvar))
        allocate(Us_store(2)%U(nx2(1), nx2(2), nvar))

        ! Each 'domain' will probably have multiple nesting communicators. Make
        ! 4 here, although there could be more or less.
        allocate(Us(1)%two_way_nesting_comms(4), Us(2)%two_way_nesting_comms(4))

        ! Make up some data for the outer domain
        do k = 1, nvar
            do j = 1, nx1(2)
                do i = 1, nx1(1)
                    x = (ll1(1) + (i-0.5_dp)*dx1(1))
                    y = (ll1(2) + (j-0.5_dp)*dx1(2))
                    z = k*k
                    Us(1)%U(i,j,k) = x + y + z
                end do
            end do
        end do

        ! Make up some data for the inner domain, which 
        ! should agree with the outer domain (with linear interpolation)
        do k = 1, nvar
            do j = 1, nx2(2)
                do i = 1, nx2(1)
                    x = (ll2(1) + (i-0.5_dp)*dx2(1))
                    y = (ll2(2) + (j-0.5_dp)*dx2(2))
                    z = k*k
                    Us(2)%U(i,j,k) = x + y + z
                end do
            end do
        end do

        ! Copy for comparison later
        Us_store(1)%U = Us(1)%U
        Us_store(2)%U = Us(2)%U

        ! Left side comms -- child to parent
        call Us(2)%two_way_nesting_comms(child_comms_index)%initialise( &
            my_dx = dx2, &
            neighbour_dx = dx1, &
            ! Star halo to send
            ijk_to_send = [ nest_ratio + [1, nest_ratio*min_exchange_halo] , &
                [nest_ratio+1, nx2(2) - nest_ratio], [1, nvar]], &
            ! Receive all cells on the outer edge of the child domain
            ijk_to_recv = [[1, nest_ratio], [1, nx2(2)], [1, nvar]], &
            neighbour_domain_index = 1, &
            my_domain_index = 2, &
            neighbour_domain_comms_index=parent_comms_index, &
            my_domain_comms_index=child_comms_index, &
            neighbour_domain_image_index = parent_image_index, &
            my_domain_image_index = this_image_local,&
            use_wetdry_limiting = .false.)

        ! x/y index of parent domain cell inside lower left corner of child
        ! domain
        xi = nint((ll2(1) - ll1(1))/dx1(1)) + 1 
        yi = nint((ll2(2) - ll1(2))/dx1(2)) + 1 

        ! Left side comms -- parent to child
        call Us(1)%two_way_nesting_comms(parent_comms_index)%initialise( &
            my_dx = dx1, &
            neighbour_dx = dx2, &
            ! Send all cells on the outer edge of the child domain
            ijk_to_send = [ xi + [0,0], yi + [0, nx2(2)/nest_ratio - 1], [1, nvar]], &
            ! Star halo recv
            ijk_to_recv = [ xi + [1,min_exchange_halo], &
                yi + [1, nx2(2)/nest_ratio - 1 - 1], [1, nvar]], &
            neighbour_domain_index = 2, &
            my_domain_index = 1, &
            neighbour_domain_comms_index=child_comms_index, &
            my_domain_comms_index=parent_comms_index, &
            neighbour_domain_image_index = child_image_index, &
            my_domain_image_index = this_image_local,&
            use_wetdry_limiting = .false.)

        call allocate_p2p_comms

        !
        ! Basic tests of cell ratios
        !
        if(any(abs(Us(2)%two_way_nesting_comms(child_comms_index)%cell_ratios - 1.0_dp/nest_ratio) > EPS)) then
            write(log_output_unit,*) 'FAIL: Problems with cell ratios', __LINE__,&
                 __FILE__
            call generic_stop()
        else
            write(log_output_unit,*) 'PASS'
        end if

        if(any(abs(Us(1)%two_way_nesting_comms(parent_comms_index)%cell_ratios - nest_ratio) > EPS)) then
            write(log_output_unit,*) 'FAIL: Problems with cell ratios', __LINE__,&
                __FILE__
            call generic_stop()
        else
            write(log_output_unit,*) 'PASS'
        end if

        !
        ! Send/recv buffers have corresponding sizes
        !
        if( size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer, kind=ip) /= &
            size(Us(2)%two_way_nesting_comms(child_comms_index)%recv_buffer, kind=ip)) then

            write(log_output_unit,*) 'FAIL: Problems with send/recv buffer sizes', __LINE__,&
                __FILE__
            write(log_output_unit,*) size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer, kind=ip)
            write(log_output_unit,*) size(Us(2)%two_way_nesting_comms(child_comms_index)%recv_buffer, kind=ip)
            write(log_output_unit,*) Us(1)%two_way_nesting_comms(parent_comms_index)%nsend
            write(log_output_unit,*) Us(1)%two_way_nesting_comms(parent_comms_index)%nsend_interior
            write(log_output_unit,*) Us(2)%two_way_nesting_comms(child_comms_index)%nrecv
            write(log_output_unit,*) Us(2)%two_way_nesting_comms(child_comms_index)%nrecv_interior
            call generic_stop()
        else
            write(log_output_unit,*) 'PASS'
        end if

        if(size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer, kind=ip) /= &
            size(Us(1)%two_way_nesting_comms(parent_comms_index)%recv_buffer, kind=ip)) then

            write(log_output_unit,*) 'FAIL: Problems with send/recv buffer sizes', __LINE__,&
                __FILE__
            write(log_output_unit,*) size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer, kind=ip)
            write(log_output_unit,*) size(Us(1)%two_way_nesting_comms(parent_comms_index)%recv_buffer, kind=ip)
            call generic_stop()
        else
            write(log_output_unit,*) 'PASS'
        end if

        !
        ! my_domain_is_finer is correctly set
        !
        if(nest_ratio > 1) then
            if(Us(1)%two_way_nesting_comms(parent_comms_index)%my_domain_is_finer) then
                write(log_output_unit,*) 'FAIL: Problem with fine/coarse classification', __LINE__, &
                    __FILE__
                call generic_stop()
            else
                write(log_output_unit,*) 'PASS'
            end if        

            if(.NOT.Us(2)%two_way_nesting_comms(child_comms_index)%my_domain_is_finer) then
                write(log_output_unit,*) 'FAIL: Problem with fine/coarse classification', __LINE__, &
                    __FILE__
                call generic_stop()
            else
                write(log_output_unit,*) 'PASS'
            end if
        end if

        ! Fill the send buffer
        call Us(1)%two_way_nesting_comms(parent_comms_index)%process_data_to_send(Us(1)%U)
        call Us(2)%two_way_nesting_comms(child_comms_index)%process_data_to_send(Us(2)%U)

        ! Check for errors in send buffer
        if(dp == force_double) then
            err_tol = 1.0e-14_dp
        else
            err_tol = 1.0e-5_dp
        end if

        ! The tests below were written before we increased the send buffer to store flux data
        ! Fix that by only checking the data 'before the flux data'
        b4_flux_data_p = Us(1)%two_way_nesting_comms(parent_comms_index)%nsend_interior
        b4_flux_data_c = Us(2)%two_way_nesting_comms(child_comms_index)%nsend_interior

        test_max = maxval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer(1:b4_flux_data_p:strd))
        test_min = minval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer(1:b4_flux_data_p:strd))
        if( (test_max/max_parent_send_val > ONE_dp + err_tol) .OR. &
            (test_min/min_parent_send_val < ONE_dp - err_tol )) then

            write(log_output_unit,*) 'FAIL: Send buffer (1) range is incorrect', __LINE__, &
                __FILE__
            write(log_output_unit,*) size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer, kind=ip)
            write(log_output_unit,*) test_max
            write(log_output_unit,*) test_min
            write(log_output_unit,*) Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer(1:b4_flux_data_p)
            !call generic_stop()

        end if
        !print*, 'Us(1): ', test_max, test_min


        ! If the nest ratio is even, and we use pointwise coarse to fine sends, then include this adjustment
        ! to account for the fact that the sent value is not in the exact centre of the parent cell
        max_parent_send_val_offset = merge(dx2(1), 0.0_dp, &
            ((mod(nest_ratio, 2) == 0).and.(.not. use_averaging_for_fine_to_coarse_data_sends)) )

        test_max = maxval(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer(1:b4_flux_data_c))
        test_min = minval(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer(1:b4_flux_data_c))
        if( (test_max/(max_parent_send_val+max_parent_send_val_offset) > ONE_dp + err_tol) .OR. &
            (test_min/(min_parent_send_val-max_parent_send_val_offset) < ONE_dp - err_tol)) then
    
            write(log_output_unit,*) 'FAIL: Send buffer (2) range is incorrect', __LINE__, &
                __FILE__
            write(log_output_unit,*) Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer(1:b4_flux_data_c)
            write(log_output_unit,*) size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer, kind=ip)
            !call generic_stop()

        end if
        !print*, 'Us(2): ', test_max, test_min

#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)        
        call Us(1)%two_way_nesting_comms(parent_comms_index)%send_data(send_to_recv_buffer=.FALSE.)
        call Us(2)%two_way_nesting_comms(child_comms_index)%send_data(send_to_recv_buffer=.FALSE.)
        call communicate_p2p
#else
        ! Communicate
        call Us(1)%two_way_nesting_comms(parent_comms_index)%send_data()
        call Us(2)%two_way_nesting_comms(child_comms_index)%send_data()
#endif
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        sync images(linked_p2p_images)
#endif

        ! Copy the recv buffer to the main array 
        call Us(1)%two_way_nesting_comms(parent_comms_index)%process_received_data(Us(1)%U)
        call Us(2)%two_way_nesting_comms(child_comms_index)%process_received_data(Us(2)%U)

        if(use_averaging_for_fine_to_coarse_data_sends .or. mod(nest_ratio, 2_ip) == 1) then
            ! Check U has not changed, up to numerical precision
            ! This will hold in 2 cases, because of the linearity of the send function:
            !     a) Whenever we use averaging for fine-to-coarse sends
            !     b) If we use pointwise sends with an odd number nest ratio [so the send value is in the middle]
            ! We use the if statement above to filter out other cases
            do i = 1, 2
                if( maxval(abs(Us(i)%U - Us_store(i)%U)/Us_store(i)%U) > err_tol) then
                    write(log_output_unit,*) 'FAIL: Error in nesting communication, case ', i, __LINE__, &
                        __FILE__
                    write(log_output_unit,*) maxval(abs(Us(i)%U - Us_store(i)%U))
                    write(log_output_unit,*) maxval(abs(Us(i)%U - Us_store(i)%U)/Us_store(i)%U)
                    !call generic_stop()
                else
                    write(log_output_unit,*) 'PASS'
                end if
            end do
        end if

        call deallocate_p2p_comms

    end subroutine

    !
    ! Main unit-test routine
    !
    subroutine test_nested_grid_comms_mod

        ! Run a range of nesting ratios -- in development, some cases
        ! exposed bugs which were not apparent with other nesting ratios
        write(log_output_unit,*) '   nest ratio = 1'
        call test_nested_grid_comms_mod_workhorse(1_ip)
        write(log_output_unit,*) '   nest ratio = 2'
        call test_nested_grid_comms_mod_workhorse(2_ip)
        write(log_output_unit,*) '   nest ratio = 3'
        call test_nested_grid_comms_mod_workhorse(3_ip)
        write(log_output_unit,*) '   nest ratio = 4'
        call test_nested_grid_comms_mod_workhorse(4_ip)
        write(log_output_unit,*) '   nest ratio = 5'
        call test_nested_grid_comms_mod_workhorse(5_ip)
        write(log_output_unit,*) '   nest ratio = 7'
        call test_nested_grid_comms_mod_workhorse(7_ip)
        write(log_output_unit,*) '   nest ratio = 101'
        call test_nested_grid_comms_mod_workhorse(101_ip)

    end subroutine
end module
