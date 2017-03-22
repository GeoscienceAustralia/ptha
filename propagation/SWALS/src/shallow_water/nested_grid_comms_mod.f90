! Module for two-way nesting communication
!

module nested_grid_comms_mod

    use global_mod, only: dp, ip, charlen
    use stop_mod, only: generic_stop
    ! Note: coarray_point2point_comms_mod works in serial or parallel 
    use coarray_point2point_comms_mod, only: include_in_p2p_send_buffer, &
        send_to_p2p_comms, recv_from_p2p_comms, linked_p2p_images, &
        allocate_p2p_comms, deallocate_p2p_comms
    use reshape_array_mod, only: repack_rank1_array
    !use mpi

    use iso_c_binding, only: C_DOUBLE

    implicit none

    private
    public two_way_nesting_comms_type, test_nested_grid_comms_mod

    !
    ! Useful constants
    !
    real(dp), parameter :: ONE_dp = 1.0_dp, ZERO_dp = 0.0_dp, HALF_dp = 0.5_dp
    real(dp), parameter :: EPS = 1.0e-06_dp
    ! Dimensions of the problem -- we solve the 2D shallow water equations
    integer(ip), parameter :: SPATIAL_DIM = 2

    !
    ! The main type of the module
    !
    type :: two_way_nesting_comms_type
        !
        ! The type for two-way nesting has:
        ! -- Info on the cell indices on our grid containing values which need
        ! to be sent to another grid
        ! -- A function which copies the data [perhaps with some averaging] to
        ! a buffer for sending. 
        ! -- A buffer which can send the data. This can be a coarray, and we
        ! might not need to fill it completely, so we will probably need info
        ! on how much data to put in.
        ! -- The index of the domain that we send to, and the image_index we
        ! send to (in the coarray case)
        ! 
        ! -- Info on the cell indices on our grid which need to receive values
        !  from another grid
        ! -- A buffer which can receive data from another grid. This is
        !  probably a coarray, and we might not use all of it, so we will
        !  probably also need to know how much data to use
        ! -- A function which can convert the received data into a form
        !  suitable for our grid [e.g. if we received data from a coarser grid,
        ! then some interpolation would be required]. We can skip this if we
        ! always send the correct data, although that means more communication
        ! than is optimal [CURRENTLY SKIPPED]
        ! -- The index of the domain we receive from, and the image_index we
        ! receive from (in the coarray case) [NOTE THIS WILL BE THE SAME AS
        ! ABOVE, for 2-way nesting]
        !

        real(dp) :: cell_ratios(SPATIAL_DIM)
        !    = max( nint(my_dx/neighbour_dx), nint(neighbour_dx/my_dx) )

        logical :: my_domain_is_finer

        !
        ! Variables used to send data
        !

        !
        ! send_inds =  min_i, min_j, min_k,
        !              max_i, max_j, max_k
        integer(ip) :: send_inds(2, SPATIAL_DIM + 1) 
        real(dp), allocatable :: send_buffer(:)
        integer(ip) :: nsend ! We only need to send send_buffer(1:nsend)

        !
        ! Variables used to receive data
        !


        !
        ! recv_inds = min_i, min_j, min_k,
        !             max_i, max_j, max_k
        integer(ip) :: recv_inds(2, SPATIAL_DIM + 1) 
        real(dp), allocatable :: recv_buffer(:)
        integer(ip) :: nrecv ! We only need to use recv_buffer(1:nrecv)

        !
        ! Communication info -- the index of communicating 'domain' in the vector of 'domains',
        ! and the index of the two_way_nesting_comms_type inside the latter
        ! We are communicating with domain(neighbour_domain_index)[neighbour_domain_image_index],
        !  specifically with the 'neighbour_domain_comms_index' index of its two_way_nesting_comms array
        !
        integer(ip) :: neighbour_domain_index = -1_ip
        integer(ip) :: neighbour_domain_comms_index = -1_ip
        integer(ip) :: my_domain_index = -1_ip
        integer(ip) :: my_domain_comms_index = -1_ip

        ! ... more communication info for coarrays
#if defined(COARRAY) 
        integer(ip) :: neighbour_domain_image_index = -1_ip
        integer(ip) :: my_domain_image_index = -1_ip
#else
        integer(ip) :: neighbour_domain_image_index = 1_ip
        integer(ip) :: my_domain_image_index = 1_ip
#endif

        character(len=charlen) :: send_ID, recv_ID

        contains

        procedure :: initialise => initialise_two_way_nesting_comms
        procedure :: process_data_to_send => process_data_to_send
        procedure :: send_data => send_data
        procedure :: process_received_data => process_received_data

    end type

    contains

    ! Set up to way nesting communicator type
    !
    ! @param two_way_nesting_comms The two-way-nesting-comms-type to set up
    ! @param my_dx [dx,dy] of my domain
    ! @param neighbour_dx [dx, dy] of neighbour domain
    ! @param ijk_to_send array with send_inds = [[min_i, min_j, min_k],
    !                                            [max_i, max_j, max_k]]
    !    giving the ranges of indices from which we extract data to send
    ! @param ijk_to_recv array with recv_inds = [[min_i, min_j, min_k],
    !                                            [max_i, max_j, max_k]]
    !    giving the ranges of indices from which we will put received data
    ! @param neighbour_domain_index integer giving the index of the neighbour
    !    domain (assuming all domains are in a vector)
    ! @param my_domain_index integer giving the index of the my domain
    !    (assuming all domains are in a vector)
    ! @param neighbour_domain_image_index optional image_index holding
    !    neighbour domain (for coarrays)
    ! @param my_domain_image_index optional image_index holding my domain (for
    !    coarrays)
    !
    subroutine initialise_two_way_nesting_comms(&
        two_way_nesting_comms, &
        my_dx, neighbour_dx, &
        ijk_to_send, ijk_to_recv, &
        neighbour_domain_index, my_domain_index, &
        neighbour_domain_comms_index, my_domain_comms_index,&
        neighbour_domain_image_index, my_domain_image_index)

        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
        real(dp), intent(in) :: my_dx(SPATIAL_DIM), neighbour_dx(SPATIAL_DIM)
        integer(ip), intent(in) :: ijk_to_send(2, SPATIAL_DIM+1), &
            ijk_to_recv(2, SPATIAL_DIM+1)
        integer(ip), intent(in) :: neighbour_domain_index, my_domain_index
        integer(ip), intent(in) :: neighbour_domain_comms_index, my_domain_comms_index
        integer(ip), optional, intent(in) :: neighbour_domain_image_index, &
            my_domain_image_index

        ! Local variables
        real(dp) :: cell_ratios(SPATIAL_DIM)
        integer(ip):: iL, iU, jL, jU , kL, kU, sbs, rbs
        character(charlen) :: n_char1, n_char2, n_char3, char1, char2, char3
      
        ! Find the ratios of the cell sizes 
        cell_ratios = my_dx/neighbour_dx

        if( cell_ratios(1) >= ONE_dp ) then
            ! The current domain is 'coarser' than the other, or the same size
            two_way_nesting_comms%my_domain_is_finer = .FALSE.

            ! Perform various sanity checks

            ! Both dimensions should be coarser or equal
            if(cell_ratios(2) < ONE_dp - EPS) then
                print*, 'Nesting cell dimensions not consistently larger'
                call generic_stop()
            end if
           
            ! cell_ratios should be an integer (up to floating point) 
            if(any(abs(cell_ratios - nint(cell_ratios)) > EPS)) then
                print*, 'Apparent non-integer cell ratios'
                call generic_stop()
            end if

        else
            ! The current domain is 'finer' than the other
            two_way_nesting_comms%my_domain_is_finer = .TRUE.

            ! Both dimensions should be finer or equal
            if(cell_ratios(2) > ONE_dp + EPS) then
                print*, 'Nesting cell dimensions not consistently smaller'
                call generic_stop()
            end if

            ! 1/cell_ratios  should be an integer (up to floating point)
            if(any(abs(ONE_dp/cell_ratios - nint(ONE_dp/cell_ratios)) > EPS)) then
                print*, 'Apparent non-integer cell ratios'
                call generic_stop()
            end if

        end if

        two_way_nesting_comms%cell_ratios = cell_ratios

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
        write(n_char1, '(I0)') two_way_nesting_comms%neighbour_domain_index
        write(n_char2, '(I0)') two_way_nesting_comms%neighbour_domain_comms_index
        write(n_char3, '(I0)') two_way_nesting_comms%neighbour_domain_image_index
        write(char1, '(I0)') two_way_nesting_comms%my_domain_index
        write(char2, '(I0)') two_way_nesting_comms%my_domain_comms_index
        write(char3, '(I0)') two_way_nesting_comms%my_domain_image_index

        ! Set up ijk to send/recv
        two_way_nesting_comms%send_inds = ijk_to_send
        two_way_nesting_comms%recv_inds = ijk_to_recv

        ! Shorthand
        iL = ijk_to_send(1,1)
        iU = ijk_to_send(2,1)
        jL = ijk_to_send(1,2)
        jU = ijk_to_send(2,2)
        kL = ijk_to_send(1,3)
        kU = ijk_to_send(2,3)

        ! Sanity check the send/recv indices
        if(two_way_nesting_comms%my_domain_is_finer) then 
            ! Finer grid should be perfectly contained in coarser grid
            ! (i.e. iU - iL + 1 should be a multiple of the relative cell sizes)
            if(mod(iU - iL + 1, int(nint(ONE_dp/cell_ratios(1)))) /= 0) then
                print*, 'Imperfectly nested i'
                call generic_stop()
            end if

            if(mod(jU - jL + 1, int(nint(ONE_dp/cell_ratios(2)))) /= 0) then
                print*, 'Imperfectly nested j'
                call generic_stop()
            end if
        end if

        if( (iL > iU).OR.(jL > jU).OR.(kL > kU) ) then
            print*, 'indices to send has min > max'
            call generic_stop()
        end if
        
        ! Figure out how much data we need to send 
        ! Make this the same length as the region we will send to
        if(two_way_nesting_comms%my_domain_is_finer) then
            sbs = (iU - iL + 1) * (jU - jL + 1) * (kU - kL + 1) / &
                int(nint(product(ONE_dp/cell_ratios)))
        else
            sbs = (iU - iL + 1) * (jU - jL + 1) * (kU - kL + 1) * &
                int(nint(product(cell_ratios)))
        end if
        two_way_nesting_comms%nsend = sbs

        ! Figure out how much data we need to receive
        rbs = product(ijk_to_recv(2,:) - ijk_to_recv(1,:) + 1)
        two_way_nesting_comms%nrecv = rbs

        ! Allocate buffers    
        allocate(two_way_nesting_comms%send_buffer(sbs))
        allocate(two_way_nesting_comms%recv_buffer(rbs))

        !
        ! Prepare communication
        !
        two_way_nesting_comms%send_ID = 'nesting_from_' // (trim(char1)) // '_' // &
            (trim(char2)) // '_' // (trim(char3)) // '_to_' // &
            (trim(n_char1)) // '_' // (trim(n_char2)) // '_' // &
            (trim(n_char3))
        
        two_way_nesting_comms%recv_ID = 'nesting_from_' // (trim(n_char1)) // '_' // &
            (trim(n_char2)) // '_' // (trim(n_char3)) // '_to_' // &
            (trim(char1)) // '_' // (trim(char2)) // '_' // &
            (trim(char3))

        ! Make sure the point2point coarray communicator has space for this variable
        call include_in_p2p_send_buffer(two_way_nesting_comms%send_buffer, &
            buffer_label = two_way_nesting_comms%send_ID, &
            receiver_image = two_way_nesting_comms%neighbour_domain_image_index)

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
        integer(ip) :: i,j,k
        integer(ip) :: iL, iU, jL, jU, kL, kU
        integer(ip) :: iCounter, jCounter, kCounter, ijkCounter
        integer(ip) :: iR, jR, kR
        integer(ip) :: inv_cell_ratios_ip(2), cell_ratios_ip(2)
        integer(ip) :: m, n, send_ind
        real(dp) :: product_cell_ratios, inv_product_cell_ratios
        real(dp) :: inv_cell_ratios(2), cell_ratios(2)
        real(dp) :: dU_di, dU_dj, send_val

        ! **Clear** the send buffer
        two_way_nesting_comms%send_buffer = ZERO_dp

        ! Define useful constants
        cell_ratios = two_way_nesting_comms%cell_ratios
        inv_cell_ratios = ONE_dp / cell_ratios 
        inv_cell_ratios_ip = int(nint(inv_cell_ratios))
        cell_ratios_ip = int(nint(cell_ratios))
        product_cell_ratios = product(cell_ratios)
        inv_product_cell_ratios = ONE_dp/product_cell_ratios

        ! Lower/upper spatial indices involved in nesting
        iL = two_way_nesting_comms%send_inds(1, 1)
        iU = two_way_nesting_comms%send_inds(2, 1)
        jL = two_way_nesting_comms%send_inds(1, 2)
        jU = two_way_nesting_comms%send_inds(2, 2)

        ! Lower/upper variable indices involved in nesting.
        ! (NOTE, the code logic currently assumes a sequence of variables
        !  is sent, so e.g. if variable 2,4 are sent, then variable 3 should be
        !  as well) 
        kL = two_way_nesting_comms%send_inds(1, 3)
        kU = two_way_nesting_comms%send_inds(2, 3)

        ! Number of variables that are sent
        kR = (kU - kL + 1)

        ! Number of i/j cells in communication zone
        if(two_way_nesting_comms%my_domain_is_finer) then
            ! iR will be less than (iU - iL + 1) by some integer factor
            iR = int( nint((iU - iL + 1) * cell_ratios(1)) )
            ! Same for jR
            jR = int( nint((jU - jL + 1) * cell_ratios(2)) )
        else
            ! iR will be greater than (iU - iL + 1) by some integer factor
            iR = int( nint((iU - iL + 1) * cell_ratios(1)) )
            ! Same for jR
            jR = int( nint((jU - jL + 1) * cell_ratios(2)) )
        end if

        ! We use iCounter, jCounter, kCounter to track the index of the send
        ! buffer that needs modification. 
        iCounter = 0
        jCounter = 0
        kCounter = 0
        ! Loop over variables (third index of U)
        do k = kL, kU

            kCounter = (k - kL) + 1

            ! Loop over spatial dimensions, with operations depending on
            ! whether we do a 'fine-to-coarse' or a 'coarse-to-fine' op

            if(two_way_nesting_comms%my_domain_is_finer) then
                ! Fine-to-coarse

                jCounter = 0
                do j = jL, jU

                    ! Every time j-jL grows by 1/cell_ratios(2), we are
                    ! modifying another 'jCounter' index in the send buffer
                    if( mod(j - jL, inv_cell_ratios_ip(2)) == 0 ) then
                        jCounter = jCounter + 1
                    end if

                    iCounter = 0
                    do i = iL, iU

                        ! Every time i-iL grows by 1/cell_ratios(1), we are
                        ! modifying another 'iCounter' index in the send buffer
                        if(mod(i - iL, inv_cell_ratios_ip(1)) == 0) then
                            iCounter = iCounter + 1
                        end if

                        ! Since the send_buffer is a rank 1 array, we have to
                        ! pack the data into it
                        ijkCounter = iCounter + (jCounter - 1) * iR + &
                            (kCounter - 1) * iR * jR 

                        ! Average all 'fine' cells inside the coarse cell
                        two_way_nesting_comms%send_buffer(ijkCounter) = &
                            two_way_nesting_comms%send_buffer(ijkCounter) + &
                            product_cell_ratios * U(i,j,k)

                    end do
                end do

            else
                ! Coarse-to-fine

                do j = jL, jU
                    do i = iL, iU
                        ! Local derivarives. Note this assumes that the send
                        ! region is not on a boundary of U.
                        dU_di = HALF_dp * ( U(i+1 , j, k) - U(i-1, j, k) )
                        dU_dj = HALF_dp * ( U(i , j+1, k) - U(i, j-1, k) )
                        ! FIXME: The definition above could cause issues at
                        ! wet-dry fronts, or shocks.

                        ! Extend the coarse cell values to a value for each
                        ! fine cell, using linear interpolation
                        do m = 1, cell_ratios_ip(2)
                            do n = 1, cell_ratios_ip(1)
                                ! Compute the value we should send
                                send_val = U(i,j,k) + &
                                    dU_di * ((n - HALF_dp) * inv_cell_ratios(2) - HALF_dp) + &
                                    dU_dj * ((m - HALF_dp) * inv_cell_ratios(1) - HALF_dp)

                                ! Compute the index of the send_buffer that we should send it to
                                send_ind = &
                                    ! Number of cells for previous variables
                                    (iU - iL + 1) * (jU - jL + 1) * (k - kL) * int(nint(product_cell_ratios)) + &
                                    ! Number of cells for previous j, with this k
                                    (iU - iL + 1) * (j - jL) * int(nint(product_cell_ratios)) + &
                                    ! Number of cells for previous m
                                    (iU - iL + 1) * (m - 1) * cell_ratios_ip(1) + &
                                    ! Cells for i and n
                                    (i - iL) * cell_ratios_ip(1) + n

                                ! Debugging
                                !if(cell_ratios_ip(1) == 7_ip) then
                                !    print*, k, j, i, m, n, send_ind, send_val, &
                                !        (iU - iL + 1) * (jU - jL + 1) * (k - kL) * int(nint(product_cell_ratios)), &
                                !        (iU - iL + 1) * (j - jL) * int(nint(product_cell_ratios)), &
                                !        (iU - iL + 1) * (m - 1) * cell_ratios_ip(1),  &
                                !        (i - iL) * cell_ratios_ip(1) + n
                                !end if

                                two_way_nesting_comms%send_buffer(send_ind) = send_val
                            end do !n
                        end do !m
                    end do !i
                end do !j
            end if !my_domain_is_finer
        end do !k

    end subroutine

    ! Copy send_buffer to the right recv_buffer. This should be called AFTER process_data_to_send.
    !
    subroutine send_data(two_way_nesting_comms)
        class(two_way_nesting_comms_type), intent(inout) :: two_way_nesting_comms
  
        ! Send to point2point coarray communication buffer 
        call send_to_p2p_comms(two_way_nesting_comms%send_buffer, &
            buffer_label=two_way_nesting_comms%send_ID, &
            put_in_recv_buffer=.TRUE.)
     
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
  
        integer(ip) :: iL, iU, jL, jU, kL, kU, nrecv 

        ! Shorthand
        iL = two_way_nesting_comms%recv_inds(1,1)
        iU = two_way_nesting_comms%recv_inds(2,1)
        jL = two_way_nesting_comms%recv_inds(1,2)
        jU = two_way_nesting_comms%recv_inds(2,2)
        kL = two_way_nesting_comms%recv_inds(1,3)
        kU = two_way_nesting_comms%recv_inds(2,3)

        nrecv = two_way_nesting_comms%nrecv

        ! Copy from the point2point coarray communicator
        !call recv_from_p2p_comms(U(iL:iU, jL:jU, kL:kU), & 
        !    buffer_label = two_way_nesting_comms%recv_ID)

        ! FIXME: We can eliminate recv_buffer here. But the tests rely on it.
        !two_way_nesting_comms%recv_buffer = reshape(U(iL:iU, jL:jU, kL:kU), [nrecv])

        call recv_from_p2p_comms(two_way_nesting_comms%recv_buffer(:), &
            buffer_label = two_way_nesting_comms%recv_ID)

        !U(iL:iU, jL:jU, kL:kU) = reshape(two_way_nesting_comms%recv_buffer, &
        !    [iU - iL + 1, jU - jL + 1, kU - kL + 1])
        call repack_rank1_array(two_way_nesting_comms%recv_buffer, &
            U(iL:iU, jL:jU, kL:kU))

        
     
    end subroutine

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
        integer(ip) :: nx2_on_nest_ratio(2) = [9_ip, 8_ip]
        integer(ip) :: nx2(2) 
        real(dp) :: dx2(2)

        integer(ip), parameter :: parent_comms_index = 3, child_comms_index = 2

        ! Type which serves the role of the array holding the domains, and something
        ! to copy the above into, to ensure no changes due to communication
        type array_of_arrays_r3_dp_type
            real(dp), allocatable :: U(:,:,:)
            type(two_way_nesting_comms_type), allocatable :: two_way_nesting_comms(:)
        end type

        type(array_of_arrays_r3_dp_type) :: Us(2), Us_store(2)

        ! Local variables
        integer(ip) :: i, j, k, xi, yi, this_image_local, num_images_local, ierr
        integer(ip) :: parent_image_index, child_image_index
        real(dp) :: x, y, z, err_tol

#ifdef COARRAY
        this_image_local = this_image()
        num_images_local = num_images()
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
            my_domain_image_index = this_image_local)

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
            my_domain_image_index = this_image_local)

        call allocate_p2p_comms

        !
        ! Basic tests of cell ratios
        !
        if(any(abs(Us(2)%two_way_nesting_comms(child_comms_index)%cell_ratios - 1.0_dp/nest_ratio) > EPS)) then
            print*, 'FAIL: Problems with cell ratios'
            call generic_stop()
        else
            print*, 'PASS'
        end if

        if(any(abs(Us(1)%two_way_nesting_comms(parent_comms_index)%cell_ratios - nest_ratio) > EPS)) then
            print*, 'FAIL: Problems with cell ratios'
            call generic_stop()
        else
            print*, 'PASS'
        end if

        !
        ! Send/recv buffers have corresponding sizes
        !
        if(size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer) /= &
            size(Us(2)%two_way_nesting_comms(child_comms_index)%recv_buffer)) then

            print*, 'FAIL: Problems with send/recv buffer sizes'
            print*, size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)
            print*, size(Us(2)%two_way_nesting_comms(child_comms_index)%recv_buffer)
            call generic_stop()
        else
            print*, 'PASS'
        end if

        if(size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer) /= &
            size(Us(1)%two_way_nesting_comms(parent_comms_index)%recv_buffer)) then

            print*, 'FAIL: Problems with send/recv buffer sizes'
            print*, size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer)
            print*, size(Us(1)%two_way_nesting_comms(parent_comms_index)%recv_buffer)
            call generic_stop()
        else
            print*, 'PASS'
        end if

        !
        ! my_domain_is_finer is correctly set
        !
        if(nest_ratio > 1) then
            if(Us(1)%two_way_nesting_comms(parent_comms_index)%my_domain_is_finer) then
                print*, 'FAIL: Problem with fine/coarse classification'
                call generic_stop()
            else
                print*, 'PASS'
            end if        

            if(.NOT.Us(2)%two_way_nesting_comms(child_comms_index)%my_domain_is_finer) then
                print*, 'FAIL: Problem with fine/coarse classification'
                call generic_stop()
            else
                print*, 'PASS'
            end if
        end if

        ! Fill the send buffer
        call Us(1)%two_way_nesting_comms(parent_comms_index)%process_data_to_send(Us(1)%U)
        call Us(2)%two_way_nesting_comms(child_comms_index)%process_data_to_send(Us(2)%U)

        ! Check for errors in send buffer
        if(dp == C_DOUBLE) then
            err_tol = 1.0e-14_dp
        else
            err_tol = 1.0e-5_dp
        end if
        
        if( (maxval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)/32.0_dp > ONE_dp + err_tol) .OR. &
            (minval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)/8.0_dp < ONE_dp - err_tol )) then

            print*, 'FAIL: Send buffer (1) range is incorrect'
            print*, size(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)
            print*, maxval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)
            print*, minval(Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer)
            print*, Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer
            call generic_stop()

        end if
        if( (maxval(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer)/32.0_dp  > ONE_dp + err_tol) .OR. &
            (minval(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer)/8.0_dp < ONE_dp - err_tol)) then
    
            print*, 'FAIL: Send buffer (2) range is incorrect'
            print*, Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer
            print*, size(Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer)
            call generic_stop()

        end if

        ! Communicate
        !Us(1)%two_way_nesting_comms(parent_comms_index)%recv_buffer = Us(2)%two_way_nesting_comms(child_comms_index)%send_buffer
        !Us(2)%two_way_nesting_comms(child_comms_index)%recv_buffer = Us(1)%two_way_nesting_comms(parent_comms_index)%send_buffer
        call Us(1)%two_way_nesting_comms(parent_comms_index)%send_data()
        call Us(2)%two_way_nesting_comms(child_comms_index)%send_data()
#ifdef COARRAY
        sync images(linked_p2p_images)
#endif

        ! Copy the recv buffer to the main array 
        call Us(1)%two_way_nesting_comms(parent_comms_index)%process_received_data(Us(1)%U)
        call Us(2)%two_way_nesting_comms(child_comms_index)%process_received_data(Us(2)%U)

        ! Should not have changed, up to floating point
        ! Check U has not changed
        do i = 1, 2
            if( maxval(abs(Us(i)%U - Us_store(i)%U)/Us_store(i)%U) > err_tol) then
                print*, 'FAIL: Error in nesting communication, case ', i
                print*, maxval(abs(Us(i)%U - Us_store(i)%U))
                call generic_stop()
            else
                print*, 'PASS'
            end if
        end do

        call deallocate_p2p_comms

    end subroutine

    !
    ! Main unit-test routine
    !
    subroutine test_nested_grid_comms_mod

        ! Run a range of nesting ratios -- in development, some cases
        ! exposed bugs which were not apparent with other nesting ratios
        print*, '   nest ratio = 1'
        call test_nested_grid_comms_mod_workhorse(1_ip)
        print*, '   nest ratio = 2'
        call test_nested_grid_comms_mod_workhorse(2_ip)
        print*, '   nest ratio = 3'
        call test_nested_grid_comms_mod_workhorse(3_ip)
        print*, '   nest ratio = 4'
        call test_nested_grid_comms_mod_workhorse(4_ip)
        print*, '   nest ratio = 5'
        call test_nested_grid_comms_mod_workhorse(5_ip)
        print*, '   nest ratio = 7'
        call test_nested_grid_comms_mod_workhorse(7_ip)
        print*, '   nest ratio = 101'
        call test_nested_grid_comms_mod_workhorse(101_ip)

    end subroutine
end module
