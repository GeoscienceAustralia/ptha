module coarray_point2point_comms_mod
!
! This module was originally written to use coarrays, but it can also be compiled to only require mpi.
! Many comments refer only to the coarray case.
!


!
!! Module for doing coarray 'point-to-point' communications of real
!! arrays (rank from 1 to 4), in a 'multiple-program, multiple data' style.
!
!! Send an arbitrary (and varying) number of arrays between any pairs of images.
!! The array dimensions do not have to be consistent among images.
!!
!! Background
!! ------------
!!
!! Coarrays by definition must be the same size on each image. It is not always
!! straightforward to apply coarrays to problems where we need to: 
!!    A) send different numbers of variables between pairs of communicating images, or; 
!!    B) send different sizes of variables between pairs of communicating images.
!!
!! For example, consider a 2D grid-based PDE solver, where a single program might
!! contain multiple nested grids (to allow high resolution in target areas). Each
!! of these grids can be partitioned across a number of images, which communicate 
!! boundary data to each other (i.e. halo exchange). The nesting
!! boundaries will force an irregular communication pattern on the code, depending on
!! the geometric placement of each grid, as well as on how they are partitioned.
!!
!! An approach to the 'irregular communication problem', implemented here, is to
!! use one coarray to do all send/recv communications. By definition, 
!! this gives a separate chunk of communicable memory (of the same size) on each
!! image. For each array 'x' that we want to communicate, we use a distinct
!! contiguous slice of the coarray for communication. 
!!
!! This module takes care of the details of communication, and provides the user
!! with a 'hopefully' simple interface.
!!
!! Note that memory will necessarily be wasted if the amount of data to
!! receive is unequal among images. In the target applications, the amount
!! of data to send is quite small compared with the overall memory usage, and in
!! such cases this wasted memory is insignificant.
!!
!! Usage
!! -----
!! 
!!     ! The general outline of the usage is below, ignoring details which are
!!     ! not relevant to the application (such as allocating variables, treating
!!     ! cases with 'this_image - 1 == zero', etc). 
!!     ! For a 'real' example, see the unit-test subroutine
!!     ! 'test_coarray_point2point_comms_mod'.
!!
!!
!!     ! import the module
!!     use coarray_point2point_comms_mod
!!
!!     ! Communicate arrays with possibly varying shapes on each image
!!     real(dp) :: x_send(10), y_send(this_image())
!!     real(dp) :: x_recv(10), y_recv(this_image() - 1) 
!!
!!     ! .... define x_send, y_send ....
!!
!!     ! Let the module know which arrays will be communicated, and where they
!!     ! will be received. Each set of array data is given a label which will be
!!     ! used in sending and receiving. Note we do not have to communicate the
!!     ! same arrays on all images
!!
!!     ! associate the label 'y_comms1' with communication of the 'y_' variables to
!!     ! image 'this_image() + 1' (for simplicity of exposition, ignore the special 
!!     ! treatment when (this_image()+1) > num_images() )
!!     call include_in_p2p_send_buffer(y_send, buffer_label='y_comms1', receiver_image=this_image() + 1_ip)
!!     if(this_image() == 1) then
!!         ! If we are on image 1, then associate the label 'x_comms1' with
!!         ! communication of the 'x_' variables to image 2
!!         call include_in_p2p_send_buffer(x_send, buffer_label='x_comms1', receiver_image=2_ip)
!!     end if
!!     
!!     ! Once we have identified which arrays are communicated using
!!     ! 'include_in_p2p_send_buffer', we can allocate communication buffers.
!!     ! At this stage, we should not call 'include_in_p2p_send_buffer' anymore.
!!     ! All images should call this subroutine, even if they are not sending or receiving,
!!     ! because the coarray must be allocated on all images at once.
!!     call allocate_p2p_comms
!!
!!     ! Now the main computational work begins
!!     do while (.... main iteration loop ... )
!!         
!!         !... update x, y, based on values at previous iteration
!!         
!!         ! Copy x,y to the receiver image 
!!         !
!!         ! send y_send to image (this_image()+1) -- note the receive image was
!!         ! defined above
!!         call send_to_p2p_comms(y_send, buffer_label='y_comms1')
!!         if(this_image() == 1) then
!!             ! send x_send to image 2
!!             call send_to_p2p_comms(x_send, buffer_label='x_comms1')
!!         end if
!!
!!         ! Make sure images we receive from have sent their data.
!!         ! This needs to be called on all communicating images
!!         if(size(linked_p2p_images) > 0) sync images(linked_p2p_images)
!!         ! linked_p2p_images is an array containing all image_indexes that
!!         ! we send to or receive from
!!
!!         ! Copy from recv buffer back to computational array
!!         call recv_from_p2p_comms(y_recv, buffer_label='y_comms1')
!!         if(this_image() == 2) then
!!             call recv_from_p2p_comms(x_recv, buffer_label='x_comms1')
!!         end if
!!       
!!         ! Do something with the received data, update x_, y_, more computations ....
!!
!!     end do

    ! kinds for double, integer, and default character length
    use global_mod, only: dp, ip, charlen, real_bytes
    ! Use real64 as a send buffer for the buffer_label character id's
    ! Use int32 for ints that opencoarrays sends [doesn't yet support e.g. int64]
    use iso_fortran_env, only: real64, int32
    ! routines to efficiently convert between rank1 and rankN arrays (n=1,2,3,4)
    use reshape_array_mod, only: flatten_array, repack_rank1_array
    use logging_mod, only: log_output_unit
    use qsort_mod, only: sort_index
    use iso_c_binding, only: c_int ! For call to sort_index, which uses C
#if defined(COARRAY_PROVIDE_CO_ROUTINES)
    use coarray_intrinsic_alternatives, only: co_broadcast, co_max, co_sum, sync_all_generic, this_image2, num_images2
#elif defined(COARRAY)
    use coarray_intrinsic_alternatives, only: sync_all_generic, this_image2, num_images2
#endif
#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    ! The 'inner loop' communication will use mpi rather than coarrays
    use mpi
#endif
    implicit none


    private

    ! Open-coarrays integer precision -- use this because opencoarrays cannot send
    ! any integer precision (without effort)
    integer, parameter :: ocaIP = int32

    public :: test_coarray_point2point_comms_mod
        ! Unit test subroutine

    ! Methods we need to use
    public :: include_in_p2p_send_buffer
        !! Append another point2point send to the data structures
    public :: allocate_p2p_comms
        !! Allocate arrays to do the point2point sends
    public :: deallocate_p2p_comms
        !! Clean up arrays to do the point2point sends
    public :: send_to_p2p_comms, recv_from_p2p_comms
        !! Send and unpack the buffers
    public :: communicate_p2p 
        !! Do all sends
    public :: print_p2p_comms
        !! Print info about the point2point comms data structures
    public :: size_of_send_recv_buffers
        !! Compute the size of the send/receive buffers

    ! This must be public to allow control of syncs
    public :: linked_p2p_images
    protected :: linked_p2p_images
        !! Array with the coarray indices of images we communicate with

    ! Useful to have the send/recv buffers publically visible for debugging
    ! But it kills the parallel efficiency on NCI for some reason! At 256
    ! cores, the run-time nearly doubled, with a huge increase in the time spent
    ! in comms
    !public :: send_buffer, recv_buffer
    !protected :: send_buffer, recv_buffer

    public:: integer_to_id
        !! This is useful, if we need to make distinct buffer_labels based on
        !! integers



    !
    ! Key module data below here.
    !


    ! len of character used for buffer_label
    integer(ocaIP), parameter :: p2p_id_len = charlen

    logical, parameter :: reorder_send_data_by_image = .true.
        !! Do we order the send data by image? If true, then 'communicate_p2p' can
        !! involve less individual sends. There does not seem to be any reason
        !! not to do this.

    !
    ! Main send buffer + metadata about what we send
    ! These are private to the module. 
    ! Ideally they would live inside a derived-type (which would permit a program with >1 multidomain).
    ! But is the compiler support up to this for coarrays?
    real(dp), allocatable :: send_buffer(:)
    integer(ip), allocatable :: send_start_index(:)
    integer(ip), allocatable :: send_size(:)
    integer(ocaIP), allocatable :: sendto_image_index(:)
    integer(ip), allocatable :: sendto_start_index(:)
    character(len=p2p_id_len), allocatable :: send_buffer_label(:)
    
    ! Note: the 'ith' send_array is sent from:
    !     send_buffer(send_start_index(i) + (0:(send_size(i) - 1) )  )
    ! to:
    !     recv_buffer(sendto_start_index(i) + (0:(send_size(i)-1)) )[sendto_image_index(i)]

    !
    ! Main receive buffer + metadata about what we receive
    ! These are private to the module
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
    real(dp), allocatable :: recv_buffer(:)[:]
#else
    real(dp), allocatable :: recv_buffer(:)
#endif
    integer(ip), allocatable :: recv_start_index(:)
    integer(ip), allocatable :: recv_size(:)
    integer(ocaIP), allocatable :: recvfrom_image_index(:)
    character(len=p2p_id_len), allocatable :: recv_buffer_label(:)

    !
    ! Note: the 'ith' recv_array comes from some part of the
    ! send_buffer on [recvfrom_image_index(i)]. We don't store the exact
    ! slice that it originates from, since we use 'put' communication here

#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
    ! We need to communicate various integer arrays when initialising
    integer(ocaIP), allocatable :: work_coarray(:, :)[:]
    ! This is used to communicate a character of up to length p2p_id_len, using
    ! transfer
    real(real64), allocatable  :: real64_coarray(:,:)[:]
#else
    integer(ocaIP), allocatable :: work_coarray(:,:)
    real(real64), allocatable :: real64_coarray(:,:)
#endif

    ! store indices of all images we send to or receive from
    ! Make this public for the user to control sync's
    integer(ocaIP), allocatable :: linked_p2p_images(:)

    ! Useful to ensure we allocate/deallocate as required
    logical :: have_allocated_p2p_comms = .FALSE.

    ! Use these variables instead of calls to this_image(), num_images().
    ! Trick to generalise the code to work with or without coarrays
    ! (in the 'without' case, we are in serial). 
    integer(ocaIP) :: this_image_local
    integer(ocaIP) :: num_images_local

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    ! Variables to use with MPI_IAllToAllv. There are requirements on kind of
    ! integer, etc.
    integer, allocatable:: mympi_send_counts(:), mympi_recv_counts(:) 
    integer, allocatable:: mympi_send_displacements(:), mympi_recv_displacements(:)
    ! Variables required for MPI_ISEND and MPI_IRECV. This is an alternative to using alltoallv.
    integer, allocatable :: mpi_recv_requests(:), mpi_send_requests(:)

    logical, parameter :: mpi_timestep_loop_use_alltoallv = .false.
        ! If TRUE then use mpi_alltoallv in the inner communication loop
        ! Otherwise use mpi_isend/irecv

#ifdef REALFLOAT
    integer :: mympi_dp = MPI_REAL
#else
    integer :: mympi_dp = MPI_DOUBLE_PRECISION
#endif
    integer :: mympi_int = MPI_INTEGER
    integer :: mympi_real64 = MPI_DOUBLE_PRECISION
#endif

    interface include_in_p2p_send_buffer
        !! Allow 'include_in_p2p_send_buffer' to apply to input arrays with rank from 1 to 4
        module procedure include_in_p2p_send_buffer_rank1, include_in_p2p_send_buffer_rank2, &
            include_in_p2p_send_buffer_rank3, include_in_p2p_send_buffer_rank4
    end interface

    interface send_to_p2p_comms
        !! Allow 'send_to_p2p_comms' to apply to input arrays with rank from 1 to 4
        module procedure send_to_p2p_comms_rank1, send_to_p2p_comms_rank2, &
            send_to_p2p_comms_rank3, send_to_p2p_comms_rank4
    end interface

    ! Allow 'recv_from_p2p_comms' to apply to output arrays with rank from 1 to 4
    interface recv_from_p2p_comms
        module procedure recv_from_p2p_comms_rank1, recv_from_p2p_comms_rank2, &
            recv_from_p2p_comms_rank3, recv_from_p2p_comms_rank4
    end interface
   
    contains

    subroutine include_in_p2p_send_buffer_generic(send_array_size, buffer_label, &
        receiver_image)
        !!
        !! Define the need for space to communicate 'send_array' to image 'receiver_image'.
        !!
        !! Identify the communication with a string 'buffer_label', which can also
        !! be used to receive the sent data.
        !!
        !! Note the actual allocations happen later (once we know all the arrays
        !! we'd like to send).

        integer(ip), intent(in) :: send_array_size !! A one dimensional real array with kind dp
        character(len=p2p_id_len), intent(in) :: buffer_label
            !! A character string (len=p2p_id_len) used to
            !! identify the communication data in both the send and recv buffers
        integer(ip), intent(in) :: receiver_image !! image index which will receive the data

        if ( allocated(send_buffer) ) then
            write(log_output_unit,*) 'send_buffer is already allocated. Cannot create ', & 
                'more send_buffer space after allocation'
            call local_stop
        end if

        !write(log_output_unit,*) '    DEBUG:', trim(buffer_label), ' ', send_array_size, receiver_image, allocated(send_size), size(send_size)

        !
        ! Append space for the send_array to the metadata describing the
        ! send_buffer
        !
        if(allocated(send_size)) then
            ! ! Array will go in: 
            ! send_buffer( &
            !    send_start_index(buffer_label) + &
            !    (0:(send_size(send_array_coomms_id) - 1)) &
            !    )
            send_start_index = [send_start_index, &
                send_start_index(size(send_start_index, kind=ip)) + &
                    send_size(size(send_size, kind=ip))]
            send_size = [send_size, int(send_array_size, ip)]
            sendto_image_index = [sendto_image_index, int(receiver_image, ocaIP)]
            if(.not. any(linked_p2p_images == receiver_image)) then
                linked_p2p_images = [linked_p2p_images, int(receiver_image, ocaIP)]
            end if

            send_buffer_label = [send_buffer_label, buffer_label]
        else
            ! Array will go in send_buffer(1:size(send_array))
            send_size = [int(send_array_size, ip)]
            send_start_index = [1]
            sendto_image_index = [int(receiver_image,ocaIP)]
            linked_p2p_images = sendto_image_index
            send_buffer_label = [buffer_label]
        end if

    end subroutine

    ! For rank1 send_array's
    subroutine include_in_p2p_send_buffer_rank1(send_array, buffer_label, &
        receiver_image)
        !! Make space to send 'send_array' to the image 'receiver_image'. The send is associated
        !! with the buffer_label, which will be used to hide the book-keeping.
        real(dp), intent(in) :: send_array(:) !! A rank-1 array with the data to send
        character(len=p2p_id_len), intent(in) :: buffer_label !! The label associated with the sent data
        integer(ip), intent(in) :: receiver_image !! The image that the data should be sent to

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank2 send_array's
    subroutine include_in_p2p_send_buffer_rank2(send_array, buffer_label, &
        receiver_image)
        real(dp), intent(in) :: send_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank3 send_array's
    subroutine include_in_p2p_send_buffer_rank3(send_array, buffer_label, &
        receiver_image)
        real(dp), intent(in) :: send_array(:,:, :)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank4 send_array's
    subroutine include_in_p2p_send_buffer_rank4(send_array, buffer_label, &
        receiver_image)
        real(dp), intent(in) :: send_array(:,:, :, :)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(n, buffer_label, &
            receiver_image)
    end subroutine

    subroutine allocate_p2p_comms
        !!
        !! Allocate the send_buffer, assuming all calls to
        !! include_in_p2p_send_buffer have already been made
        !!

        integer(ocaIP) :: desired_size, desired_size_local, i, j, n, k, n2
        character(p2p_id_len) :: charlabel

        ! Variables for call to sort_index
        integer(c_int), allocatable :: send_data_order(:), send_data_sort_criterion(:)
        integer(c_int) :: n1
#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS        
        ! MPI variables
        integer :: ierr
        integer, allocatable :: wc3_all(:,:)
#endif

        if(have_allocated_p2p_comms) then
            ! Send message to a few places to make it more likely we see it
            print*, 'Error: trying to allocate p2p comms when it is already allocated'
            write(log_output_unit,*) 'Error: trying to allocate p2p comms when it is already allocated'
            error stop
        end if

        have_allocated_p2p_comms = .TRUE.

#ifdef COARRAY
        this_image_local = this_image2()
        num_images_local = num_images2()
        ! Approach using MPI
        !call mpi_comm_rank(MPI_COMM_WORLD, this_image_local, ierr)
        !this_image_local = this_image_local + 1
        !call mpi_comm_size(MPI_COMM_WORLD, num_images_local, ierr)
#else
        this_image_local = 1
        num_images_local = 1
#endif

        ! Allocate the send buffer, to have the maximum required size.
        if(allocated(send_size)) then
            desired_size = sum(send_size)
        else
            desired_size = 0
        end if
        allocate(send_buffer(desired_size))
        !call co_max(desired_size)
        !allocate(send_buffer(desired_size)[*])

        ! Fill with a value which is suggestive of problems, in case of
        ! out-of-bounds mistakes
        if(size(send_buffer, kind=ip) > 0) send_buffer = HUGE(1.0_dp)

        if(allocated(send_size) .OR. allocated(sendto_image_index)) then
            if(size(send_size, kind=ip) /= size(sendto_image_index, kind=ip)) then
                stop 'BUG: send_size should have equal length to sendto_image_index'
            end if
        end if

        !
        ! Here, we optionally re-order the send metadata
        ! (i.e. change the order that the send arrays are packed). This
        ! can allow sends to occur 'in groups' which may potentially have
        ! speed benefits
        !
        if(reorder_send_data_by_image .and. size(send_size, kind=ip) > 0 .and. num_images_local > 1) then
        
            ! Get the index of the send-to data
            n1 = size(send_size, kind=ip)
            allocate(send_data_order(n1), send_data_sort_criterion(n1))
            send_data_order = [(i, i=1, n1)]
            ! Make the order so that images above the current image are near the start 
            send_data_sort_criterion = modulo((sendto_image_index - this_image_local), int(num_images_local, c_int))
          
            !print*, 'n11: ', n1 
            !print*, 'send_data_order1: ', send_data_order 
            !print*, 'send_data_sort_criterion1: ', send_data_sort_criterion
            !print*, 'sendto_image_index1: ', sendto_image_index
            !print*, 'send_size1: ', send_size
            !print*, 'send_start_index1: ', send_start_index
          
            if(maxval(send_data_sort_criterion) /= minval(send_data_sort_criterion)) then 
                call sort_index(send_data_order, send_data_sort_criterion, n1)

                ! Reorder the key data
                ! -- sendto_image_index
                ! -- send_size
                ! -- send_buffer_label
                ! -- send_start_index
                sendto_image_index = sendto_image_index(send_data_order)
                send_size = send_size(send_data_order)
                send_buffer_label = send_buffer_label(send_data_order)
                send_start_index(1) = 1
                do i = 2, n1
                    send_start_index(i) = send_start_index(i-1) + send_size(i-1)
                end do
            end if

            !print*, 'n1: ', n1
            !print*, 'send_data_order: ', send_data_order
            !print*, 'sendto_image_index: ', sendto_image_index
            !print*, 'send_size: ', send_size
            !print*, 'send_start_index: ', send_start_index

            ! Clean up
            deallocate(send_data_order, send_data_sort_criterion)
        end if

        !
        ! Allocate the receive buffer to have the maximum required size
        !
        ! First we need to find how big it should be
        desired_size_local = 0
        do i = 1, num_images_local
            if(allocated(send_size)) then
                desired_size = sum(send_size, mask=(sendto_image_index == i))
            else
                desired_size = 0
            end if
#ifdef COARRAY
            call co_sum(desired_size)
#endif
            if(i == this_image_local) then
                desired_size_local = desired_size
            end if
        end do
        desired_size = desired_size_local

        ! Make the recv_buffer -- used for communication every timestep
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call co_max(desired_size)
        allocate(recv_buffer(desired_size))
        call sync_all_generic
#elif defined(COARRAY)
        call co_max(desired_size)
        allocate(recv_buffer(desired_size)[*])
#else
        allocate(recv_buffer(desired_size))
#endif
        ! Fill with a value which is suggestive of problems, in case of
        ! out-of-bounds mistakes
        recv_buffer = -HUGE(1.0_dp)

        !
        ! Next find out how many communications occur. We do this by
        ! broadcasting the send_size and sendto_image_index for each image, then
        ! adding the information to the recv_size / recv_start_index /
        ! recvfrom_image_index
        !
        if(allocated(send_size)) then
            desired_size = size(send_size, kind=ip)
        else
            desired_size = 0
        end if

        ! Make 'work_coarray' which is used to send/receive data between images
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call co_max(desired_size)
        allocate( work_coarray(desired_size, 3))
        allocate(wc3_all(size(work_coarray(:,3), kind=ip), num_images_local))
        call sync_all_generic
#elif defined(COARRAY)
        call co_max(desired_size)
        allocate( work_coarray(desired_size, 3)[*] )
#else
        allocate( work_coarray(desired_size, 3))
#endif

        ! Make real64_coarary to communicate the send_buffer_label (can send length
        ! p2p_id_len characters using 'transfer')
        !n = (p2p_id_len ) / real64 + 1
        if(modulo(storage_size(charlabel), storage_size(real(1.0, real64))) == 0) then
            n = storage_size(charlabel)/storage_size(real(1.0, real64))
        else
            write(log_output_unit, *) "storage_size(charlabel) cannot be evenly divided into real64's "
            flush(log_output_unit)
            call local_stop
        end if
        ! Determine number of empty characters ' ' required to fill a real64
        if(modulo(storage_size(real(1.0, real64)), storage_size(" ")) == 0) then
            n2 = storage_size(real(1.0, real64))/storage_size(" ")
        else
            write(log_output_unit, *) "storage_size(real 64) cannot be evenly divided into empty space characters "
            flush(log_output_unit)
            call local_stop
        end if

#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        allocate( real64_coarray(n, desired_size))
        call sync_all_generic
#elif defined(COARRAY) 
        allocate( real64_coarray(n, desired_size)[*] )
#else
        allocate( real64_coarray(n, desired_size))
#endif

        do i = 1, num_images_local
            ! Broadcast the send metadata from image i to all images
            if( this_image_local == i) then
                work_coarray = 0 
                real64_coarray = transfer(repeat(" ", n2), real(1.0, real64))

                if(allocated(send_size)) then
                    n = size(sendto_image_index, kind=ip)
                    work_coarray(1:n, 1) = sendto_image_index
                    work_coarray(1:n, 2) = send_size
                    ! To broadcast the send_buffer_label, convert to a real
                    do j = 1, n
                        real64_coarray(:, j) = transfer(send_buffer_label(j), &
                            real64_coarray(:,1), size(real64_coarray(:,1), kind=ip))
                    end do
                else
                    n = 0
                end if
            end if

#if defined(COARRAY) 
            call co_broadcast(work_coarray, source_image = i)
            call co_broadcast(real64_coarray, source_image = i)
#endif

            do j = 1, size(work_coarray(:,1), kind=ip)        
                ! If the current image receives data from image i, note that
                if(work_coarray(j,1) == this_image_local) then
   
                    ! Record the receive-from image index
                    if(allocated(recvfrom_image_index)) then
                        recvfrom_image_index = [recvfrom_image_index, i ]
                    else
                        recvfrom_image_index = [i]
                    end if 

                    ! updated the 'linked images' with image i
                    if(allocated(linked_p2p_images)) then
                        if(.not. any(linked_p2p_images == i)) then
                            linked_p2p_images = [linked_p2p_images, i]
                        end if
                    else
                        linked_p2p_images = [i] 
                    end if

                    if(allocated(recv_size)) then
                        recv_size = [recv_size, work_coarray(j, 2) * 1_ip ]
                    else
                        recv_size = [work_coarray(j,2) * 1_ip]
                    end if

                    if(allocated(recv_buffer_label)) then
                        recv_buffer_label = [recv_buffer_label, &
                            transfer(real64_coarray(:, j), charlabel)]
                    else
                        recv_buffer_label = [transfer(real64_coarray(:, j), charlabel)]
                    end if

                    ! Compute the recv_start_index
                    if(allocated(recv_start_index)) then
                        n = size(recv_start_index, kind=ip)
                        recv_start_index = [recv_start_index, &
                            ! Sum of previous start index + previous recv_size
                            recv_start_index(n) + recv_size(n)]
                    else
                        recv_start_index = [1_ip]
                    end if

                    ! Put the recv_start_index into the work_coarray, so
                    ! we can send it back to the recvfrom_image_index
                    n = size(recv_start_index, kind=ip)
                    work_coarray(j,3) = recv_start_index(n)
                end if
            end do
#ifdef COARRAY  
            call sync_all_generic 
#endif


            ! Finally, broadcast the recv_start_indices back to the sendto_image_start_index 
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            ! Broadcast the 3rd column of work_coarray to image i
            call mpi_gather(work_coarray(:,3), size(work_coarray(:,3)), mympi_int, wc3_all, &
                size(work_coarray(:,3)), mympi_int, i-1, MPI_COMM_WORLD, ierr)
#endif
            if(this_image_local == i) then
                if(allocated(send_size)) then
                    allocate(sendto_start_index(size(send_size, kind=ip)))
                    if(size(send_size, kind=ip) > 0) then
                        do j = 1, size(send_size, kind=ip)
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS) 
                            sendto_start_index(j) = wc3_all(j,sendto_image_index(j))
#elif defined(COARRAY)
                            sendto_start_index(j) = work_coarray(j,3)[sendto_image_index(j)]
#else
                            sendto_start_index(j) = work_coarray(j,3)
#endif
                        end do
                    end if
                end if
            end if
#ifdef COARRAY
            call sync_all_generic
#endif
        end do

        deallocate(work_coarray)
        deallocate(real64_coarray)
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS) 
        call sync_all_generic
        deallocate(wc3_all)
#endif


        ! Check that send_buffer_label does not have repeated values
        if(allocated(send_buffer_label)) then
            do i = 1, (size(send_buffer_label, kind=ip) - 1)
                do j = i+1, size(send_buffer_label, kind=ip)
                    if(send_buffer_label(i) == send_buffer_label(j)) then
                        write(log_output_unit,*) 'Error: repeated send_buffer_labels: ', &
                            send_buffer_label(i)
                    end if
                end do
            end do
        end if

        ! Check that recv_buffer_label does not have repeated values
        if(allocated(recv_buffer_label)) then
            do i = 1, (size(recv_buffer_label, kind=ip) - 1)
                do j = i+1, size(recv_buffer_label, kind=ip)
                    if(recv_buffer_label(i) == recv_buffer_label(j)) then
                        write(log_output_unit,*) 'Error: repeated recv_buffer_labels: ', &
                            recv_buffer_label(i)
                    end if
                end do
            end do
        end if

        if(allocated(recv_size)) then
            if(size(recv_start_index, kind=ip) /= size(recv_size, kind=ip)) then
                write(log_output_unit,*) 'Error:  size(recv_start_index) /= size(recv_size) ', &
                    size(recv_start_index, kind=ip), size(recv_size, kind=ip)
                call local_stop 
            end if
        end if

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS

        ! Allocate key input variables for mpi_comms
        if(.not. reorder_send_data_by_image) then
            write(log_output_unit, *) 'ERROR: If COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined,'
            write(log_output_unit, *) 'then must have reorder_send_data_by_image=.true, so that'
            write(log_output_unit, *) 'there is only one send from image A to image B'
            call local_stop
        end if

        if(mpi_timestep_loop_use_alltoallv) then
            ! Define variables needed for mpi_alltoallv

            allocate(mympi_send_counts(num_images_local), mympi_send_displacements(num_images_local))
            allocate(mympi_recv_counts(num_images_local), mympi_recv_displacements(num_images_local))
            !print*, 'BEGIN MPI STUFF'

            mympi_send_counts = 0
            mympi_recv_counts = 0
            mympi_send_displacements = 0
            mympi_recv_displacements = 0

            ! Set the values for MPI_alltoallv, based on the metadata above
            do i = 1, num_images_local
                if(allocated(sendto_image_index)) then
                    if(any(sendto_image_index == i)) then
                        ! Set the send information 

                        ! Find the first index with image == i
                        do j = 1, size(sendto_image_index, kind=ip)
                            if(sendto_image_index(j) == i) then
                                n = j 
                                exit
                            end if
                        end do
                        ! If the start index is 'p', the MPI displacement is 'p-1'
                        mympi_send_displacements(i) = send_start_index(n) - 1
                        if(size(send_size, kind=ip) > 0) then
                            mympi_send_counts(i) = sum(send_size, mask=(sendto_image_index == i))
                        else
                            mympi_send_counts(i) = 0 
                        end if
                    end if
                end if
                ! As above, for recv information
                if(allocated(recvfrom_image_index)) then
                    if(any(recvfrom_image_index == i)) then
                        ! Set the recv information 

                        ! Find the first index with image == i
                        do j = 1, size(recvfrom_image_index, kind=ip)
                            if(recvfrom_image_index(j) == i) then
                                n = j 
                                exit
                            end if
                        end do
                        ! If the start index is 'p', the MPI displacement is 'p-1'
                        mympi_recv_displacements(i) = recv_start_index(n) - 1
                        if(size(recv_size, kind=ip) > 0) then
                            mympi_recv_counts(i) = sum(recv_size, mask=(recvfrom_image_index == i))
                        else
                            mympi_recv_counts(i) = 0
                        end if
                    end if
                end if
            end do
        else
            ! Define variables for mpi_isend/irecv - an alternative to the alltoallv approach
            allocate(mpi_recv_requests(size(recv_start_index, kind=ip)))
            mpi_recv_requests = MPI_REQUEST_NULL
            allocate(mpi_send_requests(size(send_start_index, kind=ip)))
            mpi_send_requests = MPI_REQUEST_NULL
        end if
#endif

    end subroutine

    subroutine deallocate_p2p_comms
        !!
        !! Clear all allocatable arrays in this module
        !!

        !deallocate(send_buffer, send_start_index, send_size, sendto_image_index, &
        !    sendto_start_index, send_buffer_label, recv_buffer, recv_start_index, &
        !    recv_size, recvfrom_image_index, recv_buffer_label, linked_p2p_images)

        if(allocated(send_buffer)) deallocate(send_buffer)
        if(allocated(send_start_index)) deallocate(send_start_index)
        if(allocated(send_size)) deallocate(send_size)
        if(allocated(sendto_image_index)) deallocate(sendto_image_index)
        if(allocated(sendto_start_index)) deallocate(sendto_start_index)
        if(allocated(send_buffer_label)) deallocate(send_buffer_label)

        if(allocated(recv_buffer)) deallocate(recv_buffer)
        if(allocated(recv_start_index)) deallocate(recv_start_index)
        if(allocated(recv_size)) deallocate(recv_size)
        if(allocated(recvfrom_image_index)) deallocate(recvfrom_image_index)
        if(allocated(recv_buffer_label)) deallocate(recv_buffer_label)

        if(allocated(linked_p2p_images)) deallocate(linked_p2p_images)

        if(allocated(work_coarray)) deallocate(work_coarray)
        if(allocated(real64_coarray)) deallocate(real64_coarray)

        have_allocated_p2p_comms = .FALSE.

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
        if(allocated(mympi_send_counts)) deallocate(mympi_send_counts)
        if(allocated(mympi_recv_counts)) deallocate(mympi_recv_counts)
        if(allocated(mympi_send_displacements)) deallocate(mympi_send_displacements)
        if(allocated(mympi_recv_displacements)) deallocate(mympi_recv_displacements)
        if(allocated(mpi_recv_requests)) deallocate(mpi_recv_requests)
        if(allocated(mpi_send_requests)) deallocate(mpi_send_requests)
        call sync_all_generic
#endif

    end subroutine

    subroutine find_send_buffer_label_index(buffer_label, buffer_label_int)
        !!
        !! Find the index of send_buffer_label which matches 'buffer_label'.
        !!
        !! Uses a naive scan of all buffer_labels, but should be efficient enough
        !! for typical cases with a small number of 'buffer_label' values.
        !!
        character(len=p2p_id_len), intent(in) :: buffer_label !! The buffer label to match
        integer(ip), intent(out) :: buffer_label_int !! The corresponding index in send_buffer_label
       
        integer(ip) :: i 

        !write(log_output_unit,*) 'DEBUG: ', send_buffer_label

        ! Find the integer index in the metadata corresponding to
        ! buffer_label, by finding a match with the send_buffer_label's
        buffer_label_int = -1
        do i = 1, size(send_buffer_label, kind=ip)
            if(buffer_label == send_buffer_label(i)) then
                buffer_label_int = i 
                exit
            end if
        end do

        if(buffer_label_int < 1) then
            write(log_output_unit,*) 'unrecognized buffer_label ', buffer_label
            call local_stop
        end if

    end subroutine

    subroutine put_on_recv_buffer(buffer_label_int, si, ei)
        !! 
        !! Put the send_buffer(si:ei) on the recv_buffer associated
        !! with sendto_start_index(buffer_label_int)
        !!
        !! Convenience routine for something we often have to do in
        !! send_ routines for all ranks
        !!
        integer(ip), intent(in) :: buffer_label_int, si, ei

        integer(ip):: i, recv_start_index_local, recv_end_index, recv_image

        i = buffer_label_int
        recv_start_index_local = sendto_start_index(i)
        recv_end_index = recv_start_index_local + send_size(i) - 1
        recv_image = sendto_image_index(i)

#if defined(COARRAY) && defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        write(log_output_unit,*) 'ERROR: Cannot call put_on_recv_buffer with COARRAY_USE_MPI_FOR_INTENSIVE_COMMS'
        flush(log_output_unit)
        call local_stop()
#endif
        ! put communication
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
            send_buffer(si:ei)
#else
        recv_buffer(recv_start_index_local:recv_end_index) = &
            send_buffer(si:ei)
#endif

    end subroutine

    subroutine send_to_p2p_comms_rank1(send_array, buffer_label, &
        put_in_recv_buffer)
        !!
        !! put array 'send_array' in the send buffer, in preparation for
        !! sending to the recv buffer
        !!
        ! @param send_array a rank 1 array with kind dp
        ! @param buffer_label a character string (same as was used to define
        !    the communication in other routines)
        ! @param put_in_recv_buffer logical, optional (default .TRUE.). If .TRUE.
        !    then we do the parallel 'put' here, otherwise we do not, and
        !    'communicate_p2p' must be called later
        ! 

        real(dp), intent(in) :: send_array(:) !! a rank 1 array with kind dp
        character(len=p2p_id_len), intent(in) :: buffer_label 
            !! a character string (same as was used to define the communication in other routines)
        logical, intent(in), optional :: put_in_recv_buffer
            !! If .TRUE. then we do the parallel 'put' here, otherwise we do not, and
            !! 'communicate_p2p' must be called later

        ! Generic code 
        include 'point2point_include_send_p2p.f90'
    end subroutine

    ! rank2 version of above
    subroutine send_to_p2p_comms_rank2(send_array, buffer_label, &
        put_in_recv_buffer)

        real(dp), intent(in) :: send_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code 
        include 'point2point_include_send_p2p.f90'
    end subroutine

    ! rank3 version of above
    subroutine send_to_p2p_comms_rank3(send_array, buffer_label, &
        put_in_recv_buffer)

        real(dp), intent(in) :: send_array(:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code 
        include 'point2point_include_send_p2p.f90'
    end subroutine

    ! rank4 version of above
    subroutine send_to_p2p_comms_rank4(send_array, buffer_label, &
        put_in_recv_buffer)

        real(dp), intent(in) :: send_array(:,:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code 
        include 'point2point_include_send_p2p.f90'
    end subroutine

#ifndef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS

    subroutine communicate_p2p
        !!
        !! Use coarray parallel communication to send ALL data to the recv buffer.
        !!
        !! This uses a 'put' model of communication.
        !!
        !! Note this routine is not needed if the puts are done with send_to_p2p_comms.
        !! But, this version is optimized by merging the puts (e.g. merging sends of all
        !! send arrays corresponding to the same recv_image). This can be switched off
        !! with the parameter reorder_send_data_by_image
        !!

        integer(ocaIP) :: i, start_index, end_index, recv_image, &
            recv_start_index_local, recv_end_index

        ! If there is nothing to send, exit
        if(.not. allocated(send_size)) then
            return
        end if

        if(.not. reorder_send_data_by_image) then
            ! Send ALL of the send buffers to the recv buffers, one by one
            do i = 1, size(send_size, kind=ip)
                start_index = send_start_index(i)
                end_index = start_index + send_size(i) - 1
                recv_start_index_local = sendto_start_index(i)
                recv_end_index = recv_start_index_local + send_size(i) - 1
                recv_image = sendto_image_index(i)
                ! put communication
#ifdef COARRAY
                recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                    send_buffer(start_index:end_index)
#else
                recv_buffer(recv_start_index_local:recv_end_index) = &
                    send_buffer(start_index:end_index)
#endif
            end do
        else

            ! Here, we have sorted the sends so that sends to each image are clustered together
            ! Thus we can make less sends

            ! This initialisation is a trick to get the loop to work
            start_index = 1
            end_index = start_index - 1
            recv_image = sendto_image_index(1)
            recv_start_index_local = sendto_start_index(1)
            recv_end_index = recv_start_index_local  - 1

            do i = 1, size(send_size, kind=ip)

                if(sendto_image_index(i) == recv_image) then
                    ! The 'previous' sendto image is the same as the current image
                    ! For now we just need to update the end indices
                    end_index = end_index + send_size(i)
                    recv_end_index = recv_end_index + send_size(i)
                end if

                ! If the sendto_image_index has changed, then send the data accumulated previously
                if(sendto_image_index(i) /= recv_image) then
                    ! Send the previous data
#ifdef COARRAY
                    recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                        send_buffer(start_index:end_index)
#else
                    recv_buffer(recv_start_index_local:recv_end_index) = &
                        send_buffer(start_index:end_index)
#endif

                    ! Redefine the start_index, and recv_image
                    start_index = end_index + 1
                    end_index = start_index + send_size(i) - 1
                    recv_image = sendto_image_index(i)
                    recv_start_index_local = sendto_start_index(i)
                    recv_end_index = recv_start_index_local + send_size(i) - 1
                end if

                ! If we are on the final send, we definitely need to send
                if(i == size(send_size, kind=ip)) then
#ifdef COARRAY
                    recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                        send_buffer(start_index:end_index)
#else
                    recv_buffer(recv_start_index_local:recv_end_index) = &
                        send_buffer(start_index:end_index)
#endif
                end if

            end do

        end if
    end subroutine
#endif

    !
    ! Note there are two different versions of communicate_p2p, which 
    ! are applied depending on whether COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined 
    !

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    subroutine communicate_p2p
        !!
        !!  If COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined, then use MPI for the most intensive point2point comms routine
        !!  Note this requires that the appropriate preprocessor variable has been defined
        !!
   
        integer:: mympi_my_ierr
        integer :: i, mpi_count, mpi_source, mpi_dest, mpi_ierr, mpi_tag

        if (.not. reorder_send_data_by_image) then
           write(log_output_unit, *) 'ERROR: communicate_p2p with MPI assumes reorder_send_data_by_image=.true'
           call local_stop
        end if
   
        if(mpi_timestep_loop_use_alltoallv) then
            !
            ! Version using mpi_alltoallv
            !

            !NOTE -- send_size, recv_size, send_start_index, recv_start_index, ierr, need to be of type integer
            !        also send_size, etc need to be the of size num_images()
            call mpi_alltoallv(send_buffer , mympi_send_counts, mympi_send_displacements, mympi_dp, &
                               recv_buffer , mympi_recv_counts, mympi_recv_displacements, mympi_dp, &
                               MPI_COMM_WORLD, mympi_my_ierr)
      
            if(mympi_my_ierr /= 0) then 
                write(log_output_unit,*) 'FAIL, MPI_alltoallv error: mympi_my_ierr= ', mympi_my_ierr
                write(log_output_unit,*) __LINE__,&
                    __FILE__
                call local_stop
            end if
        else
            !
            ! Version using repeated isend/irecv calls. 
            !

            if(allocated(recv_size)) then 
                ! Open up receives         
                do i = 1, size(recv_start_index, kind=ip)
                    mpi_count = recv_size(i)
                    mpi_source = recvfrom_image_index(i) - 1
                    
                    call mpi_irecv(recv_buffer(recv_start_index(i):(recv_start_index(i) + recv_size(i) - 1)), &
                        mpi_count, mympi_dp, mpi_source, MPI_ANY_TAG, MPI_COMM_WORLD, mpi_recv_requests(i), mpi_ierr)

                end do
            end if 
           
            if(allocated(send_size)) then 
                ! Do sends
                do i = 1, size(send_start_index, kind=ip)
                    mpi_count = send_size(i)
                    mpi_dest = sendto_image_index(i) - 1
                    mpi_tag = i
                    call mpi_isend(send_buffer(send_start_index(i):(send_start_index(i) + send_size(i) - 1)), &
                        mpi_count, mympi_dp, mpi_dest, mpi_tag, MPI_COMM_WORLD, mpi_send_requests(i), mpi_ierr)
                end do
            end if

            ! Ensure completion -- more strategic location of these calls would be possible
            if(allocated(send_size)) then
                mpi_count = size(send_start_index, kind=ip)
                call mpi_waitall(mpi_count, mpi_send_requests, MPI_STATUSES_IGNORE, mpi_ierr)
            end if

            if(allocated(recv_size)) then
                mpi_count = size(mpi_recv_requests, kind=ip)
                call mpi_waitall(mpi_count, mpi_recv_requests, MPI_STATUSES_IGNORE, mpi_ierr)
            end if
        end if 

    end subroutine
#endif

    !
    ! Copy data from the recv buffer to 'recv_array'
    ! 
    ! @param recv_array rank 1 array of kind dp, into which we wish to copy data
    !     associated with the label
    ! @param buffer_label character string giving a label to the communication
    !     (same as mentioned earlier)
    pure subroutine recv_from_p2p_comms_rank1(recv_array, buffer_label)

        real(dp), intent(inout) :: recv_array(:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank2 version of above
    pure subroutine recv_from_p2p_comms_rank2(recv_array, buffer_label)

        real(dp), intent(inout) :: recv_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank3 version of above
    pure subroutine recv_from_p2p_comms_rank3(recv_array, buffer_label)

        real(dp), intent(inout) :: recv_array(:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank4 version of above
    pure subroutine recv_from_p2p_comms_rank4(recv_array, buffer_label)

        real(dp), intent(inout) :: recv_array(:,:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    elemental function integer_to_id(myint) result(mychar)
        !! 
        !! Convert an integer to a 2-character ID. This can be useful in generating
        !! buffer_labels. The labels seem unique for intgers ranging from 1 to 60000,
        !! beyond that at some stage repeats occur
        !!
        integer(ip), intent(in) :: myint
        character(len=2) :: mychar

        mychar = achar(myint) // achar(myint/256)

    end function

    function size_of_send_recv_buffers() result(mysize)
        !! 
        !! Useful to know how big the send/recv buffers are
        !!
        integer(ip) :: mysize
        mysize = 0
        if(allocated(send_buffer)) then
            mysize = mysize + size(send_buffer, kind=ip)*real_bytes
        end if
        if(allocated(recv_buffer)) then
            mysize = mysize + size(recv_buffer, kind=ip)*real_bytes
        end if
    end function

    subroutine print_p2p_comms(verbose)
        !!
        !! Utility printing routine, useful for debugging
        !!
        logical, optional, intent(in) :: verbose

        integer :: i, si, ei
        logical :: verbose_in

        if(present(verbose)) then
            verbose_in = verbose
        else
            verbose_in = .TRUE.
        end if
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        critical 
#endif
        write(log_output_unit,*) ' '
        write(log_output_unit,*) '########################################'
        write(log_output_unit,*) 'image: ', this_image_local, '/', num_images_local
        write(log_output_unit,*) '########################################'
        write(log_output_unit,*) 'send_start_index: ', send_start_index
        write(log_output_unit,*) 'send_size: ', send_size
        write(log_output_unit,*) 'sendto_image_index: ', sendto_image_index
        write(log_output_unit,*) 'sendto_start_index: ', sendto_start_index
        write(log_output_unit,*) 'size(send_buffer): ', size(send_buffer, kind=ip)
        do i = 1, size(send_start_index, kind=ip)
            si = send_start_index(i)
            ei = si + send_size(i) - 1
            write(log_output_unit,*) '    ----'
            write(log_output_unit,*) '    send_buffer_label: ', trim(send_buffer_label(i))
            write(log_output_unit,*) '        send_array number: ', i
            write(log_output_unit,*) '        size:', send_size(i)
            write(log_output_unit,*) '        sendto_image: ', sendto_image_index(i)
            write(log_output_unit,*) '        sendto_start_index: ', sendto_start_index(i)
            if(verbose_in) write(log_output_unit,*) '        send_buffer: ', send_buffer(si:ei)
        end do

        !write(log_output_unit,*) 'recv_buffer: ', recv_buffer
        write(log_output_unit,*) 'size(recv_buffer): ', size(recv_buffer, kind=ip)
        write(log_output_unit,*) 'recv_start_index: ', recv_start_index
        write(log_output_unit,*) 'recv_size: ', recv_size
        write(log_output_unit,*) 'recvfrom_image_index: ', recvfrom_image_index
        do i = 1, size(recv_size, kind=ip)
            si = recv_start_index(i)
            ei = si + recv_size(i) - 1
            write(log_output_unit,*) '    ----'
            write(log_output_unit,*) '    recv_buffer_label: ', trim(recv_buffer_label(i))
            write(log_output_unit,*) '        recv_array number: ', i
            write(log_output_unit,*) '        size: ', recv_size(i)
            write(log_output_unit,*) '        recvfrom_image_index: ', recvfrom_image_index(i)
            if(verbose_in) write(log_output_unit,*) '        recv_buffer: ', recv_buffer(si:ei)
        end do

#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        end critical
#endif
    end subroutine

    subroutine local_stop
        !!
        !! Call either 'error stop' or just 'stop'. FIXME: Is this needed?
        !!
#ifdef COARRAY
        error stop  
#else
        stop
#endif 
    end subroutine

    subroutine test_coarray_point2point_comms_mod
        !!
        !! Unit tests. Do some point-2-point communication of arrays, and check that results are
        !! as expected
        !!

        implicit none

        real(dp), allocatable :: x_send(:,:,:), y_send(:), z_send(:,:)
        real(dp), allocatable :: x_recv(:,:,:), y_recv(:), z_recv(:,:)
        integer(ip) :: k, ti, ni, sendto_image, y_recv_size
        character(len=p2p_id_len) :: x_label, y_label, z_label
        real(dp) :: expected_xrecv, expected_yrecv, expected_zrecv

        !character(len=2) :: local_chars(60000)

        logical :: local_puts

#ifdef COARRAY
        ti = this_image2()
        ni = num_images2()
#else
        ti = 1
        ni = 1
#endif

        !
        ! Make some data to send/recv. Array ranks can vary, and the dimensions
        ! do not have to be consistent among images
        !
        allocate(x_send(1000, 2, 1), y_send(2 + ti), z_send(2000,1))
        allocate(x_recv(1000, 2, 1),                 z_recv(2000,1))
       
        ! Get the size of y_recv 
        y_recv_size = ti + 1 
        if(y_recv_size > ni) y_recv_size = 1
        y_recv_size = y_recv_size + 2
        allocate(y_recv(y_recv_size))

        ! Make up some data
        x_send = ti * 1.0_dp
        y_send = ti * (-1.0_dp)
        z_send = 55.5_dp ![55.5_dp, 44.4_dp]

        ! Decide which image to send 'x' to
        sendto_image = ti + 1
        if(sendto_image > ni) sendto_image = 1
        x_label = 'x_comms' 
        call include_in_p2p_send_buffer(x_send, buffer_label=x_label, &
            receiver_image=sendto_image)

        ! Decide which image to send 'y' to
        sendto_image = ti - 1
        if(sendto_image == 0) sendto_image = ni
        y_label = 'y_comms'// repeat('*', p2p_id_len-7) ! Check we can 'fill' the character length
        call include_in_p2p_send_buffer(y_send, buffer_label=y_label, &
            receiver_image=sendto_image)

        ! Only send 'z' from image 2 to image 1
        z_label = 'z_comms'
        if(mod(ti - 2, ni) == 0) then
            call include_in_p2p_send_buffer(z_send, buffer_label=z_label, &
                receiver_image = 1_ip)
        end if

        ! Allocate send/recv buffers. From now on, we cannot further call
        ! include_in_p2p_send_buffer
        call allocate_p2p_comms

        ! Run the send's twice -- once doing all puts inside 'communicate_p2p'
        do k = 1, 2
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS) 
            local_puts = .FALSE.
#else
            ! Try doing the 'put' communication in different ways
            local_puts = (k == 1)
#endif


            ! Modify the data we send
            x_send = x_send + 1 
            y_send = y_send + 1
            z_send = z_send + 1

            ! put the arrays in the send buffer
            call send_to_p2p_comms(x_send, buffer_label=x_label, &
                put_in_recv_buffer=local_puts)
            call send_to_p2p_comms(y_send, buffer_label=y_label, &
                put_in_recv_buffer=local_puts)
            if(mod(ti - 2, ni) == 0) then
                call send_to_p2p_comms(z_send, buffer_label=z_label, &
                    put_in_recv_buffer=local_puts)
            end if

            if(.not.local_puts) then
                ! Do the parallel put
                call communicate_p2p()
            end if
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            sync images(linked_p2p_images)
#endif


            ! Copy the sent data to receive buffers
            call recv_from_p2p_comms(x_recv, buffer_label=x_label)
            call recv_from_p2p_comms(y_recv, buffer_label=y_label)
            if((ti == 1)) then
                call recv_from_p2p_comms(z_recv, buffer_label=z_label)
            end if

            ! sync to prevent the 'k' loop moving ahead before we have received
            ! data
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            sync images(linked_p2p_images)
#endif

            ! test the x data
            expected_xrecv = ti - 1.0_dp
            if(expected_xrecv == 0.0_dp) then
                expected_xrecv = ni*1.0_dp
            end if
            expected_xrecv = expected_xrecv + k
            
            if(all(x_recv == expected_xrecv)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,&
                    __FILE__
                call local_stop
            end if

            ! test the y data
            expected_yrecv = ti + 1.0_dp
            if(expected_yrecv > ni*1.0_dp) then
                expected_yrecv = 1.0_dp
            end if
            expected_yrecv = -1.0_dp * expected_yrecv
            expected_yrecv = expected_yrecv + k
            if(all(y_recv == expected_yrecv)) then
                write(log_output_unit,*) 'PASS'
            else
                write(log_output_unit,*) 'FAIL', __LINE__,&
                    __FILE__
                call local_stop
            end if

            ! test the 'z' data
            if((ti == 1).and.(ni > 1)) then

                expected_zrecv = 55.5_dp + k
                if(all(z_recv == expected_zrecv)) then
                    write(log_output_unit,*) 'PASS'
                else
                    write(log_output_unit,*) 'FAIL', __LINE__,&
                        __FILE__
                    call local_stop
                end if
                
            end if

            ! Print everything
            !call print_p2p_comms()
        end do

        ! Clean up
        call deallocate_p2p_comms

        if(allocated(linked_p2p_images)) then
            write(log_output_unit,*) 'FAIL: p2p deallocation did not work', __LINE__,&
                __FILE__
        else
            write(log_output_unit,*) 'PASS'
        end if

    end subroutine

end module
