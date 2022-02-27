module coarray_point2point_comms_mod
!! Module for doing coarray 'point-to-point' communications of real
!! arrays (rank from 1 to 4), in a 'multiple-program, multiple data' style.
!!
!! Send an arbitrary (and varying) number of arrays between any pairs of images.
!! The array dimensions do not have to be consistent among images.
!! 
!! This module was originally written to use coarrays, but later extended to 
!! alternatively use MPI, either via Isend/Irecv, or via All2Allv. MANY OF
!! THE COMMENTS ONLY REFER TO THE COARRAY CASE, but the concepts translate to MPI.
!! In particular, 'image i' for corrarys corresponds to 'rank i-1' in MPI.
!!
!! Background
!! ------------
!!
!! Coarrays by definition must be the same size on each image (MPI does not have this restriction). 
!! It is not always straightforward to apply coarrays to problems where we need to:
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
!! with a simple interface.
!!
!! With coarrays, memory will necessarily be wasted if the amount of data to
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
!!     type(p2p_comms_type) :: p2p
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
!!     call include_in_p2p_send_buffer(p2p, y_send, buffer_label='y_comms1', receiver_image=this_image() + 1_ip)
!!     if(this_image() == 1) then
!!         ! If we are on image 1, then associate the label 'x_comms1' with
!!         ! communication of the 'x_' variables to image 2
!!         call include_in_p2p_send_buffer(p2p, x_send, buffer_label='x_comms1', receiver_image=2_ip)
!!     end if
!!
!!     ! Once we have identified which arrays are communicated using
!!     ! 'include_in_p2p_send_buffer', we can allocate communication buffers.
!!     ! At this stage, we should not call 'include_in_p2p_send_buffer' anymore.
!!     ! All images should call this subroutine, even if they are not sending or receiving,
!!     ! because the coarray must be allocated on all images at once.
!!     call allocate_p2p_comms(p2p)
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
!!         call send_to_p2p_comms(p2p, y_send, buffer_label='y_comms1')
!!         if(this_image() == 1) then
!!             ! send x_send to image 2
!!             call send_to_p2p_comms(p2p, x_send, buffer_label='x_comms1')
!!         end if
!!
!!         ! Make sure images we receive from have sent their data.
!!         ! This needs to be called on all communicating images
!!         if(size(p2p%linked_p2p_images) > 0) sync images(p2p%linked_p2p_images)
!!         ! p2p%linked_p2p_images is an array containing all image_indexes that
!!         ! we send to or receive from
!!
!!         ! Copy from recv buffer back to computational array
!!         call recv_from_p2p_comms(p2p, y_recv, buffer_label='y_comms1')
!!         if(this_image() == 2) then
!!             call recv_from_p2p_comms(p2p, x_recv, buffer_label='x_comms1')
!!         end if
!!
!!         ! Do something with the received data, update x_, y_, more computations ....
!!
!!     end do

    use global_mod, only: dp, ip, charlen, real_bytes
        !! kinds for double, integer, default character length, and bytes in a real
    use iso_fortran_env, only: real64, int32
        !! Use real64 as a send buffer for the buffer_label character id's
        !! Use int32 for ints that opencoarrays sends [doesn't yet support e.g. int64]
    use reshape_array_mod, only: flatten_array, repack_rank1_array
        !! routines to efficiently convert between rank1 and rankN arrays (n=1,2,3,4)
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use qsort_mod, only: sort_index
    use iso_c_binding, only: c_int
#if defined(COARRAY_PROVIDE_CO_ROUTINES)
    use coarray_intrinsic_alternatives, only: co_broadcast, co_max, co_sum, &
        sync_all_generic, this_image2, num_images2
#elif defined(COARRAY)
    use coarray_intrinsic_alternatives, only: sync_all_generic, this_image2, &
        num_images2
#endif
#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    ! The 'inner loop' communication will use mpi rather than coarrays
    use mpi
#endif
    implicit none


    private

    public :: p2p_comms_type 
        !! Main derived type that holds send/receive metadata
    public :: test_coarray_point2point_comms_mod
        !! Unit test subroutine
    public :: include_in_p2p_send_buffer
        !! Append another point2point send to the data structures
    public :: allocate_p2p_comms
        !! Allocate arrays to do the point2point sends
    public :: deallocate_p2p_comms
        !! Clean up arrays to do the point2point sends
    public :: send_to_p2p_comms, recv_from_p2p_comms
        !! Send and unpack the buffers
    public :: communicate_p2p
        !! Do all sends at once
    public :: print_p2p_comms
        !! Print info about the p2p_comms_type data 
    public :: size_of_send_recv_buffers
        !! Compute the size of the send/receive buffers

    integer, parameter :: ocaIP = int32
        !! Open-coarrays integer precision. At the time of writing 
        !! opencoarrays could not send any integer precision (without effort)

    integer(ocaIP), parameter :: p2p_id_len = charlen
        !! len of character used for buffer_label

    logical, parameter :: reorder_send_data_by_image = .true.
        !! Do we order the send data by image? If true, then 'communicate_p2p' 
        !! can sometimes combine multiple sends in one communication. There 
        !! does not seem to be any reason not to do this.

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
    !
    ! MPI specific vars
    !
#ifdef COARRAY_USE_MPI_ALLTOALLV
    logical, parameter :: mpi_timestep_loop_use_alltoallv = .true.
        !! The inner mpi communication can either use mpi_alltoallv, or mpi_isend/irecv
#else
    logical, parameter :: mpi_timestep_loop_use_alltoallv = .false.
        !! The inner mpi communication can either use mpi_alltoallv, or mpi_isend/irecv
#endif

#ifdef REALFLOAT
    integer, parameter :: mympi_dp = MPI_REAL
#else
    integer, parameter :: mympi_dp = MPI_DOUBLE_PRECISION
#endif
    integer, parameter :: mympi_int = MPI_INTEGER
    integer, parameter :: mympi_real64 = MPI_DOUBLE_PRECISION
    !
    ! End of MPI specific vars
    !
#endif

    integer(ocaIP) :: this_image_local, num_images_local
        !! Use these variables instead of calls to this_image(), num_images().
        !! This generalises the code to work with or without coarrays
        !! (in the 'without' case, we are in serial).

    type :: p2p_comms_type
        !! Main derived type, containing data and metadata to do multiple 
        !! point2point sends/receives

        logical :: have_allocated_p2p_comms = .FALSE.
            !! Record whether we have allocated the data in this type.
        integer(ocaIP), allocatable :: linked_p2p_images(:)
            !! store indices of all images we send to or receive from

        real(dp), allocatable :: send_buffer(:)
            !! Main send buffer
        integer(ip), allocatable :: send_start_index(:)
            !! The 'ith' send_array is sent from:
            !!     p2p%send_buffer(p2p%send_start_index(i) + (0:(p2p%send_size(i) - 1) )  )
        integer(ip), allocatable :: send_size(:)
            !! The 'ith' send_array is sent from:
            !!     p2p%send_buffer(p2p%send_start_index(i) + (0:(p2p%send_size(i) - 1) )  )
        integer(ocaIP), allocatable :: sendto_image_index(:)
            !! The 'ith' send_array is sent to:
            !!     p2p%recv_buffer(p2p%sendto_start_index(i) + (0:(p2p%send_size(i)-1)) )[p2p%sendto_image_index(i)]
        integer(ip), allocatable :: sendto_start_index(:)
            !! The 'ith' send_array is sent to:
            !!     p2p%recv_buffer(p2p%sendto_start_index(i) + (0:(p2p%send_size(i)-1)) )[p2p%sendto_image_index(i)]
        character(len=p2p_id_len), allocatable :: send_buffer_label(:)
            !! Labels help to keep track of matching sends and receives

#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        real(dp), allocatable :: recv_buffer(:)[:]
            !! Main receive buffer -- consider making this a pointer to a section of a 'top-level' recv buffer.
            !! That design could work-around memory wastage (because fortran coarrays are the same size on every image).
#else
        real(dp), allocatable :: recv_buffer(:)
            !! Main receive buffer
#endif
        integer(ip), allocatable :: recv_start_index(:)
            !! The ith receive is placed in p2p%recv_buffer(p2p%recv_start_index(i) + (0:(p2p%recv_size(i)-1)))
        integer(ip), allocatable :: recv_size(:)
            !! The ith receive is placed in p2p%recv_buffer(p2p%recv_start_index(i) + (0:(p2p%recv_size(i)-1)))
        integer(ocaIP), allocatable :: recvfrom_image_index(:)
            !! The 'ith' recv_array comes from some part of:
            !!     p2p%send_buffer on image [p2p%recvfrom_image_index(i)]. 
            !! We don't store the exact slice that it originates from, since 
            !! we use 'put' communication here.
        character(len=p2p_id_len), allocatable :: recv_buffer_label(:)
            !! Labels help to keep track of matching sends and receives

#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        integer(ocaIP), allocatable :: work_coarray(:,:)[:]
            !! We need to communicate various integer arrays when initialising
        real(real64), allocatable  :: real64_coarray(:,:)[:]
            !! Used to communicate a character of up to length p2p_id_len, using
            !! "transfer"
#else
        integer(ocaIP), allocatable :: work_coarray(:,:)
            !! We need to communicate various integer arrays when initialising
        real(real64), allocatable :: real64_coarray(:,:)
            !! Used to communicate a character of up to length p2p_id_len, using
            !! "transfer"
#endif

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
        integer, allocatable:: mympi_send_counts(:), mympi_recv_counts(:)
            !! Variables to use with MPI_IAllToAllv or MPI_Isend/Irecv
        integer, allocatable:: mympi_send_displacements(:), mympi_recv_displacements(:)
            !! Variables to use with MPI_IAllToAllv or MPI_Isend/Irecv
        integer, allocatable :: mpi_recv_requests(:), mpi_send_requests(:)
            !! Variables required for MPI_ISEND and MPI_IRECV.
#endif

    end type

    interface include_in_p2p_send_buffer
        !! Allow 'include_in_p2p_send_buffer' to apply to input arrays with rank from 1 to 4
        module procedure include_in_p2p_send_buffer_rank1, &
            include_in_p2p_send_buffer_rank2, &
            include_in_p2p_send_buffer_rank3, include_in_p2p_send_buffer_rank4
    end interface

    interface send_to_p2p_comms
        !! Allow 'send_to_p2p_comms' to apply to input arrays with rank from 1 to 4
        module procedure send_to_p2p_comms_rank1, send_to_p2p_comms_rank2, &
            send_to_p2p_comms_rank3, send_to_p2p_comms_rank4
    end interface

    interface recv_from_p2p_comms
        !! Allow 'recv_from_p2p_comms' to apply to output arrays with rank from 1 to 4
        module procedure recv_from_p2p_comms_rank1, recv_from_p2p_comms_rank2, &
            recv_from_p2p_comms_rank3, recv_from_p2p_comms_rank4
    end interface

    contains

    subroutine include_in_p2p_send_buffer_generic(p2p, send_array_size, &
            buffer_label, receiver_image)
        !!
        !! Record that p2p should to communicate 'send_array' to image 'receiver_image'.
        !!
        !! Associate the communication with a string 'buffer_label', which can also
        !! be used to receive the sent data.
        !!
        !! Note the actual allocations happen later (once we know all the arrays
        !! we'd like to send).

        type(p2p_comms_type), intent(inout) :: p2p
            !! The p2p_comms_type that is to do this communication
        integer(ip), intent(in) :: send_array_size 
            !! A one dimensional real array with kind dp
        character(len=p2p_id_len), intent(in) :: buffer_label
            !! A character string used to identify the communication data in 
            !! both the send and recv buffers
        integer(ip), intent(in) :: receiver_image 
            !! image index which will receive the data

        if ( allocated(p2p%send_buffer) ) then
            write(log_output_unit,*) 'p2p%send_buffer is already allocated. ', &
                'Cannot create more p2p%send_buffer space after allocation'
            call generic_stop
        end if

        !write(log_output_unit,*) '    DEBUG:', trim(buffer_label), ' ', send_array_size, receiver_image, &
        !    allocated(p2p%send_size), size(p2p%send_size)

        !
        ! Append space for the send_array to the metadata describing the
        ! p2p%send_buffer
        !
        if(allocated(p2p%send_size)) then
            ! Array will go in:
            !     p2p%send_buffer( p2p%send_start_index(buffer_label) + (0:(p2p%send_size(send_array_coomms_id) - 1)))
            p2p%send_start_index = [&
                p2p%send_start_index, &
                p2p%send_start_index(size(p2p%send_start_index, kind=ip)) + &
                    p2p%send_size(size(p2p%send_size, kind=ip))]
            p2p%send_size = [p2p%send_size, int(send_array_size, ip)]
            p2p%sendto_image_index = [p2p%sendto_image_index, &
                                      int(receiver_image, ocaIP)]
            if(.not. any(p2p%linked_p2p_images == receiver_image)) then
                p2p%linked_p2p_images = [p2p%linked_p2p_images, &
                                         int(receiver_image, ocaIP)]
            end if

            p2p%send_buffer_label = [p2p%send_buffer_label, buffer_label]
        else
            ! Array will go in p2p%send_buffer(1:size(send_array))
            p2p%send_size = [int(send_array_size, ip)]
            p2p%send_start_index = [1]
            p2p%sendto_image_index = [int(receiver_image,ocaIP)]
            p2p%linked_p2p_images = p2p%sendto_image_index
            p2p%send_buffer_label = [buffer_label]
        end if

    end subroutine

    ! For rank1 send_array's
    subroutine include_in_p2p_send_buffer_rank1(p2p, send_array, buffer_label, &
        receiver_image)
        !! Record that p2p will send 'send_array' to the image 'receiver_image'. 
        !! The send is associated with the buffer_label, which will be used to 
        !! hide the book-keeping.
        type(p2p_comms_type), intent(inout) :: p2p
            !! The p2p_comms_type that will do the send
        real(dp), intent(in) :: send_array(:) 
            !! A rank-1 array with the data to send
        character(len=p2p_id_len), intent(in) :: buffer_label 
            !! The label associated with the sent data
        integer(ip), intent(in) :: receiver_image 
            !! The image that the data should be sent to

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(p2p, n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank2 send_array's
    subroutine include_in_p2p_send_buffer_rank2(p2p, send_array, buffer_label, &
        receiver_image)
        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(p2p, n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank3 send_array's
    subroutine include_in_p2p_send_buffer_rank3(p2p, send_array, buffer_label, &
        receiver_image)
        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:, :)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(p2p, n, buffer_label, &
            receiver_image)
    end subroutine

    ! For rank4 send_array's
    subroutine include_in_p2p_send_buffer_rank4(p2p, send_array, buffer_label, &
        receiver_image)
        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:, :, :)
        character(len=p2p_id_len), intent(in) :: buffer_label
        integer(ip), intent(in) :: receiver_image

        integer(ip) :: n

        n = size(send_array, kind=ip)
        call include_in_p2p_send_buffer_generic(p2p, n, buffer_label, &
            receiver_image)
    end subroutine

    subroutine allocate_p2p_comms(p2p)
        !!
        !! Allocate the send_buffer, assuming all calls to
        !! include_in_p2p_send_buffer have already been made for p2p
        !!
        type(p2p_comms_type), intent(inout) :: p2p

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

        if(p2p%have_allocated_p2p_comms) then
            write(log_output_unit,*) &
                'Error: trying to allocate p2p comms when it is already allocated'
            call generic_stop
        end if

        p2p%have_allocated_p2p_comms = .TRUE.

#ifdef COARRAY
        this_image_local = this_image2()
        num_images_local = num_images2()
#else
        this_image_local = 1
        num_images_local = 1
#endif

        ! Allocate the send buffer
        if(allocated(p2p%send_size)) then
            desired_size = sum(p2p%send_size)
        else
            desired_size = 0
        end if
        allocate(p2p%send_buffer(desired_size))

        ! Fill with a value which is suggestive of problems, in case of
        ! out-of-bounds errors
        if(size(p2p%send_buffer, kind=ip) > 0) p2p%send_buffer = HUGE(1.0_dp)

        if(allocated(p2p%send_size) .OR. allocated(p2p%sendto_image_index)) then
            if(size(p2p%send_size, kind=ip) /= size(p2p%sendto_image_index, kind=ip)) then
                stop 'BUG: p2p%send_size should have equal length to p2p%sendto_image_index'
            end if
        end if

        !
        ! Here, we optionally re-order the send metadata
        ! (i.e. change the order that the send arrays are packed). This
        ! can allow multiple sends to be combined, which may potentially have
        ! speed benefits
        !
        if(reorder_send_data_by_image .and. (size(p2p%send_size, kind=ip) > 0) &
            .and. (num_images_local > 1)) then

            ! Get the index of the send-to data
            n1 = size(p2p%send_size, kind=ip)
            allocate(send_data_order(n1), send_data_sort_criterion(n1))
            send_data_order = [(i, i=1, n1)]
            ! Make the order so that images above the current image are near the start.
            ! In some experiments this was found better than always sending to image 1
            ! first (perhaps because that will tend to make all images send to image 1, then to image 2, then ...,
            ! so might overload image 1, then image 2, etc.)
            send_data_sort_criterion = modulo(&
                (p2p%sendto_image_index - this_image_local), &
                int(num_images_local, c_int))

            !print*, 'n11: ', n1
            !print*, 'send_data_order1: ', send_data_order
            !print*, 'send_data_sort_criterion1: ', send_data_sort_criterion
            !print*, 'p2p%sendto_image_index1: ', p2p%sendto_image_index
            !print*, 'p2p%send_size1: ', p2p%send_size
            !print*, 'p2p%send_start_index1: ', p2p%send_start_index

            if(maxval(send_data_sort_criterion) /= minval(send_data_sort_criterion)) then

                call sort_index(send_data_order, send_data_sort_criterion, n1)

                ! Reorder the key data
                ! -- p2p%sendto_image_index
                ! -- p2p%send_size
                ! -- p2p%send_buffer_label
                ! -- p2p%send_start_index
                p2p%sendto_image_index = p2p%sendto_image_index(send_data_order)
                p2p%send_size = p2p%send_size(send_data_order)
                p2p%send_buffer_label = p2p%send_buffer_label(send_data_order)
                p2p%send_start_index(1) = 1
                do i = 2, n1
                    p2p%send_start_index(i) = &
                        p2p%send_start_index(i-1) + p2p%send_size(i-1)
                end do
            end if

            !print*, 'n1: ', n1
            !print*, 'send_data_order: ', send_data_order
            !print*, 'p2p%sendto_image_index: ', p2p%sendto_image_index
            !print*, 'p2p%send_size: ', p2p%send_size
            !print*, 'p2p%send_start_index: ', p2p%send_start_index

            ! Clean up
            deallocate(send_data_order, send_data_sort_criterion)
        end if

        !
        ! Allocate the receive buffer
        !
        ! First we need to find how big it should be
        desired_size_local = 0
        do i = 1, num_images_local
            if(allocated(p2p%send_size)) then
                desired_size = sum(p2p%send_size, mask=(p2p%sendto_image_index == i))
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

        ! Make the p2p%recv_buffer -- used for communication every timestep
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        !call co_max(desired_size) ! FIXME: Do we need this co-max? For coarray yes, but not for MPI?
        allocate(p2p%recv_buffer(desired_size))
        call sync_all_generic
#elif defined(COARRAY)
        call co_max(desired_size) ! Coarrays are the same size on every image. Not required for MPI. FIXME: Consider making
                                  ! "p2p%recv_buffer" point to a section of a top-level coarray
        allocate(p2p%recv_buffer(desired_size)[*])
#else
        allocate(p2p%recv_buffer(desired_size))
#endif
        ! Fill with a value which is suggestive of problems, in case of
        ! out-of-bounds mistakes
        p2p%recv_buffer = -HUGE(1.0_dp)

        !
        ! Next find out how many communications occur. We do this by
        ! broadcasting the p2p%send_size and p2p%sendto_image_index for each image, then
        ! adding the information to the p2p%recv_size / p2p%recv_start_index /
        ! p2p%recvfrom_image_index
        !
        if(allocated(p2p%send_size)) then
            desired_size = size(p2p%send_size, kind=ip)
        else
            desired_size = 0
        end if

        ! Make 'p2p%work_coarray' which is used to send/receive data between images
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call co_max(desired_size)
        allocate(p2p%work_coarray(desired_size, 3))
        allocate(wc3_all(size(p2p%work_coarray(:,3), kind=ip), num_images_local))
        call sync_all_generic
#elif defined(COARRAY)
        call co_max(desired_size)
        allocate( p2p%work_coarray(desired_size, 3)[*] )
#else
        allocate( p2p%work_coarray(desired_size, 3))
#endif

        ! Make real64_coarary to communicate the p2p%send_buffer_label (can send length
        ! p2p_id_len characters using 'transfer')
        if(modulo(storage_size(charlabel), storage_size(real(1.0, real64))) == 0) then
            n = storage_size(charlabel)/storage_size(real(1.0, real64))
        else
            write(log_output_unit, *) "storage_size(charlabel) cannot be evenly divided into real64's "
            call generic_stop
        end if
        ! Determine number of empty characters ' ' required to fill a real64
        if(modulo(storage_size(real(1.0, real64)), storage_size(" ")) == 0) then
            n2 = storage_size(real(1.0, real64))/storage_size(" ")
        else
            write(log_output_unit, *) "storage_size(real 64) cannot be evenly divided into empty space characters "
            call generic_stop
        end if

#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        allocate( p2p%real64_coarray(n, desired_size))
        call sync_all_generic
#elif defined(COARRAY)
        allocate( p2p%real64_coarray(n, desired_size)[*] )
#else
        allocate( p2p%real64_coarray(n, desired_size))
#endif

        do i = 1, num_images_local
            ! Broadcast the send metadata from image i to all images
            if( this_image_local == i) then
                p2p%work_coarray = 0
                p2p%real64_coarray = transfer(repeat(" ", n2), real(1.0, real64))

                if(allocated(p2p%send_size)) then
                    n = size(p2p%sendto_image_index, kind=ip)
                    p2p%work_coarray(1:n, 1) = p2p%sendto_image_index
                    p2p%work_coarray(1:n, 2) = p2p%send_size
                    ! To broadcast the p2p%send_buffer_label, convert to a real
                    do j = 1, n
                        p2p%real64_coarray(:, j) = transfer(&
                            p2p%send_buffer_label(j), &
                            p2p%real64_coarray(:,1), &
                            size(p2p%real64_coarray(:,1), kind=ip))
                    end do
                else
                    n = 0
                end if
            end if

#if defined(COARRAY)
            call co_broadcast(p2p%work_coarray, source_image = i)
            call co_broadcast(p2p%real64_coarray, source_image = i)
#endif

            do j = 1, size(p2p%work_coarray(:,1), kind=ip)
                ! If the current image receives data from image i, note that
                if(p2p%work_coarray(j,1) == this_image_local) then

                    ! Record the receive-from image index
                    if(allocated(p2p%recvfrom_image_index)) then
                        p2p%recvfrom_image_index = [p2p%recvfrom_image_index, i ]
                    else
                        p2p%recvfrom_image_index = [i]
                    end if

                    ! updated the 'linked images' with image i
                    if(allocated(p2p%linked_p2p_images)) then
                        if(.not. any(p2p%linked_p2p_images == i)) then
                            p2p%linked_p2p_images = [p2p%linked_p2p_images, i]
                        end if
                    else
                        p2p%linked_p2p_images = [i]
                    end if

                    if(allocated(p2p%recv_size)) then
                        p2p%recv_size = [p2p%recv_size, p2p%work_coarray(j, 2) * 1_ip ]
                    else
                        p2p%recv_size = [p2p%work_coarray(j,2) * 1_ip]
                    end if

                    if(allocated(p2p%recv_buffer_label)) then
                        p2p%recv_buffer_label = [p2p%recv_buffer_label, &
                            transfer(p2p%real64_coarray(:, j), charlabel)]
                    else
                        p2p%recv_buffer_label = &
                            [transfer(p2p%real64_coarray(:, j), charlabel)]
                    end if

                    ! Compute the p2p%recv_start_index
                    if(allocated(p2p%recv_start_index)) then
                        n = size(p2p%recv_start_index, kind=ip)
                        p2p%recv_start_index = [p2p%recv_start_index, &
                            ! Sum of previous start index + previous p2p%recv_size
                            p2p%recv_start_index(n) + p2p%recv_size(n)]
                    else
                        p2p%recv_start_index = [1_ip]
                    end if

                    ! Put the p2p%recv_start_index into the p2p%work_coarray, so
                    ! we can send it back to the p2p%recvfrom_image_index
                    n = size(p2p%recv_start_index, kind=ip)
                    p2p%work_coarray(j,3) = p2p%recv_start_index(n)
                end if
            end do
#ifdef COARRAY
            call sync_all_generic
#endif

            ! Finally, broadcast the recv_start_indices back to the sendto_image_start_index
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            ! Broadcast the 3rd column of p2p%work_coarray to image i
            call mpi_gather(p2p%work_coarray(:,3), size(p2p%work_coarray(:,3)), &
                mympi_int, wc3_all, size(p2p%work_coarray(:,3)), mympi_int, &
                i-1, MPI_COMM_WORLD, ierr)
#endif
            if(this_image_local == i) then
                if(allocated(p2p%send_size)) then
                    allocate(p2p%sendto_start_index(size(p2p%send_size, kind=ip)))
                    if(size(p2p%send_size, kind=ip) > 0) then
                        do j = 1, size(p2p%send_size, kind=ip)
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
                            p2p%sendto_start_index(j) = &
                                wc3_all(j,p2p%sendto_image_index(j))
#elif defined(COARRAY)
                            p2p%sendto_start_index(j) = &
                                p2p%work_coarray(j,3)[p2p%sendto_image_index(j)]
#else
                            p2p%sendto_start_index(j) = &
                                p2p%work_coarray(j,3)
#endif
                        end do
                    end if
                end if
            end if
#ifdef COARRAY
            call sync_all_generic
#endif
        end do

        deallocate(p2p%work_coarray)
        deallocate(p2p%real64_coarray)
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call sync_all_generic
        deallocate(wc3_all)
#endif

        ! Check that p2p%send_buffer_label does not have repeated values
        if(allocated(p2p%send_buffer_label)) then
            do i = 1, (size(p2p%send_buffer_label, kind=ip) - 1)
                do j = i+1, size(p2p%send_buffer_label, kind=ip)
                    if(p2p%send_buffer_label(i) == p2p%send_buffer_label(j)) then
                        write(log_output_unit,*) &
                            'Error: repeated p2p%send_buffer_labels: ', &
                            p2p%send_buffer_label(i)
                    end if
                end do
            end do
        end if

        ! Check that p2p%recv_buffer_label does not have repeated values
        if(allocated(p2p%recv_buffer_label)) then
            do i = 1, (size(p2p%recv_buffer_label, kind=ip) - 1)
                do j = i+1, size(p2p%recv_buffer_label, kind=ip)
                    if(p2p%recv_buffer_label(i) == p2p%recv_buffer_label(j)) then
                        write(log_output_unit,*) &
                            'Error: repeated p2p%recv_buffer_labels: ', &
                            p2p%recv_buffer_label(i)
                    end if
                end do
            end do
        end if

        if(allocated(p2p%recv_size)) then
            if(size(p2p%recv_start_index, kind=ip) /= size(p2p%recv_size, kind=ip)) then
                write(log_output_unit,*) &
                    'Error:  size(p2p%recv_start_index) /= size(p2p%recv_size) ', &
                    size(p2p%recv_start_index, kind=ip), &
                    size(p2p%recv_size, kind=ip)
                call generic_stop
            end if
        end if

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS

        ! Allocate key input variables for mpi_comms
        if(.not. reorder_send_data_by_image) then
            write(log_output_unit, *) 'ERROR: If COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined,'
            write(log_output_unit, *) 'then must have reorder_send_data_by_image=.true, so that'
            write(log_output_unit, *) 'there is only one send from image A to image B'
            call generic_stop
        end if

        if(mpi_timestep_loop_use_alltoallv) then
            ! Define variables needed for mpi_alltoallv

            allocate(p2p%mympi_send_counts(num_images_local), &
                p2p%mympi_send_displacements(num_images_local))
            allocate(p2p%mympi_recv_counts(num_images_local), &
                p2p%mympi_recv_displacements(num_images_local))

            p2p%mympi_send_counts = 0
            p2p%mympi_recv_counts = 0
            p2p%mympi_send_displacements = 0
            p2p%mympi_recv_displacements = 0

            ! Set the values for MPI_alltoallv, based on the metadata above
            do i = 1, num_images_local
                if(allocated(p2p%sendto_image_index)) then
                    if(any(p2p%sendto_image_index == i)) then
                        ! Set the send information

                        ! Find the first index with image == i
                        do j = 1, size(p2p%sendto_image_index, kind=ip)
                            if(p2p%sendto_image_index(j) == i) then
                                n = j
                                exit
                            end if
                        end do
                        ! If the start index is 'p', the MPI displacement is 'p-1'
                        p2p%mympi_send_displacements(i) = p2p%send_start_index(n) - 1
                        if(size(p2p%send_size, kind=ip) > 0) then
                            p2p%mympi_send_counts(i) = &
                                sum(p2p%send_size, mask=(p2p%sendto_image_index == i))
                        else
                            p2p%mympi_send_counts(i) = 0
                        end if
                    end if
                end if
                ! As above, for recv information
                if(allocated(p2p%recvfrom_image_index)) then
                    if(any(p2p%recvfrom_image_index == i)) then
                        ! Set the recv information

                        ! Find the first index with image == i
                        do j = 1, size(p2p%recvfrom_image_index, kind=ip)
                            if(p2p%recvfrom_image_index(j) == i) then
                                n = j
                                exit
                            end if
                        end do
                        ! If the start index is 'p', the MPI displacement is 'p-1'
                        p2p%mympi_recv_displacements(i) = p2p%recv_start_index(n) - 1
                        if(size(p2p%recv_size, kind=ip) > 0) then
                            p2p%mympi_recv_counts(i) = &
                                sum(p2p%recv_size, mask=(p2p%recvfrom_image_index == i))
                        else
                            p2p%mympi_recv_counts(i) = 0
                        end if
                    end if
                end if
            end do
        else
            ! Define variables for mpi_isend/irecv - an alternative to the alltoallv approach
            allocate(p2p%mpi_recv_requests(size(p2p%recv_start_index, kind=ip)))
            p2p%mpi_recv_requests = MPI_REQUEST_NULL
            allocate(p2p%mpi_send_requests(size(p2p%send_start_index, kind=ip)))
            p2p%mpi_send_requests = MPI_REQUEST_NULL
        end if
#endif

    end subroutine

    subroutine deallocate_p2p_comms(p2p)
        !!
        !! Clear all allocatable arrays in this module
        !!
        type(p2p_comms_type), intent(inout) :: p2p

        if(allocated(p2p%send_buffer)) deallocate(p2p%send_buffer)
        if(allocated(p2p%send_start_index)) deallocate(p2p%send_start_index)
        if(allocated(p2p%send_size)) deallocate(p2p%send_size)
        if(allocated(p2p%sendto_image_index)) deallocate(p2p%sendto_image_index)
        if(allocated(p2p%sendto_start_index)) deallocate(p2p%sendto_start_index)
        if(allocated(p2p%send_buffer_label)) deallocate(p2p%send_buffer_label)

        if(allocated(p2p%recv_buffer)) deallocate(p2p%recv_buffer)
        if(allocated(p2p%recv_start_index)) deallocate(p2p%recv_start_index)
        if(allocated(p2p%recv_size)) deallocate(p2p%recv_size)
        if(allocated(p2p%recvfrom_image_index)) deallocate(p2p%recvfrom_image_index)
        if(allocated(p2p%recv_buffer_label)) deallocate(p2p%recv_buffer_label)

        if(allocated(p2p%linked_p2p_images)) deallocate(p2p%linked_p2p_images)

        if(allocated(p2p%work_coarray)) deallocate(p2p%work_coarray)
        if(allocated(p2p%real64_coarray)) deallocate(p2p%real64_coarray)

        p2p%have_allocated_p2p_comms = .FALSE.

#ifdef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS
        if(allocated(p2p%mympi_send_counts)) deallocate(p2p%mympi_send_counts)
        if(allocated(p2p%mympi_recv_counts)) deallocate(p2p%mympi_recv_counts)
        if(allocated(p2p%mympi_send_displacements)) deallocate(p2p%mympi_send_displacements)
        if(allocated(p2p%mympi_recv_displacements)) deallocate(p2p%mympi_recv_displacements)
        if(allocated(p2p%mpi_recv_requests)) deallocate(p2p%mpi_recv_requests)
        if(allocated(p2p%mpi_send_requests)) deallocate(p2p%mpi_send_requests)
        call sync_all_generic
#endif

    end subroutine

    subroutine find_send_buffer_label_index(p2p, buffer_label, buffer_label_int)
        !!
        !! Find the index of p2p%send_buffer_label which matches 'buffer_label'.
        !!
        !! Uses a naive scan of all buffer_labels, but should be efficient enough
        !! for typical cases with a small number of 'buffer_label' values.
        !!
        type(p2p_comms_type), intent(in) :: p2p
        character(len=p2p_id_len), intent(in) :: buffer_label !! The buffer label to match
        integer(ip), intent(out) :: buffer_label_int !! The corresponding index in p2p%send_buffer_label

        integer(ip) :: i

        !write(log_output_unit,*) 'DEBUG: ', p2p%send_buffer_label

        ! Find the integer index in the metadata corresponding to
        ! buffer_label, by finding a match with the p2p%send_buffer_label's
        buffer_label_int = -1
        do i = 1, size(p2p%send_buffer_label, kind=ip)
            if(buffer_label == p2p%send_buffer_label(i)) then
                buffer_label_int = i
                exit
            end if
        end do

        if(buffer_label_int < 1) then
            write(log_output_unit,*) 'unrecognized buffer_label ', buffer_label
            call generic_stop
        end if

    end subroutine

    subroutine put_on_recv_buffer(p2p, buffer_label_int, si, ei)
        !!
        !! Put the p2p%send_buffer(si:ei) on the p2p%recv_buffer associated
        !! with p2p%sendto_start_index(buffer_label_int)
        !!
        !! Convenience routine for something we often have to do in
        !! send_ routines for all ranks
        !!
        type(p2p_comms_type), intent(inout) :: p2p
        integer(ip), intent(in) :: buffer_label_int, si, ei

        integer(ip):: i, recv_start_index_local, recv_end_index, recv_image

        i = buffer_label_int
        recv_start_index_local = p2p%sendto_start_index(i)
        recv_end_index = recv_start_index_local + p2p%send_size(i) - 1
        recv_image = p2p%sendto_image_index(i)

#if defined(COARRAY) && defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        write(log_output_unit,*) 'ERROR: Cannot call put_on_recv_buffer with COARRAY_USE_MPI_FOR_INTENSIVE_COMMS'
        flush(log_output_unit)
        call generic_stop
#endif
        ! put communication
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        p2p%recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
            p2p%send_buffer(si:ei)
#else
        p2p%recv_buffer(recv_start_index_local:recv_end_index) = &
            p2p%send_buffer(si:ei)
#endif

    end subroutine

    subroutine send_to_p2p_comms_rank1(p2p, send_array, buffer_label, &
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
        type(p2p_comms_type), intent(inout) :: p2p
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
    subroutine send_to_p2p_comms_rank2(p2p, send_array, buffer_label, &
        put_in_recv_buffer)
        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code
        include 'point2point_include_send_p2p.f90'
    end subroutine

    ! rank3 version of above
    subroutine send_to_p2p_comms_rank3(p2p, send_array, buffer_label, &
        put_in_recv_buffer)
        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code
        include 'point2point_include_send_p2p.f90'
    end subroutine

    ! rank4 version of above
    subroutine send_to_p2p_comms_rank4(p2p, send_array, buffer_label, &
        put_in_recv_buffer)

        type(p2p_comms_type), intent(inout) :: p2p
        real(dp), intent(in) :: send_array(:,:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label
        logical, intent(in), optional :: put_in_recv_buffer

        ! Generic code
        include 'point2point_include_send_p2p.f90'
    end subroutine

#ifndef COARRAY_USE_MPI_FOR_INTENSIVE_COMMS

    subroutine communicate_p2p(p2p)
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
        type(p2p_comms_type), intent(inout) :: p2p

        integer(ocaIP) :: i, start_index, end_index, recv_image, &
            recv_start_index_local, recv_end_index

        ! If there is nothing to send, exit
        if(.not. allocated(p2p%send_size)) then
            return
        end if

        if(.not. reorder_send_data_by_image) then
            ! Send ALL of the send buffers to the recv buffers, one by one
            do i = 1, size(p2p%send_size, kind=ip)
                start_index = p2p%send_start_index(i)
                end_index = start_index + p2p%send_size(i) - 1
                recv_start_index_local = p2p%sendto_start_index(i)
                recv_end_index = recv_start_index_local + p2p%send_size(i) - 1
                recv_image = p2p%sendto_image_index(i)
                ! put communication
#ifdef COARRAY
                p2p%recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                    p2p%send_buffer(start_index:end_index)
#else
                p2p%recv_buffer(recv_start_index_local:recv_end_index) = &
                    p2p%send_buffer(start_index:end_index)
#endif
            end do
        else

            ! Here, we have sorted the sends so that sends to each image are clustered together
            ! Thus we can make less sends

            ! This initialisation is a trick to get the loop to work
            start_index = 1
            end_index = start_index - 1
            recv_image = p2p%sendto_image_index(1)
            recv_start_index_local = p2p%sendto_start_index(1)
            recv_end_index = recv_start_index_local  - 1

            do i = 1, size(p2p%send_size, kind=ip)

                if(p2p%sendto_image_index(i) == recv_image) then
                    ! The 'previous' sendto image is the same as the current image
                    ! For now we just need to update the end indices
                    end_index = end_index + p2p%send_size(i)
                    recv_end_index = recv_end_index + p2p%send_size(i)
                end if

                ! If the p2p%sendto_image_index has changed, then send the data accumulated previously
                if(p2p%sendto_image_index(i) /= recv_image) then
                    ! Send the previous data
#ifdef COARRAY
                    p2p%recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                        p2p%send_buffer(start_index:end_index)
#else
                    p2p%recv_buffer(recv_start_index_local:recv_end_index) = &
                        p2p%send_buffer(start_index:end_index)
#endif

                    ! Redefine the start_index, and recv_image
                    start_index = end_index + 1
                    end_index = start_index + p2p%send_size(i) - 1
                    recv_image = p2p%sendto_image_index(i)
                    recv_start_index_local = p2p%sendto_start_index(i)
                    recv_end_index = recv_start_index_local + p2p%send_size(i) - 1
                end if

                ! If we are on the final send, we definitely need to send
                if(i == size(p2p%send_size, kind=ip)) then
#ifdef COARRAY
                    p2p%recv_buffer(recv_start_index_local:recv_end_index)[recv_image] = &
                        p2p%send_buffer(start_index:end_index)
#else
                    p2p%recv_buffer(recv_start_index_local:recv_end_index) = &
                        p2p%send_buffer(start_index:end_index)
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
    subroutine communicate_p2p(p2p)
        !!
        !!  If COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined, then use MPI for the most intensive point2point comms routine
        !!  Note this requires that the appropriate preprocessor variable has been defined
        !!
        type(p2p_comms_type), intent(inout) :: p2p

        integer:: mympi_my_ierr
        integer :: i, mpi_count, mpi_source, mpi_dest, mpi_ierr, mpi_tag

        if (.not. reorder_send_data_by_image) then
           write(log_output_unit, *) &
               'ERROR: communicate_p2p with MPI assumes reorder_send_data_by_image=.true'
           call generic_stop
        end if

        if(mpi_timestep_loop_use_alltoallv) then
            !
            ! Version using mpi_alltoallv
            !

            !NOTE -- p2p%send_size, p2p%recv_size, p2p%send_start_index, 
            !        p2p%recv_start_index, ierr, need to be of type integer
            !        also p2p%send_size, etc need to be the of size num_images()
            call mpi_alltoallv(&
                p2p%send_buffer , p2p%mympi_send_counts, p2p%mympi_send_displacements, mympi_dp, &
                p2p%recv_buffer , p2p%mympi_recv_counts, p2p%mympi_recv_displacements, mympi_dp, &
                MPI_COMM_WORLD, mympi_my_ierr)

            if(mympi_my_ierr /= 0) then
                write(log_output_unit,*) &
                    'FAIL, MPI_alltoallv error: mympi_my_ierr= ', mympi_my_ierr
                write(log_output_unit,*) __LINE__,&
                    __FILE__
                call generic_stop
            end if
        else
            !
            ! Version using repeated isend/irecv calls.
            !

            if(allocated(p2p%recv_size)) then
                ! Open up receives
                do i = 1, size(p2p%recv_start_index, kind=ip)
                    mpi_count = p2p%recv_size(i)
                    mpi_source = p2p%recvfrom_image_index(i) - 1

                    call mpi_irecv(&
                        p2p%recv_buffer(p2p%recv_start_index(i):(p2p%recv_start_index(i) + p2p%recv_size(i) - 1)), &
                        mpi_count, mympi_dp, mpi_source, MPI_ANY_TAG, &
                        MPI_COMM_WORLD, p2p%mpi_recv_requests(i), mpi_ierr)

                end do
            end if

            if(allocated(p2p%send_size)) then
                ! Do sends
                do i = 1, size(p2p%send_start_index, kind=ip)
                    mpi_count = p2p%send_size(i)
                    mpi_dest = p2p%sendto_image_index(i) - 1
                    mpi_tag = i
                    call mpi_isend(&
                        p2p%send_buffer(p2p%send_start_index(i):(p2p%send_start_index(i) + p2p%send_size(i) - 1)), &
                        mpi_count, mympi_dp, mpi_dest, mpi_tag, MPI_COMM_WORLD, &
                        p2p%mpi_send_requests(i), mpi_ierr)
                end do
            end if

            ! Ensure completion -- more strategic location of these calls would be possible
            if(allocated(p2p%send_size)) then
                mpi_count = size(p2p%send_start_index, kind=ip)
                call mpi_waitall(mpi_count, p2p%mpi_send_requests, MPI_STATUSES_IGNORE, mpi_ierr)
            end if

            if(allocated(p2p%recv_size)) then
                mpi_count = size(p2p%mpi_recv_requests, kind=ip)
                call mpi_waitall(mpi_count, p2p%mpi_recv_requests, MPI_STATUSES_IGNORE, mpi_ierr)
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
    pure subroutine recv_from_p2p_comms_rank1(p2p, recv_array, buffer_label)

        type(p2p_comms_type), intent(in) :: p2p
        real(dp), intent(inout) :: recv_array(:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank2 version of above
    pure subroutine recv_from_p2p_comms_rank2(p2p, recv_array, buffer_label)

        type(p2p_comms_type), intent(in) :: p2p
        real(dp), intent(inout) :: recv_array(:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank3 version of above
    pure subroutine recv_from_p2p_comms_rank3(p2p, recv_array, buffer_label)

        type(p2p_comms_type), intent(in) :: p2p
        real(dp), intent(inout) :: recv_array(:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    ! rank4 version of above
    pure subroutine recv_from_p2p_comms_rank4(p2p, recv_array, buffer_label)

        type(p2p_comms_type), intent(in) :: p2p
        real(dp), intent(inout) :: recv_array(:,:,:,:)
        character(len=p2p_id_len), intent(in) :: buffer_label

        include 'point2point_include_recv_p2p.f90'

    end subroutine

    function size_of_send_recv_buffers(p2p) result(mysize)
        !!
        !! Useful to know how big the send/recv buffers are
        !!
        type(p2p_comms_type), intent(in) :: p2p
        integer(ip) :: mysize

        mysize = 0
        if(allocated(p2p%send_buffer)) then
            mysize = mysize + size(p2p%send_buffer, kind=ip)*real_bytes
        end if
        if(allocated(p2p%recv_buffer)) then
            mysize = mysize + size(p2p%recv_buffer, kind=ip)*real_bytes
        end if
    end function

    subroutine print_p2p_comms(p2p, verbose)
        !!
        !! Utility printing routine, useful for debugging
        !!
        type(p2p_comms_type), intent(in) :: p2p
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
        write(log_output_unit,*) 'p2p%send_start_index: ', p2p%send_start_index
        write(log_output_unit,*) 'p2p%send_size: ', p2p%send_size
        write(log_output_unit,*) 'p2p%sendto_image_index: ', p2p%sendto_image_index
        write(log_output_unit,*) 'p2p%sendto_start_index: ', p2p%sendto_start_index
        write(log_output_unit,*) 'size(p2p%send_buffer): ', size(p2p%send_buffer, kind=ip)
        do i = 1, size(p2p%send_start_index, kind=ip)
            si = p2p%send_start_index(i)
            ei = si + p2p%send_size(i) - 1
            write(log_output_unit,*) '    ----'
            write(log_output_unit,*) '    p2p%send_buffer_label: ', trim(p2p%send_buffer_label(i))
            write(log_output_unit,*) '        send_array number: ', i
            write(log_output_unit,*) '        size:', p2p%send_size(i)
            write(log_output_unit,*) '        sendto_image: ', p2p%sendto_image_index(i)
            write(log_output_unit,*) '        p2p%sendto_start_index: ', p2p%sendto_start_index(i)
            if(verbose_in) write(log_output_unit,*) '        p2p%send_buffer: ', p2p%send_buffer(si:ei)
        end do

        !write(log_output_unit,*) 'p2p%recv_buffer: ', p2p%recv_buffer
        write(log_output_unit,*) 'size(p2p%recv_buffer): ', size(p2p%recv_buffer, kind=ip)
        write(log_output_unit,*) 'p2p%recv_start_index: ', p2p%recv_start_index
        write(log_output_unit,*) 'p2p%recv_size: ', p2p%recv_size
        write(log_output_unit,*) 'p2p%recvfrom_image_index: ', p2p%recvfrom_image_index
        do i = 1, size(p2p%recv_size, kind=ip)
            si = p2p%recv_start_index(i)
            ei = si + p2p%recv_size(i) - 1
            write(log_output_unit,*) '    ----'
            write(log_output_unit,*) '    p2p%recv_buffer_label: ', trim(p2p%recv_buffer_label(i))
            write(log_output_unit,*) '        recv_array number: ', i
            write(log_output_unit,*) '        size: ', p2p%recv_size(i)
            write(log_output_unit,*) '        p2p%recvfrom_image_index: ', p2p%recvfrom_image_index(i)
            if(verbose_in) write(log_output_unit,*) '        p2p%recv_buffer: ', p2p%recv_buffer(si:ei)
        end do

#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        end critical
#endif
    end subroutine

    subroutine test_coarray_point2point_comms_mod
        !!
        !! Unit tests. Do some point-2-point communication of arrays, and check that results are
        !! as expected
        !!

        implicit none

        type(p2p_comms_type) :: p2p
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
        call include_in_p2p_send_buffer(p2p, x_send, buffer_label=x_label, &
            receiver_image=sendto_image)

        ! Decide which image to send 'y' to
        sendto_image = ti - 1
        if(sendto_image == 0) sendto_image = ni
        y_label = 'y_comms'// repeat('*', p2p_id_len-7) ! Check we can 'fill' the character length
        call include_in_p2p_send_buffer(p2p, y_send, buffer_label=y_label, &
            receiver_image=sendto_image)

        ! Only send 'z' from image 2 to image 1
        z_label = 'z_comms'
        if(mod(ti - 2, ni) == 0) then
            call include_in_p2p_send_buffer(p2p, z_send, buffer_label=z_label, &
                receiver_image = 1_ip)
        end if

        ! Allocate send/recv buffers. From now on, we cannot further call
        ! include_in_p2p_send_buffer
        call allocate_p2p_comms(p2p)

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
            call send_to_p2p_comms(p2p, x_send, buffer_label=x_label, &
                put_in_recv_buffer=local_puts)
            call send_to_p2p_comms(p2p, y_send, buffer_label=y_label, &
                put_in_recv_buffer=local_puts)
            if(mod(ti - 2, ni) == 0) then
                call send_to_p2p_comms(p2p, z_send, buffer_label=z_label, &
                    put_in_recv_buffer=local_puts)
            end if

            if(.not.local_puts) then
                ! Do the parallel put
                call communicate_p2p(p2p)
            end if
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            sync images(p2p%linked_p2p_images)
#endif


            ! Copy the sent data to receive buffers
            call recv_from_p2p_comms(p2p, x_recv, buffer_label=x_label)
            call recv_from_p2p_comms(p2p, y_recv, buffer_label=y_label)
            if(ti == 1) then
                call recv_from_p2p_comms(p2p, z_recv, buffer_label=z_label)
            end if

            ! sync to prevent the 'k' loop moving ahead before we have received
            ! data
#if defined(COARRAY) && !defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
            sync images(p2p%linked_p2p_images)
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
                call generic_stop
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
                call generic_stop
            end if

            ! test the 'z' data
            if((ti == 1).and.(ni > 1)) then

                expected_zrecv = 55.5_dp + k
                if(all(z_recv == expected_zrecv)) then
                    write(log_output_unit,*) 'PASS'
                else
                    write(log_output_unit,*) 'FAIL', __LINE__,&
                        __FILE__
                    call generic_stop
                end if

            end if

            ! Print everything
            !call print_p2p_comms(p2p)
        end do

        ! Clean up
        call deallocate_p2p_comms(p2p)

        if(allocated(p2p%linked_p2p_images)) then
            write(log_output_unit,*) 'FAIL: p2p deallocation did not work', __LINE__,&
                __FILE__
        else
            write(log_output_unit,*) 'PASS'
        end if

    end subroutine

end module
