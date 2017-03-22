!
! Various coarray based parallel communication routines
!
#define POINT2POINT

! Compile with -DTIMERLOCAL to add timing to the code 
#ifdef TIMERLOCAL
#   define TIMER_START(tname) call timer%timer_start(tname)
#   define TIMER_STOP(tname)  call timer%timer_end(tname)
#else
#   define TIMER_START(tname)
#   define TIMER_STOP(tname)
#endif

module coarray_utilities_mod
    use global_mod, only: dp, ip, charlen
#ifdef COARRAY
    ! Need this if events are to be used. But as of 4/11/2016 it does not
    ! seem opencoarrays supports events inside a derived type
    !USE iso_fortran_env, only: event_type
#endif

#ifdef POINT2POINT
    ! Communicate using the point2point communications mod, rather than by
    ! using coarrays directly
    use coarray_point2point_comms_mod
    use reshape_array_mod, only: flatten_array
#endif

#ifdef TIMERLOCAL
    use timer_mod
#endif

    implicit none

    private
    public:: partitioned_domain_NESW_comms_type, test_coarray_utilities_mod

#ifdef TIMERLOCAL
    type(timer_type) :: timer
#endif

#ifdef COARRAY
    !! comms%put_complete records the event that the data has been put into our
    !! buffer by a neighbouring image.
    !! FIXME: Put into derived type when opencoarray supports this
    !TYPE(event_type), ALLOCATABLE :: comms_event_put_complete(:)[:]
    !! comms%get_complete records the event that the data from each boundary has been
    !! received by our neighbour image [i.e. copied from its buffer to the main array]
    !TYPE(event_type), ALLOCATABLE :: comms_event_get_complete(:)[:]

     !real(dp), allocatable :: dummy_buffer(:)[:,:]
#endif


    !
    ! Type for managing communications in a domain that is partitioned into a 'grid'
    ! of sub-domains.
    !
    type :: partitioned_domain_NESW_comms_type

        ! We want to be able to compile the code with non-coarray compilers.
        ! The preprocessing flag should ensure that the type 'does nothing gracefully'
        ! in case we don't have coarrays
#ifdef COARRAY

#ifndef POINT2POINT
        ! Main communication buffers
        ! Note: currently not all these would be needed -- either send or recv
        ! buffers alone could be used
        ! 
        real(dp), allocatable:: north_send_buffer(:,:,:)[:,:], &
            north_recv_buffer(:,:,:)[:,:]
        real(dp), allocatable:: east_send_buffer(:,:,:)[:,:], &
            east_recv_buffer(:,:,:)[:,:]
        real(dp), allocatable:: south_send_buffer(:,:,:)[:,:], &
            south_recv_buffer(:,:,:)[:,:]
        real(dp), allocatable:: west_send_buffer(:,:,:)[:,:], &
            west_recv_buffer(:,:,:)[:,:]

#endif   
        ! dummy buffer, used to identify co-subscripts
        !real(dp), allocatable :: dummy_buffer(:)[:,:]
#endif

        ! 2D co-subscripts of this_image() that are associated with the buffers
        integer(4):: ti_xy(2)
        ! The domain is split into co_size = [nx, ny] sub-domains
        integer(4):: co_size_xy(2)

        ! Dimensions of the parent grid that we are communicating from/to
        ! This is typically domain%U. It INCLUDES halo's
        integer(4):: array_shape(3) ! [nx, ny, nvar]
        ! x/y dimensions of grid on this_image(), EXCLUDING halo
        integer(4) :: interior_nx(2) 

        ! Co-subscripts for N,E,S,W neighbours. 
        !     neighbour_im_xy(1,1:2) -- North
        !     neighbour_im_xy(2,1:2) -- East
        !     neighbour_im_xy(3,1:2) -- South
        !     neighbour_im_xy(4,1:2) -- West
        ! Note that the function image_index(coarray, co-subscripts) requires
        ! 'co-subscripts' to be of rank 1 -- and this seems strict -- so e.g. we cannot
        ! pass neighbour_im_xy(1,1:2) directly to the function, we have to flatten it.
        integer(4):: neighbour_im_xy(4,2)
        !
        ! image indices of N,E,S,W neighbours (index 1 = North, then head
        ! clockwise to index 4 = West). Zero means we are on a boundary
        integer(4):: neighbour_images(4) = [0,0,0,0]
        !
        ! Hold non-zero indices of neighbour_images
        integer(4), allocatable:: neighbour_images_keep(:)
        !
        ! Hold the number of halo cells on the N,E,S,W boundary that are
        ! communicated from other domains
        integer(ip):: halo_pad(4) = [0, 0, 0, 0]

        contains

        procedure :: initialise => setup_partitioned_domain_comms
        procedure :: communicate => communicate_halos
        procedure :: print => print_comms
        ! These allowed more control over comms, but the user has to control sync's
        procedure :: copy_halos_to_buffer => copy_halos_to_buffer
        procedure :: copy_halos_from_buffer => copy_halos_from_buffer

    end type

    contains

    !
    ! Convenience printing function
    !
    subroutine print_comms(comms)
        class(partitioned_domain_NESW_comms_type), intent(in):: comms
        integer(ip):: i, j, tmp(2)

    ! If COARRAY is not defined, then remove all code from this subroutine
#ifdef COARRAY

        do j = 1, num_images()
            if(this_image() == j) then
                print*, 'Image: ', this_image()
                print*, 'ti_xy: ', comms%ti_xy
                print*, 'neighbour_im_xy: '
                do i = 1, 4
                    tmp = comms%neighbour_im_xy(i,1:2)
                    print*, '   (coindex:)', tmp
                end do
                print*, 'neighbour_images_full: ', comms%neighbour_images
                print*, 'neighbour_images: ', &
                    comms%neighbour_images(comms%neighbour_images_keep)
                print*, 'halo_pad: ', comms%halo_pad
            end if
            sync all
        end do
#endif

    end subroutine


    ! Alternative to get co-subscripts corresponding to an image index.
    !
    ! If U is a coarray with co-rank = 2 from [1:dims(1), 1:dims(2)], with
    ! dims(1)*dims(2) = num_images(), then this routine is an alternative
    ! to this_image(U). 
    !
    ! It's useful if we want to know how coarrays would be arranged, without
    ! actually allocating one
    !
    subroutine corank_from_image(dims, image_ind, corank)
        integer(ip), intent(in) :: dims(2), image_ind
        integer(4), intent(out) :: corank(2)

        ! Deliberate integer division here
        corank(2) = (image_ind - 1)/dims(1) + 1
        corank(1) = image_ind - (corank(2) - 1)*dims(1)


    end subroutine 

    ! Alternative to get image index from co-subscripts .
    ! Assumes co-rank is 2, and lower co bounds are all 1.
    !
    ! If U is a coarray with co-rank = 2 from [1:dims(1), 1:dims(2)], with
    ! dims(1)*dims(2) = num_images(), then this routine is an alternative to
    ! image_index(U, corank)
    !
    subroutine image_from_corank(dims, corank, image_ind)
        integer(ip), intent(in) :: dims(2)
        integer(4), intent(in) :: corank(2)
        integer(4), intent(out) :: image_ind

        if( (any(corank > dims)).OR.(any(corank < 1))) then
            image_ind = 0
        else
            image_ind = (corank(2) - 1) * dims(1) + corank(1)
        end if

    end subroutine 

    !
    ! Compute new lower_left, lw, nx for partitioned domains, and allocate
    ! communication buffers. Note that the domain would generally not have its
    ! main grids allocated at this stage (since the new lower_left, lw, and nx
    ! would be used to do that).
    ! 
    ! @param comms partitioned_domain_NESW_comms_type which will take care of
    !    coarray communication for the partitioned domain
    ! @param co_size_xy integer array of length 2 [nx_ca,ny_ca]. Divide the
    !    domain into [nx,ny] sub-domains, each with the same number of cells
    !    (excluding edge regions where communication happens)
    ! @param global_ll original domain lower left [x,y]
    ! @param global_lw original domain [length,width]
    ! @param global_nx original domain [nx, ny]
    ! @param local_ll new domain lower left [x,y]
    ! @param local_lw new domain [length, width]
    ! @param local_nx new domain [nx, ny]
    !
    subroutine setup_partitioned_domain_comms(comms, co_size_xy, global_ll, &
        global_lw, global_nx, local_ll, local_lw, local_nx)

        class(partitioned_domain_NESW_comms_type), intent(inout):: comms
        real(dp), intent(in):: global_ll(2), global_lw(2)
        integer(ip), intent(in):: global_nx(2), co_size_xy(2)
        real(dp), intent(out):: local_ll(2), local_lw(2)
        integer(ip), intent(out):: local_nx(2)
        
        integer(4) :: offset(2), i, counter, tmp(2), nvar, halo_width, ti
        integer(ip) :: interior_nx(2)
        real(dp):: dx(2)
        real(dp), allocatable :: dummy_array(:)
        character(len=charlen) :: buffer_label


#ifndef COARRAY

        ! If COARRAY is not defined, then just set the local_ll/lw/nx to the
        ! global values
        local_ll = global_ll
        local_lw = global_lw
        local_nx = global_nx

#else

        !ALLOCATE(comms_event_put_complete(4)[*])
        !ALLOCATE(comms_event_get_complete(4)[*])

        ! Number of variables and halo width (hard coded for now, might need to
        ! change later)
        nvar = 4
        halo_width = 2
    
        if(product(co_size_xy) /= num_images()) then
            print*, 'Error: Need to split domain into number of coarray images'
            sync all
            error stop
        end if

        dx = global_lw/(1.0_dp*global_nx)

        ! Define local domain width/nx -- these will be redefined below, to add
        ! in halo regions.
        local_lw = global_lw/(1.0_dp*co_size_xy)
        interior_nx = global_nx/co_size_xy 
        local_nx = interior_nx

        comms%interior_nx = interior_nx

        if(all(interior_nx*co_size_xy == global_nx)) then
#ifndef POINT2POINT
            ! Notice that here, interior_nx refers to the dimensions that the coarray OWNS
            allocate(&
                comms%south_send_buffer(interior_nx(1), halo_width, nvar)&
                    [1:co_size_xy(1), *], &
                comms%south_recv_buffer(interior_nx(1), halo_width, nvar)&
                    [1:co_size_xy(1), *], &
                comms%north_send_buffer(interior_nx(1), halo_width, nvar)&
                    [1:co_size_xy(1), *], &
                comms%north_recv_buffer(interior_nx(1), halo_width, nvar)&
                    [1:co_size_xy(1), *], &
                comms%east_send_buffer(halo_width, interior_nx(2), nvar)&
                    [1:co_size_xy(1), *], &
                comms%east_recv_buffer(halo_width, interior_nx(2), nvar)&
                    [1:co_size_xy(1), *], &
                comms%west_send_buffer(halo_width, interior_nx(2), nvar)&
                    [1:co_size_xy(1), *], &
                comms%west_recv_buffer(halo_width, interior_nx(2), nvar)&
                    [1:co_size_xy(1), *] )

            !allocate(dummy_buffer(1)[1:co_size_xy(1), *])
#else
            !! Use this to get neighbour co-subscripts if we don't allocate the
            !! above buffers, and for 1D comms later
            !allocate(dummy_buffer(&
            !    maxval(comms%interior_nx) * nvar * halo_width)&
            !    [1:co_size_xy(1), *])
#endif

                
        else
            print*, 'Error: domain shape does not evenly divide into coarray partition'
            error stop
        end if

        !!! Get coarray images of neighbours
        !comms%ti_xy = this_image(dummy_buffer) !
        !ti = this_image()
        !!!print*, 'before: ', ti, comms%ti_xy , '     ', co_size_xy
        !! Alternative
        call corank_from_image(co_size_xy, this_image(), comms%ti_xy)
        !print*, '    after: ', ti, comms%ti_xy
        
        ! Lower-left of this sub-domain, if there is no halo
        ! This is corrected to account for the halo below
        local_ll = global_ll + (comms%ti_xy - [1,1])*local_lw

        ! Find images of N,E,S,W neighbours
        do i = 1, 4
            select case(i)
            case(1) ! N
                offset = [0,1]
            case(2) ! E
                offset = [1,0]
            case(3) ! S
                offset = [0,-1]
            case(4) ! W
                offset = [-1,0]
            end select
            comms%neighbour_im_xy(i,1:2) = comms%ti_xy + offset 
            ! This 'copy' is required to avoid the rank being 2 in the second argument of 'image_index',
            ! which causes an error
            tmp = comms%neighbour_im_xy(i,1:2)
            !comms%neighbour_images(i) = image_index(dummy_buffer, tmp)
            !print*, 'before: ', this_image(), comms%neighbour_images(i), '    ', co_size_xy
            call image_from_corank(co_size_xy, tmp, comms%neighbour_images(i))
            !print*, 'after: ', this_image(), comms%neighbour_images(i)
        end do
    
        ! Array with the indices of neighbour images we want to keep
        allocate(comms%neighbour_images_keep(count(comms%neighbour_images > 0)))
        counter = 0
        do i = 1, 4
            if(comms%neighbour_images(i) > 0) then
                counter = counter+1
                comms%neighbour_images_keep(counter) = i
            end if
        end do

        ! Add in boundary for communication, if we need it
        if(comms%neighbour_images(1) > 0) then
            ! Boundary on north
            local_nx = local_nx + [0,halo_width]
            local_lw = local_lw + [0,halo_width]*dx
            comms%halo_pad(1) = halo_width
            
#ifdef POINT2POINT
            ! Allocate the size of the north buffer
            allocate(dummy_array(halo_width * comms%interior_nx(1) * nvar))
            buffer_label = 'N2S_NESW'
            call include_in_p2p_send_buffer(dummy_array, buffer_label, receiver_image=comms%neighbour_images(1))
            deallocate(dummy_array)
#endif            

        end if
        if(comms%neighbour_images(2) > 0) then
            ! Boundary on east
            local_nx = local_nx + [halo_width,0]
            local_lw = local_lw + [halo_width,0]*dx
            comms%halo_pad(2) = halo_width
#ifdef POINT2POINT
            ! Allocate the size of the east buffer
            allocate(dummy_array(halo_width * comms%interior_nx(2) * nvar))
            buffer_label='E2W_NESW'
            call include_in_p2p_send_buffer(dummy_array, buffer_label, receiver_image=comms%neighbour_images(2))
            deallocate(dummy_array)
#endif            
        end if
        if(comms%neighbour_images(3) > 0) then
            ! Boundary on south (so lower-left must also be changed)
            local_nx = local_nx + [0,halo_width]
            local_lw = local_lw + [0,halo_width]*dx
            local_ll = local_ll - [0,halo_width] * dx
            comms%halo_pad(3) = halo_width
#ifdef POINT2POINT
            ! Allocate the size of the south buffer
            allocate(dummy_array(halo_width * comms%interior_nx(1) * nvar))
            buffer_label='S2N_NESW'
            call include_in_p2p_send_buffer(dummy_array, buffer_label, receiver_image=comms%neighbour_images(3))
            deallocate(dummy_array)
#endif            
        end if
        if(comms%neighbour_images(4) > 0) then
            ! Boundary on west (so lower-left must also be changed)
            local_nx = local_nx + [halo_width,0]
            local_lw = local_lw + [halo_width,0]*dx
            local_ll = local_ll - [halo_width, 0]*dx
            comms%halo_pad(4) = halo_width
#ifdef POINT2POINT
            ! Allocate the size of the west buffer
            allocate(dummy_array(halo_width * comms%interior_nx(2) * nvar))
            buffer_label='W2E_NESW'
            call include_in_p2p_send_buffer(dummy_array, buffer_label, receiver_image=comms%neighbour_images(4))
            deallocate(dummy_array)
#endif            
        end if

        ! Useful to have the dimensions of the array we communicate from
        comms%array_shape = [local_nx(1), local_nx(2), nvar]

        ! Pre-initialise the get
        !event post(comms_event_get_complete(1)[this_image()])
        !event post(comms_event_get_complete(2)[this_image()])
        !event post(comms_event_get_complete(3)[this_image()])
        !event post(comms_event_get_complete(4)[this_image()])

        !deallocate(dummy_buffer)
!#ifndef POINT2POINT
!        sync all
!#endif

#endif

    end subroutine

    ! Copy the halo data from U to the send/recv buffer
    ! 
    ! @param comms communicator
    ! @param U the main array
    !
    subroutine copy_halos_to_buffer(comms, U)
        class(partitioned_domain_NESW_comms_type), intent(inout):: comms
        !class(domain_type), intent(inout):: domain
        real(dp), intent(inout):: U(:,:,:)

        integer(ip):: Np, Ep, Sp, Wp, ca(2), nx, ny, nvar, arr_size
        character(len=charlen) :: buffer_label

#ifdef COARRAY
        ! Get the number of 'padded' cells on NESW sides, over which communication will occur
        Np = comms%halo_pad(1)
        Ep = comms%halo_pad(2)
        Sp = comms%halo_pad(3)
        Wp = comms%halo_pad(4)

        ! Dimension of U = [nx, ny, nvar]
        nx = comms%array_shape(1)
        ny = comms%array_shape(2)
        nvar = comms%array_shape(3) 

        ! Note that the coarray size is (nx - Ep - Wp), (ny - Np - Sp), and this is constant
        ! on all images irrespective of the halo_pad, because that is how nx,ny were defined

        if(comms%neighbour_images(1) /= 0) then
#ifndef POINT2POINT
            comms%north_send_buffer = U((1+Wp):(nx-Ep), (ny-Np - 1):(ny-Np), 1:nvar)
            ca = comms%neighbour_im_xy(1,1:2)
            ! Ensure that the neighbour has already got the previous data
            !event wait(comms_event_get_complete(1))
            comms%south_recv_buffer(:,:,:)[ca(1), ca(2)] = comms%north_send_buffer(:,:,:)
            !event post(comms_event_put_complete(3)[comms%neighbour_images(1)])
#else
            buffer_label = 'N2S_NESW'
            call send_to_p2p_comms(U((1+Wp):(nx-Ep), (ny-Np - 1):(ny-Np), 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(2) /= 0) then
#ifndef POINT2POINT
            comms%east_send_buffer = U((nx-Ep-1):(nx-Ep), (1+Sp):(ny-Np), 1:nvar)
            ca = comms%neighbour_im_xy(2,1:2)
            ! Ensure that the neighbour has already got the previous data
            !event wait(comms_event_get_complete(2))
            comms%west_recv_buffer(:,:,:)[ca(1), ca(2)] = comms%east_send_buffer(:,:,:)
            !event post(comms_event_put_complete(4)[comms%neighbour_images(2)])
#else
            buffer_label = 'E2W_NESW'
            call send_to_p2p_comms(U((nx-Ep-1):(nx-Ep), (1+Sp):(ny-Np), 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(3) /= 0) then
#ifndef POINT2POINT
            comms%south_send_buffer = U((1+Wp):(nx-Ep), (1+Sp):(2+Sp), 1:nvar) 
            ca = comms%neighbour_im_xy(3,1:2)
            ! Ensure that the neighbour has already got the previous data
            !event wait(comms_event_get_complete(3))
            comms%north_recv_buffer(:,:,:)[ca(1), ca(2)] = comms%south_send_buffer(:,:,:)
            !event post(comms_event_put_complete(1)[comms%neighbour_images(3)])
#else
            buffer_label = 'S2N_NESW'
            call send_to_p2p_comms(U((1+Wp):(nx-Ep), (1+Sp):(2+Sp), 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(4) /= 0) then
#ifndef POINT2POINT
            comms%west_send_buffer = U((1+Wp):(2+Wp), (1+Sp):(ny-Np), 1:nvar)
            ca = comms%neighbour_im_xy(4,1:2)
            !event wait(comms_event_get_complete(4))
            comms%east_recv_buffer(:,:,:)[ca(1), ca(2)] = comms%west_send_buffer(:,:,:)
            !event post(comms_event_put_complete(2)[comms%neighbour_images(4)])
#else
            buffer_label = 'W2E_NESW'
            call send_to_p2p_comms(U((1+Wp):(2+Wp), (1+Sp):(ny-Np), 1:nvar), buffer_label)
#endif
        end if

#endif
    end subroutine

    ! Assuming recv buffer has been updated, copy the data to U
    ! 
    ! @param comms communicator
    ! @param U the main array
    !
    subroutine copy_halos_from_buffer(comms, U)
        class(partitioned_domain_NESW_comms_type), intent(inout):: comms
        !class(domain_type), intent(inout):: domain
        real(dp), intent(inout):: U(:,:,:)

        integer(ip):: Np, Ep, Sp, Wp, ca(2), nx, ny, nvar, arr_size
        character(len=charlen) :: buffer_label

#ifdef COARRAY
        ! Get the number of 'padded' cells on NESW sides, over which communication will occur
        Np = comms%halo_pad(1)
        Ep = comms%halo_pad(2)
        Sp = comms%halo_pad(3)
        Wp = comms%halo_pad(4)

        ! Dimension of U = [nx, ny, nvar]
        nx = comms%array_shape(1)
        ny = comms%array_shape(2)
        nvar = comms%array_shape(3) 

        !! Copy from buffer to U
        if(comms%neighbour_images(1) /= 0) then
#ifndef POINT2POINT
            !event wait(comms_event_put_complete(1))
            U((1+Wp):(nx-Ep), (ny - Np+1):ny, 1:nvar) = comms%north_recv_buffer
            !event post(comms_event_get_complete(3)[comms%neighbour_images(1)])
#else
            buffer_label = 'S2N_NESW'
            call recv_from_p2p_comms(U((1+Wp):(nx-Ep), (ny - Np+1):ny, 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(2) /= 0) then
#ifndef POINT2POINT
            !event wait(comms_event_put_complete(2))
            U((nx-Ep+1):nx, (1+Sp):(ny-Np), 1:nvar) = comms%east_recv_buffer
            !event post(comms_event_get_complete(4)[comms%neighbour_images(2)])
#else
            buffer_label = 'W2E_NESW'
            call recv_from_p2p_comms(U((nx-Ep+1):nx, (1+Sp):(ny-Np), 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(3) /= 0) then
#ifndef POINT2POINT
            !event wait(comms_event_put_complete(3))
            U((1+Wp):(nx-Ep), 1:Sp, 1:nvar) = comms%south_recv_buffer
            !event post(comms_event_get_complete(1)[comms%neighbour_images(3)])
#else
            buffer_label = 'N2S_NESW'
            call recv_from_p2p_comms(U((1+Wp):(nx-Ep), 1:Sp, 1:nvar), buffer_label)
#endif
        end if
        if(comms%neighbour_images(4) /= 0) then
#ifndef POINT2POINT
            !event wait(comms_event_put_complete(4))
            U(1:Wp, (1+Sp):(ny-Np), 1:nvar) = comms%west_recv_buffer
            !event post(comms_event_get_complete(2)[comms%neighbour_images(4)])
#else
            buffer_label = 'E2W_NESW'
            call recv_from_p2p_comms(U(1:Wp, (1+Sp):(ny-Np), 1:nvar), buffer_label)
#endif
        end if

#endif
    end subroutine

    !
    ! Communicate 'halo' information between sub-domains
    !
    ! @param comms the partitioned_domain_NESW_comms_type 
    ! @param domain the related domain object
    !
    subroutine communicate_halos(comms, U)
        class(partitioned_domain_NESW_comms_type), intent(inout):: comms
        !class(domain_type), intent(inout):: domain
        real(dp), intent(inout):: U(:,:,:)

        integer(ip):: Np, Ep, Sp, Wp, ca(2), nx, ny, nvar

        ! If COARRAY is not defined, then do nothing!
#ifdef COARRAY

        TIMER_START('to_buffer')
        ! Step 1
        ! Copy from U to buffers
        call copy_halos_to_buffer(comms, U)

        TIMER_STOP('to_buffer')
        ! Need to have finshed the communication before we can procced to copy
        ! buffers into U
        sync images(comms%neighbour_images(comms%neighbour_images_keep))
        !sync all

        TIMER_START('from_buffer')
        call copy_halos_from_buffer(comms, U)
        TIMER_STOP('from_buffer')

        TIMER_START('sync')
        ! Need to sync again, to prevent the possibility that
        ! before the above buffer copy is complete, another image
        ! could have completed the copy, proceeded and updated the recv buffer again.
        !
        ! Probably this could be more efficiently done with events
        sync images(comms%neighbour_images(comms%neighbour_images_keep))
        !sync all
        !sync memory
        TIMER_STOP('sync')
#endif
    end subroutine

    !
    ! Main test
    !
    subroutine test_coarray_utilities_mod()

        real(dp), allocatable :: U(:,:,:), xs(:), ys(:), storage(:)
        type(partitioned_domain_NESW_comms_type) :: comms
        real(dp):: local_ll(2), local_dx(2), local_lw(2)
        real(dp):: global_ll(2), global_dx(2), global_lw(2)
        integer(ip):: local_nx(2), global_nx(2)
        integer(4):: coarray_size(2) 
        integer(ip):: i, j, ti, ni, ncomms
        real(dp):: error_tol

#ifndef COARRAY
        print*, '    (skipped as code was not compiled with -DCOARRAY)'
#else

        ! Choose coarray layout. Might not be the most time efficient, but good to test more
        if(mod(num_images(),2) == 0) then
            if(mod(num_images(), 4) == 0) then
                coarray_size = [num_images()/4, 4]
            else
                coarray_size = [num_images()/2, 2]
            end if
        else
            coarray_size = [1, num_images()]
        end if

        ! Make up some values
        global_lw = [3600.0_dp, 1800.0_dp]
        global_nx = [384_ip, 192_ip]
        global_dx = global_lw/global_nx
        global_ll = [0._dp, 0._dp]

        ! Get the domain info
        call comms%initialise(coarray_size, global_ll, global_lw, global_nx, &
            local_ll, local_lw, local_nx)
#ifdef POINT2POINT
        call allocate_p2p_comms
#endif        
        local_dx = global_dx

        ! Allocate the 'local patch' of U
        allocate(U(local_nx(1), local_nx(2), 4), xs(local_nx(1)), ys(local_nx(2)), storage(maxval(local_nx)))
        ! Compute the 'local patch' values for xs, ys
        xs = (/ ( local_ll(1) + (i-0.5_dp)*local_lw(1)/local_nx(1), i=1, local_nx(1) )  /)
        ys = (/ ( local_ll(2) + (i-0.5_dp)*local_lw(2)/local_nx(2), i=1, local_nx(2) )  /)

        ! Communicate many times
        do ncomms = 1, 100

            ! Make some large numbers
            do j = 1, local_nx(2)
                do i = 1, local_nx(1)
                    U(i,j,1) = xs(i) + 0*ys(j) + 1 + ncomms 
                    U(i,j,2) = 0*xs(i) + ys(j) + 2 + ncomms
                    U(i,j,3) = xs(i) + 0*ys(j) + 3 + ncomms
                    U(i,j,4) = 0*xs(i) + ys(j) + 4 + ncomms
                end do
            end do
        
            ! Do parallel exchange, which should cause no changes to U given how it was defined
            call comms%communicate(U)

            do j = 1, local_nx(2)
                do i = 1, local_nx(1)
                    U(i,j,1) = U(i,j,1) - (xs(i) + 0*ys(j) + 1 + ncomms) 
                    U(i,j,2) = U(i,j,2) - (0*xs(i) + ys(j) + 2 + ncomms) 
                    U(i,j,3) = U(i,j,3) - (xs(i) + 0*ys(j) + 3 + ncomms) 
                    U(i,j,4) = U(i,j,4) - (0*xs(i) + ys(j) + 4 + ncomms) 
                end do
            end do

            ! We may tolerate tiny errors caused by the different ways xs/ys would
            ! be calculated (although not if the dimensions are sufficiently regular)
            error_tol = 0.0_dp !maxval(local_dx)*EPSILON(max(maxval(xs), maxval(ys)))*3000
            ! Might need to increase this to e.g. 2.0e-05 when using single precision with high core counts

            if(all(abs(U) <= error_tol)) then
                !print*, 'PASS image', this_image(), maxval(abs(U)), local_ll + comms%halo_pad([4,3])*local_dx, &
                !    local_ll + local_lw - comms%halo_pad([2,1])*local_dx, ncomms
            else
                print*, 'FAIL image', this_image(), maxval(abs(U)), local_ll + comms%halo_pad([4,3])*local_dx, &
                    local_ll + local_lw - comms%halo_pad([2,1])*local_dx, ncomms
                error stop
            end if
        end do

        print*, 'PASS'

#ifdef POINT2POINT
        call deallocate_p2p_comms
#endif

#ifdef TIMERLOCAL
        if(this_image() == 1) then
            call timer%print()
        end if
#endif
        
#endif

    end subroutine

end module

