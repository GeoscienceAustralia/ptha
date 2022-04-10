! Generic code used in receiving arrays of all ranks (using "include"), see
! coarray_point2point_comms_mod.f90, particularly subroutine
! 'recv_from_p2p_comms'
!
!

        integer(ip) :: i, buffer_label_int, si, ei

        if(.not. p2p%have_allocated_p2p_comms) then
            error stop 'Need to call allocate_p2p_comms before trying to receive'
        end if

        ! If there is nothing to receive, exit
        if(.not. allocated(p2p%recv_buffer_label)) then
            return
        end if

        ! Find the integer recv index corresponding to the label
        buffer_label_int = -1
        do i = 1, size(p2p%recv_buffer_label, kind=ip)
            if(buffer_label == p2p%recv_buffer_label(i)) then
                buffer_label_int = i
                exit
            end if
        end do

        ! Ensure there was a match
        if(buffer_label_int < 1) then
            !print*, 'unrecognized buffer_label :', buffer_label
            error stop 'unrecognized buffer_label '
        end if

        ! Check the array sizes agree
        if(size(recv_array, kind=ip) /= p2p%recv_size(buffer_label_int)) then
            !print*, 'size(recv_array) = ', size(recv_array), &
            !    ' does not match buffer size (', &
            !    recv_size(buffer_label_int), ') for ', &
            !    ' comms ', buffer_label, ' on image ', this_image_local
            error stop 'size of recv_array does not match size of p2p%send_array'
        end if

        !Start/end indices
        si = p2p%recv_start_index(buffer_label_int)
        ei = si + size(recv_array, kind=ip) - 1

        ! Main copy. Note this routine is generic, with recv_array having rank 1
        ! to 4
        call repack_rank1_array(p2p%recv_buffer(si:ei), recv_array)
