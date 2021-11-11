! Generic code used in sending arrays of all ranks, see
! coarray_point2point_comms_mod.f90, particularly subroutine
! 'send_to_p2p_comms
!
! I don't know how to code this without repetition, except for using 'include'
!
        integer(ip) :: si, ei, buffer_label_int
        logical :: put_in_recv_buffer_local

        if(.not. have_allocated_p2p_comms) then
            error stop 'Need to call allocate_p2p_comms before trying to send'
        end if

        ! If there is nothing to send, exit
        if(.not. allocated(send_buffer)) then
            return
        end if

        call find_send_buffer_label_index(buffer_label, buffer_label_int)

        ! Copy the send array to the right part of the send buffer
        si = send_start_index(buffer_label_int)
        ei = si + size(send_array, kind=ip) - 1
        ! The following routine is generic, for send_array having rank 1 to 4.
        call flatten_array(send_array, send_buffer(si:ei))
        ! Intrinsic alternative (but would this make a copy?
        !send_buffer(si:ei) = pack(send_array, .true.)

        !
        ! Communicate to the recv buffer on another image
        !
        if(present(put_in_recv_buffer)) then
            put_in_recv_buffer_local = put_in_recv_buffer
        else
            put_in_recv_buffer_local = .TRUE.
        end if

        if(put_in_recv_buffer_local) then
            call put_on_recv_buffer(buffer_label_int, si, ei)
        end if
