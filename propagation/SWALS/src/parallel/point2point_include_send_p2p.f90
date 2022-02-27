! Generic code used in sending arrays of all ranks (using "include"), see
! coarray_point2point_comms_mod.f90, particularly subroutine
! 'send_to_p2p_comms'
!
        integer(ip) :: si, ei, buffer_label_int
        logical :: put_in_recv_buffer_local

        if(.not. p2p%have_allocated_p2p_comms) then
            error stop 'Need to call allocate_p2p_comms before trying to send'
        end if

        ! If there is nothing to send, exit
        if(.not. allocated(p2p%send_buffer)) then
            return
        end if

        call find_send_buffer_label_index(p2p, buffer_label, buffer_label_int)

        ! Copy the send array to the right part of the send buffer
        si = p2p%send_start_index(buffer_label_int)
        ei = si + size(send_array, kind=ip) - 1
        ! The following routine is generic, for send_array having rank 1 to 4.
        call flatten_array(send_array, p2p%send_buffer(si:ei))

        !
        ! Communicate to the recv buffer on another image
        !
        if(present(put_in_recv_buffer)) then
            put_in_recv_buffer_local = put_in_recv_buffer
        else
            put_in_recv_buffer_local = .TRUE.
        end if

        if(put_in_recv_buffer_local) then
            call put_on_recv_buffer(p2p, buffer_label_int, si, ei)
        end if
