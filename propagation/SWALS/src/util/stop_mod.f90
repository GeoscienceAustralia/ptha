module stop_mod
    !! Provide a 'stop' command that works with both CAF and non-CAF code
    use logging_mod, only: log_output_unit
    implicit none

    contains

    subroutine generic_stop()
        !! Alternative to 'stop' and 'error stop' that tries to flush log files before halting.
        !! This makes it more likely that a useful error message is recorded.
        integer :: ierr, i
        logical :: is_open_file_unit

        print*, 'FAIL -- call to generic stop'

        ! Search for file units in the first 100 files, and flush if open
        ! A better but more complex approach would be to track open file units
        do i = 1, 10000
            inquire(unit=i, opened=is_open_file_unit)
            if(is_open_file_unit) call flush(i)
        end do

        ! Make sure we do not miss the log output unit
        inquire(unit=log_output_unit, opened=is_open_file_unit)
        if(is_open_file_unit) call flush(log_output_unit)

        ! Call the most relevant 'stop' command
#ifdef COARRAY
        error stop
#else
        stop
#endif
    end subroutine

end module
