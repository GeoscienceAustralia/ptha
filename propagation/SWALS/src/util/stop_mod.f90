! Provide a 'stop' command that works with both CAF and non-CAF code
module stop_mod
#ifdef MPI
    use mpi
#endif
    implicit none

    contains

    subroutine generic_stop()
        integer :: ierr

        call flush()
        ! Call the most relevant 'stop' command
#ifdef COARRAY
        error stop
#elif defined(MPI)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#else
        stop
#endif
    end subroutine

end module
