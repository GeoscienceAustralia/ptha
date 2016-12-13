! Provide a 'stop' command that works with both CAF and non-CAF code
MODULE stop_mod
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE generic_stop()
#ifdef COARRAY
    call flush()
    error stop
#else
    call flush()
    stop
#endif
    END SUBROUTINE

END MODULE
