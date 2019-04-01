!
! Module to help have an 'image specific' log for the coarray case
! We 'use' log_output_unit to output generic print type statements.
! By default it just goes to stdout, but it can also be sent to a file, for
! easier interpretation of parallel outputs
!
module logging_mod

    use global_mod, only: charlen
    use iso_fortran_env, only: output_unit
    implicit none

    public :: log_output_unit, send_log_output_to_file

    integer :: log_output_unit = output_unit

    contains

    ! Call this to send the log to a file (image specific, if coarrays are used)
    subroutine send_log_output_to_file(filename_prefix)
        character(len=*), intent(in) :: filename_prefix

        character(len=charlen) :: log_filename, ti_char

        log_filename = trim(filename_prefix) // '.log'

#ifdef COARRAY
        write(ti_char, '(I0.20)') this_image()
        log_filename = trim(filename_prefix) // '_image_' // trim(ti_char) // '.log'
#endif

        ! Update the output unit
        open(newunit=log_output_unit, file=log_filename)

    end subroutine

end module
