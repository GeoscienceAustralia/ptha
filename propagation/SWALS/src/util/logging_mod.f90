module logging_mod
    !!
    !! Module to make an 'image specific' file log when using coarrays. 
    !! We 'use' log_output_unit to output generic print type statements to the log.
    !! By default log_output_unit=stdout. But it can also be directed to a file.
    !!

    use global_mod, only: charlen
    use coarray_intrinsic_alternatives
    use iso_fortran_env, only: output_unit
    implicit none

    public :: log_output_unit, send_log_output_to_file

    integer :: log_output_unit = output_unit !! File unit for output

    contains

    subroutine send_log_output_to_file(filename_prefix)
        !! Direct the log_output_unit to a file, with name like (filename_prefix + '.log'). 
        !! If coarrays are used, then the name will be like (filename_prefix + '_image_' + my_image + '.log')
        character(len=*), intent(in) :: filename_prefix

        character(len=charlen) :: log_filename, ti_char

        log_filename = trim(filename_prefix) // '.log'

#ifdef COARRAY
        write(ti_char, '(I20.20)') this_image2()
        log_filename = trim(filename_prefix) // '_image_' // trim(ti_char) // '.log'
#endif

        ! Update the output unit
        open(newunit=log_output_unit, file=log_filename)

    end subroutine

end module
