module timer_mod
    !
    ! Module to time parts of code
    !
    use global_mod, only: dp, charlen, ip
#ifdef COARRAY
    ! Assumes that mpi is available with coarrays, which should
    ! be the case e.g. for opencoarrays
    use mpi, only: mpi_wtime
#endif
#ifndef NOOPENMP
    use omp_lib
#endif
    use iso_fortran_env, only: output_unit
    use iso_c_binding, only: C_DOUBLE

    implicit none

    integer, parameter, private:: max_timers = 100

    private
    public:: timer_type

    type:: timer_type
        character(len=charlen) :: names(max_timers)
        real(C_DOUBLE):: start(max_timers) = 0.0_dp
        real(C_DOUBLE):: total(max_timers) = 0.0_dp
        integer:: ntimers = 0
        integer:: last_index = 1

        contains

        procedure:: timer_start => timer_start
        procedure:: timer_end => timer_end
        procedure:: print => timer_print
        procedure:: find_index => find_index

    end type

contains

    subroutine find_index(timer, tname, tname_index)
        class(timer_type), intent(inout):: timer
        character(len=*), intent(in):: tname
        integer, intent(out):: tname_index

        integer:: i

        tname_index = 0

        ! Search for tname in the timer names
        if(timer%ntimers > 0) then

            do i = 1, timer%ntimers
                if(tname == timer%names(i)) then
                    tname_index = i
                    timer%last_index = i
                    return
                end if
            end do

        end if

    end subroutine


    subroutine timer_start(timer, tname)
        class(timer_type), intent(inout):: timer
        character(len=*), intent(in):: tname

        integer:: tname_index, i

        call find_index(timer, tname, tname_index)

        ! If not found, add it
        if(tname_index == 0) then
            tname_index = timer%ntimers + 1
            timer%ntimers = tname_index
            timer%names(tname_index) = tname
        end if

        ! Set the start time
#ifdef COARRAY
        timer%start(tname_index) = mpi_wtime()
#elif defined(NOOPENMP)
        call cpu_time(timer%start(tname_index))
#else
        timer%start(tname_index) = omp_get_wtime()
#endif
       
    end subroutine


    subroutine timer_end(timer, tname)
        class(timer_type), intent(inout):: timer
        character(len=*), intent(in):: tname

        integer tname_index, i
        real(C_DOUBLE) current_time

        ! Get the current time 
#ifdef COARRAY
        current_time = mpi_wtime()
#elif defined(NOOPENMP)
        call cpu_time(current_time)
#else
        current_time =  omp_get_wtime()
#endif

        call find_index(timer, tname, tname_index)

        ! If not found, error
        if(tname_index == 0) then
            print*, tname
            stop 'No timer with this name'
        end if

        if(timer%start(tname_index) == 0.0_C_DOUBLE) then
            stop 'Timer not started'
        end if
       
        timer%total(tname_index) = timer%total(tname_index) + (current_time - timer%start(tname_index))
        timer%start(tname_index) = 0.0_C_DOUBLE

    end subroutine

    !
    ! @param timer
    ! @param output_file_unit Optional unit of file to write the information to. If not provide, uses output_unit
    ! from iso_fortran_env, which prints to screen
    !
    subroutine timer_print(timer, output_file_unit)
        class(timer_type), intent(in):: timer
        integer(ip), optional:: output_file_unit

        integer:: i, out_unit
        real(dp):: timer_sum

        if(present(output_file_unit)) then
            out_unit = output_file_unit
        else
            out_unit = output_unit
        end if

        timer_sum = sum(timer%total)

        write(out_unit, *) '############################################'
        write(out_unit, *) 'WALLCLOCK TIME SUMMARY'
        write(out_unit, *) '--------------------------------------------'
        write(out_unit, *) 'Details'
        if(timer%ntimers > 0) then
            do i = 1, timer%ntimers
                write(out_unit, *) '    ', TRIM(timer%names(i)), ': ', timer%total(i), &
                    ' : ', timer%total(i)/timer_sum*100.0_dp, TRIM('%')
            end do
        end if
        write(out_unit, *) '---------------'
        write(out_unit, *) 'Total WALLCLOCK time: ', sum(timer%total)
        write(out_unit, *) '############################################'
    end subroutine

end module
