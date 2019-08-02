module timer_mod
    !
    ! Module to time parts of code
    !
    ! Notice we do not use dp/ip precision variables -- instead use C_INT, C_DOUBLE

#ifdef COARRAY
    ! Assumes that mpi is available with coarrays, which should
    ! be the case e.g. for opencoarrays
    use mpi, only: mpi_wtime
#endif
#ifndef NOOPENMP
    use omp_lib
#endif
    use iso_fortran_env, only: output_unit
    use iso_c_binding, only: C_DOUBLE, C_INT

    implicit none

    integer, parameter, private:: max_timers = 100
    integer, parameter, private:: charlen_timer_names=128
    real(C_DOUBLE), parameter :: unstarted = -HUGE(1.0_C_DOUBLE)

    private
    public:: timer_type

    type:: timer_type
        character(len=charlen_timer_names) :: names(max_timers)
        real(C_DOUBLE):: start(max_timers) = unstarted
        real(C_DOUBLE):: total(max_timers) = 0.0_C_DOUBLE
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

        integer:: tname_index

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

        integer:: tname_index
        real(C_DOUBLE):: current_time

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

        if(timer%start(tname_index) == unstarted) then
            stop 'Timer not started'
        end if
       
        timer%total(tname_index) = timer%total(tname_index) + (current_time - timer%start(tname_index))
        timer%start(tname_index) = unstarted

    end subroutine

    !
    ! @param timer
    ! @param output_file_unit Optional unit of file to write the information to. If not provide, uses output_unit
    ! from iso_fortran_env, which prints to screen
    !
    subroutine timer_print(timer, output_file_unit)
        class(timer_type), intent(in):: timer
        integer(C_INT), optional:: output_file_unit

        integer:: i, out_unit
        real(C_DOUBLE):: timer_sum

        if(present(output_file_unit)) then
            out_unit = output_file_unit
        else
            out_unit = output_unit
        end if

        timer_sum = sum(timer%total)

        write(out_unit, "(A)") '############################################'
        write(out_unit, "(A)") 'WALLCLOCK TIME SUMMARY'
        write(out_unit, "(A)") '--------------------------------------------'
        write(out_unit, "(A)") 'Details'
        if(timer%ntimers > 0) then
            do i = 1, timer%ntimers
                write(out_unit, "(A5, A50, A2, F20.8, A, F20.16, A2)") '    ', TRIM(timer%names(i)), ': ', timer%total(i), &
                    ' : ', timer%total(i)/timer_sum*100.0_C_DOUBLE, TRIM(' %')
            end do
        end if
        write(out_unit, "(A)") '---------------'
        write(out_unit, "(A, F20.8)") 'Total WALLCLOCK time: ', sum(timer%total)
        write(out_unit, "(A)") '############################################'
    end subroutine

end module
