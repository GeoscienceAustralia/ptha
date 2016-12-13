MODULE timer_mod
    !
    ! Module to time parts of code
    !
    USE global_mod, only: dp, charlen, ip
    USE omp_lib
    USE iso_fortran_env, only: output_unit
    IMPLICIT NONE

    INTEGER, PARAMETER, PRIVATE:: max_timers = 100

    PRIVATE
    PUBLIC:: timer_type

    TYPE:: timer_type
        CHARACTER(len=charlen) :: names(max_timers)
        REAL(dp):: start(max_timers) = 0.0_dp
        REAL(dp):: total(max_timers) = 0.0_dp
        INTEGER:: ntimers = 0
        INTEGER:: last_index = 1

        CONTAINS

        PROCEDURE:: timer_start => timer_start
        PROCEDURE:: timer_end => timer_end
        PROCEDURE:: print => timer_print
        PROCEDURE:: find_index => find_index

    END TYPE

CONTAINS

    SUBROUTINE find_index(timer, tname, tname_index)
        CLASS(timer_type), INTENT(INOUT):: timer
        CHARACTER(len=*), INTENT(IN):: tname
        INTEGER, INTENT(OUT):: tname_index

        INTEGER:: i

        tname_index = 0

        ! Search for tname in the timer names
        IF(timer%ntimers > 0) THEN

            ! Loop
            DO i = 1, timer%ntimers
                IF(tname == timer%names(i)) THEN
                    tname_index = i
                    timer%last_index = i
                    RETURN
                END IF
            END DO

        END IF

    END SUBROUTINE


    SUBROUTINE timer_start(timer, tname)
        CLASS(timer_type), INTENT(INOUT):: timer
        CHARACTER(len=*), INTENT(IN):: tname

        INTEGER:: tname_index, i

        CALL find_index(timer, tname, tname_index)

        ! If not found, add it
        IF(tname_index == 0) THEN
            tname_index = timer%ntimers + 1
            timer%ntimers = tname_index
            timer%names(tname_index) = tname
        END IF

        ! Set the start time
#ifdef NOOPENMP
        CALL cpu_time(timer%start(tname_index))
#else
        timer%start(tname_index) = omp_get_wtime()
#endif
       
    END SUBROUTINE


    SUBROUTINE timer_end(timer, tname)
        CLASS(timer_type), INTENT(INOUT):: timer
        CHARACTER(len=*), INTENT(IN):: tname

        INTEGER tname_index, i
        REAL(dp) current_time

        ! Get the current time 
#ifdef NOOPENMP
        call cpu_time(current_time)
#else
        current_time =  omp_get_wtime()
#endif

        call find_index(timer, tname, tname_index)

        ! If not found, error
        IF(tname_index == 0) THEN
            PRINT*, tname
            STOP('No timer with this name')
        END IF

        IF(timer%start(tname_index) == 0._dp) THEN
            STOP('Timer not started')
        END IF
       
        timer%total(tname_index) = timer%total(tname_index) + (current_time - timer%start(tname_index))
        timer%start(tname_index) = 0.0

    END SUBROUTINE

    !
    ! @param timer
    ! @param output_file_unit Optional unit of file to write the information to. If not provide, uses output_unit
    ! from iso_fortran_env, which prints to screen
    !
    SUBROUTINE timer_print(timer, output_file_unit)
        CLASS(timer_type), INTENT(IN):: timer
        INTEGER(ip), OPTIONAL:: output_file_unit

        INTEGER:: i, out_unit
        REAL(dp):: timer_sum

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
        IF(timer%ntimers > 0) THEN
            DO i = 1, timer%ntimers
                write(out_unit, *) '    ', TRIM(timer%names(i)), ': ', timer%total(i), &
                    ' : ', timer%total(i)/timer_sum*100.0_dp, TRIM('%')
            END DO
        END IF
        write(out_unit, *) '---------------'
        write(out_unit, *) 'Total WALLCLOCK time: ', sum(timer%total)
        write(out_unit, *) '############################################'
    END SUBROUTINE

END MODULE
