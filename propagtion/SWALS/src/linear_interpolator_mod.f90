MODULE linear_interpolator_mod
    !
    ! Type for linear interpolation
    !
    USE global_mod, only: dp, ip
    USE stop_mod, only: generic_stop
    IMPLICIT NONE


    TYPE linear_interpolator_type
        REAL(dp), ALLOCATABLE:: xs_local(:), ys_local(:)
        REAL(dp), POINTER :: xs(:), ys(:)
        INTEGER(ip) :: n

        CONTAINS
        ! To build the an interpolator from x,y data (with x monotonic
        !   increasing, NO REPEATED VALUES), we do
        ! CALL linear_interpolator%initialise(x,y)
        ! or to not copy x,y
        ! CALL linear_interpolator%initialise(x,y, copy_data=.FALSE.)
        !
        PROCEDURE:: initialise => initialise_linear_interpolator
        !
        ! CALL linear_interpolator%eval(xout, yout)
        ! will update yout with values interpolated at xout
        !
        PROCEDURE:: eval => eval_linear_interpolator
        !
        ! CALL linear_interpolator%finalise()
        ! to clear pointers, deallocate data, etc
        !
        PROCEDURE:: finalise => finalise_linear_interpolator

    END TYPE

    CONTAINS

    SUBROUTINE initialise_linear_interpolator(linear_interpolator, x, y, copy_data)
        CLASS(linear_interpolator_type), TARGET, INTENT(INOUT):: linear_interpolator
        REAL(dp), TARGET, INTENT(IN):: x(:), y(:)
        LOGICAL, OPTIONAL, INTENT(IN):: copy_data
        INTEGER(ip):: i, n
        LOGICAL:: use_pointers
        

        if (present(copy_data)) then
            use_pointers = merge(.FALSE., .TRUE., copy_data)
        else
            use_pointers = .FALSE.
        end if

        if(use_pointers) then
            ALLOCATE(linear_interpolator%xs_local(n))
            linear_interpolator%xs_local = x
            linear_interpolator%xs => linear_interpolator%xs_local
            ALLOCATE(linear_interpolator%ys_local(n))
            linear_interpolator%ys_local = y
            linear_interpolator%ys => linear_interpolator%ys_local
        else
            linear_interpolator%xs => x 
            linear_interpolator%ys => y 
        end if
       
        n = size(x) 
        linear_interpolator%n = n
        if(n /= size(y)) then
            print*, 'initialise_linear_interpolator error: size of x and y must be equal'
            call generic_stop()
        end if

        DO i = 1, size(linear_interpolator%xs)-1
            if(linear_interpolator%xs(i) >= linear_interpolator%xs(i+1)) then
                print*, 'initialise_linear_interpolator error: x must be monotonic increasing (no repeated values)'
                call generic_stop()
            end if
        END DO

    END SUBROUTINE

    SUBROUTINE finalise_linear_interpolator(linear_interpolator)
        CLASS(linear_interpolator_type), INTENT(INOUT):: linear_interpolator

        if(allocated(linear_interpolator%xs_local)) then        
            DEALLOCATE(linear_interpolator%xs_local, linear_interpolator%ys_local)
        end if

        linear_interpolator%xs => NULL()
        linear_interpolator%ys => NULL()

    END SUBROUTINE

    SUBROUTINE eval_linear_interpolator(linear_interpolator, output_x, output_y)
        CLASS(linear_interpolator_type), INTENT(IN):: linear_interpolator
        REAL(dp), INTENT(IN):: output_x(:)
        REAL(dp), INTENT(OUT):: output_y(:)
        INTEGER(ip):: i, n

        n = size(output_x)
        if(n /= size(output_y)) then
            print*, 'eval_linear_interpolator error: size of x and y must be equal'
            call generic_stop()
        end if

        CALL linear_interpolation(linear_interpolator%n, linear_interpolator%xs, &
            linear_interpolator%ys, n, output_x, output_y)

    END SUBROUTINE

    !
    ! Suppose x is a SORTED vector of length n with x(i) <= x(i+1).
    ! We want to find the index corresponding to the value in 'x' that is nearest y.
    ! This is a useful operation e.g. for interpolation of sorted input data
    !
    PURE SUBROUTINE nearest_index_sorted(n, x, y, output)
        INTEGER(ip), INTENT(IN) :: n
        REAL(dp), INTENT(IN) :: x(n)
        REAL(dp), INTENT(IN) :: y      
        INTEGER(ip), INTENT(OUT) :: output

        INTEGER :: i, upper, lower

        IF(y < x(1)) THEN
            output = 1
        ELSE 
            IF (y > x(n)) THEN
                output = n
            ELSE
                lower = 1
                upper = n
                i = (lower + upper)/2 !floor(0.5_dp*(lower + upper))
                DO WHILE ((x(i) > y).OR.(x(i+1) < y))
                   IF(x(i) > y) THEN
                       upper = i
                   ELSE
                       lower = i
                   END IF
                   i = (lower + upper)/2 !floor(0.5_dp*(lower + upper))
                END DO

                IF ( y - x(i) > x(i+1) - y) THEN
                    output = i+1
                ELSE
                    output = i
                END IF

            END IF
        END IF

    END SUBROUTINE

    !
    ! Suppose input_x, input_y are some data with input_x being sorted and increasing 
    ! We are given output_x, and wish to get output_y by linearly interpolating the
    ! original series.
    PURE SUBROUTINE linear_interpolation(n_input, input_x, input_y, n_output, output_x, output_y)
        INTEGER(ip), INTENT(IN):: n_input, n_output
        REAL(dp), INTENT(IN):: input_x(n_input), input_y(n_input), output_x(n_output)
        REAL(dp), INTENT(OUT):: output_y(n_output)

        INTEGER(ip):: i, j
        REAL(dp):: gradient

        DO i = 1, n_output
            ! set j = nearest index to output_x(i) in input_x
            CALL nearest_index_sorted(n_input, input_x, output_x(i), j)
            ! interpolate
            if (input_x(j) > output_x(i)) then
                if(j > 1) then
                    gradient = (input_y(j) - input_y(j-1))/(input_x(j) - input_x(j-1))
                    output_y(i) = input_y(j-1) +  gradient * (output_x(i) - input_x(j-1))
                else
                    output_y(i) = input_y(1)
                end if
            else
                if(j < n_input) then
                    gradient = (input_y(j+1) - input_y(j))/(input_x(j+1) - input_x(j))
                    output_y(i) = input_y(j) +  gradient * (output_x(i) - input_x(j))
                else
                    output_y(i) = input_y(n_input)
                end if
            end if
        END DO

    END SUBROUTINE

    SUBROUTINE test_linear_interpolator_mod
        REAL(dp):: x(5), y(5), xout(9), yout_true(9), yout(9)
        TYPE(linear_interpolator_type):: li
        INTEGER(ip):: i, j
        LOGICAL:: copy_data

        ! Data to interpolate from
        x = [1.0_dp, 4.0_dp, 4.01_dp, 10.0_dp, 11.0_dp]
        y = [-1.0_dp, -4.0_dp, 4.01_dp, 10.0_dp, 110.0_dp]

        ! X values to interpolate at
        xout      = [-10.0_dp, 1.0_dp, 2.5_dp,  4.0_dp, 4.005_dp, 6.0_dp, 10.5_dp, 11.0_dp, 13.0_dp]
        ! Y values that we should get from interpolating
        yout_true = [-1.0_dp, -1.0_dp, -2.5_dp, -4.0_dp, 0.005_dp, 6.0_dp, 60.0_dp, 110.0_dp, 110.0_dp]

        DO j = 1, 2
            ! Test cases with/without pointers
            if(j == 1) then
                copy_data = .true.
            else
                copy_data = .false.
            end if

            ! Set it up correctly
            call li%initialise(x,y, copy_data)
            if(all(abs(li%xs - x) < 1.0e-10_dp)) then
                print*, 'PASS'
            else
                print*, 'FAIL', li%xs - x
            end if
            
            if(all(abs(li%ys - y) < 1.0e-10_dp)) then
                print*, 'PASS'
            else
                print*, 'FAIL', li%ys - y
            end if

            ! Do the interpolation
            call li%eval(xout, yout)

            DO i = 1, size(xout)
                if(abs(yout(i) - yout_true(i)) < 1.0e-6_dp) then
                    print*, 'PASS'
                else
                    print*, 'FAIL', xout(i), yout(i), yout_true(i)
                end if
            END DO

            ! Clean up
            call li%finalise()

            if(allocated(li%ys_local) .OR. allocated(li%xs_local)) then
                print*, 'FAIL'
            else
                print*, 'PASS'
            end if
            if(associated(li%ys) .OR. associated(li%xs)) then
                print*, 'FAIL'
            else
                print*, 'PASS'
            end if

        END DO

    END SUBROUTINE

END MODULE

