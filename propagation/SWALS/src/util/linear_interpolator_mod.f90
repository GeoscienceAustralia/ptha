module linear_interpolator_mod
    !
    ! Type for linear interpolation
    !
    use global_mod, only: dp, ip
    use stop_mod, only: generic_stop
    implicit none


    type linear_interpolator_type
        real(dp), allocatable:: xs_local(:), ys_local(:)
        real(dp), pointer :: xs(:), ys(:)
        integer(ip) :: n

        contains
        ! To build the an interpolator from x,y data (with x monotonic
        !   increasing, NO REPEATED VALUES), we do
        ! CALL linear_interpolator%initialise(x,y)
        ! or to not copy x,y
        ! CALL linear_interpolator%initialise(x,y, copy_data=.FALSE.)
        !
        procedure:: initialise => initialise_linear_interpolator
        !
        ! CALL linear_interpolator%eval(xout, yout)
        ! will update yout with values interpolated at xout
        !
        procedure:: eval => eval_linear_interpolator
        !
        ! CALL linear_interpolator%finalise()
        ! to clear pointers, deallocate data, etc
        !
        procedure:: finalise => finalise_linear_interpolator

    end type

    contains

    subroutine initialise_linear_interpolator(linear_interpolator, x, y, copy_data)
        class(linear_interpolator_type), target, intent(inout):: linear_interpolator
        real(dp), target, intent(in):: x(:), y(:)
        logical, optional, intent(in):: copy_data
        integer(ip):: i, n
        logical:: use_pointers
        

        if (present(copy_data)) then
            use_pointers = merge(.FALSE., .TRUE., copy_data)
        else
            use_pointers = .FALSE.
        end if

        n = size(x) 
        linear_interpolator%n = n
        if(n /= size(y)) then
            print*, 'initialise_linear_interpolator error: size of x and y must be equal'
            call generic_stop()
        end if

        ! Set xs, ys, optionally making a local copy of the data
        if(.not. use_pointers) then
            allocate(linear_interpolator%xs_local(n))
            linear_interpolator%xs_local = x
            linear_interpolator%xs => linear_interpolator%xs_local
            allocate(linear_interpolator%ys_local(n))
            linear_interpolator%ys_local = y
            linear_interpolator%ys => linear_interpolator%ys_local
        else
            linear_interpolator%xs => x 
            linear_interpolator%ys => y 
        end if
      
        ! Check xs is monotonic 
        do i = 1, size(linear_interpolator%xs)-1
            if(linear_interpolator%xs(i) >= linear_interpolator%xs(i+1)) then
                print*, 'initialise_linear_interpolator error: x must be monotonic increasing (no repeated values)'
                call generic_stop()
            end if
        end do

    end subroutine

    subroutine finalise_linear_interpolator(linear_interpolator)
        class(linear_interpolator_type), intent(inout):: linear_interpolator

        if(allocated(linear_interpolator%xs_local)) then        
            deallocate(linear_interpolator%xs_local, linear_interpolator%ys_local)
        end if

        linear_interpolator%xs => NULL()
        linear_interpolator%ys => NULL()

    end subroutine

    subroutine eval_linear_interpolator(linear_interpolator, output_x, output_y)
        class(linear_interpolator_type), intent(in):: linear_interpolator
        real(dp), intent(in):: output_x(:)
        real(dp), intent(out):: output_y(:)
        integer(ip):: n

        n = size(output_x)
        if(n /= size(output_y)) then
            print*, 'eval_linear_interpolator error: size of x and y must be equal'
            call generic_stop()
        end if

        call linear_interpolation(linear_interpolator%n, linear_interpolator%xs, &
            linear_interpolator%ys, n, output_x, output_y)

    end subroutine

    !
    ! Suppose x is a SORTED vector of length n with x(i) <= x(i+1).
    ! We want to find the index corresponding to the value in 'x' that is nearest y.
    ! This is a useful operation e.g. for interpolation of sorted input data
    !
    pure subroutine nearest_index_sorted(n, x, y, output)
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y      
        integer(ip), intent(out) :: output

        integer :: i, upper, lower

        if(y < x(1)) then
            output = 1
        else 
            if (y > x(n)) then
                output = n
            else
                lower = 1
                upper = n
                ! Deliberate integer division
                i = (lower + upper)/2 !floor(0.5_dp*(lower + upper))
                do while ((x(i) > y).or.(x(i+1) < y))
                   if(x(i) > y) then
                       upper = i
                   else
                       lower = i
                   end if
                   ! Deliberate integer division
                   i = (lower + upper)/2 !floor(0.5_dp*(lower + upper))
                end do

                if ( y - x(i) > x(i+1) - y) then
                    output = i+1
                else
                    output = i
                end if

            end if
        end if

    end subroutine

    !
    ! Suppose input_x, input_y are some data with input_x being sorted and increasing 
    ! We are given output_x, and wish to get output_y by linearly interpolating the
    ! original series.
    pure subroutine linear_interpolation(n_input, input_x, input_y, n_output, output_x, output_y)
        integer(ip), intent(in):: n_input, n_output
        real(dp), intent(in):: input_x(n_input), input_y(n_input), output_x(n_output)
        real(dp), intent(out):: output_y(n_output)

        integer(ip):: i, j
        real(dp):: gradient

        do i = 1, n_output
            ! set j = nearest index to output_x(i) in input_x
            call nearest_index_sorted(n_input, input_x, output_x(i), j)
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
        end do

    end subroutine

    subroutine test_linear_interpolator_mod
        real(dp):: x(5), y(5), xout(9), yout_true(9), yout(9)
        type(linear_interpolator_type):: li
        integer(ip):: i, j
        logical:: copy_data

        ! Data to interpolate from
        x = [1.0_dp, 4.0_dp, 4.01_dp, 10.0_dp, 11.0_dp]
        y = [-1.0_dp, -4.0_dp, 4.01_dp, 10.0_dp, 110.0_dp]

        ! X values to interpolate at
        xout      = [-10.0_dp, 1.0_dp, 2.5_dp,  4.0_dp, 4.005_dp, 6.0_dp, 10.5_dp, 11.0_dp, 13.0_dp]
        ! Y values that we should get from interpolating
        yout_true = [-1.0_dp, -1.0_dp, -2.5_dp, -4.0_dp, 0.005_dp, 6.0_dp, 60.0_dp, 110.0_dp, 110.0_dp]

        do j = 1, 2
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

            ! Check that xs_local / ys_local are used or not appropriately
            if(copy_data) then
                if(allocated(li%xs_local) .and. allocated(li%ys_local)) then
                    print*, 'PASS'
                else
                    print*, 'FAIL'
                end if
            else
                ! Should not make a copy
                if(allocated(li%xs_local) .or. allocated(li%ys_local)) then
                    print*, 'FAIL'
                else
                    print*, 'PASS'
                end if
            end if

            ! Do the interpolation
            call li%eval(xout, yout)

            do i = 1, size(xout)
                if(abs(yout(i) - yout_true(i)) < 1.0e-6_dp) then
                    print*, 'PASS'
                else
                    print*, 'FAIL', xout(i), yout(i), yout_true(i)
                end if
            end do

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

        end do

    end subroutine

end module

