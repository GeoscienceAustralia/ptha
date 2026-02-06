module quadratic_extrapolation_mod
    ! Use results at 3 time levels to extrapolate a solution in time by fitting a quadratic

    use global_mod, only: dp, ip
    use stop_mod, only: generic_stop
    use logging_mod, only: log_output_unit

    implicit none

    private
    public :: quadratic_extrapolation_type, test_quadratic_extrapolation_mod

    real(dp), parameter :: missing_time = HUGE(1.0_dp)*0.987654321_dp 
        !! Time value to represent uninitialised times, unlikely to appear elsewhere

    type quadratic_extrapolation_type
        !! Given three (time, value) points with ordered times
        !!   (t*_{-2}, y_{-2}), (t*_{-1}, y_{-1}), (t*_{0}, y_{0)}
        !! we want to approximate y_{1} at another time t*_{1}
        !! using a quadratic fit to the three points
        real(dp) :: tn2 = missing_time, tn1 = missing_time, t0 = missing_time
            !! Times corresponding to yn2, yn1, y0. FIXME: Consider making this force_double
        integer(ip) :: nx = -1, ny = -1, nz = -1
            !! Dimensions of arrays that are extrapolated in time, y0(nx, ny, nz)
        real(dp), allocatable :: yn2(:,:,:), yn1(:,:,:), y0(:,:,:)
            !! Values observed at times tn2, tn1, t0, which will be used for extrapolation

        contains

        procedure, non_overridable :: setup => setup_quadratic_extrapolation_type
        procedure, non_overridable :: append_value_and_time_for_quadratic_extrapolation => &
            append_value_and_time_for_quadratic_extrapolation
        procedure, non_overridable :: extrapolate_in_time => extrapolate_in_time

    end type

    contains

    subroutine setup_quadratic_extrapolation_type(qet, nx, ny, nz)
        class(quadratic_extrapolation_type), intent(inout) :: qet
        integer(ip), intent(in) :: nx, ny, nz

        ! Make space to store historic values for quadratic extrapolation
        qet%nx = nx
        qet%ny = ny
        qet%nz = nz
        if(allocated(qet%yn2)) deallocate(qet%yn2)
        if(allocated(qet%yn1)) deallocate(qet%yn1)
        if(allocated(qet%y0)) deallocate(qet%y0)
        allocate(qet%yn2(nx, ny, nz), qet%yn1(nx, ny, nz), qet%y0(nx, ny, nz))

    end subroutine 

    subroutine append_value_and_time_for_quadratic_extrapolation(qet, new_t0, new_y0)
        !! Append value at a new time to be used for quadratic extrapolation,
        !! and replace the old 'n-2' time with the old 'n-1' time and the old 'n-1' time with the old 'n' time
        class(quadratic_extrapolation_type), intent(inout) :: qet
        real(dp), intent(in) :: new_t0, new_y0(:,:,:)

        integer(ip) :: i, j, k
  
        ! Replace the "n-2" value with the "n-1" value 
        qet%tn2 = qet%tn1 
        call swap_alloc_rank3(qet%yn1, qet%yn2)

        ! Replace the "n-1" value with the "n" value
        qet%tn1 = qet%t0
        call swap_alloc_rank3(qet%y0, qet%yn1)

        ! Replace the "n" value with the new value
        qet%t0  = new_t0
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(qet, new_y0) COLLAPSE(2)
        do k = 1, qet%nz
            do j = 1, qet%ny
                do i = 1, qet%nx
                    qet%y0(i,j,k) = new_y0(i,j,k)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
  
        contains 
        subroutine swap_alloc_rank3(x0, x1)
            !! Swap x0 and x0, both being rank 3 real(dp) arrays
            real(dp), intent(inout), allocatable :: x0(:,:,:), x1(:,:,:)
            real(dp), allocatable :: tmp(:,:,:)     

            call move_alloc(x0, tmp) ! x0 is deallocated, tmp = x0
            call move_alloc(x1, x0) ! x1 is deallocated, x0 = x1
            call move_alloc(tmp, x1) ! tmp is deallocated, x1 = tmp = x0
        end subroutine
    end subroutine

    subroutine extrapolate_in_time(qet, new_time, new_y, do_nothing_if_missing_times, did_extrapolate)
        !! Say we have 3 (time, value) points with ordered times
        !!   (t*_{-2}, y_{-2}), (t*_{-1}, y_{-1}), (t*_{0}, y_{0)}
        !! and we want to approximate y_{1} at a future time t*_{1}
        !! using a quadratic fit to the three points
        !!
        !! For simplicity define a new time that is zero at t*_{0}
        !!    t_{i} = t*_{i} - t*_{0}
        !!
        !! A polynomial fit is:
        !!    a_{0} + a_{1} t_{i} + a_{2} t_{i}^2 = y_{i}
        !! So
        !!    a_{0}                                 = y_{0}
        !!    a_{0} + a_{1} t_{-1} + a_{2} t_{-1}^2 = y_{-1}
        !!    a_{0} + a_{1} t_{-2} + a_{2} t_{-2}^2 = y_{-2}
        !!
        !! Only a_{1} and a_{2} are non-trivial unknowns, and they satisfy
        !!   [ [ t_{-1},  t_{-1}^2 ]   %*% [a_{1}  = [y{-1} - y{0},
        !!     [ t_{-2},  t_{-2}^2 ] ]      a_{2}]    y{-2} - y{0} ]
        !!
        !!           M                 %*% x       =     b  
        !!
        !! The LHS 2x2 matrix M has an inverse M^{-1}
        !! 1/(t_{-1}*t_{-2}^2 - t_{-1}^2*t_{-2}) * [ [t_{-2}^2,  -t_{-1}^2 ]
        !!                                           [-t_{-2} ,  t_{-1}    ] ]
        !!
        !! and this gives (introducing premult = 1/(t_{-1}*t_{-2}^2 - t_{-1}^2*t_{-2}) )
        !!
        !! premult * [ t_{-2}^2 * ( y_{-1} - y_{0} ) - t_{-1}^2 * ( y_{-2} - y_{0} )   = [a_{1}
        !!             -t_{-2}  * ( y_{-1} - y_{0} ) + t_{-1}   * ( y_{-2} - y_{0} ) ]    a_{2} ]
        !!
        !! So the value of a_{0} + a_{1} t_{1} + a_{2} t_{1}^2 = y1 is
        !!    y1 = y0 + 
        !!	      t1   * premult * ( tn2^2 * ( yn1 - y0 ) - tn1^2 * ( yn2 - y0 ) ) + 
        !!         t1^2 * premult * ( -tn2  * ( yn1 - y0 ) + tn1   * ( yn2 - y0 ) )
        !! which can be rearranged as
        !!    y1 = y0 * (1 + t1*premult*(-tn2^2 + tn1^2) + t1^2 * premult * (tn2 - tn1)) + 
        !!         yn1* (    t1*premult*tn2^2            + t1^2 * premult * (-tn2)) + 
        !!         yn2* (    t1 * premult * (-tn1^2)     + t1^2 * premult * tn1 )
        !!       = y0 * a + yn1 * b + yn2 * c

        class(quadratic_extrapolation_type), intent(in) :: qet
        real(dp), intent(in) :: new_time
            !! Time at which we wish to extrapolate
        real(dp), intent(inout) :: new_y(:,:,:)
            !! Used to store extrapolated values
        logical, optional, intent(in) :: do_nothing_if_missing_times
            !! If .false. (default) and this routine is called before three time levels
            !! have been appended to qet, then thow an error. Otherwise return silently with new_y unchanged.
        logical, optional, intent(out) :: did_extrapolate

        logical :: do_nothing_if_invalid
        real(dp) :: tn2, tn1, t0, t1, a, b, c, premult
        integer(ip) :: i, j, k

        if(present(do_nothing_if_missing_times)) then
            do_nothing_if_invalid = do_nothing_if_missing_times
        else
            do_nothing_if_invalid = .false.
        end if

        ! It's possible that qet is not yet storing 3 time levels
        if(qet%tn2 == missing_time .or. qet%tn1 == missing_time .or. qet%t0 == missing_time) then
            if(do_nothing_if_invalid) then
                if(present(did_extrapolate)) did_extrapolate = .false.
                ! Silently do nothing 
                return
            else
                write(log_output_unit, *) 'Error: called extrapolate_in_time before initialising fully, times are ', &
                    qet%tn2, ', ', qet%tn1, ', ', qet%t0
                call generic_stop
            end if
        end if

        if(present(did_extrapolate)) did_extrapolate = .true.

        ! Recentre the times so that t0 = 0
        ! FIXME: Beware introduction of round-off errors here. Realistically not a problem with double precision
        tn2 = qet%tn2 - qet%t0
        tn1 = qet%tn1 - qet%t0
        t0  = 0.0_dp
        t1  = new_time - qet%t0

        ! Extrapolate y
        premult = 1.0_dp/(tn1*tn2**2 - tn2*tn1**2)
        a = (1.0_dp + t1*premult*(-tn2**2 + tn1**2) + t1**2 * premult * (tn2 - tn1))
        b = (    t1*premult*tn2**2          + t1**2 * premult * (-tn2))
        c = (    t1 * premult * (-tn1**2)     + t1**2 * premult * tn1 )

        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(qet, new_y, a, b, c) COLLAPSE(2)
        do k = 1, qet%nz
            do j = 1, qet%ny
                do i = 1, qet%nx
                ! y1 = y0 * a + yn1 * b + yn2 * c
                new_y(i,j,k) = qet%y0(i,j,k) * a + qet%yn1(i,j,k) * b + qet%yn2(i,j,k) * c
                end do
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

    subroutine test_quadratic_extrapolation_mod
        type(quadratic_extrapolation_type) :: qet
        integer(ip), parameter :: nx = 2, ny = 1, nz = 2
        real(dp) :: t0, tn1, tn2, t1
        real(dp) :: y0(nx, ny, nz), yn1(nx, ny, nz), yn2(nx, ny, nz), y1(nx, ny, nz), y1_approx(nx, ny, nz)
        real(dp), parameter :: errtol = max(10*spacing(10.0_dp), 1.0e-10_dp)

        call qet%setup(nx, ny, nz)

        t1 = 1.0_dp
        y1 = -2.0_dp

        ! Check if we can force 'do nothing' when qet is not properly initialised
        !call qet%extrapolate_in_time(t1, y1) !FAILS
        call qet%extrapolate_in_time(t1, y1, do_nothing_if_missing_times=.TRUE.) ! Does nothing
        if(all(y1 == -2.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Choose time values so that known solution is 3*y0 - 3*yn1 + yn2
        !
        tn2 = -2.0_dp
        tn1 = -1.0_dp
        t0 = 0.0_dp
        t1 = 1.0_dp

        yn2 = 5.1_dp
        yn1 = 3.2_dp
        y0  = 1.0_dp

        call qet%append_value_and_time_for_quadratic_extrapolation(tn2, yn2)
        call qet%append_value_and_time_for_quadratic_extrapolation(tn1, yn1)
        call qet%append_value_and_time_for_quadratic_extrapolation(t0, y0)
        call qet%extrapolate_in_time(t1, y1)

        ! Check the error
        y1 = y1 - (3.0_dp * y0 - 3.0_dp * yn1 + yn2)
        if(any(abs(y1) > errtol)) then
            print*, 'FAIL', maxval(abs(y1))
        else
            print*, 'PASS'
        end if

        !
        ! Chose values with solution known from a preliminary R implementation. This problem has uneven times
        !

        t0 = 0.4900490049_dp
        y0 = 2.46620225979_dp
        tn1 = 0.2340490049_dp
        yn1 = 2.3288972602_dp
        tn2 = -0.3291509951_dp
        yn2 = 1.99424077945_dp
        t1 = 0.89217286456_dp
        y1_approx = 2.66318916921_dp ! Answer from the other code

        call qet%append_value_and_time_for_quadratic_extrapolation(tn2, yn2)
        call qet%append_value_and_time_for_quadratic_extrapolation(tn1, yn1)
        call qet%append_value_and_time_for_quadratic_extrapolation(t0, y0)
        call qet%extrapolate_in_time(t1, y1)
        y1 = y1 - y1_approx
        if(any(abs(y1) > errtol)) then
            print*, 'FAIL', maxval(abs(y1))
        else
            print*, 'PASS'
        end if

    end subroutine
end module
