!
! Some useful routines for extrapolation and slope-limiting. These do not rely on the domain_type.
!
module extrapolation_limiting_mod

    use global_mod, only: dp, ip, charlen
    implicit none

    real(dp), parameter :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp

    contains
    !
    ! minmod function, which is used in some gradient limiters
    !
    ! @param a,b real numbers
    !
    elemental function minmod(a,b) result(minmod_ab)
        real(dp), intent(in):: a, b
        real(dp):: minmod_ab
        
        minmod_ab = merge(min(abs(a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
        !minmod_ab = sign(one_dp, a) * max(zero_dp, min(abs(a), sign(one_dp, a)*b))

    end function

    !
    ! minmod subroutine, which is used in some gradient limiters
    !
    ! @param a,b real numbers
    !
    elemental subroutine minmod_sub(a,b, minmod_ab)
        real(dp), intent(in):: a, b
        real(dp), intent(out):: minmod_ab

        if((a>0.0_dp .and. b>0.0_dp).or.(a<0.0_dp .and. b<0.0_dp)) then
            minmod_ab = min(abs(a), abs(b))*sign(ONE_dp, a)
        else
            minmod_ab = ZERO_dp
        end if

    end subroutine

    !
    ! Get the "gradient times dx" around U_local
    ! @param U_local -- variable at the cell of interest, say x
    ! @param U_lower -- variable at (x - dx)
    ! @param U_upper -- variable at (x + dx)
    ! @param theta -- limiter parameter
    ! @param gradient_dx -- hold output
    ! @param n -- length of each vector
    !
    subroutine limited_gradient_dx_vectorized(U_local, U_lower, U_upper, theta, gradient_dx, n)
        integer(ip), intent(in):: n
        real(dp), intent(in):: U_local(n), U_lower(n), U_upper(n), theta(n) 
        real(dp), intent(out) :: gradient_dx(n)

        character(len=charlen), parameter :: limiter_type = 'MC' !'Minmod2' !'Superbee_variant' ! 'MC'! 'nolimit' !

        integer(ip) :: i
        real(dp):: a, b, c, d, e, th, sa, sb, half_sasb

        if(limiter_type == 'MC') then

            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                th = theta(i)
                sa = sign(ONE_dp,a)
                sb = sign(ONE_dp,b)
                half_sasb = HALF_dp * (sa + sb)
                d = min(abs(a), abs(b)) * half_sasb * th ! Limit on local gradient
                e = HALF_dp * (a + b)
                c = merge(ZERO_dp, e, d == ZERO_dp) 
                ! NOTE: IF d /= 0, then clearly d, c have the same sign
                ! We exploit this to avoid a further minmod call (which seems
                ! expensive)
                gradient_dx(i) = merge(min(c, d), max(c, d), d > ZERO_dp)
            end do

        else if(limiter_type == "Superbee_variant") then

            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                ! Divide by 1.6 which is the default 'max theta' in the rk2 algorithms
                th = theta(i) * 2.0_dp/1.6_dp
                !d = minmod(a, th*b)
                d = merge(min(abs(a), abs(th*b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
                !e = minmod(th*a, b)
                e = merge(min(abs(th*a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
                if(abs(e) > abs(d)) then
                    b = e
                else
                    b = d
                endif
                gradient_dx(i) = b 
            end do

        else if(limiter_type == "Minmod2") then
            ! Same as 'MC' (!!!)
            !$OMP SIMD
            do i = 1, n

                a = U_upper(i) - U_local(i)
                b = U_local(i) - U_lower(i)
                th = theta(i)
                e = HALF_dp * (a + b)
                a = a * th
                b = b * th
                if(b > ZERO_dp .and. a > ZERO_dp) then
                    ! Positive slopes
                    if(b < a) a = b
                    d = min(e, a)
                else if(b < ZERO_dp .and. a < ZERO_dp) then
                    ! Negative slopes
                    if(b > a) a = b
                    d = max(e, a)
                else
                    d = ZERO_dp
                endif
                gradient_dx(i) = d
            end do

        else if(limiter_type == 'nolimit') then

            gradient_dx = HALF_dp * (U_upper - U_lower)

        else 
            gradient_dx = ZERO_dp
        end if
            

    end subroutine


    !
    ! Very limited testing of the domain routines.
    ! More important are the validation tests.
    !
    subroutine test_extrapolation_limiting_mod
    
        integer(ip), parameter :: N = 10
        real(dp) :: U(N), U_lower(N), U_upper(N)
        real(dp) :: theta(N), extrapolation_sign(N)
        real(dp) :: desired_answer(N), answer(N)
        integer(ip) :: i

        theta = 1.0_dp
        U = (/(i*1.0_dp - 5.5_dp, i = 1, N)/)
        U_lower = U - 1.0_dp
        U_upper = U + 1.0_dp

        ! Basic extrapolation tests, positive side
        extrapolation_sign=1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * (U + U_upper)
        call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side, reduced theta
        theta = 0.5_dp
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 0.5_dp * U + 0.5_dp*0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! As above, negative side, zero theta
        theta = 0.0_dp
        extrapolation_sign = 1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = 1.0_dp * U + 0.0_dp*0.5_dp * (U + U_lower)
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Local min.
        theta = 4.0_dp
        extrapolation_sign = 1.0_dp
        U_lower = U + 2.0_dp
        U_upper = U + 3.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U 
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Local max.
        theta = 4.0_dp
        U_lower = U - 2.0_dp
        U_upper = U - 3.0_dp
        extrapolation_sign = -1.0_dp
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U 
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        ! Limited gradients
        ! 
        ! Central gradient = 2, lower-gradient=1, upper-gradient=3
        U_lower = U - 1.0_dp
        U_upper = U + 3.0_dp
        extrapolation_sign = 1.0_dp
        theta = 1.0_dp ! Limited-gradient = 1
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U + 0.5_dp * 1.0_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        theta = 1.5_dp ! Limited-gradient = 1.5
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U + 0.5_dp * 1.5_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        
        extrapolation_sign = -1.0_dp ! As above, negative side
        call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
        desired_answer = U - 0.5_dp * 1.5_dp
        call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

        contains

            subroutine assert_equal_within_tol(desired_answer, answer, line)
                real(dp) :: desired_answer(N), answer(N)
                integer :: line

                real(dp) :: tol

                tol = spacing(1.0_dp) * 100

                if(all(abs(desired_answer - answer) < tol)) then
                    print*, 'PASS'
                else
                    print*, 'FAIL', line, desired_answer - answer
                end if

            end subroutine

    end subroutine

end module
