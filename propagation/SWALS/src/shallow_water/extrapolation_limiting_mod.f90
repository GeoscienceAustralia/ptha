module extrapolation_limiting_mod
    !!
    !! Some useful routines for extrapolation and slope-limiting. These do not rely on the domain_type.
    !!

    use global_mod, only: dp, ip, charlen
    implicit none

    real(dp), parameter :: HALF_dp = 0.5_dp, ZERO_dp = 0.0_dp, ONE_dp=1.0_dp

    contains

    elemental function minmod(a,b) result(minmod_ab)
        !!
        !! minmod function, which is used in some gradient limiters
        !!
        !! @param a,b real numbers
        !!
        real(dp), intent(in):: a, b
        real(dp):: minmod_ab
        
        minmod_ab = merge(min(abs(a), abs(b))*sign(ONE_dp,a), ZERO_dp, sign(ONE_dp,a) == sign(ONE_dp,b))
        !minmod_ab = sign(one_dp, a) * max(zero_dp, min(abs(a), sign(one_dp, a)*b))

    end function

    elemental subroutine minmod_sub(a,b, minmod_ab)
        !!
        !! minmod subroutine, which is used in some gradient limiters
        !!
        !! @param a,b real numbers
        !!
        real(dp), intent(in):: a, b
        real(dp), intent(out):: minmod_ab

        if((a>ZERO_dp .and. b>ZERO_dp).or.(a<ZERO_dp .and. b<ZERO_dp)) then
            minmod_ab = min(abs(a), abs(b))*sign(ONE_dp, a)
        else
            minmod_ab = ZERO_dp
        end if

    end subroutine

    subroutine limited_gradient_dx_vectorized(U_local, U_lower, U_upper, theta, gradient_dx, n)
        !!
        !! Slope limiter -- get the "gradient times dx" around U_local, given the upper and lower values of U.
        !!
        integer(ip), intent(in):: n !! Length of the input arrays
        real(dp), intent(in):: U_local(n) !! U_{i}
        real(dp), intent(in):: U_lower(n) !! U_{i-1}
        real(dp), intent(in):: U_upper(n) !! U_{i+1}
        real(dp), intent(in):: theta(n) !! Limiter parameter, often in [1-2].
        real(dp), intent(out) :: gradient_dx(n) !! Output

        character(len=charlen), parameter :: limiter_type = 'MC' !'Minmod2' !'Superbee_variant' ! 'MC'! 'nolimit' !
        !! Type of limiter

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
    ! This is the same as the above routine, but elemental so not necessarily vectorized.
    ! It allows trying a different loop structure elsewhere in the code.
    !
    pure elemental subroutine limited_gradient_dx(U_local, U_lower, U_upper, theta, gradient_dx)
        !!
        !! Slope limiter -- get the "gradient times dx" around U_local, given the upper and lower values of U.
        !!
        real(dp), intent(in):: U_local !! U_{i}
        real(dp), intent(in):: U_lower !! U_{i-1}
        real(dp), intent(in):: U_upper !! U_{i+1}
        real(dp), intent(in):: theta !! Limiter parameter, often in [1-2].
        real(dp), intent(out) :: gradient_dx !! Output

        real(dp):: a, b, c, d, e, th, sa, sb, half_sasb

        a = U_upper - U_local
        b = U_local - U_lower
        th = theta
        sa = sign(ONE_dp,a)
        sb = sign(ONE_dp,b)
        half_sasb = HALF_dp * (sa + sb)
        d = min(abs(a), abs(b)) * half_sasb * th ! Limit on local gradient
        e = HALF_dp * (a + b)
        c = merge(ZERO_dp, e, d == ZERO_dp) 
        ! NOTE: IF d /= 0, then clearly d, c have the same sign
        ! We exploit this to avoid a further minmod call (which seems
        ! expensive)
        gradient_dx = merge(min(c, d), max(c, d), d > ZERO_dp)

    end subroutine

    subroutine test_extrapolation_limiting_mod
        !! Unit tests
    
        integer(ip), parameter :: N = 10
        real(dp) :: U(N), U_lower(N), U_upper(N)
        real(dp) :: theta(N), extrapolation_sign(N)
        real(dp) :: desired_answer(N), answer(N)
        integer(ip) :: i, version

        ! Test both versions of the limited gradient extrapolation routine
        do version = 1, 2

            theta = 1.0_dp
            U = (/(i*1.0_dp - 5.5_dp, i = 1, N)/)
            U_lower = U - 1.0_dp
            U_upper = U + 1.0_dp

            ! Basic extrapolation tests, positive side
            extrapolation_sign=1.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = 0.5_dp * (U + U_upper)
            call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

            ! As above, negative side
            extrapolation_sign = -1.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = 0.5_dp * (U + U_lower)
            call assert_equal_within_tol(desired_answer, U + HALF_dp*extrapolation_sign*answer, __LINE__)

            ! As above, negative side, reduced theta
            theta = 0.5_dp
            extrapolation_sign = -1.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = 0.5_dp * U + 0.5_dp*0.5_dp * (U + U_lower)
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            ! As above, negative side, zero theta
            theta = 0.0_dp
            extrapolation_sign = 1.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = 1.0_dp * U + 0.0_dp*0.5_dp * (U + U_lower)
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            ! Local min.
            theta = 4.0_dp
            extrapolation_sign = 1.0_dp
            U_lower = U + 2.0_dp
            U_upper = U + 3.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = U 
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            ! Local max.
            theta = 4.0_dp
            U_lower = U - 2.0_dp
            U_upper = U - 3.0_dp
            extrapolation_sign = -1.0_dp
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = U 
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            ! Limited gradients
            ! 
            ! Central gradient = 2, lower-gradient=1, upper-gradient=3
            U_lower = U - 1.0_dp
            U_upper = U + 3.0_dp
            extrapolation_sign = 1.0_dp
            theta = 1.0_dp ! Limited-gradient = 1
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = U + 0.5_dp * 1.0_dp
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            theta = 1.5_dp ! Limited-gradient = 1.5
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = U + 0.5_dp * 1.5_dp
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)

            
            extrapolation_sign = -1.0_dp ! As above, negative side
            call test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
            desired_answer = U - 0.5_dp * 1.5_dp
            call assert_equal_within_tol(desired_answer, U+HALF_dp*extrapolation_sign*answer, __LINE__)
        end do

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

            ! Use this to test one or the other variant of the extrapolation routine
            subroutine test_both_limited_gradient_dx_vectorized(version, U, U_lower, U_upper, theta, answer, N)
                integer(ip), intent(in) :: N, version
                real(dp), intent(in) :: U(N), U_lower(N), U_upper(N)
                real(dp), intent(in) :: theta(N)
                real(dp), intent(inout) :: answer(N)

                if(version == 1) then
                    call limited_gradient_dx_vectorized(U, U_lower, U_upper, theta, answer, N)
                else if (version == 2) then
                    call limited_gradient_dx(U, U_lower, U_upper, theta, answer)
                else
                    print*, 'unknown version of limited_gradient_dx in test suite'
                    stop
                end if

            end subroutine

    end subroutine

end module
