module fortran_stage_vs_rate_curve

    use iso_c_binding, only: C_DOUBLE, C_INT
    implicit none

    ! These should match R's integer and double
    integer, parameter :: ip = C_INT
    integer, parameter :: dp = C_DOUBLE

    contains

    !
    ! This fortran subroutine performs a particular 'bottleneck' calculation,
    ! where we compute a stage-vs-exceedance-rate curve from a tabulated Mw-vs-exceedance-rate curve.
    !
    ! It's very specific. But the R interface can be more intuitive
    !
    subroutine make_stage_vs_rate_curve_site( rate_curve_Mw, rate_curve_exrate, N1, dMw, &
            unique_Mw, N2, event_Mw_to_unique_Mw_match, event_conditional_prob, N3, &
            sorted_event_stages, sorted_stage_inds, &
            output_stages, output_stage_exrates, N4) bind(C, name='make_stage_vs_rate_curve_site')

        ! Mw-vs-exceedance rate, tabulated
        real(dp), intent(in) :: rate_curve_Mw(N1), rate_curve_exrate(N1)
        integer(ip), intent(in) :: N1
        !
        ! The constant spacing between unique Mw values
        real(dp) :: dMw(1)
        !
        ! Unique Mw values in the scenarios, sorted (e.g. often 7.2, 7.3, ... 9.8)
        real(dp), intent(in) :: unique_Mw(N2)
        integer(ip), intent(in) :: N2
        !
        ! An index such that unique_Mw(event_Mw_to_unique_Mw_match(j)) will give the
        ! magnitude of scenario j
        integer(ip), intent(in) :: event_Mw_to_unique_Mw_match(N3)
        !
        ! The conditional probability of scenario j, given that an event with the same magnitude occurred
        real(dp), intent(in) :: event_conditional_prob(N3)
        integer(ip), intent(in) :: N3
        !
        ! The event stages, ALREADY SORTED in decreasing order. The array must be pre-sorted 
        ! for computational efficiency -- as typically this routine is called many times with 
        ! the same sorted_event_stages.
        ! In practice we would take the "event_stages" ordered by events, then in R, do:
        !     sorted_event_stages_and_indices = sort(event_stages, index.return=TRUE, decreasing=TRUE)
        !     sorted_event_stages = sorted_event_stages_and_indices$x
        !     sorted_stage_inds = sorted_event_stages_and_indices$ix
        real(dp), intent(in) :: sorted_event_stages(N3)
        ! As above, an index so that event_conditional_prob(sorted_stage_inds(j)) gives the 
        ! conditional probability of the scenario corresponding to sorted_event_stages(j). 
        integer(ip), intent(in) :: sorted_stage_inds(N3)
        !
        ! The stages at which we want to compute the exceedance rates
        real(dp), intent(in) :: output_stages(N4)
        ! An array to return the exceedance rates
        real(dp), intent(inout) :: output_stage_exrates(N4)
        integer(ip), intent(in) :: N4
        
        ! Local variables
        real(dp) :: max_mw_max, unique_rates(N2), rate_curve_approx(N2+1), rate_curve_approx_Mws(N2+1)
        real(dp) :: scenario_rates(N3), sorted_exrates_with_cap(N3+2), sorted_stages_with_cap(N3+2)
        integer(ip) :: i

        !
        ! Below in comments, we show the R code that was translated into Fortran.
        !

        !# In the PTHA18 'all rate curves' data structure, this holds
        !max_mw_max = max(rate_curve_mw)
        max_mw_max = rate_curve_mw(N1)

        !dMw = 0.1 # Size of magnitude bins
        !# Individual rates for each scenario = (conditional-probability) *
        !# (incremental-rate)
        !#
        !# This calculation method (where we zero rates when Mw > Mw-max) is
        !# identical to the treatment in 'rptha'. This circumvents any
        !# interpolation issues -- ( Mw > Mw-max) always has zero rate.
        !f1 = approxfun(rate_curve_mw, rate_curve_exrate, rule=2, ties='ordered')
        !unique_rates = 
        !    (f1(unique_Mw - dMw/2) * (unique_Mw - dMw/2 <= max_mw_max) - 
        !     f1(unique_Mw + dMw/2) * (unique_Mw + dMw/2 <= max_mw_max) )
        rate_curve_approx_Mws(1:N2) = unique_Mw(1:N2) - dMw(1)*0.5_dp
        rate_curve_approx_Mws(N2+1) = unique_Mw(N2) + dMw(1)*0.5_dp
        call linear_interpolation(N1, rate_curve_Mw, rate_curve_exrate, N2+1, rate_curve_approx_Mws, rate_curve_approx)
        where(rate_curve_approx_Mws > max_mw_max) rate_curve_approx = 0.0_dp
        unique_rates(1:N2) = rate_curve_approx(1:N2) - rate_curve_approx(2:(N2+1))

        !scenario_rates = p2_event_cp * unique_rates[event_Mw_to_unique_Mw_match]
        do i = 1, N3 
            scenario_rates(i) = event_conditional_prob(i) * unique_rates(event_Mw_to_unique_Mw_match(i))
        end do

        ! The stages were already sorted. Just add the 'cap' so we can neatly treat zero exceedance-rates
        sorted_stages_with_cap(1) = HUGE(1.0_dp)
        sorted_stages_with_cap(2) = sorted_event_stages(1) + 1.0e-03_dp ! Same 1.0e-03 constant used in the R code
        sorted_stages_with_cap(3:(N3+2)) = sorted_event_stages

        ! sorted_exrates_with_cap = cumsum( c(0, 0, scenario_rates[sorted_stages_and_inds$ix]) )
        sorted_exrates_with_cap(1:2) = 0.0_dp
        do i = 1, N3
            sorted_exrates_with_cap(i+2) = sorted_exrates_with_cap(i+1) + scenario_rates(sorted_stage_inds(i))
        end do

        ! The fortran interpolator requires that sorted_stages_with_cap is in increasing order
        sorted_stages_with_cap = sorted_stages_with_cap((N3+2):1:-1)
        sorted_exrates_with_cap = sorted_exrates_with_cap((N3+2):1:-1)

        ! If we have any equal stages, then make their exrates equal as well (corresponding to the largest value).
        ! This is equivalent to the R treatment (which would employ 'ties = max' in 'approx')
        do i = 2, N3 ! No need to touch indices (N3+1, N3+2), which are "the cap"
            if(sorted_stages_with_cap(i) == sorted_stages_with_cap(i-1)) then
                sorted_exrates_with_cap(i) = sorted_exrates_with_cap(i-1)
            end if
        end do

        !stage_exceedance_rates = approx(sorted_stages_with_cap, sorted_exrates_with_cap, xout=stages, rule=2:1, ties=max)$y
        call linear_interpolation(N3+2, sorted_stages_with_cap, sorted_exrates_with_cap, N4, &
            output_stages, output_stage_exrates)

    end subroutine

    !
    ! Below are a couple of routines which support linear interpolation -- copies from
    ! ptha/propagation/SWALS/src/util/linear_interpolation_mod.f90
    !

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
        real(dp), intent(inout):: output_y(n_output)

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

end module
