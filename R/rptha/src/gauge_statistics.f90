module gauge_statistics_mod
    !
    use iso_c_binding
    implicit none

    integer, parameter:: dp = C_DOUBLE, ip = C_INT
    integer, parameter:: num_gauge_statistics = 5

    contains

    ! Compute multi-gauge summary statistics
    !
    subroutine gauge_statistics(gauge_times, flow_var, &
        stage_threshold_for_arrival_time, &
        ngauges, ntimes, nflow_var, output) BIND(C, name='gauge_statistics')

        real(dp), intent(in) :: gauge_times(ntimes), stage_threshold_for_arrival_time
        integer(ip), intent(in) :: ngauges, ntimes, nflow_var
        real(dp), intent(in) :: flow_var(ngauges, ntimes, nflow_var)
        real(dp), intent(inout) :: output(ngauges, num_gauge_statistics)

        integer(ip) :: i, j, k
        !real(dp) :: stage(ntimes) ! Small chance of stack overflow
        real(dp), allocatable :: stage(:), work(:)
        allocate(stage(ntimes), work(ntimes))

        do i = 1, ngauges
            ! It is assumed that flow_var(:,:,1) contains the stage time_series
            ! at every gauge
            stage = flow_var(i, :, 1)

            ! Maximum stage
            output(i,1) = maxval(stage)
        
            ! Zero-crossing period
            output(i,2) = zero_crossing_period(gauge_times, stage, work, ntimes)

            ! Peak to trough
            output(i,3) = output(i,1) - minval(stage)

            ! Arrival time
            output(i,4) = -1.0 ! 'Later R will treat this as missing value'
            do j = 1, ntimes
                if(abs(stage(j)) > stage_threshold_for_arrival_time) then
                    output(i,4) = gauge_times(j)
                    ! Break out of the loop
                    exit
                end if
            end do

            ! Initial stage
            output(i,5) = stage(1)

        end do
    end subroutine

    ! Compute the zero crossing period of vector 'x' measured at times 't'
    !
    ! Returns -1 if there are not enough crossings to compute this [we
    ! need at least 2 up crossings and 2 down crossings]
    !
    ! @param t times -- must be monotonic increasing
    ! @param x values at time t
    ! @param work a real work array (to avoid repeated allocation)
    ! @param ntimes length of t, x, and work
    ! @return average of up-crossing and down-crossing periods
    !
    function zero_crossing_period(t, x, work, ntimes) result(zcp)
        integer(ip), intent(in) :: ntimes
        real(dp), intent(in) :: t(ntimes), x(ntimes)
        real(dp), intent(inout) :: work(ntimes)
        real(dp) :: zcp

        integer(ip) :: i, nup, ndown
        real(dp) :: first_up_crossing, first_down_crossing
        real(dp) :: last_up_crossing, last_down_crossing
        logical :: found_first_up_crossing, found_first_down_crossing
        real(dp) :: s01, t_when_x_equals_zero 


        ! Sign of x
        work(1) = 0.0_dp 
        do i = 2, ntimes
            ! positive --> up crossing
            ! negative --> down crossing
            ! 0 --> no sign change
            !work(i) = sign(1.0_dp, x(i)) - sign(1.0_dp, x(i-1))
            if(x(i) > 0.0_dp .and. x(i-1) <= 0.0_dp) then
                work(i) = 1.0_dp
            else
                if(x(i) < 0.0_dp .and. x(i-1) >= 0.0_dp) then
                    work(i) = -1.0_dp
                else
                    work(i) = 0.0_dp
                end if
            end if
        end do

        found_first_up_crossing = .FALSE.
        found_first_down_crossing = .FALSE.
        nup = 0
        ndown = 0
        do i = 2, ntimes
            if(work(i) == 0.0_dp) cycle
           
            ! A number in [0-1] giving the position of the zero between
            ! t(i-1) and t(i) 
            s01 = -x(i-1)/(x(i)-x(i-1))
            t_when_x_equals_zero = t(i-1) + s01 * (t(i)-t(i-1))

            if(work(i) > 0.0_dp) then
                ! Up crossing
                nup = nup + 1
                last_up_crossing = t_when_x_equals_zero
                if(.not. found_first_up_crossing) then
                    first_up_crossing = last_up_crossing
                    found_first_up_crossing = .TRUE.
                endif
            else
                ! down crossing
                ndown = ndown + 1
                last_down_crossing = t_when_x_equals_zero 
                if(.not. found_first_down_crossing) then
                    first_down_crossing = last_down_crossing
                    found_first_down_crossing = .TRUE.
                endif
            end if
                
        end do

        !! DEBUG
        !print*, 'dc: ', ndown, last_down_crossing, first_down_crossing
        !print*, 'uc: ', nup, last_up_crossing, first_up_crossing

        if( (ndown < 2).or.(nup < 2) ) then
            zcp = -1.0_dp
        else
            zcp = 0.5_dp * ( &
                (last_down_crossing - first_down_crossing)/(ndown-1) + &
                (last_up_crossing - first_up_crossing)/(nup-1) )
        end if

    end function

end module
