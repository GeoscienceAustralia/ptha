module tsunami_arrival_time_mod
    !!
    !! Computing tsunami arrival times, maxima, and exceedance-rates for 
    !! tsunami scenarios with "maxima >= wave_threshold" and "arrival_time <= time_threshold"
    !!
    use iso_c_binding, only: C_DOUBLE, C_INT
    implicit none

    integer, parameter :: dp = C_DOUBLE, ip = C_INT
        !! Control precision of reals and integers by using real(dp) and integer(ip)

    contains

    subroutine tsunami_arrival_time_and_maxima_f(n, stage, time, msl, arrival_fraction_of_maxima, &
        tsunami_maxima, tsunami_arrival_time) bind(C, name='tsunami_arrival_time_and_maxima_fortran')
        !!
        !! Given a tsunami time-series, and a definition of the mean-sea-level and arrival time 
        !! (determined as a fraction of the tsunami maxima), compute the tsunami arrival time and maxima.
        !!
        integer(ip), intent(in) :: n
            !! Length of stage(:) and time(:)
        real(dp), intent(in) :: stage(n)
            !! Array of tsunami stage values (in m). The tsunami waveform is taken to be (stage - msl).
        real(dp), intent(in) :: time(n)
            !! Array of times (seconds post-earthquake) with length equal to stage
        real(dp), intent(in) ::  msl(1)
            !! Value of mean sea level (typically 0.0_dp), used to define the tsunami waveform as (stage - msl)
        real(dp), intent(in) ::  arrival_fraction_of_maxima(1)
            !! The arrival time is defined as the first time that (stage - msl) >= (arrival_fraction_of_maxima * tsunami_maxima)
        real(dp), intent(out) :: tsunami_maxima(1)
            !! Output tsunami maxima ( maxval(stage) - msl )
        real(dp), intent(out) :: tsunami_arrival_time(1)
            !! Output the first time when the (stage - msl) >= (arrival_fraction_of_maxima * tsunami_maxima) 

        real(dp) :: arrival_stage_threshold
        integer(ip) :: i, arrival_index
        
        if((arrival_fraction_of_maxima(1) > 1.0_dp) .or. (arrival_fraction_of_maxima(1) < 0.0_dp)) then
            !stop 'arrival_fraction_of_maxima must be within [0.0, 1.0]'
            tsunami_maxima(1) = -HUGE(1.0_dp)
            tsunami_arrival_time(1) = -HUGE(1.0_dp)
            return
        end if
        
        tsunami_maxima(1) = maxval(stage) - msl(1)

        ! The tsunami arrival is defined by the stage exceedng this value
        arrival_stage_threshold = arrival_fraction_of_maxima(1) * tsunami_maxima(1) + msl(1)

        ! Find the index when the tsunami arrives
        do i = 1, size(stage)
            if( stage(i) >= arrival_stage_threshold) then
                arrival_index = i
                exit
            end if
        end do
       
        !if(size(time) == 1) then
        !    ! Case where time(1) gives the timestep
        !    tsunami_arrival_time = (arrival_index - 1)*time(1)
        !else if(size(time) == size(stage)) then
            ! Case where time(:) is an array of times 
            tsunami_arrival_time(1) = time(arrival_index)
        !else
        !    stop ' "time" must either be an array of times with length equal to "stage", ' // &
        !         'or an array of length=1 giving the timestep' 
        !endif

    end subroutine


    subroutine exceedance_rate_given_maxima_and_arrival_time_f(ns, scenario_maxima, scenario_arrival_times, scenario_rates, &
            nmax, tsunami_maxima_in_output, narrival, arrival_times_in_output, output_exceedance_rate_by_maxima_time) &
            bind(C, name='exceedance_rate_given_maxima_and_arrival_time_fortran')
            !!
            !! Suppose we have a set of scenarios with known maxima, arrival_times, and occurrance rates.
            !! This routine computes a matrix of exceedance_rates, for every combination of "tsunami_maxima_in_output"
            !! and "arrival_times_in_output".
            !!
        integer(ip), intent(in) :: ns
            !! Number of scenarios, size of scenario_maxima, scenario_arrival_times, scenario_rates
        real(dp), intent(in) :: scenario_maxima(ns)
            !! The tsunami maxima for each scenario 
        real(dp), intent(in) :: scenario_arrival_times(ns)
            !! The arrival time for each scenario
        real(dp), intent(in) :: scenario_rates(ns)
            !! The occurrance rate of each individual scenario (events/year). So sum(scenario_rates) gives the rate of any scenario,
            !! and (scenario_rates)/sum(scenario_rates) gives the scenario conditional probabilitiy
        integer(ip), intent(in) :: nmax
            !! size of the tsunami_maxima_in_output array
        real(dp), intent(in) :: tsunami_maxima_in_output(nmax)
            !! The tsunami maxima for which we wish to compute exceedance-rates for the given arrival times
        integer(ip), intent(in) :: narrival
            !! size of the arrival_times_in_output array 
        real(dp), intent(in) :: arrival_times_in_output(narrival)
            !! The arrival times for which we wish to compute nonexceedance-rates for the given tsunami maxima
        real(dp), intent(out) :: output_exceedance_rate_by_maxima_time(nmax, narrival)
            !! On output, the i,j entry of this matrix will contain 
            !!    sum(scenario_rates, &
            !!        mask=( (scenario_maxima >= tsunami_maxima_in_output(i)) &
            !!               (scenario_arrival_times <= arrival_times_in_output(j)) )
            !! NOTE the use of ">=" for scenario_maxima, but "<=" for arrival times. The idea is that large maxima, and small
            !! arrivals, are both "bad". Then this sign convention might be more useful (?). 

        integer(ip) :: i, j
       
        !if(any(shape(output_exceedance_rate_by_maxima_time) /= &
        !     [size(tsunami_maxima_in_output), size(arrival_times_in_output)])) then
        !     stop ' Wrong shape for output exceedance_rate matrix '
        ! end if

        !if((size(scenario_maxima) /= size(scenario_arrival_times)) .or. &
        !   (size(scenario_rates)  /= size(scenario_arrival_times)) ) then
        !    stop 'scenario_maxima, scenario_rates, and scenario_arrival_times must all have the same size'
        !end if

        ! Compute the exceedance-rate for every pair of tsunami_maxima_in_output, arrival_times_in_output
        do i = 1, size(tsunami_maxima_in_output)
            do j = 1, size(arrival_times_in_output)
                ! This calculation might be sped-up with clever pre-sorting of data. But I doubt it will
                ! matter in the applications of interest.
                 output_exceedance_rate_by_maxima_time(i,j) = sum(scenario_rates, &
                     mask = ( (scenario_maxima >= tsunami_maxima_in_output(i)) .and. &
                              (scenario_arrival_times <= arrival_times_in_output(j))))
            end do
        end do


    end subroutine

end module
