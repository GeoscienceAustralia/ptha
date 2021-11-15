module date_to_numeric_mod
    !!
    !! Convert dates and times to/from Julian days
    !!

    use global_mod, only: ip
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use iso_c_binding, only: C_DOUBLE

    implicit none

    ! Here we always use double-precision reals
    integer, parameter, private :: rp = C_DOUBLE

    integer, parameter :: julian_day_1970_01_01 = 2440588
        !! the date '1970-01-01' as a julian day
    integer, parameter :: julian_day_1992_01_01 = 8035 + julian_day_1970_01_01
        !! the date '1992-01-01' as a julian day

    integer, parameter :: days_in_month_regular_year(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    integer, parameter :: days_in_month_leap_year(12)    = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    contains

    function is_leap_year(yyyy)
        !! Is year yyyy a leap-year?
        integer(ip), intent(in) :: yyyy
        logical :: is_leap_year
        ! Leap years -- note centuries are only leap years if they divide 400, e.g.
        ! 2100, 2200, 2300 are not leap-years, but 2400 is.
        is_leap_year = (mod(yyyy, 4) == 0) .and. ((mod(yyyy, 100) /= 0) .or. (mod(yyyy, 400) == 0))
    end function

    subroutine check_datetime(yyyy, mm, dd, hh, mi, ss)
            !! Check that a date of the form yyyy/mm/dd/hh/mi/ss is valid
            !! The hh/mi/ss are optional.
            !! Because this can throw errors, it is not elemental
        integer(ip), intent(in) :: yyyy, mm, dd
            !! Year/month/day
        integer(ip), optional, intent(in) :: hh, mi, ss
            !! Hour/minute/second

        integer :: max_days

            if(mm < 1 .or. mm > 12) then
                write(log_output_unit,*) 'Month must be from 1-12, but got ', mm
                call generic_stop
            end if

           ! Figure out how many days in the month
            if(is_leap_year(yyyy)) then
                max_days = days_in_month_leap_year(mm)
            else
                max_days = days_in_month_regular_year(mm)
            end if

            if(dd < 1 .or. dd > max_days) then
                write(log_output_unit, *) 'The number of days in month ', mm , &
                    ' should be between 1 and ', max_days, ', but I got ', dd
                call generic_stop
            end if

            if(present(hh)) then
                if(hh < 0 .or. hh > 23) then
                    write(log_output_unit, *) 'The hour should range in 0 to 23, but I got ', hh
                    call generic_stop
                end if
            end if

            if(present(mi)) then
                if(mi < 0 .or. mi > 59) then
                    write(log_output_unit, *) 'The minute should range in 0 to 59, but I got ', mi
                    call generic_stop
                end if
            end if

            if(present(ss)) then
                if(ss < 0 .or. ss > 59) then
                    write(log_output_unit, *) 'The second should range in 0 to 59, but I got ', ss
                    call generic_stop
                end if
            end if

    end subroutine

    elemental subroutine ymd_to_julian_day(yyyy, mm, dd, julian)
        !! Converts a calendar date to a Julian date (i.e. integer number of days since start of Julian calendar)
        !! Fliegel and Van Flandern (1968) CACM 11(10):657
        !! Various discussions of this online, e.g.
        !! https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/271933
        !! https://stackoverflow.com/questions/17354100/date-time-difference-in-fortran
        integer, intent(in) :: yyyy, mm, dd
            !! year, month, day
        integer, intent(out) :: julian
            !! Number of days since start of julian calendar

        julian = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + 367*(mm-2-((mm-14)/12)*12)/12- &
            3*((yyyy + 4900 + (mm - 14)/12)/100)/4

    end subroutine

    elemental subroutine julian_day_to_ymd (jd, yyyy, mm, dd)
        !! Convert a julian date (days since start of julian calendar) into year/month/day
        !! Fliegel and Van Flandern (1968) CACM 11(10):657
        !! See discussion at:
        !! https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/271933
        integer,intent(in) :: jd
            !! Number of days since start of julian calendar
        integer,intent(out) :: yyyy,mm,dd
            !! year, month, day

        integer :: l,n

        l = jd + 68569
        n = 4*l/146097
        l = l - (146097*n + 3)/4
        yyyy = 4000*(l + 1)/1461001
        l = l - 1461*yyyy/4 + 31
        mm = 80*l/2447
        dd = l - 2447*mm/80
        l = mm/11
        mm = mm + 2 - 12*l
        yyyy = 100*(n - 49) + yyyy + l

    end subroutine


    elemental subroutine datetime_to_days(yyyy, mm, dd, hh, mi, ss, days, julian_origin)
        !! Get the "real number of days" since an origin date (default origin = 1970-01-01)
        !!
        integer(ip), intent(in) :: yyyy, mm, dd, hh, mi, ss
            !! Date in yyyy/mm/dd/hh/mm/ss format, where the hours is in 24hr format
        real(rp), intent(out) :: days
            !! Real "number of days since the julian_origin"
        integer(ip), optional, intent(in) :: julian_origin
            !! Optional origin for 'days' (as an integer number of days). Defaults to the
            !! julian calendar date of 1970-01-01 (=2440588) -- like R's julian() function

        integer(ip):: yyyy_l, mm_l, dd_l, hh_l, mi_l, ss_l, jd_l, jo_l

        if(present(julian_origin)) then
            jo_l = julian_origin
        else
            jo_l = julian_day_1970_01_01
        end if

        call ymd_to_julian_day(yyyy, mm, dd, jd_l)

        days = (jd_l - jo_l)

        days = days + (1.0_rp/24.0_rp)*hh + (1.0_rp/(24.0_rp*60.0_rp))*mi + &
            (1.0_rp/(24.0_rp*3600.0_rp))*ss

    end subroutine

    elemental subroutine days_to_datetime(days, yyyy, mm, dd, hh, mi, ss, julian_origin)
        !! Get the date (yyyy/mm/dd HH:MI:SS) for a given 'real number of days' from an origin date (default origin = 1970-01-01)
        !!
        real(rp), intent(in) :: days
            !! Real "number of days since the julian_origin"
        integer(ip), intent(out) :: yyyy, mm, dd, hh, mi, ss
            !! Date in yyyy/mm/dd/hh/mm/ss format, where the hours is in 24hr format
        integer(ip), optional, intent(in) :: julian_origin
            !! Optional origin for 'days' (as an integer number of days). Defaults to the
            !! julian calendar date of 1970-01-01 (=2440588) -- like R's julian() function

        integer(ip):: jo_l
        real(rp) :: fractional_days
        integer(ip) :: integer_days

        if(present(julian_origin)) then
            jo_l = julian_origin
        else
            jo_l = julian_day_1970_01_01
        end if

        integer_days = floor(days) + jo_l
        fractional_days = days - floor(days)

        if(1.0_rp - fractional_days <= 0.5_rp/(3600.0_rp*24.0_rp)) then
            ! We are within 0.5 seconds of having an extra day. This will round up
            integer_days = integer_days + 1
            fractional_days = 0.0_rp
        end if

        call julian_day_to_ymd(integer_days, yyyy, mm, dd)

        hh = floor(fractional_days*24.0_rp)
        mi = floor( (fractional_days*24.0_rp - hh)*60.0_rp )
        ss = nint(fractional_days*24.0_rp*60.0_rp*60.0_rp - hh*60.0_rp*60.0_rp - mi*60.0_rp)

        ! Deal with rounding
        if(ss == 60) then
            ss = 0
            mi = mi + 1
        end if

        if(mi == 60) then
            mi = 0
            hh = hh+1
        end if

    end subroutine



    subroutine test_date_to_numeric_mod

        integer(ip) :: yyyy, mm, dd, jd
        integer(ip) :: yyyy2(2), mm2(2), dd2(2), jd2(2)  ! Test vectorized versions
        integer(ip) :: hh2(2), mi2(2), ss2(2)
        integer(ip) :: yyyy2b(2), mm2b(2), dd2b(2), jd2b(2)  ! Test vectorized versions
        integer(ip) :: hh2b(2), mi2b(2), ss2b(2)
        real(rp) :: days(2), expected_days(2)

        ! Check leap years -- noting centuries that do not divide 400 are not leap years
        if(is_leap_year(1992) .and. is_leap_year(2000) .and. &
            (.not. is_leap_year(2100)) .and. is_leap_year(2400)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Check start of 1970
        yyyy = 1970; mm = 01; dd = 01
        call ymd_to_julian_day(yyyy, mm, dd, jd)

        if(jd == julian_day_1970_01_01) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__
        end if

        ! Check start of 1992
        yyyy = 1992; mm = 01; dd = 01
        call ymd_to_julian_day(yyyy, mm, dd, jd)

        if(jd == julian_day_1992_01_01) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__
        end if

        ! Check vectorized versions
        yyyy2 = [1970, 1992]; mm2 = 01; dd2= 01
        call ymd_to_julian_day(yyyy2, mm2, dd2, jd2)
        if(all(jd2 == [julian_day_1970_01_01, julian_day_1992_01_01])) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__
        end if

        ! Check inverse
        call julian_day_to_ymd(jd2, yyyy2, mm2, dd2)
        if(all(yyyy2 == [1970, 1992] .and. mm2 == 01 .and. dd2 == 01)) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__
        end if


        ! Check some ad-hoc julian days, with comparison values from R's julian
        yyyy2 = [2020, 1900]
        mm2 = [07, 02]
        dd2 = [20, 15]
        hh2 = [16, 23]
        mi2 = [20, 59]
        ss2 = [31, 59]
        call datetime_to_days(yyyy2, mm2, dd2, hh2, mi2, ss2, days)
        expected_days = [18463.680914351851243_rp, -25521.000011574073142_rp] ! From R's julian()
        ! Should agree to nearest second
        if(all(abs(days - expected_days) < 0.5_rp/(24.0_rp*3600.0_rp))) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__, days, (days - expected_days)
        end if

        ! Check the inverse
        call days_to_datetime(days, yyyy2b, mm2b, dd2b, hh2b, mi2b, ss2b)
        if(all(yyyy2b == yyyy2 .and. mm2b == mm2 .and. dd2b == dd2 .and. hh2b == hh2 .and. mi2b == mi2 .and. ss2b == ss2)) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__, &
                yyyy2b, mm2b, dd2b, hh2b, mi2b, ss2b
        end if

        ! Check the dates
        call check_datetime(yyyy2(1), mm2(1), dd2(1), hh2(1), mi2(1), ss2(1))
        call check_datetime(yyyy2(1), mm2(1), dd2(1))
        call check_datetime(yyyy2(2), mm2(2), dd2(2), hh2(2), mi2(2), ss2(2))
        call check_datetime(yyyy2(2), mm2(2), dd2(2))

        ! More ad-hoc checks -- large/small dates
        yyyy2 = [3000, 0001]
        call datetime_to_days(yyyy2, mm2, dd2, hh2, mi2, ss2, days)
        expected_days = [376400.68091435183305_rp, -719116.00001157412771_rp]

        ! Should agree to nearest second
        if(all(abs(days - expected_days) < 0.5_rp/(24.0_rp*3600.0_rp))) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__, days, (days - expected_days)
        end if

        ! Check rounding -- take a day within < 1/2 second of 1970-01-01.
        days = julian_day_1970_01_01 - 0.4_rp/(3600.0_rp * 24.0_rp)
        call days_to_datetime(days, yyyy2, mm2, dd2, hh2, mi2, ss2, julian_origin = [0, 0])
        if(all(yyyy2 == 1970 .and. dd2 == 01 .and. mm2 == 01 .and. hh2 == 00 .and. mi2 == 00 .and. ss2 == 00)) then
            print*, 'PASS'
        else
            print*, 'FAIL', __LINE__, yyyy2, dd2, mm2, hh2, mi2, ss2
        end if

    end subroutine

end module
