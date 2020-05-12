module forcing_mod
    !!
    !! SWALS allows the user to define a subroutine "domain%forcing_subroutine" to 
    !! apply user-specified forcing terms. This module contains some subroutines that
    !! might be useful to use there.
    !!

    use global_mod, only: dp, ip
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    implicit none

    contains

    subroutine add_forcing_work_to_U_over_time(domain, time, dt, start_time, end_time)
        !!
        !! Add the field domain%forcing_work to domain%U. Do this smoothly over the time interval between start_time and end_time,
        !! such that the sum of the contributions over all time-steps in this interval adds up to domain%forcing_work. 
        !! A typical use-case is to perturb the stage and elevation with the co-seismic vertical deformation from an earthquake over
        !! some time-interval (or instantaneously if start_time = end_time)
        !!
        type(domain_type), intent(inout) :: domain        
            !! The domain to be updated. Note that it MUST have domain%forcing_work allocated to be the same
            !! shape as domain%U. Furthermore, domain%forcing_work must be populated with the total forcing 
            !! to be applied between start_time and end_time
        real(dp), intent(in) :: time
            !! Update the domain for a time-step from (time-dt, time). Note "time" is the END of the timestep
        real(dp), intent(in) :: dt
            !! Update the domain for a time-step from (time-dt, time). Note "time" is the END of the timestep
        real(dp), intent(in) :: start_time, end_time
            !! Apply the entire forcing at a steady rate between start_time and end_time

        real(dp) :: update_fraction, local_dt
        integer(ip) :: i, j, k

        if( ((time - dt) <= end_time) .and. (time >= start_time)) then
            ! The forcing is only applied between start_time and end_time - so
            ! we need to make sure the current time-step overlapped that period.
            ! (Note we are at the end of the timestep, and domain%time was already updated).

            ! Start with a few error checks
            if(.not. allocated(domain%forcing_work)) then
                write(log_output_unit, *) &
                    'Error: domain%forcing work must be allocated to use add_forcing_work_to_U_over_time'
                call generic_stop
            end if

            if(.not. all(shape(domain%forcing_work) == shape(domain%U))) then
                write(log_output_unit, *) &
                    'Error: domain%forcing work must have the same shape as domain%U'
                call generic_stop
            end if

            if(start_time > end_time) then
                write(log_output_unit, *) &
                    'Error: start_time must be > end_time in add_forcing_work_to_U_over_time'
                call generic_stop
            end if

            ! Determine what fraction of domain%forcing_work should be added to domain%U
            if(start_time == end_time) then
                update_fraction = 1.0_dp
            else
                ! The time-step might not be completely between start_time and end_time.
                ! This is the fraction of the forcing we actually need to apply 
                local_dt = min(time, end_time) - max((time-dt), start_time)
                update_fraction = local_dt/(end_time - start_time)
            end if
            !print*, update_fraction


            ! Apply the forcing
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, update_fraction)
            !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do k = 1, size(domain%U,3)
                do j = 1, domain%nx(2)
                    do i = 1, domain%nx(1)
                        domain%U(i,j,k) = domain%U(i,j,k) + update_fraction * domain%forcing_work(i,j,k)
                    end do
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        end if

    end subroutine


    subroutine test_forcing_mod
        !! Test the forcing subroutines by setting up a domain, and checking the forcing works

        real(dp), parameter :: global_lw(2) = [50.0_dp, 30.0_dp]
        real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp]
        integer(ip), parameter :: global_nx(2) = [50, 30]

        type(domain_type) :: domain
        real(dp), allocatable :: backup_U(:,:,:)
        integer(ip) :: i, j
        real(dp) :: err_tol, err, time, dt

        err_tol = spacing(10.0_dp)*10

        call domain%allocate_quantities(global_lw, global_nx, global_ll, create_output_files=.FALSE.)

        ! Make up domain%U
        do j = 1, domain%nx(2)
            do i = 1, domain%nx(1)
                domain%U(i,j,STG) = domain%x(i) + domain%y(j)
                domain%U(i,j,ELV) = domain%x(i) + domain%y(j) - 10.0_dp
                    ! bed is below stage
                domain%U(i,j,UH) = domain%x(i) 
                domain%U(i,j,VH) = domain%y(j) 
            end do
        end do

        ! Backup the initial condition for later usage
        allocate(backup_U(domain%nx(1), domain%nx(2), 4))
        backup_U = domain%U

        allocate(domain%forcing_work(domain%nx(1), domain%nx(2), 4))

        ! Make up the forcing
        domain%forcing_work(:,:,1) = 10.0_dp
        domain%forcing_work(:,:,2) = 5.0_dp
        domain%forcing_work(:,:,3) = 1.0_dp
        domain%forcing_work(:,:,4) = -1.0_dp

        ! Suppose we just evolved from 0.0 to 10.0
        time = 10.0_dp
        dt = 10.0_dp

        !
        ! Tests of add_forcing_work_to_U_over_time
        !

        ! Test 1 -- time-step exactly equals range of forcing
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = 0.0_dp, end_time = 10.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work)))
        if(err < err_tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__ 
        end if
        domain%U = backup_U
            ! Reset U for the next test

        ! Test 2 -- time-step starts before range of forcing
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = 4.0_dp, end_time = 10.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work)))
        if(err < err_tol ) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__
        end if
        domain%U = backup_U

        ! Test 3 -- time-step starts after range of forcing
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = -4.0_dp, end_time = 10.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work * 10.0_dp/14.0_dp)))
        if(err < err_tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__, &
               minval(domain%U - (backup_U + domain%forcing_work * 10.0_dp/14.0_dp)), &
               maxval(domain%U - (backup_U + domain%forcing_work * 10.0_dp/14.0_dp))
        end if
        domain%U = backup_U

        ! Test 4 -- time-step fully inside the forcing time
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = -15.0_dp, end_time = 15.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work * 10.0_dp/30.0_dp)))
        if(err < err_tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__, &
               minval(domain%U - (backup_U + domain%forcing_work * 10.0_dp/30.0_dp)), &
               maxval(domain%U - (backup_U + domain%forcing_work * 10.0_dp/30.0_dp))

        end if
        domain%U = backup_U

        ! Test 5 -- time-step overshoots the forcing end time
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = -15.0_dp, end_time = 5.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work * 5.0_dp/20.0_dp)))
        if(err < err_tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__, &
               minval(domain%U - (backup_U + domain%forcing_work * 5.0_dp/20.0_dp)), &
               maxval(domain%U - (backup_U + domain%forcing_work * 5.0_dp/20.0_dp))

        end if
        domain%U = backup_U

        ! Test 6 -- time-step doesn't overlap with the forcing 
        call add_forcing_work_to_U_over_time(domain, time, dt, start_time = -15.0_dp, end_time = -5.0_dp)
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work * 0.0_dp)))
        if(err == 0.0_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__, &
               minval(domain%U - (backup_U + domain%forcing_work * 0.0_dp)), &
               maxval(domain%U - (backup_U + domain%forcing_work * 0.0_dp))
        end if
        domain%U = backup_U

        ! Test 7 -- a timestepping type test, where we should end up adding the full contribution of forcing_work
        ! Use some erratic numbers -- it shouldn't matter so long as we fully step through start_time to end_time
        time = -10.01242_dp
        dt = 1.4444_dp
        do j = 1, 30
            ! When the loop is complete we will have covered the time
            call add_forcing_work_to_U_over_time(domain, time, dt, start_time = -1.332_dp, end_time = 5.15430_dp)
            time = time + dt
        end do
        err = maxval(abs(domain%U - (backup_U + domain%forcing_work * 1.0_dp)))
        if(err < err_tol) then
            print*, 'PASS'
        else
            print*, 'FAIL', err, err_tol, __LINE__, &
               minval(domain%U - (backup_U + domain%forcing_work)), &
               maxval(domain%U - (backup_U + domain%forcing_work))

        end if
        domain%U = backup_U


        deallocate(backup_U)

    end subroutine

end module
