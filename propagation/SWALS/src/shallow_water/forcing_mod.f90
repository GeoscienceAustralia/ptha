module forcing_mod
    !!
    !! SWALS allows the user to define a subroutine "domain%forcing_subroutine" to 
    !! apply user-specified forcing terms. In general that can be defined by the user in their 
    !! program. However, there are some common forcing cases, e.g., often we want to add a spatially variable source
    !! over some specified time-interval. This module contains some subroutines and
    !! derived_types to help with such cases.
    !!

    use global_mod, only: dp, ip, charlen
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use iso_c_binding, only: c_loc, c_f_pointer
    implicit none

    private
    public add_forcing_work_to_U_over_time, &
        forcing_patch_type, apply_forcing_patch, &
        forcing_patch_array_type, apply_forcing_patch_array, &
        test_forcing_mod


    type forcing_patch_type
        !! Suppose we want to apply a forcing over some rectangular subset of domain%U,
        !! say domain%U(i0:i1, j0:j1, k0:k1). The forcing is applied between some specified start_time and end_time. 
        !! This data structure can help with that
        real(dp), allocatable :: forcing_work(:,:,:)
            !! This is added to domain%U(i0:i1, j0:j1, k0:k1) over the specified time interval
        integer(ip) :: i0 = HUGE(1_ip), i1=-HUGE(1_ip), j0=HUGE(1_ip), j1=-HUGE(1_ip), k0=STG, k1=ELV
            !! An index denoting where the forcing is applied (i.e. to domain%U(i0:i1, j0:j1, k0:k1). 
            !! The default value will cause an error if used in apply_forcing_patch
        real(dp) :: start_time = HUGE(1.0_dp), end_time = -HUGE(1.0_dp)
            !! Controls the time-interval over which the forcing is applied. Must have start_time >= end_time.
            !! The default values will throw an error if used in apply_forcing_patch
        contains

        procedure:: setup => setup_forcing_patch
            !! Allocate forcing_patch%forcing_work and set start_time, end_time, and spatial indices
        procedure:: finalise => deallocate_forcing_patch
            !! Deallocate forcing_patch%forcing_work and revert to default start_time, end_time, and spatial indices
    end type

    type forcing_patch_array_type
        !! Type to hold an array of forcing patches. This could be used e.g. to simulate an earthquake
        !! where each unit-source had a different start-time and rise time.
        type(forcing_patch_type), allocatable :: forcing_patches(:)
            !! Hold multiple forcing_patches. 
    end type


    contains

    subroutine add_forcing_work_to_U_over_time(U, forcing_work, time, dt, start_time, end_time)
        !!
        !! Add the field forcing_work(:,:,:) to U(:,:,:) (the latter is typically domain%U). Do this smoothly over the time interval
        !! between start_time and end_time, such that the sum of the contributions over all time-steps in this interval adds up to
        !! domain%forcing_work.  A typical use-case is to perturb the stage and elevation with the co-seismic vertical deformation
        !! from an earthquake over some time-interval (or instantaneously if start_time = end_time)
        !!
        real(dp), intent(inout) :: U(:,:,:)
            !! The field to be updated, typically domain%U
        real(dp), intent(in) :: forcing_work(:,:,:)
            !! The forcing -- it MUST have the same shape as U. Furthermore, forcing_work must be populated 
            !! with the total forcing to be applied between start_time and end_time
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
            ! (Note we are at the end of the timestep, and time was already updated).

            if(.not. all(shape(forcing_work) == shape(U))) then
                write(log_output_unit, *) &
                    'Error: forcing work must have the same shape as U'
                call generic_stop
            end if

            ! !This if-statement must always be .FALSE. given the outer if
            !if(start_time > end_time) then
            !    write(log_output_unit, *) &
            !        'Error: start_time must be > end_time in add_forcing_work_to_U_over_time'
            !    call generic_stop
            !end if

            ! Determine what fraction of forcing_work should be added to U
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
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(U, forcing_work, update_fraction)
            !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do k = 1, size(U,3)
                do j = 1, size(U, 2)
                    do i = 1, size(U, 1)
                        U(i,j,k) = U(i,j,k) + update_fraction * forcing_work(i,j,k)
                    end do
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        end if

    end subroutine

    subroutine apply_forcing_patch_base(domain, dt, fp)
        !!
        !! Apply a forcing patch to the domain over the given timestep.
        !! This is called by various other subroutines that use the forcing_patch_type or forcing_patch_array_type,
        !! so is separated here. Ultimately it calls add_forcing_work_to_U_over_time, which does most of the work.
        !!
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt
        type(forcing_patch_type), intent(in) :: fp

        integer(ip) :: i0, i1, j0, j1, k0, k1

        ! Unpack the indices where the forcing applies
        i0 = fp%i0
        i1 = fp%i1
        j0 = fp%j0
        j1 = fp%j1
        k0 = fp%k0
        k1 = fp%k1

        ! Do a few sanity checks
        if(i0 > i1 .or. j0 > j1) then
            write(log_output_unit, *) 'Error in apply_forcing_patch: the forcing_patch does not have indices set correctly'
            call generic_stop
        end if
        if(fp%start_time > fp%end_time) then
            write(log_output_unit, *) 'Error in apply_forcing_patch: the forcing_patch start/end times are not set correctly'
            call generic_stop
        end if

        if(.not. ( (k0 <= ELV) .and. (k0 >= STG) .and. (k1 <= ELV) .and. (k0 <= k1))) then
            write(log_output_unit, *) 'Error in apply_forcing_patch_base: the forcing_patch%forcing_work does not'
            write(log_output_unit, *) 'have 3rd dimension compatible with domain%U', k0, k1
            call generic_stop
        end if

        if(k1 >= ELV .and. k0 <= ELV .and. (.not. domain%support_elevation_forcing)) then
            write(log_output_unit, *) 'Error in apply_forcing_patch_base: It appears the forcing_patch is setup to'
            write(log_output_unit, *) 'change the elevation. However, this is not allowed because'
            write(log_output_unit, *) 'the current timestepping method has domain%support_elevation_forcing=.FALSE.'
            write(log_output_unit, *) 'You can try manually setting the domain%support_elevation_forcing=.TRUE. '
            write(log_output_unit, *) 'before allocating the domain. This will work IF the numerical scheme is able'
            write(log_output_unit, *) 'to optionally support elevation forcing (by doing extra work).'
            write(log_output_unit, *) 'domain%myid: ', domain%myid
            write(log_output_unit, *) 'domain%timestepping_method: ', trim(domain%timestepping_method)
            call generic_stop
        end if

        ! All seems well so apply the forcing. 
        call add_forcing_work_to_U_over_time(domain%U(i0:i1, j0:j1, k0:k1), fp%forcing_work, domain%time, &
            dt, fp%start_time, fp%end_time)

    end subroutine

    subroutine apply_forcing_patch(domain, dt)
        !!
        !! Convenience subroutine. If domain%forcing_context_cptr contains a forcing_patch, then this subroutine can be
        !! passed to domain%forcing_subroutine to apply the forcing. Ultimately it will call add_forcing_work_to_U_over_time,
        !! taking care of the spatial subsetting, with some error checking.
        !!
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt

        type(forcing_patch_type), pointer :: fp

        call c_f_pointer(domain%forcing_context_cptr, fp)

        call apply_forcing_patch_base(domain, dt, fp)

        fp => NULL()

    end subroutine

    subroutine apply_forcing_patch_array(domain, dt)
        !!
        !! Convenience subroutine. If domain%forcing_context_cptr contains a forcing_patch_array_type, then this subroutine can be
        !! passed to domain%forcing_subroutine to apply the forcing for every forcing patch it contains. Ultimately it will
        !! repeatedly call add_forcing_work_to_U_over_time, once for each forcing_patch, taking care of the spatial subsetting, 
        !! with some error checking.
        !!
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt

        type(forcing_patch_array_type), pointer :: fpat
        integer(ip) :: i

        call c_f_pointer(domain%forcing_context_cptr, fpat)

        do i = 1, size(fpat%forcing_patches)
            call apply_forcing_patch_base(domain, dt, fpat%forcing_patches(i))
        end do

        fpat => NULL()

    end subroutine

    subroutine setup_forcing_patch(forcing_patch, start_time, end_time, i0, i1, j0, j1, k0, k1)
        !! Convenience routine to setup a forcing_patch_type. This sets the start and end times, and 
        !! allocates forcing_patch%forcing_work(i0:i1, j0:j1, k0:k1). If k0,k1 are not provided then
        !! they are set to STG,ELV as used in the domain object.
        class(forcing_patch_type), intent(inout) :: forcing_patch
            !! The input forcing patch
        real(dp), intent(in) :: start_time, end_time
            !! The time range over which forcing_patch%forcing_work will be added. Beware the values
            !! of forcing_patch%forcing_work will have to be specified separately
        integer(ip), intent(in) :: i0, i1, j0, j1
            !! Range for first dimension is i0:i1, range for second dimension is j0:j1
        integer(ip), optional, intent(in) :: k0, k1
            !! Range for third dimension is k0:k1, by default that is STG:ELV

        forcing_patch%start_time = start_time
        forcing_patch%end_time   = end_time

        forcing_patch%i0 = i0
        forcing_patch%i1 = i1
        forcing_patch%j0 = j0
        forcing_patch%j1 = j1

        ! Allocate the forcing_work component
        if( (.not. present(k0)) .and. (.not. present(k1)) )  then
            ! In this case we have not specified the bounds of the 3rd dimension.
            ! Just use STG:ELV

            forcing_patch%k0 = STG
            forcing_patch%k1 = ELV

            ! Make the same dimensions as domain%U
            allocate(forcing_patch%forcing_work(i0:i1, j0:j1, STG:ELV))
        else
            ! In this case we have specified the bounds of the 3rd dimension.
            ! In that case both should be provided
            if(present(k0) .and. present(k1)) then
                ! User-specified k0:k1
                forcing_patch%k0 = k0
                forcing_patch%k1 = k1

                ! The 3rd dimension indices are specified
                allocate(forcing_patch%forcing_work(i0:i1, j0:j1, k0:k1))

            else
                ! Not permitted

                write(log_output_unit, *) 'Error in calling forcing_patch%setup: Must either provide BOTH k0 and k1, or neither'
                call generic_stop
            end if
        end if

        ! Initialise it to something.
        forcing_patch%forcing_work = 0.0_dp

    end subroutine

    subroutine deallocate_forcing_patch(forcing_patch)
        !! Deallocate forcing_patch%forcing_work and revert to the default start_time, end_time, and spatial indices
        class(forcing_patch_type), intent(inout) :: forcing_patch
            !! A forcing patch to deallocate

        deallocate(forcing_patch%forcing_work)

        ! Set the other parameters to the typical defaults
        forcing_patch%i0 = HUGE(1_ip)
        forcing_patch%i1 = -HUGE(1_ip)
        forcing_patch%j0 = HUGE(1_ip)
        forcing_patch%j1 = -HUGE(1_ip)
        forcing_patch%start_time = HUGE(1.0_dp)
        forcing_patch%end_time = -HUGE(1.0_dp)

        forcing_patch%k0 = STG
        forcing_patch%k1 = ELV

    end subroutine

    subroutine setup_domain_for_test(domain, global_lw, global_nx, global_ll)
        !! Routine to setup the domain for the test, and make the forcing_context

        type(domain_type), intent(inout) :: domain
        real(dp) :: global_lw(2), global_ll(2)
        integer(ip) :: global_nx(2)
        type(forcing_patch_type), pointer :: init_context

        integer(ip) :: i, j

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

        ! Make up the forcing -- first allocate it, so that there is a unique one for the domain
        ! When the routine exits this memory will not be freed because init_context is a pointer (not an allocatable)
        allocate(init_context)
        call init_context%setup(start_time = 0.0_dp, end_time = 0.0_dp, i0=1, i1=domain%nx(1), j0=1, j1=domain%nx(2))
            ! Later we will manually change start_time, end_time 

        init_context%forcing_work(:,:,1) = 10.0_dp
        init_context%forcing_work(:,:,2) = 5.0_dp
        init_context%forcing_work(:,:,3) = 1.0_dp
        init_context%forcing_work(:,:,4) = -1.0_dp
        ! Now store the result in the domain -- this is what will be used later to do stuff
        domain%forcing_context_cptr = c_loc(init_context)


    end subroutine

    subroutine test_forcing_mod
        !! Test the forcing subroutines by setting up 2 domains, and checking the forcing works

        real(dp), parameter :: global_lw(2) = [50.0_dp, 30.0_dp], local_lw(2) = [10.0_dp, 10.0_dp]
        real(dp), parameter :: global_ll(2) = [0.0_dp, 0.0_dp], local_ll(2) = [1.0_dp, 1.0_dp]
        integer(ip), parameter :: global_nx(2) = [50, 30], local_nx(2) = [23, 90]

        type(domain_type), allocatable :: domains(:)
        integer(ip) :: i, j
        real(dp) :: err_tol, err, time, dt
        type(forcing_patch_type), pointer :: domain_forcing_context
        integer(ip), parameter :: nd = 2

        err_tol = spacing(10.0_dp)*10

        allocate(domains(nd))

        ! Setup the main test domain
        call setup_domain_for_test(domains(1), global_lw, global_nx, global_ll)
        ! Backup the initial condition for later usage
        domains(1)%backup_U = domains(1)%U

        ! Setup another test domain
        call setup_domain_for_test(domains(2), local_lw, local_nx, local_ll)
        domains(2)%backup_U = domains(2)%U

        ! Suppose we just evolved from 0.0 to 10.0
        time = 10.0_dp
        dt = 10.0_dp

        !
        ! Tests of add_forcing_work_to_U_over_time
        !

        ! Test 1 -- time-step exactly equals range of forcing
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
                ! unpack the forcing_context, which is stored as a c_ptr in the domain
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                time, dt, start_time = 0.0_dp, end_time = 10.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work)))
            if(err < err_tol) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__ 
            end if
            domains(j)%U = domains(j)%backup_U
        end do


        ! Test 2 -- time-step starts before range of forcing
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                time, dt, start_time = 4.0_dp, end_time = 10.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work)))
            if(err < err_tol ) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        ! Test 3 -- time-step starts after range of forcing
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                time, dt, start_time = -4.0_dp, end_time = 10.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/14.0_dp)))
            if(err < err_tol) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__, &
                   minval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/14.0_dp)), &
                   maxval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/14.0_dp))
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        ! Test 4 -- time-step fully inside the forcing time
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                time, dt, start_time = -15.0_dp, end_time = 15.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/30.0_dp)))
            if(err < err_tol) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__, &
                   minval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/30.0_dp)), &
                   maxval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 10.0_dp/30.0_dp))
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        ! Test 5 -- time-step overshoots the forcing end time
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work,&
                time, dt, start_time = -15.0_dp, end_time = 5.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 5.0_dp/20.0_dp)))
            if(err < err_tol) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__, &
                   minval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 5.0_dp/20.0_dp)), &
                   maxval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 5.0_dp/20.0_dp))
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        ! Test 6 -- time-step doesn't overlap with the forcing 
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                time, dt, start_time = -15.0_dp, end_time = -5.0_dp)
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 0.0_dp)))
            if(err == 0.0_dp) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__, &
                   minval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 0.0_dp)), &
                   maxval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 0.0_dp))
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        ! Test 7 -- a timestepping type test, where we should end up adding the full contribution of forcing_work
        ! Use some erratic numbers -- it shouldn't matter so long as we fully step through start_time to end_time
        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            time = -10.01242_dp
            dt = 1.4444_dp
            do i = 1, 30
                ! When the loop is complete we will have covered the time
                call add_forcing_work_to_U_over_time(domains(j)%U, domain_forcing_context%forcing_work, &
                    time, dt, start_time = -1.332_dp, end_time = 5.15430_dp)
                time = time + dt
            end do
            err = maxval(abs(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work * 1.0_dp)))
            if(err < err_tol) then
                print*, 'PASS'
            else
                print*, 'FAIL', err, err_tol, __LINE__, &
                   minval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work)), &
                   maxval(domains(j)%U - (domains(j)%backup_U + domain_forcing_context%forcing_work))
            end if
            domains(j)%U = domains(j)%backup_U
        end do

        do j = 1, nd
            call c_f_pointer(domains(j)%forcing_context_cptr, domain_forcing_context)
            deallocate(domain_forcing_context%forcing_work)
            deallocate(domain_forcing_context)
        end do

    end subroutine

end module
