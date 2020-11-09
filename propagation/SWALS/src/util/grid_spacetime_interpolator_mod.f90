module grid_spacetime_interpolator_mod
    !! Interpolation of multi-gridded time-series.

    use global_mod, only : dp, ip, charlen
    use linear_interpolator_mod, only: nearest_index_sorted
    use logging_mod, only: log_output_unit
    use stop_mod, only: generic_stop
    use iso_c_binding, only: c_ptr
    implicit none

    real(dp), parameter :: cx = 2.7_dp, cy = 3.05_dp
        !! Parameters for test case

    type grid_spacetime_interpolator_type
        !! Interpolation of multi-gridded time-series.
        !! Bilinear inteprolation in space, linear interpolation in time.

        integer(ip) :: nx(3)
            !! Dimensions of grids used for interpolating
        real(dp) :: g1_time, g2_time
            !! Times for grids g1/g2
        real(dp) :: dt
            !! The time between g1_time and g2_time. If we try to interpolate past g2_time,
            !! then g1 is replaced with g2, and get_grid_at_time is called to update g2.
        real(dp), allocatable :: x(:), y(:)
            !! x/y coordinates of the fields
        real(dp), allocatable :: g1(:,:,:), g2(:,:,:)
            !! Fields used to interpolate at the 'first time' and 'second time'
            !! The first dimension has size = size(x)
            !! The second dimension has size = size(y)
            !! The third dimension has size = "number of grids". Often this will be 1,
            !! but sometimes we want several grids (e.g. east/north velocity)
        type(c_ptr) :: interp_dataptr
            !! Used to store any data required by the get_grid_at_time function

       
        procedure(get_grid_at_time), pointer, nopass :: get_grid_at_time => NULL()
            !! User-provided subroutine to update the g1 and g2 values 

        contains

        procedure :: initialise => initialise_grid_spacetime_interpolator
        procedure :: setup_to_read_at_time => setup_to_read_at_time
        procedure :: interpolate => interpolate_from_grid_spacetime
        procedure :: interpolate_grid_in_time => interpolate_full_grids_time_only
        procedure :: finalise => finalise_grid_spacetime_interpolator
            
    end type

    interface
        subroutine get_grid_at_time(time, x, y, grid)
            !! A subroutine with this interface must be provided to the grid_spacetime_interpolator
            !! It populates the grid given the time, and the x/y coordinates associated with the grid.
            import dp
            implicit none
            real(dp), intent(in) :: time
            real(dp), intent(in) :: x(:), y(:)
            real(dp), intent(inout) :: grid(:,:,:)
        end subroutine
    end interface

    contains

    subroutine initialise_grid_spacetime_interpolator(gsi, x, y, start_time, dt, ng)
        !! Setup the grid_spacetime_interpolator 
        !! Allocate memory, set gsi%g1, gsi%g2, etc.

        class(grid_spacetime_interpolator_type), intent(inout) :: gsi
        real(dp), intent(in) :: x(:), y(:)
            !! x/y values for the grid on which interpolation occurs 
        real(dp), intent(in) :: start_time, dt
            !! the time associated with g1, and the time between interpolation grids g1, g2
        integer(ip), intent(in), optional :: ng
            !! Size of the third dimension of gsi%g1, gsi%g2

        integer(ip) :: i, ng_local

        ! By default assume third dimension of gsi%g1, gsi%g2 has size 1
        if(present(ng)) then
            ng_local = ng
        else
            ng_local = 1_ip
        end if

        do i = 1, (size(x, kind=ip)-1)
            if(x(i+1) <= x(i)) then
                write(log_output_unit, *) "x must be increasing", x
                call generic_stop
            end if
        end do
        do i = 1, (size(y, kind=ip)-1)
            if(y(i+1) <= y(i)) then
                write(log_output_unit, *) "y must be increasing", y
                call generic_stop 
            end if
        end do

        ! Grid size
        gsi%nx = [size(x, kind=ip), size(y, kind=ip), ng_local]

        ! Grid coordinates
        allocate(gsi%x(size(x, kind=ip)), gsi%y(size(y, kind=ip)))
        gsi%x = x
        gsi%y = y

        ! Grids
        allocate(gsi%g1( gsi%nx(1), gsi%nx(2), gsi%nx(3)))
        allocate(gsi%g2( gsi%nx(1), gsi%nx(2), gsi%nx(3)))

        ! Time associated with grids
        gsi%dt = dt
        gsi%g1_time = start_time
        gsi%g2_time = start_time + dt

        ! Read the grids
        call gsi%get_grid_at_time(gsi%g1_time, x, y, gsi%g1)
        call gsi%get_grid_at_time(gsi%g2_time, x, y, gsi%g2)

    end subroutine

    subroutine setup_to_read_at_time(gsi, time)
        !! Ensure that gsi contains grids to interpolate at the given time. 
        !! Although we can integrate this check into the interpolation routine, that is
        !! potentially problematic in parallel with threads (because many threads might try
        !! to update g1/g2 from a file). To avoid that, we can call this routine separately on a single thread,
        !! before doing interpolation in parallel.
        class(grid_spacetime_interpolator_type), intent(inout) :: gsi
        real(dp), intent(in) :: time
            !! Time at which we will want to interpolate

        if(time < gsi%g1_time) then
            write(log_output_unit, *) "Currently cannot interpolate at time < gsi%g1_time=", gsi%g1_time
            call generic_stop
        end if

        if(time > gsi%g2_time) then
            ! Need to update the gsi grids, because we asked for a time outside the range of g1/g2

            if(time > gsi%g2_time + gsi%dt) then
                write(log_output_unit,*) "Currently cannot take take steps larger than gsi%dt"
                call generic_stop
            end if

            ! g1 ==> g2
            gsi%g1_time = gsi%g2_time
            gsi%g1 = gsi%g2
            
            ! Update g2
            gsi%g2_time = gsi%g2_time + gsi%dt
            call gsi%get_grid_at_time(gsi%g2_time, gsi%x, gsi%y, gsi%g2)

        end if

        if(.not. (time >= gsi%g1_time .and. time <= gsi%g2_time)) then
            write(log_output_unit, *) "Logical error in grid_spacetime setup_to_read_at_time"
            call generic_stop
        end if

    end subroutine

    subroutine interpolate_full_grids_time_only(gsi, time, values, suppress_checks)
        !!
        !! Interpolate linearly in time between gsi%g1 and gsi%g2. Return an array with 
        !! the same dimensions as gsi%g1 and gsi%g2
        !!
        class(grid_spacetime_interpolator_type), intent(inout) :: gsi
        real(dp), intent(in) :: time
            !! Time at which interpolated grid values are required
        real(dp), intent(inout) :: values(:,:,:)
            !! Stores output
        logical, optional :: suppress_checks
            !! If .TRUE., do not check that time is inside the current 'gsi grids time range', or the dimensions of the inputs. 
            !! This check will provoke an automatic update of the grids if required, so if you suppress it, ensure you know that no
            !! update is needed.

        logical :: check
        integer(ip) :: i, j, k
        real(dp) :: w1_t, w2_t

        if(present(suppress_checks)) then
            check = .not.(suppress_checks)
        else
            check = .true.
        end if

        
        if(check) then
            ! Ensure that time is within [gsi%g1_time, gsi%g2_time]
            ! If not, update gsi
            call gsi%setup_to_read_at_time(time)

        end if

        if(.not.(all(shape(gsi%g1) == shape(values)))) then
            write(log_output_unit, *) "Error in grid_spacetime_interpolator_type: shape of values must equal shape of gsi grids"
            call generic_stop
        end if


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(gsi, values, time)

        ! Time-interpolation weights
        w1_t = (gsi%g2_time - time)/(gsi%g2_time - gsi%g1_time)
        w2_t = 1.0_dp - w1_t

        !$OMP DO COLLAPSE(3)
        do k = 1, size(values, 3, kind=ip)
            do j = 1, size(values, 2, kind=ip)
                do i = 1, size(values, 1, kind=ip)
                    values(i, j, k) = w1_t * gsi%g1(i,j,k) + w2_t * gsi%g2(i,j,k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine

    subroutine interpolate_from_grid_spacetime(gsi, time, xs, ys, values, suppress_checks)
        !! Get 'values' given coordinates xs, ys at the given time. 
        !! Use bilinear interpolation in space, and linear interpolation in time.
        !! Note that xs, ys can be completely unstructured -- however, this means the lookup is relatively slow (finding the
        !! nearest cell involves 2 bisection searches for each point). 

        class(grid_spacetime_interpolator_type), intent(inout) :: gsi
        real(dp), intent(in) :: time
            !! Time when output is desired
        real(dp), intent(in) :: xs(:), ys(:)
            !! Coordinates where output is desired
        real(dp), intent(inout) :: values(:,:)
            !! Holds output
        logical, optional, intent(in) :: suppress_checks
            !! If .TRUE., do not check that time is inside the current 'gsi grids time range', or the dimensions of the inputs. 
            !! This check will provoke an automatic update  if required, so if you suppress it, ensure you know that no update is
            !! needed.

        real(dp) :: w1_t, w2_t, w0_x, w1_x, w0_y, w1_y, v0, v1, g1_val, g2_val
        integer(ip) :: i, i0, i1, j0, j1, k
        logical :: check

        if(present(suppress_checks)) then
            check = .not.(suppress_checks)
        else
            check = .true.
        end if

        if(check) then

            call gsi%setup_to_read_at_time(time)

            if(.not. (size(xs, kind=ip) == size(ys, kind=ip) .and. size(xs, kind=ip) == size(values, 1, kind=ip))) then
                write(log_output_unit, *) "Sizes of xs, ys and the first-dim of values should be equal"
                call generic_stop
            end if

            if(size(values, 2, kind=ip) /= size(gsi%g1, 3, kind=ip)) then
                write(log_output_unit, *) "Error: The second dimension of values should have"
                write(log_output_unit, *) "       the same size as the third dimension of gsi%g1"
            end if

        end if

        ! Time-interpolation weights
        w1_t = (gsi%g2_time - time)/(gsi%g2_time - gsi%g1_time)
        w2_t = 1.0_dp - w1_t
        !print*, 'Time weights: ', w1_t, w2_t

        ! Do the interpolation
        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(gsi, xs, ys, values, time)
        do i = 1, size(xs, kind=ip)

            ! Get x-indices for bilinear interpolation: i0 < i1
            call nearest_index_sorted(gsi%nx(1), gsi%x, xs(i), i0)
            if(xs(i) > gsi%x(i0)) then
                i1 = min(i0+1, gsi%nx(1))
            else
                i1 = i0
                i0 = max(i0-1, 1)
            end if
            ! Get weights for x bilinear interpolation
            if(i1 == i0) then
                w0_x = 0.5_dp
            else
                w0_x = (gsi%x(i1) - xs(i))/(gsi%x(i1) - gsi%x(i0))
            end if
            w1_x = 1.0_dp - w0_x

            !print*, xs(i), gsi%x(i0), gsi%x(i1), i0, i1, w0_x, w1_x

            ! Get y-indices for bilinear interpolation: j0 < j1
            call nearest_index_sorted(gsi%nx(2), gsi%y, ys(i), j0)
            if(ys(i) > gsi%y(j0)) then
                j1 = min(j0+1, gsi%nx(2))            
            else
                j1 = j0
                j0 = max(j0-1, 1)
            end if
            ! Get weights for y bilinear interpolation
            if(j1 == j0) then
                w0_y = 0.5_dp
            else
                w0_y = (gsi%y(j1) - ys(i))/(gsi%y(j1) - gsi%y(j0))
            end if
            w1_y = 1.0_dp - w0_y
            
            !print*, ys(i), gsi%y(j0), gsi%y(j1), j0, j1, w0_y, w1_y

            do k = 1, size(gsi%g1, 3, kind=ip)
                !! Bilinear interpolation on g1
                ! Pure j0 value
                v0 = w0_x * gsi%g1(i0, j0, k) + w1_x * gsi%g1(i1, j0, k)
                ! Pure j1 value
                v1 = w0_x * gsi%g1(i0, j1, k) + w1_x * gsi%g1(i1, j1, k)
                ! Interpolated in y direction
                g1_val = w0_y * v0 + w1_y * v1
                !print*, v0, v1, g1_val

                !! Bilinear interpolation on g2
                ! Pure j0 value
                v0 = w0_x * gsi%g2(i0, j0, k) + w1_x * gsi%g2(i1, j0, k)
                ! Pure j1 value
                v1 = w0_x * gsi%g2(i0, j1, k) + w1_x * gsi%g2(i1, j1, k)
                ! Interpolated in y direction
                g2_val = w0_y * v0 + w1_y * v1
                !print*, v0, v1, g2_val

                values(i,k) = w1_t * g1_val + w2_t * g2_val
            end do

        end do
        !$OMP END PARALLEL DO

    end subroutine

    subroutine finalise_grid_spacetime_interpolator(gsi)
        !! Clean-up the grid_spacetime_interpolator
        class(grid_spacetime_interpolator_type), intent(inout) :: gsi

        deallocate(gsi%x, gsi%y, gsi%g1, gsi%g2)
        gsi%get_grid_at_time => NULL()

        ! Crazy values to help throw errors with incorrect usage
        gsi%g1_time = HUGE(1.0_dp)
        gsi%g2_time = -HUGE(1.0_dp)
        gsi%dt = -HUGE(1.0_dp)

    end subroutine

    ! Use this for testing the grid interpolator
    subroutine for_test_get_grid_at_time(time, x, y, grid)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: x(:), y(:)
            real(dp), intent(inout) :: grid(:,:,:)

            integer(ip) :: j, k

            do k = 1, size(grid,3, kind=ip)
            do j = 1, size(grid,2, kind=ip)
                grid(:,j,k) = (cx * x + cy * y(j) + time)*k
            end do
            end do

    end subroutine

    subroutine test_grid_spacetime_interpolator_mod
        !! Basic unit-test routine
        type(grid_spacetime_interpolator_type) :: gsi
        real(dp) :: xs(10), ys(20), time, dt, start_time
        real(dp) :: dx, dy, x0, y0
        integer(ip), parameter :: ng = 2
        integer(ip) :: i, j, k
        real(dp) :: tx(5), ty(5), tv(5,ng), test_grid(10, 20, ng)
        real(dp) :: errtol 
        
        errtol = 3*spacing(100.0_dp) 

        ! x-coordinates of gsi
        dx = 0.5_dp
        x0 = -1.0_dp
        do i = 1, size(xs)
            xs(i) = x0 + (i-1)*dx
        end do
        !print*, xs

        ! y-coordinates of gsi
        dy = 0.1_dp
        y0 = 12.0_dp
        do i = 1, size(ys)
            ys(i) = y0 + (i-1)*dy
        end do
        !print*, ys

        ! Setup the gsi object
        start_time = 35.0_dp
        dt = 15.0_dp
        gsi%get_grid_at_time => for_test_get_grid_at_time
        call gsi%initialise(xs, ys, start_time, dt, ng)

        !
        ! Setup the test data -- off the grid in this case
        !
        tx = x0 - 1.0_dp
        ty = y0 - 1.0_dp
        time = start_time + 5.0_dp
        call gsi%interpolate(time, tx, ty, tv)
        do k = 1, ng
            if(all(abs(tv(:,k) - k*(time + cx * x0 + cy * y0)) < errtol)) then
                print*, 'PASS'
            else
                print*, 'FAIL', abs(tv(:,k) - k*(time + cx * x0 + cy * y0))
            end if
        end do

        !
        ! Setup the test data -- on the grid in this case
        !
        tx = x0 + dx * ([1.0 , 3.0, 6.3, 4.0, 8.0] - 1.0)
        ty = y0 + dy * ([20.0, 4.0, 5.2, 1.0, 2.0] - 1.0)
        time = start_time + 5.0_dp
        call gsi%interpolate(time, tx, ty, tv)
        do k = 1, ng
            if(all(abs(tv(:,k) - k*(time + cx * tx + cy * ty)) < errtol)) then
                print*, 'PASS'
            else
                print*, 'FAIL', abs(tv(:,k) - k*(time + cx * tx + cy * ty))
            end if
        end do

        ! As before, with time large enough to trigger a grid update
        time = start_time + 5.0_dp + dt
        call gsi%interpolate(time, tx, ty, tv) ! This will cause g1 => g2, and g2 to be updated using for_test_get_grid_at_time
        do k = 1, ng
            if(all(abs(tv(:,k) - k*(time + cx * tx + cy * ty)) < errtol)) then
                print*, 'PASS'
            else
                print*, 'FAIL', abs(tv(:,k) - k*(time + cx * tx + cy * ty))
            end if
        end do

        ! Finally check we can interpolate the grids in time
        time = 0.6_dp*gsi%g1_time + 0.4_dp*gsi%g2_time
        if(abs(time - (0.6_dp*50.0_dp + 0.4_dp*65.0_dp)) < errtol) then
            print*, 'PASS'
        else
            print*, 'FAIL', time, time - (0.6_dp*50.0_dp + 0.4_dp*65.0_dp), errtol
        end if
        call gsi%interpolate_grid_in_time(time, test_grid)
        if(all(abs(test_grid - (0.6_dp * gsi%g1 + 0.4_dp * gsi%g2)) < errtol)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        call gsi%finalise

    end subroutine

end module
