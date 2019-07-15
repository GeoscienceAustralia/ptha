module burn_into_grid_mod
    use global_mod, only: ip, dp, charlen
    use logging_mod
    use stop_mod
    use file_io_mod, only: read_csv_into_array
    implicit none

    private

    public test_burn_into_grid_mod, burn_xyz_into_grid, burn_line_into_grid, xyz_lines_type

    type rank2_allocatable_type
        real(dp), allocatable :: xyz(:,:)
    end type

    ! Type to hold a collection of xyz lines
    type xyz_lines_type
        type(rank2_allocatable_type), allocatable :: lines(:)
        contains
        procedure :: read_from_csv => read_xyz_lines_from_csv
        procedure :: burn_into_grid => burn_lines_into_grid
    end type

    contains

    subroutine read_xyz_lines_from_csv(xyz_lines, line_files, skip_header)
        class(xyz_lines_type), intent(inout) :: xyz_lines
        character(len=charlen), intent(in) :: line_files(:)
        integer(ip), optional, intent(in) :: skip_header
    
        integer(ip) :: i, nl, skip_h

        if(present(skip_header)) then
            skip_h = skip_header
        else
            skip_h = 0
        end if

        nl = size(line_files)

        allocate(xyz_lines%lines(nl))

        do i = 1, nl
            call read_csv_into_array(xyz_lines%lines(i)%xyz, line_files(i), skip_header=skip_h)
            if(size(xyz_lines%lines(i)%xyz, 1) /= 3) then
                write(log_output_unit, *) "Error in read_xyz_lines_from_csv: file ", line_files(i)
                write(log_output_unit, *) "does not seem to be a 3 column csv file"
                call generic_stop
            end if
        end do

    end subroutine

    ! Given a set of points x,y,z, and a grid of values with prescribed lower_left and upper_right,
    ! burn the z values into the grid at cells containing x,y
    ! @param x, y, z rank-1 arrays with point xyz coordinates
    ! @param grid rank-2 array with grid values
    ! @param lower_left, upper_right coordinates defining the grid extent
    ! @param burn_type character controlling when we burn -- either 'point_value' (default, always burn), or 'max' (only burn if z >
    ! grid value) or 'min' (only burn if z < grid_value)
    subroutine burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right, burn_type)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(inout) :: grid(:,:)
        real(dp), intent(in) :: lower_left(2), upper_right(2)
        character(len=*), optional :: burn_type

        real(dp) :: dx(2)
        integer(ip) :: i, i0, j0
        character(len=charlen) :: burnt

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if
        

        if(size(x) /= size(y) .or. size(x) /= size(z)) then
            write(log_output_unit, *) "Error in burn_xyz_into_grid: x,y, and z must have the same length"
            call generic_stop
        end if     

        if(.not. all(upper_right > lower_left)) then
            write(log_output_unit, *) "Error in burn_xyz_into_grid: upper_right must be > lower_left"
            call generic_stop
        end if

        ! Cell size
        dx = (upper_right - lower_left)/(1.0_dp * shape(grid))

        do i = 1, size(x)
            ! Indices of x,y
            i0 = floor((x(i)-lower_left(1))/dx(1)) + 1
            j0 = floor((y(i)-lower_left(2))/dx(2)) + 1
            ! Burn z if it is inside the grid
            if(i0 > 0 .and. j0 > 0 .and. i0 <= size(grid, 1) .and. j0 <= size(grid, 2)) then
                select case (burnt)
                case('point_value')
                    grid(i0, j0) = z(i)
                case('max')
                    grid(i0, j0) = max(grid(i0, j0), z(i))
                case('min')
                    grid(i0, j0) = min(grid(i0, j0), z(i))
                case default
                    write(log_output_unit, *) "Error in burn_xyz_into_grid: unknown value of burn_type"
                    call generic_stop
                end select
            end if
        end do

    end subroutine

    ! As above, but interpolate along xyz as required to ensure we don't miss cells even if the
    ! line segments have length >> dx.
    ! Note that the interpolation is a bit hap-hazard (in terms of exactly which interpolation point is hit)
    subroutine burn_line_into_grid(x, y, z, grid, lower_left, upper_right, burn_type)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(inout) :: grid(:,:)
        real(dp), intent(in) :: lower_left(2), upper_right(2)
        character(len=*), optional :: burn_type

        real(dp) :: dx(2), xi, yi, zi, xip1, yip1, zip1, xj(1), yj(1), zj(1)
        integer(ip) :: i, i0, j0, np, j
        character(len=charlen) :: burnt

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if

        if(size(x) /= size(y) .or. size(x) /= size(z)) then
            write(log_output_unit, *) "Error in burn_line_into_grid: x,y, and z must have the same length"
            call generic_stop
        end if     

        if(.not. all(upper_right > lower_left)) then
            write(log_output_unit, *) "Error in burn_line_into_grid: upper_right must be > lower_left"
            call generic_stop
        end if

        ! Cell size
        dx = (upper_right - lower_left)/(1.0_dp * shape(grid))

        do i = 1, size(x)
            if(i < size(x)) then
                xi = x(i)
                yi = y(i)
                zi = z(i)
                xip1 = x(i+1)
                yip1 = y(i+1)
                zip1 = z(i+1)

                ! Interpolate np+1 points along the line segment between points i and (i+1)
                np = ceiling(max((xip1-xi)/dx(1), (yip1-yi)/dx(2))) + 1
                do j = 0, np
                    xj = xi + (j*(xip1-xi))/np
                    yj = yi + (j*(yip1-yi))/np
                    zj = zi + (j*(zip1-zi))/np
                    call burn_xyz_into_grid(xj, yj, zj, grid, lower_left, upper_right, burnt)
                end do
            else
                call burn_xyz_into_grid(x(i:i), y(i:i), z(i:i), grid, lower_left, upper_right, burnt)
            end if
        end do

    end subroutine

    ! As above, working with our list of lines
    subroutine burn_lines_into_grid(xyz_lines, grid, lower_left, upper_right, burn_type)
        class(xyz_lines_type), intent(in) :: xyz_lines
        real(dp), intent(inout) :: grid(:,:)
        real(dp), intent(in) :: lower_left(2), upper_right(2)
        character(len=*), optional :: burn_type

        character(len=charlen) :: burnt
        integer(ip) :: i

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if

        do i = 1, size(xyz_lines%lines)
            call burn_line_into_grid(&
                xyz_lines%lines(i)%xyz(1,:), xyz_lines%lines(i)%xyz(2,:), xyz_lines%lines(i)%xyz(3,:), &
                grid, lower_left, upper_right, burnt)
        end do

    end subroutine

    subroutine test_burn_into_grid_mod()
        real(dp) :: grid(10,10)
        real(dp) :: lower_left(2), upper_right(2)
        real(dp) :: x(5), y(5), z(5)
        real(dp) :: expected_grid(10,10)
        real(dp) :: lx(3), ly(3), lz(3)
        integer(ip) :: i

        lower_left = [0.0, 0.0]
        upper_right = [10.0, 10.0]

        x = [1.6_dp, 2.2_dp, 3.6_dp, 8.7_dp, 9.9_dp]
        y = [0.5_dp, 1.3_dp, 2.6_dp, 3.7_dp, 4.8_dp]
        z = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]

        ! Regular case
        grid = 0.0_dp
        call burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right)

        ! The solution
        expected_grid = 0.0_dp
        do i = 1, 5
            expected_grid(int(x(i))+1, int(y(i))+1) = 1.0_dp*i
        end do

        if(all(grid == expected_grid)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Alternative
        !
        grid = 10.0_dp
        call burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right, burn_type = 'max')
        if(all(grid == 10.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Alternative
        grid = -10.0_dp
        call burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right, burn_type = 'min')
        if(all(grid == -10.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Nontrivial use of max
        grid = -10.0_dp
        call burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right, burn_type = 'max')
        where(expected_grid == 0.0_dp) expected_grid = -10.0_dp
        if(all(grid == expected_grid)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Line version
        !
        lx = [-1.0_dp ,  5.5_dp,  5.5_dp]
        ly = [3.5_dp  ,  3.5_dp, 11.5_dp]
        lz = lx + ly
        grid = 0.0_dp
        call burn_line_into_grid(lx, ly, lz, grid, lower_left, upper_right)
        expected_grid = 0.0
        expected_grid(1:6,4) = [4.12500000_dp, 4.93750000_dp,  5.75000000_dp, 7.37500000_dp, 8.18750000_dp, 9.00000000_dp]
        expected_grid(6,5:10) = [9.88888931_dp, 10.7777777_dp, 11.6666670_dp, 13.4444447_dp, 14.3333330_dp, 15.2222223_dp]
        if(all(abs(grid - expected_grid) < 1.0e-05_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

    end subroutine

end module
