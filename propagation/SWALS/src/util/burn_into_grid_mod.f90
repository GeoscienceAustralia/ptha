module burn_into_grid_mod
    !!
    !! Allows definition of 'xyz lines' which can denote the elevation of linear features, and
    !! provides routines to burn these into grids. For example, this can be used to ensure that 
    !! an elevation grid includes a continuous representation of a breakwall or riverwall.
    !!
    use global_mod, only: ip, dp, charlen, force_double
    use logging_mod
    use stop_mod
    use file_io_mod, only: read_csv_into_array
    implicit none

    private

    public test_burn_into_grid_mod, xyz_lines_type

    type rank2_allocatable_type
        !! Store the xyz values of the linear features we want to burn into the grid
        real(dp), allocatable :: xyz(:,:)
    end type

    type xyz_lines_type
        !! Type to hold a collection of xyz lines, and enable 'burning' the z values into a grid (i.e. setting the grid values along
        !! the lines) 
        type(rank2_allocatable_type), allocatable :: lines(:)
        contains
        procedure :: read_from_csv => read_xyz_lines_from_csv
        procedure :: burn_into_grid => burn_lines_into_grid
    end type

    contains

    subroutine read_xyz_lines_from_csv(xyz_lines, line_files, skip_header)
        !! Read a set of line_files (each containing x,y,z data for a xyz-line in csv format)
        !! into an xyz_lines_type object. 
        class(xyz_lines_type), intent(inout) :: xyz_lines !! The xyz_lines to be set from the line_files.
        character(len=charlen), intent(in) :: line_files(:) !! Array with csv file names containing xyz lines data.
        integer(ip), optional, intent(in) :: skip_header !! How many header rows to skip when reading the files (default 0)?
    
        integer(ip) :: i, nl, skip_h

        if(present(skip_header)) then
            skip_h = skip_header
        else
            skip_h = 0
        end if

        nl = size(line_files, kind=ip)

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

    subroutine burn_xyz_into_grid(x, y, z, grid, lower_left, upper_right, burn_type)
        !! Given a set of points x,y,z, and a grid of values with prescribed lower_left and upper_right,
        !! burn the z values into the grid at cells containing x,y
        real(dp), intent(in) :: x(:), y(:), z(:) !! rank-1 arrays with point xyz coordinates
        real(dp), intent(inout) :: grid(:,:) !! rank-2 array with grid values
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! coordinates defining the grid extent
        character(len=*), optional :: burn_type
        !! character controlling when we burn -- either 'point_value' (default, always burn), or 'max' (only burn if z >
        !! grid value) or 'min' (only burn if z < grid_value)

        real(dp) :: dx(2)
        integer(ip) :: i, i0, j0
        character(len=charlen) :: burnt

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if
        

        if(size(x, kind=ip) /= size(y, kind=ip) .or. size(x, kind=ip) /= size(z, kind=ip)) then
            write(log_output_unit, *) "Error in burn_xyz_into_grid: x,y, and z must have the same length"
            call generic_stop
        end if     

        if(.not. all(upper_right > lower_left)) then
            write(log_output_unit, *) "Error in burn_xyz_into_grid: upper_right must be > lower_left"
            call generic_stop
        end if

        ! Cell size
        dx = (upper_right - lower_left)/(1.0_dp * shape(grid))

        do i = 1, size(x, kind=ip)
            ! Indices of x,y
            i0 = floor((x(i)-lower_left(1))/dx(1)) + 1
            j0 = floor((y(i)-lower_left(2))/dx(2)) + 1
            ! Burn z if it is inside the grid
            if(i0 > 0 .and. j0 > 0 .and. i0 <= size(grid, 1, kind=ip) .and. j0 <= size(grid, 2, kind=ip)) then
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

    subroutine burn_line_into_grid(x, y, z, grid, lower_left, upper_right, burn_type)
        !!
        !! Similar to burn_xyz_into_grid, but interpolates along xyz as required to ensure we don't miss cells even if the
        !! line segments have length >> dx.
        !! Note that the interpolation is a bit hap-hazard (in terms of exactly which interpolation point is hit).
        !!
        real(dp), intent(in) :: x(:), y(:), z(:) !! Coordinates of the 3D line
        real(dp), intent(inout) :: grid(:,:) !! Grid into which we burn the elevations
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! Coordinates of the grid
        character(len=*), optional :: burn_type 
        !! character controlling when we burn -- either 'point_value' (default, always burn), or 'max' (only burn if z >
        !! grid value) or 'min' (only burn if z < grid_value)

        real(dp) :: dx(2), xi, yi, zi, xip1, yip1, zip1, xj(1), yj(1), zj(1)
        integer(ip) :: i, i0, j0, np, j
        real(force_double) :: n1, n2
        character(len=charlen) :: burnt

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if

        if(size(x, kind=ip) /= size(y, kind=ip) .or. size(x, kind=ip) /= size(z, kind=ip)) then
            write(log_output_unit, *) "Error in burn_line_into_grid: x,y, and z must have the same length"
            call generic_stop
        end if     

        if(.not. all(upper_right > lower_left)) then
            write(log_output_unit, *) "Error in burn_line_into_grid: upper_right must be > lower_left"
            call generic_stop
        end if

        ! Cell size
        dx = (upper_right - lower_left)/(1.0_dp * shape(grid))

        do i = 1, size(x, kind=ip)
            if(i < size(x, kind=ip)) then
                xi = x(i)
                yi = y(i)
                zi = z(i)
                xip1 = x(i+1)
                yip1 = y(i+1)
                zip1 = z(i+1)

                ! Interpolate np+1 points along the line segment between points i and (i+1)
                n1 = abs(real(xip1, force_double)- real(xi, force_double))/real(dx(1), force_double)
                n2 = abs(real(yip1, force_double)- real(yi, force_double))/real(dx(2), force_double)
                np = ceiling(max(n1, n2)) + 1
                do j = 0, np
                    xj = xi + (j*(xip1-xi))/np
                    yj = yi + (j*(yip1-yi))/np
                    zj = zi + (j*(zip1-zi))/np
                    !print*, np, j, xj, yj, zj
                    call burn_xyz_into_grid(xj, yj, zj, grid, lower_left, upper_right, burnt)
                end do
            else
                call burn_xyz_into_grid(x(i:i), y(i:i), z(i:i), grid, lower_left, upper_right, burnt)
            end if
        end do

    end subroutine

    subroutine burn_lines_into_grid(xyz_lines, grid, lower_left, upper_right, burn_type)
        !!
        !! Burn a set of xyz lines into a grid, interpolating along the lines as required.
        !!
        class(xyz_lines_type), intent(in) :: xyz_lines !! The xyz lines to burn into the grid
        real(dp), intent(inout) :: grid(:,:) !! The grid
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! The grid extent
        character(len=*), optional :: burn_type
        !! character controlling when we burn -- either 'point_value' (default, always burn the z value into the grid), or 'max' 
        !! (only burn if z > grid value) or 'min' (only burn if z < grid_value)

        character(len=charlen) :: burnt
        integer(ip) :: i

        if(present(burn_type)) then
            burnt = burn_type
        else
            burnt = 'point_value'
        end if

        do i = 1, size(xyz_lines%lines, kind=ip)
            call burn_line_into_grid(&
                xyz_lines%lines(i)%xyz(1,:), xyz_lines%lines(i)%xyz(2,:), xyz_lines%lines(i)%xyz(3,:), &
                grid, lower_left, upper_right, burnt)
        end do

    end subroutine

    subroutine test_burn_into_grid_mod()
        !! Unit tests
        real(dp) :: grid(10,10)
        real(dp) :: lower_left(2), upper_right(2)
        real(dp) :: x(5), y(5), z(5)
        real(dp) :: expected_grid(10,10)
        real(dp) :: lx(3), ly(3), lz(3)
        integer(ip) :: i
        real(dp), parameter :: eps = 1.0e-03_dp

        lower_left = [0.0_dp, 0.0_dp]
        upper_right = [10.0_dp, 10.0_dp]

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
        ! (the 'eps' in the line definition helps the test to pass in both single and double precision,
        !  otherwise round-off can shift the results).
        !
        lx = [-1.0_dp ,  5.5_dp,  5.5_dp] + eps
        ly = [3.5_dp  ,  3.5_dp, 11.5_dp] - eps
        lz = lx + ly
        grid = 0.0_dp
        call burn_line_into_grid(lx, ly, lz, grid, lower_left, upper_right)
        expected_grid = 0.0_dp
        expected_grid(1:6,4) = [4.12500000_dp, 4.93750000_dp,  5.75000000_dp, 7.37500000_dp, 8.18750000_dp, 9.00000000_dp]
        expected_grid(6,5:10) = [9.88888931_dp, 10.7777777_dp, 11.6666670_dp, 13.4444447_dp, 14.3333330_dp, 15.2222223_dp]
        if(all(abs(grid - expected_grid) < 1.0e-05_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !print*, 'grid(1:6, 4): ', grid(1:6, 4)
        !print*, 'expected_grid(1:6, 4): ', expected_grid(1:6, 4) 
        !print*, 'grid(6, 5:10): ', grid(6, 5:10)
        !print*, 'expected_grid(6, 5:10): ', expected_grid(6, 5:10) 

    end subroutine

end module
