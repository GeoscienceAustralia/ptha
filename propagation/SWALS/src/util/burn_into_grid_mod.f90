module burn_into_grid_mod
    !!
    !! Define 'xyz_lines_type' for burning linear features into grids, and
    !! 'polygons_values_type' for burning values inside polygons into a grid.
    !! In the latter case there is also a simple subroutine for the
    !! single-polygon case (burn_polygon_value_into_grid).
    !!
    !! These routines are often useful, for example, to make elevation grids
    !! have a continuous representation of a breakwall, or set different
    !! initial stage values in different parts of the domain.
    !!
    use global_mod, only: ip, dp, charlen, force_double
    use logging_mod
    use stop_mod
    use file_io_mod, only: read_csv_into_array, read_character_file
    use points_in_poly_mod, only: point_in_poly
    implicit none

    private

    public test_burn_into_grid_mod, xyz_lines_type, polygons_values_type, burn_polygon_value_into_grid

    type rank2_allocatable_type
        !! Store the xyz values of the linear features we want to burn into the grid
        real(dp), allocatable :: xyz(:,:)
    end type

    type xyz_lines_type
        !! Type to hold a collection of xyz lines, and enable 'burning' the z values into a grid (i.e. setting the grid values along
        !! the lines)
        type(rank2_allocatable_type), allocatable :: lines(:)
        contains
        procedure, non_overridable :: read_from_csv => read_xyz_lines_from_csv
        procedure, non_overridable :: burn_into_grid => burn_lines_into_grid
    end type

    type polyvalue_type
        !! Store a polygon and a value
        character(len=charlen) :: poly_file = ""
        real(dp), allocatable :: vertices(:,:) !! vertices
        real(dp) :: poly_value = -HUGE(1.0_dp) !! Value to be burned
    end type

    type polygons_values_type
        !! Hold a collection of polygon/value pairs, and enable burning
        !! the values into a grid at points inside the respective polygons
        type(polyvalue_type), allocatable :: polyvalue(:)
        contains
        procedure, non_overridable :: setup => setup_polygons_values_type_from_csv_and_value
        procedure, non_overridable :: setup_from_csv_with_file_value_columns => &
            setup_pvt_from_csv_with_file_value_columns
        procedure, non_overridable :: burn_into_grid => burn_polygons_values_type_into_grid
    end type

    contains

    subroutine setup_pvt_from_csv_with_file_value_columns(pvt, file_value_csv, &
        skip_header_file_value_csv, skip_header_polygons, verbose)
        !! Use a two-column csv file containing
        !!    "polygon_filename", polygon_value
        !! to setup a polygons_values_type
        class(polygons_values_type), intent(inout) :: pvt !! The polygons_values_type to be initialised
        character(len=*), intent(in) :: file_value_csv
            !! A csv file with two columns ("file", "value"). 
            !! In each row, the "file" entry is a filename with the polygon geometry (2 column csv file)
            !! while the "value" entry is a real number (value to set)
        integer, optional, intent(in) :: skip_header_file_value_csv
            !! Number of header rows in "file_value_csv" that should be skipped before reading
        integer, optional, intent(in) :: skip_header_polygons
            !! Number of header rows in the files containing the polygon geometries
        logical, optional, intent(in) :: verbose
            !! Print more info on the files

        integer(ip) :: skip_header_f, skip_header_p, fid, n, nfiles, cma, i
        character(len=charlen) :: buffer
        character(len=charlen), allocatable :: polygon_files(:), file_lines(:)
        real(dp), allocatable :: polygon_values(:)
        logical :: verbose_messages

        skip_header_f = 0_ip
        if(present(skip_header_file_value_csv)) skip_header_f = skip_header_file_value_csv
        
        skip_header_p = 0_ip
        if(present(skip_header_polygons)) skip_header_p = skip_header_polygons

        verbose_messages = .FALSE.
        if(present(verbose)) verbose_messages = verbose

        ! Read the input file as an array of characters
        call read_character_file(file_value_csv, file_lines, '(A)')
        n = size(file_lines)
        nfiles = n - skip_header_f ! Number of files
        allocate(polygon_values(nfiles), polygon_files(nfiles))

        if(verbose_messages) then
            write(log_output_unit, "(A)") 'Setting up polygons_values_type using ' // trim(file_value_csv)
            write(log_output_unit, "(A, I8)") '    Number of files =', nfiles
            write(log_output_unit, "(A)") '    Files, values'
            flush(log_output_unit)
        end if

        do i = 1, nfiles
            ! Read each part separately
            buffer = file_lines(i + skip_header_f)
            cma = index(buffer, ",") ! Comma index
            if(cma <= 0) then
                write(log_output_unit, *) 'Error reading ', trim(file_value_csv)
                write(log_output_unit, *) 'No comma on line ', i+skip_header_f
                flush(log_output_unit)
                call generic_stop
            end if
            polygon_files(i) = adjustl(buffer(:cma-1)) ! Set the polygon file
            read(buffer(cma+1:), *) polygon_values(i) ! Set the value
            if(verbose_messages) then 
                write(log_output_unit, "(A, ES25.12E3)") '    ' // trim(polygon_files(i)) // ',',  polygon_values(i)
            end if
        end do

        call pvt%setup(polygon_files, polygon_values, skip_header=skip_header_p)
        
    end subroutine

    subroutine setup_polygons_values_type_from_csv_and_value(pvt, polygon_files, polygon_values, skip_header)
        !! Use a set of polygon_files (each containing x,y vertices in csv format) and associated values
        !! to setup a polygons_values_type
        class(polygons_values_type), intent(inout) :: pvt !! The polygons_values_type to be initialised
        character(len=charlen), intent(in) :: polygon_files(:) !! Array with csv file names containing polygon vertices
        real(dp), intent(in) :: polygon_values(:) !! Array of the same length as polygon_files, containing their values
        integer(ip), optional, intent(in) :: skip_header !! How many header rows to skip when reading the polygon files (default 0)?

        integer(ip) :: i, nl, skip_h

        skip_h = 0
        if(present(skip_header)) skip_h = skip_header

        ! Setup
        nl = size(polygon_files, kind=ip)
        if(size(polygon_values) /= nl) then
            write(log_output_unit, *) 'Error: polygon_files must have the same length as polygon_files'
            call generic_stop
        end if
        allocate(pvt%polyvalue(nl))

        do i = 1, nl
            ! Setup each polygon/value pair
            ! Vertices
            pvt%polyvalue(i)%poly_file = polygon_files(i)
            call read_csv_into_array(pvt%polyvalue(i)%vertices, polygon_files(i), skip_header=skip_h)
            ! Value
            pvt%polyvalue(i)%poly_value = polygon_values(i)

            if(size(pvt%polyvalue(i)%vertices, 1) /= 2) then
                write(log_output_unit, *) "Error in setup_polygons_values_type_from_csv_and_value: file ", trim(polygon_files(i))
                write(log_output_unit, *) "does not seem to be a 2 column csv file"
                call generic_stop
            end if
        end do

    end subroutine

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
        character(len=*), optional, intent(in) :: burn_type
            !! character controlling when we burn -- either 
            !! 'point_value' (default, always burn poly_value into the grid), or 
            !! 'max' (only burn if poly_value > grid value), or 
            !! 'min' (only burn if poly_value < grid_value), or
            !! 'add' (add poly_value to grid_value)

        real(dp) :: dx(2)
        integer(ip) :: i, i0, j0
        character(len=charlen) :: burnt

        burnt = 'point_value'
        if(present(burn_type)) burnt = burn_type

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
                case('add')
                    grid(i0, j0) = grid(i0, j0) + z(i)
                case default
                    write(log_output_unit, *) "Error in burn_xyz_into_grid: unknown value of burn_type ", trim(burnt)
                    call generic_stop
                end select
            end if
        end do

    end subroutine

    subroutine burn_line_into_grid(x, y, z, grid, lower_left, upper_right, burn_type)
        !! Similar to burn_xyz_into_grid, but interpolates along xyz as required to ensure we don't miss cells even if the
        !! line segments have length >> dx.
        !! Note that the interpolation is a bit hap-hazard (in terms of exactly which interpolation point is hit).
        real(dp), intent(in) :: x(:), y(:), z(:) !! Coordinates of the 3D line
        real(dp), intent(inout) :: grid(:,:) !! Grid into which we burn the elevations
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! Coordinates of the grid
        character(len=*), optional, intent(in) :: burn_type
            !! character controlling when we burn -- either 
            !! 'point_value' (default, always burn poly_value into the grid), or 
            !! 'max' (only burn if poly_value > grid value), or 
            !! 'min' (only burn if poly_value < grid_value), or
            !! 'add' (add poly_value to grid_value)

        real(dp) :: dx(2), xi, yi, zi, xip1, yip1, zip1, xj(1), yj(1), zj(1)
        integer(ip) :: i, i0, j0, np, j
        real(force_double) :: n1, n2
        character(len=charlen) :: burnt

        burnt = 'point_value'
        if(present(burn_type))  burnt = burn_type

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
        !! Burn a set of xyz lines into a grid, interpolating along the lines as required.
        class(xyz_lines_type), intent(in) :: xyz_lines !! The xyz lines to burn into the grid
        real(dp), intent(inout) :: grid(:,:) !! The grid
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! The grid extent
        character(len=*), optional, intent(in) :: burn_type
            !! character controlling when we burn -- either 
            !! 'point_value' (default, always burn poly_value into the grid), or 
            !! 'max' (only burn if poly_value > grid value), or 
            !! 'min' (only burn if poly_value < grid_value), or
            !! 'add' (add poly_value to grid_value)

        character(len=charlen) :: burnt
        integer(ip) :: i

        burnt = 'point_value'
        if(present(burn_type)) burnt = burn_type

        do i = 1, size(xyz_lines%lines, kind=ip)
            call burn_line_into_grid(&
                xyz_lines%lines(i)%xyz(1,:), xyz_lines%lines(i)%xyz(2,:), xyz_lines%lines(i)%xyz(3,:), &
                grid, lower_left, upper_right, burnt)
        end do

    end subroutine

    subroutine burn_polygon_value_into_grid(poly_x, poly_y, poly_value, grid, lower_left, upper_right, burn_type)
        !! Set grid points inside the polygon to the poly_value.
        !! This can be used standalone, but for more convenience
        !! in the multi-polygon case see polygons_values_type%burn_into_grid
        real(dp), intent(in) :: poly_x(:), poly_y(:), poly_value !! Polygon vertex coordinates and value to burn
        real(dp), intent(inout) :: grid(:,:) !! The grid -- first dimension ranges from lower_left(1) to upper_right(1).
        real(dp), intent(in) :: lower_left(2), upper_right(2) !! The grid extent (coordinate range for each dimension).
        character(len=*), optional, intent(in) :: burn_type
            !! character controlling when we burn -- either 
            !! 'point_value' (default, always burn poly_value into the grid), or 
            !! 'max' (only burn if poly_value > grid value), or 
            !! 'min' (only burn if poly_value < grid_value), or
            !! 'add' (add poly_value to grid_value)

        character(len=charlen) :: burnt
        real(dp) :: x, y, min_x, max_x, min_y, max_y
        integer :: i, j
        logical :: is_in_poly

        burnt = 'point_value'
        if(present(burn_type)) burnt = burn_type

        if(size(poly_x) /= size(poly_y)) then
            write(log_output_unit, *) "Error in burn_polygon_value_into_grid: unequal lengths of poly_x, poly_y"
            call generic_stop
        end if

        if(.not. all(upper_right > lower_left)) then
            write(log_output_unit, *) "Error in burn_polyvalue_into_grid: upper_right must be > lower_left"
            call generic_stop
        end if

        min_x = minval(poly_x)
        max_x = maxval(poly_x)
        min_y = minval(poly_y)
        max_y = maxval(poly_y)

        ! Quick exit if clearly no overlap
        if(upper_right(2) < min_y .or. lower_left(2) > max_y .or. &
           upper_right(1) < min_x .or. lower_left(1) > max_x) return

        ! Search each grid point
        !$OMP PARALLEL DO DEFAULT(PRIVATE) &
        !$OMP SHARED(poly_x, poly_y, poly_value, grid, upper_right, lower_left, burnt, min_y, min_x, max_y, max_x, log_output_unit)
        do j = 1, size(grid, 2)
            ! Get y coord
            y = lower_left(2) + ((j - 0.5_dp)/size(grid, 2))*(upper_right(2) - lower_left(2))
            ! Quick exit if possible
            if(y < min_y .or. y > max_y) cycle

            do i = 1, size(grid, 1)
                ! Get x coord
                x = lower_left(1) + ((i - 0.5_dp)/size(grid, 1))*(upper_right(1) - lower_left(1))
                ! Quick exit if possible
                if(x < min_x .or. x > max_x) cycle

                ! Determine if x/y is inside
                call point_in_poly(size(poly_x), poly_x, poly_y, x, y, is_in_poly)

                if(is_in_poly) then
                    ! Set the value
                    select case (burnt)
                    case('point_value')
                        grid(i,j) = poly_value
                    case('max')
                        grid(i, j) = max(grid(i, j), poly_value)
                    case('min')
                        grid(i, j) = min(grid(i, j), poly_value)
                    case('add')
                        grid(i, j) = grid(i, j) + poly_value
                    case default
                        write(log_output_unit, *) &
                            "Error in burn_polyvalue_into_grid: unknown value of burn_type ", trim(burnt)
                        call generic_stop
                    end select
                end if
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

    subroutine burn_polygons_values_type_into_grid(pvt, grid, lower_left, upper_right, burn_type)
        !!
        !! For all polygon/value pairs in a polygons_values_type, burn the
        !! value into the grid at points inside the polygon.
        !!
        class(polygons_values_type), intent(in) :: pvt
            !! Holds all polygons and their values
        real(dp), intent(inout) :: grid(:,:)
            !! The grid -- first dimension ranges from lower_left(1) to upper_right(1).
        real(dp), intent(in) :: lower_left(2), upper_right(2)
            !! The grid extent, i.e. coordinate range for the first and second dimensions.
        character(len=*), optional, intent(in) :: burn_type
            !! character controlling when we burn -- either 
            !! 'point_value' (default, always burn poly_value into the grid), or 
            !! 'max' (only burn if poly_value > grid value), or 
            !! 'min' (only burn if poly_value < grid_value), or
            !! 'add' (add poly_value to grid_value)

        integer :: k

        if(.not. allocated(pvt%polyvalue)) then
            write(log_output_unit, *) 'Error: unallocated polyvalue array '
            call generic_stop
        end if

        do k = 1, size(pvt%polyvalue)
            if(.not. allocated(pvt%polyvalue(k)%vertices)) then
                write(log_output_unit, *) 'Error: unallocated polygon vertices at pv%polyvalue(', k, ')'
                call generic_stop
            end if
            call burn_polygon_value_into_grid(pvt%polyvalue(k)%vertices(1,:), pvt%polyvalue(k)%vertices(2,:), &
                pvt%polyvalue(k)%poly_value, &
                grid, lower_left, upper_right, burn_type)
        end do

    end subroutine

    subroutine test_burn_into_grid_mod()
        !! Unit tests
        real(dp) :: grid(10,10), grid_store(10,10)
        real(dp) :: lower_left(2), upper_right(2)
        real(dp) :: x(5), y(5), z(5)
        real(dp) :: expected_grid(10,10)
        real(dp) :: lx(3), ly(3), lz(3)
        real(dp) :: px(4), py(4), pv
        integer(ip) :: i, fid
        real(dp), parameter :: eps = 1.0e-03_dp
        type(polygons_values_type) :: pvt, pvt2

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

        !
        ! Polygon version, test 1
        !
        grid = 0.0_dp
        px = [3.0_dp, 3.0_dp, 6.0_dp, 4.0_dp]
        py = [2.0_dp, 5.0_dp, 5.0_dp, 2.0_dp]
        pv = 1.0_dp
        call burn_polygon_value_into_grid(px, py, pv, grid, lower_left, upper_right, 'point_value')

        ! Plot this situation to see where the tests below come from.
        if(abs(sum(grid) - 6.0_dp) < 3*spacing(6.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        if(all(grid(4:4, 3:3) == pv) .and. all(grid(4:5, 4) == pv) .and. all(grid(4:6, 5) == pv)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Write the above poly to a file for later testing
        open(newunit=fid, file='polygon_test_file_1.csv', form='formatted', action='write')
        write(fid, *) 'lon, lat'
        do i = 1, size(px)
            write(fid, *) px(i), ', ', py(i)
        end do
        close(fid)

        !
        ! Polygon version, test 2 -- swap 'px' and 'py'
        !
        grid = 0.0_dp
        py = [3.0_dp, 3.0_dp, 6.0_dp, 4.0_dp]
        px = [2.0_dp, 5.0_dp, 5.0_dp, 2.0_dp]
        pv = 1.0_dp
        call burn_polygon_value_into_grid(px, py, pv, grid, lower_left, upper_right, 'point_value')

        ! Plot this situation to see where the tests below come from.
        if(abs(sum(grid) - 6.0_dp) < 3*spacing(6.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        if(all(grid(3:3, 4:4) == pv) .and. all(grid(4, 4:5) == pv) .and. all(grid(5, 4:6) == pv)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if
        
        ! Write the above poly to a file for later testing
        open(newunit=fid, file='polygon_test_file_2.csv', form='formatted', action='write')
        write(fid, *) 'lon, lat'
        do i = 1, size(px)
            write(fid, *) px(i), ', ', py(i)
        end do
        close(fid)
       
        !
        ! Polygon version, test 3 -- use the polygons_values_type with 2 polygons
        !
        ! Provide the polygon files and values in 1D arrays
        call pvt%setup([character(len=charlen) :: 'polygon_test_file_1.csv', 'polygon_test_file_2.csv'], &
                       [2.0_dp, 1.0_dp], skip_header=1_ip)
        grid = 0.0_dp
        call pvt%burn_into_grid(grid, lower_left, upper_right, burn_type='max') 

        grid_store = grid ! For later comparison in test 5

        ! Plot this situation to see where the tests below come from.
        if(abs(sum(grid) - 14.0_dp) < 3*spacing(14.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Sites covered by polygon 1 should have the value 2
        if(all(grid(4:4, 3:3) == 2.0_dp) .and. all(grid(4:5, 4) == 2.0_dp) .and. all(grid(4:6, 5) == 2.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if
        ! Since there is overlap of polygon 1 and 2, only these grid values should be 1.0
        if(grid(3,4) == 1.0_dp .and. grid(5,6) == 1.0_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if


        !
        ! Polygon version, test 4 -- use 'add' burning
        !
        grid = 0.0_dp
        call pvt%burn_into_grid(grid, lower_left, upper_right, burn_type='add') 
        ! Plot this situation to see where the tests below come from.
        if(abs(sum(grid) - 18.0_dp) < 3*spacing(18.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Sites covered by both polygons should have the value 3.0
        if(all(grid(4:5, 4:5) == 3.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Remaining polygon 1 sites should have the value 2.0
        if(grid(4,3) == 2.0_dp .and. grid(6,5) == 2.0_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Remaining polygon 2 sites should have the value 1.0
        if(grid(3,4) == 1.0_dp .and. grid(5,6) == 1.0_dp) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Test 5 - a different interface for reading in the polygon value pairs
        !

        ! Provide them in a separate csv file
        open(newunit=fid, file='polygon_value_pairs_list.csv', form='formatted', action='write')
        write(fid, "(A)") 'File, value'
        write(fid, "(A)") "polygon_test_file_1.csv, 2.0"
        write(fid, "(A)") "polygon_test_file_2.csv, 1.0"
        close(fid)

        ! Result should be the same as the other version
        call pvt2%setup_from_csv_with_file_value_columns(&
            file_value_csv='polygon_value_pairs_list.csv', &
            skip_header_file_value_csv=1_ip, &
            skip_header_polygons = 1_ip, &
            verbose=.FALSE.)

        grid = 0.0_dp
        call pvt2%burn_into_grid(grid, lower_left, upper_right, burn_type='max') 

        if(all(grid == grid_store)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        end subroutine

end module
