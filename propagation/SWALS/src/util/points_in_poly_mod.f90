module points_in_poly_mod
    !
    !! Module to find whether points are inside (or on the boundary of) a polygon
    !
    !! The tests suggest being on the boundary is counted as being inside, although
    !! I would not rely on that behaviour always being true for complex polygons.
    !!
    !! Main routine is
    !!   points_in_poly(vertx, verty, px, py, is_inside)
    !! where all arguments are arrays
    !! (with size(vertx) == size(verty) & (size(px) == size(py) == size(is_inside))
    !!
    !! A discussion of the problem is here (we use a similar algorithm)
    !! https://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
    !
    use global_mod, only: ip, dp
    implicit none

    contains

    ! Find if a single point (px, py) is inside a poly defined by nvert vertices
    ! (vertx, verty).
    !
    ! @param nvert number of vertices in the polygon
    ! @param vertx x coordinates of polygon vertices, in order
    ! @param verty y coordinates of polygon vertices in order
    ! @param px x coordinates of points which might be in or outside the polygon
    ! @param py y coordinates of points which might be in or outside the polygon
    ! @param is_inside output logical array which indicates if each px,py is inside the polygon
    !
    pure subroutine point_in_poly(nvert, vertx, verty, px, py, is_inside)
        integer(ip), intent(in):: nvert
        real(dp), intent(in):: vertx(nvert), verty(nvert), px, py
        logical, intent(out):: is_inside

        integer(ip) :: i, j
        real(dp) :: g0

        is_inside = .FALSE.

        do i = 1, nvert
            ! j is an index of the next vertex
            if(i == 1) then
                j = nvert
            else
                j = i - 1
            end if

            ! Main test
            if((verty(i) > py) .neqv. (verty(j) > py)) then
                g0 = (vertx(j) - vertx(i)) * (py - verty(i)) / (verty(j) - verty(i)) + vertx(i)
                if(px < g0) then
                    is_inside = (.NOT. is_inside)
                end if
            end if

        end do
    end subroutine

    ! Generalise point_in_poly to many points
    !
    ! Use a bounding box check which could speed up some cases.
    !
    ! @param vertx x coordinates of polygon vertices, in order
    ! @param verty y coordinates of polygon vertices in order
    ! @param px x coordinates of points which might be in or outside the polygon
    ! @param py y coordinates of points which might be in or outside the polygon
    ! @param is_inside output logical array which indicates if each px,py is inside the polygon
    !
    pure subroutine points_in_poly(vertx, verty, px, py, is_inside)
        real(dp), intent(in):: vertx(:), verty(:), px(:), py(:)
        logical, intent(out):: is_inside(:)

        integer(ip):: i, nvert, np

        ! bounding box for polygon to speed up check
        real(dp):: min_vertx, max_vertx, min_verty, max_verty

        nvert = size(vertx, kind=ip)
        np = size(px, kind=ip)

        min_vertx = minval(vertx)
        max_vertx = maxval(vertx)
        min_verty = minval(verty)
        max_verty = maxval(verty)

        do i = 1, np
            if(px(i) >= min_vertx .and. px(i) <= max_vertx .and. &
               py(i) >= min_verty .and. py(i) <= max_verty) then
                call point_in_poly(nvert, vertx, verty, px(i), py(i), is_inside(i))
            else
                is_inside(i) = .FALSE.
            end if
        end do

    end subroutine


    ! Test. Call as test_points_in_poly_mod()
    !
    subroutine test_points_in_poly_mod
        integer(ip):: nvert

        ! doesn't seem to matter whether last point == first point or not
        real(dp), dimension(5) :: vertx = [0._dp, 0.0_dp, 1.0_dp, 1.0_dp, 0._dp]
        real(dp), dimension(5) :: verty = [0._dp, 1.0_dp, 1.0_dp, 0.0_dp, 0._dp]
        !real(dp), dimension(4) :: vertx = [0._dp, 0.0_dp, 1.0_dp, 1.0_dp]
        !real(dp), dimension(4) :: verty = [0._dp, 1.0_dp, 1.0_dp, 0.0_dp]

        integer(ip):: i
        logical:: is_inside_result

        type test_point_type
            real(dp):: p(2)
            logical:: is_inside
        end type

        type(test_point_type):: test_points(10)

        real(dp):: test_points_array(10, 2)
        logical:: is_inside_result_array(10), is_inside_array(10)

        nvert = size(vertx, kind=ip)

        ! Manually make test points and flag if they lie on a boundary
        !
        ! It seems the routine counts being 'on' the boundary as being inside
        ! (although that may vary in more complex cases due to round-off?)
        !
        test_points(1) = test_point_type([-1.0_dp, 0.0_dp], .FALSE.)
        test_points(2) = test_point_type([-1.0_dp, 0.5_dp], .FALSE.)
        test_points(3) = test_point_type([ 0.0_dp, 0.0_dp], .TRUE.)
        test_points(4) = test_point_type([ 0.0_dp, 0.5_dp], .TRUE.)
        test_points(5) = test_point_type([ 0.5_dp, 0.5_dp], .TRUE.)
        test_points(6) = test_point_type([ 0.1_dp, 0.9_dp], .TRUE.)
        test_points(7) = test_point_type([ 0.1_dp, 2.9_dp], .FALSE.)
        test_points(8) = test_point_type([ 0.1_dp, 0.0_dp], .TRUE.)
        test_points(9) = test_point_type([ 1.0e-10_dp, 1.0e-10_dp], .TRUE.)
        test_points(10) = test_point_type([ -1.0e-10_dp, -1.0e-10_dp], .FALSE.)

        do i = 1, 10

            call point_in_poly(size(vertx)*1_ip, vertx, verty, &
                test_points(i)%p(1), test_points(i)%p(2), is_inside_result)

            if(is_inside_result .eqv. test_points(i)%is_inside) then
                print*, 'PASS'
            else
                print*, 'FAIL on point ', i, ' : ', test_points(i)
            end if

        end do

        ! Test the array version
        do i = 1, 10
            test_points_array(i,:) = test_points(i)%p
            is_inside_array(i) = test_points(i)%is_inside
        end do

        call points_in_poly(vertx, verty, test_points_array(:,1), &
            test_points_array(:,2), is_inside_result_array)

        if(all(is_inside_result_array .eqv. is_inside_array)) then
            print*, 'PASS'
        else
            print*, 'FAIL on array test'
        end if

    end subroutine

end module


