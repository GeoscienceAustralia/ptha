module reshape_array_mod
    !! Module to efficiently reshape arrays, without calling the 'reshape' function, which seems relatively expensive.
    !!
    !! Useful for coarray_point2point_comms_mod, which uses rank1 buffers, but can be generalised
    !! to other ranks with this module

    use global_mod, only: dp, ip, charlen
    use stop_mod, only: generic_stop
    implicit none

    private
    public :: flatten_array ! Convert nd arrays to 1d [n=1,2,3,4]
    public :: repack_rank1_array ! Convert rank1 arrays to rankn [n=1,2,3,4]
    public :: test_reshape_array_mod ! Unit tests

    interface flatten_array
        !! Subroutine to copy a rank-n array into a rank-1 array.
        !! This can be faster than using the intrinsic 'reshape' function.
        ! HOWEVER, the intrinsic 'pack' function might be better?
        module procedure flatten_array_rank4, flatten_array_rank3, &
            flatten_array_rank2, flatten_array_rank1
    end interface

    interface repack_rank1_array
        !! Copy a rank-1 array into an array with some other rank but the same size.
        !! This can be faster than using the intrinsic 'reshape' function.
        module procedure repack_rank1_array_rank1, repack_rank1_array_rank2, &
            repack_rank1_array_rank3, repack_rank1_array_rank4
    end interface

    contains

    !
    ! Useful routines to repack rank-n arrays to rank1
    !
    pure subroutine flatten_array_rank4(rank4_array, output_rank1)
        real(dp), intent(in) :: rank4_array(:,:,:,:)
        real(dp), intent(inout) :: output_rank1(:)

        integer(ip) :: i, j, k, l, counter, r4shape(4)

        counter = 1
        r4shape = shape(rank4_array)
        do l = 1, r4shape(4)
            do k = 1, r4shape(3)
                do j = 1, r4shape(2)
                    do i = 1, r4shape(1)
                        output_rank1(counter) = rank4_array(i,j,k,l)
                        counter = counter + 1
                    end do
                end do
            end do
        end do

    end subroutine

    pure subroutine flatten_array_rank3(rank3_array, output_rank1)
        real(dp), intent(in) :: rank3_array(:,:,:)
        real(dp), intent(inout) :: output_rank1(:)

        integer(ip) :: i, j, k, counter, r3shape(3)

        counter = 1
        r3shape = shape(rank3_array)
        do k = 1, r3shape(3)
            do j = 1, r3shape(2)
                do i = 1, r3shape(1)
                    output_rank1(counter) = rank3_array(i,j,k)
                    counter = counter + 1
                end do
            end do
        end do

    end subroutine

    pure subroutine flatten_array_rank2(rank2_array, output_rank1)
        real(dp), intent(in) :: rank2_array(:,:)
        real(dp), intent(inout) :: output_rank1(:)

        integer(ip) :: i, j, counter, r2shape(2)

        counter = 1
        r2shape = shape(rank2_array)
        do j = 1, r2shape(2)
            do i = 1, r2shape(1)
                output_rank1(counter) = rank2_array(i,j)
                counter = counter + 1
            end do
        end do

    end subroutine

    pure subroutine flatten_array_rank1(rank1_array, output_rank1)
        real(dp), intent(in) :: rank1_array(:)
        real(dp), intent(inout) :: output_rank1(:)

        output_rank1 = rank1_array

    end subroutine


    !
    ! Routines to repack rank-1 arrays into rank-n
    !

    pure subroutine repack_rank1_array_rank1(rank1_array, output_rank1)
        real(dp), intent(in) :: rank1_array(:)
        real(dp), intent(inout) :: output_rank1(:)

        output_rank1 = rank1_array

    end subroutine

    pure subroutine repack_rank1_array_rank2(rank1_array, output_rank2)
        real(dp), intent(in) :: rank1_array(:)
        real(dp), intent(inout) :: output_rank2(:,:)

        integer(ip) :: i, j, output_shape(2), counter

        output_shape = shape(output_rank2)

        !if(product(output_shape) /= size(rank1_array)) then
        !    stop('Non-conforming dimensions in repack_rank1_array_rank2')
        !end if

        counter = 1
        do j = 1, output_shape(2)
            do i = 1, output_shape(1)
                output_rank2(i, j) = rank1_array(counter)
                counter = counter+1
            end do
        end do

    end subroutine

    pure subroutine repack_rank1_array_rank3(rank1_array, output_rank3)
        real(dp), intent(in) :: rank1_array(:)
        real(dp), intent(inout) :: output_rank3(:,:,:)

        integer(ip) :: i, j, k, output_shape(3), counter

        output_shape = shape(output_rank3)

        !if(product(output_shape) /= size(rank1_array)) then
        !    stop('Non-conforming dimensions in repack_rank1_array_rank3')
        !end if

        counter = 1
        do k = 1, output_shape(3)
            do j = 1, output_shape(2)
                do i = 1, output_shape(1)
                    output_rank3(i, j, k) = rank1_array(counter)
                    counter = counter+1
                end do
            end do
        end do
    end subroutine

    pure subroutine repack_rank1_array_rank4(rank1_array, output_rank4)
        real(dp), intent(in) :: rank1_array(:)
        real(dp), intent(inout) :: output_rank4(:,:,:,:)

        integer(ip) :: i, j, k, l, output_shape(4), counter

        output_shape = shape(output_rank4)

        !if(product(output_shape) /= size(rank1_array)) then
        !    stop('Non-conforming dimensions in repack_rank1_array_rank4')
        !end if

        counter = 1

        do l = 1, output_shape(4)
            do k = 1, output_shape(3)
                do j = 1, output_shape(2)
                    do i = 1, output_shape(1)
                        output_rank4(i, j, k, l) = rank1_array(counter)
                        counter = counter+1
                    end do
                end do
            end do
        end do

    end subroutine

    subroutine test_reshape_array_mod
        !!
        !! Unit tests
        !!
        real(dp) :: r1(10), r2(10,9), r3(10,9,8), r4(10, 9, 8, 7)
        real(dp) :: a1(10), a2(90),   a3(720),    a4(5040)
        integer(ip) :: i
        integer(ip):: newshp1(1)
        integer(ip), parameter:: newshp2(2) = [10,9], newshp3(3) = [10, 9, 8], newshp4(4) = [10, 9, 8, 7]

        ! Create some data
        a1 = (/ (i*1.0_dp, i=1, size(a1)) /)
        a2 = (/ (i*1.0_dp, i=1, size(a2)) /)
        a3 = (/ (i*1.0_dp, i=1, size(a3)) /)
        a4 = (/ (i*1.0_dp, i=1, size(a4)) /)

        ! Do the repack
        call repack_rank1_array(a1, r1)
        call repack_rank1_array(a2, r2)
        call repack_rank1_array(a3, r3)
        call repack_rank1_array(a4, r4)

        ! Compare with 'reshape', except in the case of rank1
        if(all(r1 == a1)) then
            print*, 'PASS'
        else
            print*, 'FAIL r1'
            call generic_stop()
        end if

        !newshp2 = [10, 9]
        if(all(r2 == reshape(a2, newshp2) )) then
            print*, 'PASS'
        else
            print*, 'FAIL r2'
            call generic_stop()
        end if

        !newshp3 = [10, 9, 8]
        if(all(r3 == reshape(a3, newshp3) )) then
            print*, 'PASS'
        else
            print*, 'FAIL r3'
            call generic_stop()
        end if

        !newshp4 = [10, 9, 8, 7]
        if(all(r4 == reshape(a4, newshp4) )) then
            print*, 'PASS'
        else
            print*, 'FAIL r4'
            call generic_stop()
        end if

        ! Now flatten arrays again
        call flatten_array(r1, a1)
        call flatten_array(r2, a2)
        call flatten_array(r3, a3)
        call flatten_array(r4, a4)

        ! Check it worked
        if(all(r1 == a1)) then
            print*, 'PASS'
        else
            print*, 'FAIL a1'
            call generic_stop
        end if

        newshp1 = [size(a2)]
        if(all(reshape(r2, newshp1) == a2)) then
            print*, 'PASS'
        else
            print*, 'FAIL a2'
            call generic_stop
        end if

        newshp1 = [size(a3)]
        if(all(reshape(r3, newshp1) == a3)) then
            print*, 'PASS'
        else
            print*, 'FAIL a3'
            call generic_stop
        end if

        newshp1 = [size(a4)]
        if(all(reshape(r4, newshp1) == a4)) then
            print*, 'PASS'
        else
            print*, 'FAIL a4'
            call generic_stop
        end if

    end subroutine

end module


