

module qsort_mod
    !
    !! Module for sorting, based on C's qsort
    !
    !! Contains a generic subroutine 
    !!   sort(x, size(x))
    !! which changes x to be in sorted order. 

    !! It also contains a generic subroutine
    !!   sort_index(inds, x, size(x))
    !! which puts indices in 'inds' such that x(inds) is sorted. Here it is required
    !! that "size(inds) == size(x)"

    !! It also contains a generic match subroutine
    !!   match(x1, x2, matches)
    !! which puts indices in matches such that x1(i) == x2(matches(i)), or else
    !! matches(i) == -1. 
    !
    !! The routines work with 'x' being c_int, c_float, or c_double
    !
    !! Note -- sorts of large arrays (e.g. 10^7 elements) may fail if compiled with '-Ofast',
    !! unless the stack-size is unlimited (ulimit -s unlimited).  On the other hand,
    !! I do not get this with the less aggressive compiler optimization '-O3'
    !
    use iso_c_binding

    implicit none

    private
    public :: sort, test_qsort_mod, sort_index, match

    interface

        ! Call qsort from C
        subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
          import c_ptr, c_size_t, c_funptr
          type(c_ptr), value       :: array
          integer(c_size_t), value :: elem_count
          integer(c_size_t), value :: elem_size
          type(c_funptr), value    :: compare !int(*compare)(const void *, const void *)
        end subroutine qsort !standard C library qsort

    end interface

    ! Simple 'sort', e.g. sort(x, size(x))
    interface sort
        module procedure sort_cint, sort_cfloat, sort_cdouble, sort_character
    end interface

    ! Return indices in 'sorted' order, e.g. sort_index(inds, x, size(x))
    interface sort_index
        module procedure sort_index_cint, sort_index_cfloat, sort_index_cdouble, sort_index_character
    end interface

    interface match
        module procedure match_cint, match_cfloat, match_cdouble, match_string
    end interface

    contains

    !
    ! Sort for c_int
    ! 

    pure function compare_int(i1, i2) result(compar) bind(C)
        integer(c_int), intent(in) :: i1, i2
        integer(c_int) :: compar

        if( i1 > i2 ) compar = 1
        if( i1 == i2 ) compar = 0
        if( i1 < i2 )  compar = -1

    end function

    subroutine sort_cint(array, n)
        integer, intent(in) :: n
        integer(c_int), intent(inout), target :: array(n)
       
        integer(c_size_t) :: elem_count, elem_size 

        elem_count = size(array)
        elem_size = c_sizeof(array(1))

        call qsort(c_loc(array(1)), elem_count, elem_size, c_funloc(compare_int))

    end subroutine

    !
    ! Sort for c_float
    ! 

    pure function compare_float(i1, i2) result(compar) bind(C)
        real(c_float), intent(in) :: i1, i2
        integer(c_int) :: compar

        if( i1 > i2 ) compar = 1
        if( i1 == i2 ) compar = 0
        if( i1 < i2 )  compar = -1

    end function

    subroutine sort_cfloat(array, n)
        integer, intent(in) :: n
        real(c_float), intent(inout), target :: array(n)
       
        integer(c_size_t) :: elem_count, elem_size 

        elem_count = size(array)
        elem_size = c_sizeof(array(1))

        call qsort(c_loc(array(1)), elem_count, elem_size, c_funloc(compare_float))

    end subroutine

    !
    ! Sort for c_double
    ! 

    pure function compare_double(i1, i2) result(compar) bind(C)
        real(c_double), intent(in) :: i1, i2
        integer(c_int) :: compar

        if( i1 > i2 ) compar = 1
        if( i1 == i2 ) compar = 0
        if( i1 < i2 )  compar = -1

    end function

    subroutine sort_cdouble(array, n)
        integer, intent(in) :: n
        real(c_double), intent(inout), target :: array(n)
       
        integer(c_size_t) :: elem_count, elem_size 

        elem_count = size(array)
        elem_size = c_sizeof(array(1))

        call qsort(c_loc(array(1)), elem_count, elem_size, c_funloc(compare_double))
    end subroutine

    ! For sort character, we use 'sort_index' to work around the complexity of passing fortran strings to C
    subroutine sort_character(array, n)
        integer(c_int), intent(in) :: n
        character(len=*), intent(inout) :: array(n)
       
        integer(c_int), allocatable :: index_array(:)

        allocate(index_array(n))
        call sort_index_character(index_array, array, n)
        array = array(index_array)
    
    end subroutine

    subroutine sort_index_character(inds, array, n)
        integer(c_int), intent(in) :: n
        character(len=*), intent(in) :: array(n)
        integer(c_int), intent(inout), target :: inds(n)

        integer(c_size_t) :: elem_count, elem_size 
        integer(c_int) :: i

#include "qsort_sort_index_include.inc"


    end subroutine

    !
    ! Index sort of c_int
    !
    subroutine sort_index_cint(inds, array, n)
        integer(c_int), intent(in) :: n
        integer(c_int), intent(in) :: array(n)
        integer(c_int), intent(inout), target :: inds(n)

        integer(c_size_t) :: elem_count, elem_size 
        integer(c_int) :: i

#include "qsort_sort_index_include.inc"

    end subroutine 

    !
    ! Index sort of c_float
    !
    subroutine sort_index_cfloat(inds, array, n)
        integer(c_int), intent(in) :: n
        real(c_float), intent(in) :: array(n)
        integer(c_int), intent(inout), target :: inds(n)

        integer(c_size_t) :: elem_count, elem_size 
        integer(c_int) :: i

#include "qsort_sort_index_include.inc"

    end subroutine 

    !
    ! Index sort of c_double
    !
    subroutine sort_index_cdouble(inds, array, n)
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: array(n)
        integer(c_int), intent(inout), target :: inds(n)

        integer(c_size_t) :: elem_count, elem_size 
        integer(c_int) :: i

#include "qsort_sort_index_include.inc"

    end subroutine 

    !
    ! Emulate R's 'match' routine for c_int. 
    ! 
    ! This implementation avoids an O(n*n) search by sorting both arrays first.
    !
    ! @param array1 array with values that we try to find in array2
    ! @param array2 array where we try to find matches from array1
    ! @param matches output array giving, for each array1 entry, the index of a
    !     matching entry in array2, (or else -1 for no match)
    !
    subroutine match_cint(array1, array2, matches)
        integer(c_int), intent(in) :: array1(:), array2(:)
        integer(c_int), intent(inout) :: matches(:)

        integer(c_int), allocatable :: inds1(:), inds2(:)
        integer(c_int) :: n, i, i2

#include "qsort_match_include.inc"
                
    end subroutine

    subroutine match_cfloat(array1, array2, matches)
        real(c_float), intent(in) :: array1(:), array2(:)
        integer(c_int), intent(inout) :: matches(:)

        integer(c_int), allocatable :: inds1(:), inds2(:)
        integer(c_int) :: n, i, i2

#include "qsort_match_include.inc"
                
    end subroutine

    subroutine match_cdouble(array1, array2, matches)
        real(c_double), intent(in) :: array1(:), array2(:)
        integer(c_int), intent(inout) :: matches(:)

        integer(c_int), allocatable :: inds1(:), inds2(:)
        integer(c_int) :: n, i, i2

#include "qsort_match_include.inc"
                
    end subroutine

    subroutine match_string(array1, array2, matches)
        character(len=*), intent(in) :: array1(:), array2(:)
        integer(c_int), intent(inout) :: matches(:)
        integer(c_int), allocatable :: inds1(:), inds2(:)
        integer(c_int) :: n, i, i2

#include "qsort_match_include.inc"

    end subroutine

    !
    ! Unit test
    !
    subroutine test_qsort_mod()

        integer(c_int) :: x(10), sorted_x(10), n, i, index_order_x(10), inds(10)
        real(c_float) :: y(10)
        real(c_double) :: z(10)
        integer(c_int), allocatable :: rns_inds(:), match1(:), match2(:), matches(:)
        real(c_float), allocatable :: rns(:), rns2(:)
        real(c_double), allocatable :: rns_d(:), rns2_d(:)
        integer(c_int), allocatable :: rns_int(:), rns2_int(:)
        character(len=12) :: cx(10), cx_sorted(10), cx_match(4)
        integer(c_int) :: cx_matches(10)
        logical :: failed


        ! Sort of a string
        x = [1, 3, 2, 5, 4, 7, 6, 9, 8, 1] 

        cx = ['asdf', 'fffd', 'reec', 'adfa', 'wert', 'asfa', 'ddhd', 'dghd', 'rtyu', 'rurt']
        cx_sorted = cx
        call sort(cx_sorted, 10)
        call sort_index(x, cx, 10)
        if(all(cx(x) == cx_sorted)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

      
        ! Sort of an integer 
        x = [1, 3, 2, 5, 4, 7, 6, 9, 8, 1] 

        sorted_x = [1, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        index_order_x = [1, 10, 3, 2, 5, 4, 7, 6, 9, 8]

        ! Use these variables later
        y = x * 2.0
        z = x * 3.0d0
       
        call sort_index(inds, x, 10)
        if(all(inds - index_order_x == 0)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        call sort(x, 10)
        
        if( maxval(abs(x - sorted_x)) == 0) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif


        ! Sort of c_float
        call sort_index(inds, y, 10)
        if(all(inds - index_order_x == 0)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        call sort(y, 10)

        if( maxval(abs(y - sorted_x*2.0)) == 0) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        ! Sort of c_double 
        call sort_index(inds, z, 10)
        if(all(inds - index_order_x == 0)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        call sort(z, size(z))

        if( maxval(abs(z - sorted_x*3.0)) == 0) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        !
        ! Sort a large random vector
        !
        call random_seed() ! This is deterministic
        n = 100000
        allocate(rns(n), rns2(n), rns_inds(n))
        call random_number(rns)

        ! Copy
        rns2 = rns
        call sort_index(rns_inds, rns, n)
        call sort(rns, size(rns))

        failed = .FALSE.
        do i = 1, (n-1)
            if(rns(i) > rns(i+1)) then
                failed = .TRUE.
                exit
            endif 
        end do

        if(failed) then
            print*, 'FAIL'
        else
            print*, 'PASS'
        end if
        

        if(any(rns /= rns2(rns_inds))) then
            print*, 'FAIL'
        else
            print*, 'PASS'
        endif

        !
        ! Test 'match' (simple case)
        !
        match1 = [1, -1, 1, 3, 5, 8, 9, 11, 0, 0]
        match2 = [1, 1, -2, 5, 90, 0]
        matches = 0 * match1 

        call match_cint(match1, match2, matches)
        if(all(matches == [1, -1, 1, -1, 4, -1, -1, -1, 6, 6])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Large scale test of 'match'
        ! 
        ! Make a large random vector. Then sort it and match with the original. Then 
        ! compare the match against the results from sort_index, which clearly should
        ! be the same
        !
        call random_seed() ! This is deterministic
        n = 1000
        deallocate(rns, rns2, rns_inds, matches)
        allocate(rns(n), rns2(n), rns_inds(n), matches(n))
        call random_number(rns)

        ! Copy
        rns2 = rns

        ! test using c_int inputs
        allocate(rns_int(n), rns2_int(n))
        rns_int = int(rns * 1e+08, c_int)
        rns2_int = int(rns2 * 1e+08, c_int)

        call sort_index(rns_inds, rns_int, n)
        call sort(rns_int, size(rns_int))

        ! This 'match' should be equivalent to getting the sorted indices
        call match(rns_int, rns2_int, matches)
        if(all(matches == rns_inds)) then
            print*, "PASS"
        else
            print*, 'FAIL'
        end if


        ! test as above, using c_float inputs
        call sort_index(rns_inds, rns, n)
        call sort(rns, size(rns))

        ! This 'match' should be equivalent to getting the sorted indices
        !! FIXME: rare failures? Perhaps for some particular random values?
        !! Is it because of repeated values?
        call match(rns, rns2, matches)
        if(all(matches == rns_inds)) then
            print*, "PASS"
        else
            print*, 'FAIL'
        end if

        ! test as above, using c_double inputs
        allocate(rns_d(n), rns2_d(n))
        call random_number(rns_d)
        rns2_d = rns_d
        call sort_index(rns_inds, rns_d, n)
        call sort(rns_d, size(rns))

        ! This 'match' should be equivalent to getting the sorted indices
        !! FIXME: rare failures? Perhaps for some particular random values?
        !! Is it because of repeated values?
        call match(rns_d, rns2_d, matches)
        if(all(matches == rns_inds)) then
            print*, "PASS"
        else
            print*, 'FAIL'
        end if


        ! test as above, using character inputs
        cx_match(1:3) = ['fffd', 'asdf', 'rtyu']
        cx_match(4) = 'asdfasdfaf'
        call match(cx, cx_match, cx_matches)
        if(all(cx_matches == [ 2, 1, -1, -1, -1, -1, -1, -1, 3, -1])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        cx_matches = -1
        call match(cx_match, cx, cx_matches(1:4))
        if(all(cx_matches(1:4) == [ 2, 1, 9, -1])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        endif

        ! FIXME: Check we can apply sort_index to a constant array, e.g. [2,2,2,2,2]
        ! I think it has issues

    end subroutine
end module
