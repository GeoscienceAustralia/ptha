

module qsort_mod
    !!
    !! Module for sorting.
    !!
    !! Contains a generic interface
    !!   sort(x, size(x))
    !! which changes x to be in sorted order. 
    !!
    !! It also contains a generic interface
    !!   sort_index(inds, x, size(x))
    !! which puts indices in 'inds' such that x(inds) is sorted. Here it is required
    !! that "size(inds) == size(x)"
    !!
    !! It also contains a generic match interface
    !!   match(x1, x2, matches)
    !! which puts indices in matches such that x1(i) == x2(matches(i)), or else
    !! matches(i) == -1. 
    !!
    !! The routines work with 'x' being c_int, c_float, c_double, or character(len=*)
    !!
    !! There is also an iso_c_binding interface to C's qsort, which allows sorting anything with a user-defined compare function.
    !! Originally the sorting routines were based on this, using locally defined compare functions (generally contained inside the
    !! subroutine that called qsort). However this approach was not supported with PGI the compiler, so was changed.
    use iso_c_binding

    implicit none

    private
    public :: sort, test_qsort_mod, sort_index, match

    interface

        ! Call qsort from C -- this gave me trouble with PGI compiler when used with 
        ! functions defined in contains blocks of other subroutines, so is not used below, 
        ! but is useful as one can pass arbitrary types and "compare" functions to it.
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

    subroutine sort_cint(array, n)
        integer, intent(in) :: n
        integer(c_int), intent(inout) :: array(n)

        integer(c_int), allocatable :: inds(:)

        allocate(inds(n))
        call sort_index_cint(inds, array, n)
        array = array(inds)
        deallocate(inds)       

    end subroutine

    !
    ! Sort for c_float
    ! 
    subroutine sort_cfloat(array, n)
        integer, intent(in) :: n
        real(c_float), intent(inout) :: array(n)

        integer(c_int), allocatable :: inds(:)

        allocate(inds(n))
        call sort_index_cfloat(inds, array, n)
        array = array(inds)
        deallocate(inds)       

    end subroutine

    !
    ! Sort for c_double
    ! 
    subroutine sort_cdouble(array, n)
        integer, intent(in) :: n
        real(c_double), intent(inout) :: array(n)

        integer(c_int), allocatable :: inds(:)

        allocate(inds(n))
        call sort_index_cdouble(inds, array, n)
        array = array(inds)
        deallocate(inds)       

    end subroutine


    ! For sort character, we use 'sort_index' to work around the complexity of passing fortran strings to C
    subroutine sort_character(array, n)
        integer, intent(in) :: n
        character(len=*), intent(inout) :: array(n)

        integer(c_int), allocatable :: inds(:)

        allocate(inds(n))
        call sort_index_character(inds, array, n)
        array = array(inds)
        deallocate(inds)       

    end subroutine

    subroutine sort_index_character(inds, array, n)
#define SORT_INDEX_TYPE character(len=*)
#define SORT_INDEX_TYPE2 character(len=len(array))
#include "sort_index_template2.f90"
#undef SORT_INDEX_TYPE
#undef SORT_INDEX_TYPE2

    end subroutine

    !
    ! Index sort of c_int
    !
    subroutine sort_index_cint(inds, array, n)
#define SORT_INDEX_TYPE integer(c_int)
#define SORT_INDEX_TYPE2 integer(c_int)
#include "sort_index_template2.f90"
#undef SORT_INDEX_TYPE2
#undef SORT_INDEX_TYPE

    end subroutine 

    !
    ! Index sort of c_float
    !
    subroutine sort_index_cfloat(inds, array, n)
#define SORT_INDEX_TYPE real(c_float)
#define SORT_INDEX_TYPE2 real(c_float)
#include "sort_index_template2.f90"
#undef SORT_INDEX_TYPE2
#undef SORT_INDEX_TYPE

    end subroutine 

    !
    ! Index sort of c_double
    !
    subroutine sort_index_cdouble(inds, array, n)
#define SORT_INDEX_TYPE real(c_double)
#define SORT_INDEX_TYPE2 real(c_double)
#include "sort_index_template2.f90"
#undef SORT_INDEX_TYPE2
#undef SORT_INDEX_TYPE

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
        integer(c_int) :: cx_matches(10), counter
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
        if(all(x(inds) - x(index_order_x) == 0)) then
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
        if(all(y(inds) - y(index_order_x) == 0)) then
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
        if(all(z(inds) - z(index_order_x) == 0)) then
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

        ! Check we can apply sort to a constant array, e.g. [2,2,2,2,2 ... ]
        x(1:10) = 2
        call sort(x, size(x))
        if(all(x == 2)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Check we can apply sort_index to a constant array, e.g. [2,2,2,2,2 ... ]
        call sort_index(inds, x, size(x))
        ! Test that all numbers from 1-10 appeared
        counter=0
        do i = 1, size(x)
            if(any(inds == i)) counter=counter+1
        end do
        if(counter == size(x)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if



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
        allocate(matches(size(match1)))
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
        if(all(rns2_int(matches) == rns_int)) then
            print*, "PASS"
        else
            print*, 'FAIL'
        end if


        ! test as above, using c_float inputs
        call sort_index(rns_inds, rns, n)
        call sort(rns, size(rns))

        ! This 'match' should be equivalent to getting the sorted indices
        call match(rns, rns2, matches)
        if(all(rns2(matches) == rns)) then
            print*, "PASS"
        else
            print*, 'FAIL'
        end if

        ! test as above, using c_double inputs
        allocate(rns_d(n), rns2_d(n))
        call random_number(rns_d)
        rns2_d = rns_d
        call sort_index(rns_inds, rns_d, n)
        call sort(rns_d, size(rns_d))

        ! This 'match' should be equivalent to getting the sorted indices
        call match(rns_d, rns2_d, matches)
        if(all(rns2_d(matches) == rns_d)) then
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

    end subroutine
end module
