module which_mod
    !!
    !! Emulates a few R routines: 'which', 'rle', 'cumsum', ...
    !!
    use global_mod, only: ip

    implicit none
    
    private
    public which, test_which, rle_ip, cumsum_ip, bind_arrays_ip, remove_rows_ip

    interface 
        !! When generalising rle, it is helpful for the user to pass
        !! a function defining whether two objects (indexed by integers)
        !! are equal
        function eq_fun(i1, i2) result(is_equal)
            import ip
            implicit none
            integer(ip), intent(in) :: i1, i2
            logical :: is_equal
        end function
    end interface

    contains

    pure subroutine which(logical_array, allocatable_array)
        !! Return an array with integers corresponding to indices where logical_array is .true. .
        !! Like 'which' in R.
        logical, intent(in) :: logical_array(:)
        integer(ip), allocatable, intent(out):: allocatable_array(:)
        integer(ip):: arr_len, i, counter

        if(allocated(allocatable_array)) then
            deallocate(allocatable_array)
        end if

        arr_len = count(logical_array)

        allocate(allocatable_array(arr_len))

        if(arr_len > 0) then
            counter = 0
            do i = 1, size(logical_array)
                if(logical_array(i)) then
                    counter = counter + 1
                    allocatable_array(counter) = i
                end if
            end do
        end if

    end subroutine

    subroutine rle_ip(integer_array, values, lengths, equality_function)
        !!
        !! Emulate R's function 'rle', which computes the lengths and values of runs of equal
        !! values in a vector
        integer(ip), intent(in) :: integer_array(:) !! rank 1 array with integers
        integer(ip), allocatable, intent(inout) :: values(:) !! return the values of the integer array
        integer(ip), allocatable, intent(inout) :: lengths(:) !! return the length of the 'runs' of equal values
        procedure(eq_fun), optional :: equality_function 
        !! (optional) function f(i1, i2) which returns
        !! .true. if i1 should be treated as equal to i2, and false otherwise. This
        !! can be used to generalise rle_ip, e.g. we could look for equality of rows in a
        !! matrix by having equality function return all(mymatrix(i1,:) == mymatrix(i2,:))

        integer(ip) :: n, counter, np, i
        procedure(eq_fun), pointer :: eq_fun_local

        if(present(equality_function)) then
            eq_fun_local => equality_function
        else
            eq_fun_local => equality_function_default
        end if
        
        
        if(allocated(values)) deallocate(values)
        if(allocated(lengths)) deallocate(lengths)
        
        n = size(integer_array)

        select case(n)
        case(0)
            ! If integer_array is empty, return length=0 arrays
            allocate(values(0), lengths(0))

        case(1) 
            allocate(values(1), lengths(1))
            values(1) = integer_array(1)
            lengths(1) = 1_ip 

        case default 
            ! n > 1 (typical case)
            np = 1_ip
            do i = 1, n-1
                !if(integer_array(i) /= integer_array(i+1)) np = np + 1_ip
                if(.not. eq_fun_local(i, i+1)) np = np + 1_ip
            end do
            allocate(values(np), lengths(np))

            counter = 1_ip
            np = 0_ip
            do i = 2, n

                !if(integer_array(i-1) /= integer_array(i)) then
                if(.not. eq_fun_local(i-1, i)) then 
                    np = np + 1 
                    values(np) = integer_array(i-1)
                    lengths(np) = counter
                    counter = 0_ip
                end if

                counter = counter + 1_ip

                if(i == n) then 
                    np = np + 1 
                    values(np) = integer_array(n)
                    lengths(np) = counter
                end if

            end do
        end select

        contains

            function equality_function_default(i1, i2) result(is_equal)
                integer(ip), intent(in) :: i1, i2
                logical :: is_equal

                if(integer_array(i1) == integer_array(i2)) then
                    is_equal = .true.
                else
                    is_equal = .false.
                end if
            end function

    end subroutine

    pure subroutine cumsum_ip(integer_array)
        !!
        !! Cumulative sum (for integer(ip))
        !!
        integer(ip), intent(inout) :: integer_array(:)

        integer(ip) :: sumval, n, i

        n = size(integer_array)

        if(n > 0) then
            sumval = 0_ip
            do i = 1, n
                sumval = sumval + integer_array(i)
                integer_array(i) = sumval
            end do
        end if

    end subroutine

    subroutine bind_arrays_ip(ar1, ar2, rowbind)
        !!
        !! 'Bind' rank 2 arrays 'ar1' and 'ar2', either by rows or by columns
        !!
        !! This can be used to emulate R's 'rbind' and 'cbind' functions, which
        !! combine 2-d arrays by rows or by columns
        !!
        integer(ip), intent(inout), allocatable :: ar1(:,:)
        !! allocatable array, which will be modified on output to contain rbind(ar1, ar2), or cbind(ar1, ar2)
        integer(ip), intent(in) :: ar2(:,:) !! as defined above
        logical, optional, intent(in) :: rowbind !! logical indicating whether to do a 'rbind' (TRUE) or 'cbind' (FALSE)

        integer(ip), allocatable :: tmp_ar(:,:)
        logical :: bind_rows 
        integer(ip) :: d1(2), d2(2)
       
        if(present(rowbind)) then
            bind_rows = rowbind
        else
            bind_rows = .TRUE.
        end if 

        d1 = shape(ar1)
        d2 = shape(ar2)

        if(bind_rows) then
            if(d1(2) /= d2(2)) then
                stop 'cannot bind rows of these arrays, since they have different numbers of columns'
            end if
        else
            if(d1(1) /= d2(1)) then
                stop 'cannot bind columns of these arrays, since they have different numbers of rows'
            end if
        end if

        call move_alloc(ar1, tmp_ar)
        
        if(bind_rows) then
            allocate(ar1(d1(1) + d2(1), d1(2)))
            ar1(1:d1(1),:) = tmp_ar
            ar1((d1(1) + 1):(d1(1)+d2(1)), :) = ar2
        else
            allocate(ar1(d1(1), d1(2)+d2(2)))
            ar1(:, 1:d1(2)) = tmp_ar
            ar1(:, (d1(2) + 1):(d1(2)+d2(2))) = ar2
        end if

    end subroutine 

    subroutine remove_rows_ip(ar1, indices, apply_to_rows)
        !!
        !! Remove rows from a rank-2 allocatable array ar1. Alternatively remove columns (if apply_to_rows == .FALSE.)
        !!
        integer(ip), allocatable, intent(inout) :: ar1(:,:) !! Rank 2 array
        integer(ip), intent(in) :: indices(:) !! Indices of rows (or columns) to remove
        logical, optional, intent(in) :: apply_to_rows !! If TRUE (default) remove rows, otherwise remove columns.

        logical :: rows       
        integer(ip), allocatable :: tmp_ar(:,:)
        integer(ip) :: d1(2), i, counter
        logical, allocatable :: tokeep(:)
        
        ! Apply to rows by default 
        if(present(apply_to_rows)) then
            rows = apply_to_rows 
        else
            rows = .true.
        end if

        ! Check indices are valid
        if(minval(indices) < 1) stop 'invalid indices to remove'

        d1 = shape(ar1)
        if(rows) then
            if(maxval(indices) > d1(1)) stop 'invalid row indices to remove'
            allocate(tokeep(d1(1)))
            tokeep = .true.
            do i = 1, d1(1)
                if(any(i == indices)) tokeep(i) = .false.
            end do
        else
            if(maxval(indices) > d1(2)) stop 'invalid column indices to remove'
            allocate(tokeep(d1(2)))
            tokeep = .true.
            do i = 1, d1(2)
                if(any(i == indices)) tokeep(i) = .false.
            end do
        end if

        ! Main action
        call move_alloc(ar1, tmp_ar)

        if(rows) then
            allocate(ar1(count(tokeep), d1(2)))
            ! Only copy rows with 'tokeep = true'
            counter = 1
            do i = 1, d1(1)
                if(tokeep(i)) then 
                    ar1(counter,:) = tmp_ar(i,:)
                    counter = counter + 1
                end if
            end do
        else
            allocate(ar1(d1(1), count(tokeep)))
            ! Only copy columns with 'tokeep = false'
            counter = 1
            do i = 1, d1(2)
                if(tokeep(i)) then
                    ar1(:, counter) = tmp_ar(:,i)
                    counter = counter + 1
                end if
            end do
        end if

    end subroutine
    
    subroutine test_which() 
        !! Test the module

        logical :: test_data(5)
        integer(ip), allocatable :: data_indices(:), values(:), lengths(:)
        ! Deliberately include an array with length=0 for testing
        integer(ip) :: test_data_rle(10), test_data_rle0(0), test_data_rle1(1), i

        integer(ip), allocatable :: ar1(:,:), ar2(:,:)

        !
        ! Test of 'which'
        !
        test_data = [.FALSE., .TRUE., .FALSE., .FALSE., .TRUE.]

        call which(test_data, data_indices)

        if( (data_indices(1) == 2) .AND. (data_indices(2) == 5) .AND. (size(data_indices) == 2)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if 

        !
        ! Test of 'rle'
        !

        ! Case with all unique values
        test_data_rle = [2,1,3,4,3,2,5,6,8,7]

        call rle_ip(test_data_rle, values, lengths)

        if(all(lengths == 1_ip .and. values == test_data_rle)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Case with array length = 0
        call rle_ip(test_data_rle0, values, lengths)
        if(size(lengths) == 0 .and. size(values) == 0) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Case with array length = 0
        test_data_rle1 = [3_ip]
        call rle_ip(test_data_rle1, values, lengths)
        if(size(lengths) == 1 .and. size(values) == 1 .and. lengths(1) == 1 .and. values(1) == 3) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Typical case
        test_data_rle = [2,2,3,3,3,2,5,6,8,8]
        call rle_ip(test_data_rle, values, lengths)
        if(size(lengths) == 6 .and. all(lengths == [2, 3, 1, 1, 1, 2]) .and. &
           size(values) == 6 .and. all(values == [2, 3, 2, 5, 6, 8])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if      
 
        ! Typical case
        test_data_rle = [2,2,3,3,3,2,5,6,7,8]
        call rle_ip(test_data_rle, values, lengths)
        if(size(lengths) == 7 .and. all(lengths == [2, 3, 1, 1, 1, 1, 1]) .and. &
           size(values) == 7 .and. all(values == [2, 3, 2, 5, 6, 7, 8])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if      

        ! Typical case
        test_data_rle = [2,3,3,3,3,3,3,3,3,3]
        call rle_ip(test_data_rle, values, lengths)
        if(size(lengths) == 2 .and. all(lengths == [1, 9]) .and. &
           size(values) == 2 .and. all(values == [2, 3])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if      

        ! Typical case
        test_data_rle = [2,3,3,3,3,3,3,3,3,10]
        call rle_ip(test_data_rle, values, lengths)
        if(size(lengths) == 3 .and. all(lengths == [1, 8, 1]) .and. &
           size(values) == 3 .and. all(values == [2, 3, 10])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if      

        ! Simple case with only 1 value
        test_data_rle = [2,2,2,2,2,2,2,2,2,2]
        call rle_ip(test_data_rle, values, lengths)
        if(size(lengths) == 1 .and. all(lengths == [10]) .and. &
           size(values) == 1 .and. all(values == [2])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if      

        ! Case with user provided equality function -- 'equal' means i/3 is equal
        test_data_rle = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        call rle_ip(test_data_rle, values, lengths, test_equality_function_integer_division_3)
        if(all(values/3 == [0, 1, 2, 3]) .and. all(lengths == [2, 3, 3, 2])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Test of cumsum
        !
        test_data_rle = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] 
        call cumsum_ip(test_data_rle)
        if(all(test_data_rle == [1, 3, 6, 10, 15, 21, 28, 36, 45, 55])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! test of bind_arrays, binding by column
        !
        ar1 = reshape( (/ (i, i=1, 6) /), [2,3])
        ar2 = reshape( (/ (i, i=7,10) /), [2,2] )

        call bind_arrays_ip(ar1, ar2, rowbind=.false.)

        if( all(shape(ar1) == [2,5]) .and. &
            all(ar1(1,:) == [1, 3, 5, 7, 9]) .and. &
            all(ar1(2,:) == [2, 4, 6, 8, 10]) ) then 
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if


        !
        ! test of bind_arrays, default binding (row)
        !
        ar1 = reshape( (/ (i, i=1, 6) /), [2,3])
        ar2 = reshape( (/ (i, i=7,9) /), [1,3] )
        call bind_arrays_ip(ar1, ar2)
        if(all(shape(ar1) == [3,3]) .and. &
           all(ar1(:,1) == [1,2,7]) .and. &
           all(ar1(:,2) == [3,4,8]) .and. &
           all(ar1(:,3) == [5,6,9]) ) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if
        deallocate(ar1, ar2)
    

        !
        ! Test remove_rows
        !
        ar1 = reshape( (/ (i, i=1, 9) /), [3,3])
        call remove_rows_ip(ar1, [2])
        if(all(ar1(1,:) == [1, 4, 7]) .and. &
            all(ar1(2,:) == [3, 6, 9]) .and. &
            all(shape(ar1) == [2,3]) ) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        ! Test remove_rows, with column removal
        !
        ar1 = reshape( (/ (i, i=1, 9) /), [3,3])
        call remove_rows_ip(ar1, [2], apply_to_rows=.false.)
        if( all(ar1(1,:) == [1, 7]) .and. &
            all(ar1(2,:) == [2, 8]) .and. &
            all(ar1(3,:) == [3, 9]) .and. &
            all(shape(ar1) == [3,2]) ) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Remove everything
        call remove_rows_ip(ar1, [1, 2], apply_to_rows=.false.)
        if(size(ar1) == 0 .and. all(shape(ar1) == [3,0])) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        contains
            ! Function used for a test of rle_ip
            function test_equality_function_integer_division_3(i1, i2) result(is_equal)
                integer(ip), intent(in) :: i1, i2
                logical :: is_equal
                is_equal = (i1/3 == i2/3)
            end function
    end subroutine


end module
