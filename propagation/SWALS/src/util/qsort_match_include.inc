        ! Generic code for 'match'
        ! This works irrespective of the kind of 'arary1', but due to limitations of Fortran
        ! we need to make a separate one for each kind. 
        !
        if(size(array1) /= size(matches)) then
            print*, 'Size of array1 (', size(array1), ') is not equal to size of matches (', size(matches), ')'
            stop
        end if

        allocate(inds1(size(array1)), inds2(size(array2)))

        ! Get indices giving an 'order' for array1 and array2
        n = size(inds1)
        call sort_index(inds1, array1, n)
        n = size(inds2)
        call sort_index(inds2, array2, n)

        matches = -1

        i2 = 1
        i = 1
        do while (i <= size(array1) .and. i2 <= size(array2))

            if(array1(inds1(i)) > array2(inds2(i2))) then
                ! Move forward in array2
                i2 = i2 + 1
                cycle
            end if

            if(array1(inds1(i)) == array2(inds2(i2))) then
                ! Record the match, and move onto the next entry of (sorted) array1
                matches(inds1(i)) = inds2(i2)
                i = i + 1
                cycle
            end if

            if(array1(inds1(i)) < array2(inds2(i2))) then
                ! There can be no match for the 'ith' sorted entry of array1
                i = i + 1
                cycle 
            end if
        end do

