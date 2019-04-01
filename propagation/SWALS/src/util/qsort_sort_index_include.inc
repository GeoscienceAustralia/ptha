        ! Generic code for sort_index. This works irrespective of the type
        ! of array, but due to limitations of Fortran, we need to make a 
        ! different subroutine for each kind
        inds = (/ (i, i=1, n) /)

        elem_count = size(array)
        elem_size = storage_size(inds(1))/8
      
        call qsort(c_loc(inds), elem_count, elem_size, c_funloc(compar))
    
        contains
            pure function compar(i1, i2) result(compare) bind(C)
                integer(c_int), intent(in) :: i1, i2
                integer(c_int) :: compare
        
                if(array(i1) > array(i2)) compare = 1
                if(array(i1) == array(i2)) compare = 0
                if(array(i1) < array(i2)) compare = -1

            end function

