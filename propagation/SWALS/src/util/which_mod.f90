module which_mod
    !
    ! Simple module to compute the indices which are 'TRUE' in a logical array
    !
    use global_mod, only: ip

    implicit none
    
    private
    public which, test_which


    contains

    ! Return an array with integers corresponding to indices where logical_array is .true.
    ! Like 'which' in R.
    subroutine which(logical_array, allocatable_array)
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
                    counter = counter+1
                    allocatable_array(counter) = i
                end if
            end do
        end if

    end subroutine

    ! Test the 'which' subroutine
    subroutine test_which() 

        logical :: test_data(5)
        integer(ip), allocatable :: data_indices(:)

        test_data = [.FALSE., .TRUE., .FALSE., .FALSE., .TRUE.]

        call which(test_data, data_indices)

        if( (data_indices(1) == 2) .AND. (data_indices(2) == 5) .AND. (size(data_indices) == 2)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if 

    end subroutine

end module
