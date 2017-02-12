MODULE which_mod
    !
    ! Simple module to compute the indices which are 'TRUE' in a logical array
    !
    USE global_mod, only: ip

    IMPLICIT NONE
    
    PRIVATE
    PUBLIC which, test_which


    contains

    ! Return an array with integers corresponding to indices where logical_array is .true.
    ! Like 'which' in R.
    SUBROUTINE which(logical_array, allocatable_array)
        LOGICAL, INTENT(IN) :: logical_array(:)
        INTEGER(ip), ALLOCATABLE, INTENT(OUT):: allocatable_array(:)
        INTEGER(ip):: arr_len, i, counter

        IF(allocated(allocatable_array)) THEN
            deallocate(allocatable_array)
        END IF

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

    END SUBROUTINE

    ! Test the 'which' subroutine
    SUBROUTINE test_which() 
        LOGICAL:: does_which_work

        LOGICAL :: test_data(5)
        INTEGER(ip), ALLOCATABLE :: data_indices(:)

        test_data = [.FALSE., .TRUE., .FALSE., .FALSE., .TRUE.]

        call which(test_data, data_indices)

        if( (data_indices(1) == 2) .AND. (data_indices(2) == 5) .AND. (size(data_indices) == 2)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if 

    END SUBROUTINE

END MODULE
