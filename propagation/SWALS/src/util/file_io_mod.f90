module file_io_mod
    !
    ! Various convenience routines for ascii file IO
    !

    use global_mod, only: dp, ip, charlen
    implicit none


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function count_file_lines(input_file_unit_no)
        ! Count the number of lines in input_file_unit_no
        integer(ip):: input_file_unit_no
        integer(ip):: count_file_lines

        integer(ip):: io_test

        rewind(input_file_unit_no) ! Start of file

        io_test = 0
        count_file_lines = 0

        do while (io_test >= 0)
            read(input_file_unit_no,*, iostat=io_test)
            if(io_test >= 0) count_file_lines = count_file_lines+1 
        end do

        rewind(input_file_unit_no) ! Start of file
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mkdir_p(output_folder, mkdir_status)
        ! Make a directory, recursively (i.e. like "mkdir -p output_folder")
        ! We try to work-around possible failures in parallel.
        ! It would be good to have a less hacky implementation however
        character(*), intent(in) :: output_folder
        integer(ip), optional, intent(out) :: mkdir_status

        character(len=charlen) :: mkdir_command, output_folder_name
        integer(ip) :: local_status, i

        output_folder_name = output_folder 

        mkdir_command = 'mkdir -p ' // trim(output_folder_name)

        ! When running in parallel, this can fail, but succeed if we try again
        ! See https://stackoverflow.com/questions/10627045/random-failure-of-mpi-fortran-code
        ! Therefore, we try up to 100 times in a loop
        do i = 1, 100
            !call system(trim(mkdir_command), status=local_status)
            !local_status = system(trim(mkdir_command))
            call execute_command_line(trim(mkdir_command), exitstat=local_status)
            if(local_status == 0) exit
        end do

        if(present(mkdir_status)) then
            mkdir_status = local_status
        end if

    end subroutine    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_character_file(input_file_unit_no, output_lines, format_string)
        ! Read the entire contents of input_file_unit_no into the allocatable
        ! character array 'output_lines'. Use format_string in the read
        integer(ip), intent(in):: input_file_unit_no
        character(*), intent(in):: format_string
        character(len=charlen), allocatable, intent(inout):: output_lines(:)

        integer(ip):: num_lines, i, io_test

        io_test = 0

        ! Clean out 'output_lines' in case it is already allocated
        if(allocated(output_lines)) deallocate(output_lines)

        ! Count how many lines we will need
        num_lines=count_file_lines(input_file_unit_no)
        allocate(output_lines(num_lines))

        print*, 'reading file ...' 

        do i=1,num_lines
            read(input_file_unit_no, format_string, iostat=io_test) output_lines(i) 
        end do
        rewind(input_file_unit_no)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_csv_into_array(array, csv_file, skip_header)
        ! Read the contents of csv_file into an array
        ! Allocation is automatically taken care of
        real(dp), allocatable, intent(inout):: array(:,:)
        character(len=charlen), intent(in):: csv_file
        integer(ip), optional, intent(in):: skip_header

        character(len=charlen):: local_buffer
        integer(ip) :: file_unit_no, file_columns, file_rows, i, skip_header_local

        if(present(skip_header)) then
            skip_header_local = skip_header
        else
            skip_header_local = 0
        endif
        
        ! Clean out 'array' in case it is already allocated
        if(allocated(array)) deallocate(array)

        open(newunit=file_unit_no, file=csv_file, action='read')

        file_rows = count_file_lines(file_unit_no)
        
        read(file_unit_no, '(A)') local_buffer
        rewind(file_unit_no)

        file_columns =  count(transfer(local_buffer, 'a', len(local_buffer)) == ",") + 1

        allocate(array(file_columns, file_rows - skip_header_local))
        
        do i = 1, file_rows
            if(i <= skip_header_local) then
                read(file_unit_no, *)
            else
                read(file_unit_no,*) array(:,i - skip_header_local)          
            end if
        end do

        close(file_unit_no)

    end subroutine

end module
