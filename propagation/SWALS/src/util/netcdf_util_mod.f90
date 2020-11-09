module netcdf_util
    !!
    !! Class to write domain variables to a netcdf file.
    !!

    use global_mod, only: charlen, ip, dp
    use stop_mod, only: generic_stop
    use iso_c_binding, only: C_DOUBLE
    use iso_fortran_env

#ifndef NONETCDF
    use netcdf
#endif

! Switch on/off netcdf 4 calls, which will fail if we have only compiled with netcdf 3
#ifdef NETCDF4_GRIDS
# define USING_NETCDF4_GUARD( x ) x    
#else
# define USING_NETCDF4_GUARD( x ) 
#endif

    implicit none

    private 
    public :: nc_grid_output_type, check

    ! Indices of stage, uh, vh, elevation in domain%U
    integer(ip), parameter :: STG=1, UH=2, VH=3, ELV=4
    ! Note -- the above definitions must agree with those in domain_mod.f90

    type :: nc_grid_output_type
        !!
        !! Type for saving domain data to netcdf. Holds unit numbers, etc
        !!

        character(len=charlen) :: filename
        integer :: nc_file_id
            !! Netcdf file id should be 'plain' integer

        integer :: dim_time_id, dim_x_id, dim_y_id
        integer :: var_time_id, var_x_id, var_y_id

        integer :: var_stage_id, var_uh_id, var_vh_id, var_elev_id
            !! Time varying variables
#ifdef DEBUG_ARRAY
        integer :: var_debug_array_id
            !! Time-varying array for debugging
#endif
        integer :: var_max_stage_id, var_elev0_id, var_is_priority_domain_id, &
            var_manningsq_id, var_priority_ind, var_priority_img
            !! Static variables

        integer :: flush_every_n_output_steps = 25
            !! It can be slow to flush every time we write to files. But
            !! this can be useful for interactively looking at model results. So
            !! use this variable to control the frequency of file-flushing.

        integer :: num_output_steps = 0
            ! Track number of output steps

        logical :: record_max_U
            !! Flag for whether we store max stage, etc

        logical :: time_var_store_flag(4) = [.true. , .true., .true., .true.]
            !! Options to store stage/uh/vh/elev over time.
            !! .true. or .false. for STG, UH, VH, ELV
            !! For example, if we didn't want to store ELV every time-step, we
            !! would do "nc_grid_output%time_var_store_flag(ELV) = .false."

        integer :: spatial_stride = 1
            !! Allow only storing every n'th point in x/y space.
            !! For instance, on a 1-arc-minute model, a value of 4 would
            !! store the results at 4 arc-minutes
        integer :: spatial_start(2)
        integer :: spatial_count(2)

        real(dp) :: spatial_ll_full_domain(2) = [-HUGE(1.0_dp), -HUGE(1.0_dp)]
            !! When we partition domains, this can store the lower-left of the FULL domain
            !! We can then use that to make spatial_start consistent among all domains

        contains
        procedure :: initialise => initialise
        procedure :: write_grids => write_grids
        procedure :: store_max_stage => store_max_stage
        procedure :: store_priority_domain_cells => store_priority_domain_cells
        procedure :: finalise => finalise

    end type


    contains

    !
    ! Convenience subroutine for netcdf error handling
    !
    subroutine check(error_status, line)
        integer, intent(in) :: error_status
        integer, optional, intent(in) :: line

    ! Gracefully compile in case we don't link with netcdf 
#ifndef NONETCDF

        if(error_status /= nf90_noerr) then
          print*, 'error_status: ', error_status
          print *, trim(nf90_strerror(error_status))
          if(present(line)) then
              print*, 'line ', line
          end if
          call generic_stop()
        end if

#endif

    end subroutine

    subroutine initialise(nc_grid_output, filename, output_precision, &
        record_max_U, xs, ys, attribute_names, attribute_values)
        !!
        !! Create nc file to hold gridded outputs
        !!

        class(nc_grid_output_type), intent(inout) :: nc_grid_output !!type with nc_grid_output info
        character(len=charlen), intent(in) :: filename !! Name of the output netcdf file
        integer(ip), intent(in) :: output_precision !! Fortran real kind denoting precision of output (either C_FLOAT or C_DOUBLE)
        logical, intent(in) :: record_max_U !! Do we record the max-stage ?
        real(dp), intent(in) :: xs(:), ys(:) !! array of grid x/y coordinates
        character(charlen), optional, intent(in):: attribute_names(:) !! global attribute names for netcdf file
        character(charlen), optional, intent(in):: attribute_values(:) 
            !! values for global attributes of netcdf file (same length as attribute_names)

        integer:: iNcid, output_prec, i, output_byte, output_int4, spatial_stride, spatial_start(2)
        integer:: nx, ny, first_index_relative_to_full_domain(2), output_prec_force_double
        real(dp) :: dx_local(2)
        character(len=charlen) :: local_att

        integer:: netcdf_file_type
        !logical:: using_netcdf4
       
    ! Gracefully compile in case we don't link with netcdf
#ifndef NONETCDF

        ! Setup the filetype -- currently the HDF5 format seems slower, and compression is not great
#ifdef NETCDF4_GRIDS
        netcdf_file_type = NF90_HDF5 
#else
        netcdf_file_type = NF90_CLOBBER
#endif

        nx = size(xs)
        ny = size(ys)

        ! Define the stride / count / start
        nc_grid_output%spatial_start = nc_grid_output%spatial_stride/2 + 1

        if(any(nc_grid_output%spatial_ll_full_domain /= -HUGE(1.0_dp))) then
            ! For partitioned domains we might need to adjust the spatial start for
            ! consistency with other parallel-parts of this domain -- so that we are still
            ! storing a regular sequence of points relative to the full domain

            dx_local = [(xs(nx) - xs(1))/(nx-1.0_dp), (ys(ny) - ys(1))/(ny-1.0_dp)]

            first_index_relative_to_full_domain = &
                nint( ( [xs(1), ys(1)] - nc_grid_output%spatial_ll_full_domain )/dx_local + 0.5_dp, kind=ip)

            nc_grid_output%spatial_start = modulo(&
                nc_grid_output%spatial_start - first_index_relative_to_full_domain,&
                nc_grid_output%spatial_stride) + 1_ip
        end if

        spatial_stride = nc_grid_output%spatial_stride
        spatial_start = nc_grid_output%spatial_start
        nc_grid_output%spatial_count = [ &
            size(xs(spatial_start(1):nx:spatial_stride)), &
            size(ys(spatial_start(2):ny:spatial_stride))]

        ! Usually we output single precision reals, irrespective of the model run
        ! precision [can be changed in global_mod.f90]
        ! However, for x/y/time, lets store them as double always
        if(output_precision == C_DOUBLE) then
            output_prec = NF90_REAL8
        else
            output_prec = NF90_REAL4
        end if

        output_prec_force_double = NF90_REAL8

        output_byte = NF90_BYTE
        output_int4 = NF90_INT4

        !
        ! Create output file
        !
        call check(nf90_create(filename, netcdf_file_type, nc_grid_output%nc_file_id), __LINE__)

        iNcid = nc_grid_output%nc_file_id ! Shorthand

        ! Define dimensions
        call check(nf90_def_dim(iNcid, "x", len=nc_grid_output%spatial_count(1), dimid=nc_grid_output%dim_x_id), __LINE__) 
        call check(nf90_def_dim(iNcid, "y", len=nc_grid_output%spatial_count(2), dimid=nc_grid_output%dim_y_id), __LINE__)
        call check(nf90_def_dim(iNcid, "time", NF90_UNLIMITED, nc_grid_output%dim_time_id), __LINE__)

        ! Define variables
        call check( nf90_def_var(iNcid, "x", output_prec_force_double, &
            (/ nc_grid_output%dim_x_id /), nc_grid_output%var_x_id), __LINE__)
        call check( nf90_def_var(iNcid, "y", output_prec_force_double, &
            (/ nc_grid_output%dim_y_id /), nc_grid_output%var_y_id), __LINE__ )
        call check( nf90_def_var(iNcid, "time", output_prec_force_double, &
            (/ nc_grid_output%dim_time_id /), nc_grid_output%var_time_id), __LINE__ )

        ! Stage
        if(nc_grid_output%time_var_store_flag(STG)) then
            call check( nf90_def_var(iNcid, "stage", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], &
                nc_grid_output%var_stage_id), __LINE__ )
            ! Compress
            USING_NETCDF4_GUARD( call check(nf90_def_var_deflate(iNcid, nc_grid_output%var_stage_id, 1, 1, 3), __LINE__) )
        end if
        ! UH
        if(nc_grid_output%time_var_store_flag(UH)) then
            call check( nf90_def_var(iNcid, "uh", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], &
                nc_grid_output%var_uh_id), __LINE__ )
            ! Compress
            USING_NETCDF4_GUARD( call check(nf90_def_var_deflate(iNcid, nc_grid_output%var_uh_id, 1, 1, 3), __LINE__) )
        end if
        ! VH
        if(nc_grid_output%time_var_store_flag(VH)) then
            call check( nf90_def_var(iNcid, "vh", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], &
                nc_grid_output%var_vh_id), __LINE__ )
            ! Compress
            USING_NETCDF4_GUARD( call check(nf90_def_var_deflate(iNcid, nc_grid_output%var_vh_id, 1, 1, 3), __LINE__) )
        end if
        ! Elev
        if(nc_grid_output%time_var_store_flag(ELV)) then
            call check( nf90_def_var(iNcid, "elev", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], &
                nc_grid_output%var_elev_id), __LINE__ )
            ! Compress
            USING_NETCDF4_GUARD( call check(nf90_def_var_deflate(iNcid, nc_grid_output%var_elev_id, 1, 1, 3), __LINE__) )
        end if

#ifdef DEBUG_ARRAY
        ! A time-varying rank-3 debug array, with spatial dimensions (nx, ny)
        call check(nf90_def_var(iNcid, 'debug_array', output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], &
                nc_grid_output%var_debug_array_id), __LINE__ )

#endif

        ! priority_domain -- Use 1 if the current domain is the priority domain, and zero otherwise
        call check(nf90_def_var(iNcid, "is_priority_domain", output_byte, &
            [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
            nc_grid_output%var_is_priority_domain_id), __LINE__ )
        
        ! priority_domain index
        call check(nf90_def_var(iNcid, "priority_domain_index", output_int4, &
            [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
            nc_grid_output%var_priority_ind), __LINE__ )
        ! priority domain image
        call check(nf90_def_var(iNcid, "priority_domain_image", output_int4, &
            [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
            nc_grid_output%var_priority_img), __LINE__ )

        ! Peak stage
        if(record_max_U) then

            nc_grid_output%record_max_U = .TRUE.

            call check( nf90_def_var(iNcid, "max_stage", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
                nc_grid_output%var_max_stage_id), __LINE__ )

            call check( nf90_def_var(iNcid, "elevation0", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
                nc_grid_output%var_elev0_id), __LINE__ )

            call check( nf90_def_var(iNcid, "manning_squared", output_prec, &
                [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
                nc_grid_output%var_manningsq_id), __LINE__ )
        else

            nc_grid_output%record_max_U = .FALSE.

        end if

        !
        ! Add some attributes
        !
        if(present(attribute_names) .and. present(attribute_values)) then
            do i = 1, size(attribute_names)
                call check(nf90_put_att(iNcid, nf90_global, attribute_names(i), attribute_values(i)), __LINE__)
            end do
        end if
#ifdef SRC_GIT_VERSION
        ! Add the git revision number to the file
        call check(nf90_put_att(iNcid, nf90_global, 'git_revision_number',& ! Continuation to reduce chance of > 132 char
SRC_GIT_VERSION ), &
        __LINE__)
#endif
        ! Add the command-line call
        call get_command(local_att)
        call check(nf90_put_att(iNcid, nf90_global, 'run_command', TRIM(local_att)), &
        __LINE__)
#ifndef PGI_COMPILER
        ! Add the compiler version
        call check(nf90_put_att(iNcid, nf90_global, 'compiler_version', TRIM(compiler_version())), &
        __LINE__)
        ! Add the compiler options
        call check(nf90_put_att(iNcid, nf90_global, 'compiler_options', TRIM(compiler_options())), &
        __LINE__)
#endif

        ! Finish definitions so writing can begin
        call check(nf90_enddef(iNcid), __LINE__)

        !
        ! Write a few things
        !
        call check(nf90_put_var(iNcid, nc_grid_output%var_x_id , &
            xs(spatial_start(1):nx:spatial_stride)), __LINE__ ) 
        call check(nf90_put_var(iNcid, nc_grid_output%var_y_id , &
            ys(spatial_start(2):ny:spatial_stride)), __LINE__)

#endif

    end subroutine

    subroutine write_grids(nc_grid_output, time, U, debug_array)
        !!
        !! Write variables in U to file associated with nc_grid_output
        !!

        class(nc_grid_output_type), intent(inout) :: nc_grid_output 
            !! instance of nc_grid_output_type which has been initialised
        real(dp), intent(in) :: U(:,:,:) 
            !! array with the flow data. It has size(U,3) == 4. 
            !! Also U(:,:,1) = stage, U(:,:,2) = uh, U(:,:,3) = vh, U(:,:,4) = elev
        real(dp), intent(in) :: time 
            !! Time associated with U(:,:,:)
        real(dp), optional, intent(in) :: debug_array(:,:)
            !! If compiled with -DDEBUG_ARRAY, this array will also be written to the netcdf file.
            !! Provides a useful means of debugging or reporting non-standard quantities.

        integer:: iNcid, nxy(2), spatial_start(2), spatial_stride(2)
        real(C_DOUBLE) :: current_cpu_time

        ! Gracefully compile in case we don't link with netcdf
#ifndef NONETCDF

        nc_grid_output%num_output_steps = nc_grid_output%num_output_steps + 1

        if(size(U,3) /= 4) then
            print*, 'Incorrect size of 3 rank of U when writing grids to netcdf'
            call generic_stop()
        end if

        ! Shorthand
        iNcid = nc_grid_output%nc_file_id
        spatial_start = nc_grid_output%spatial_start
        nxy(1) = size(U, dim=1)
        nxy(2) = size(U, dim=2)
        spatial_stride = nc_grid_output%spatial_stride

        ! Save the time
        call check(nf90_put_var(iNcid, nc_grid_output%var_time_id, time, &
            start=[nc_grid_output%num_output_steps]), __LINE__ )

        ! Save the stage
        if(nc_grid_output%time_var_store_flag(STG)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_stage_id, &
                U(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2),STG), &
                start=[1,1,nc_grid_output%num_output_steps]),&
                __LINE__)
        end if

        ! Save uh
        if(nc_grid_output%time_var_store_flag(UH)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_uh_id, &
                U(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2),UH), &
                start=[1,1,nc_grid_output%num_output_steps]), &
                __LINE__)
        end if

        ! Save vh
        if(nc_grid_output%time_var_store_flag(VH)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_vh_id, &
                U(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2),VH), &
                start=[1,1,nc_grid_output%num_output_steps]), &
                __LINE__)
        end if

        ! Save elev
        if(nc_grid_output%time_var_store_flag(ELV)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_elev_id, &
                U(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2),ELV), &
                start=[1,1,nc_grid_output%num_output_steps]), &
                __LINE__)
        end if

#ifdef DEBUG_ARRAY
        if(present(debug_array)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_debug_array_id, &
                debug_array(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2)), &
                start=[1,1,nc_grid_output%num_output_steps]), &
                __LINE__)
        end if
#endif

        ! It can be slow to flush everytime we write to the file (but sometimes useful for interactive work)
        if(mod(nc_grid_output%num_output_steps, nc_grid_output%flush_every_n_output_steps) == 0) call check(nf90_sync(iNcid))

#endif

    end subroutine


    subroutine store_priority_domain_cells(nc_grid_output, priority_domain_index, priority_domain_image, &
        is_priority_domain_not_periodic, my_index, my_image)
        !! Store a 1 byte integer, denoting the cells in the current domain that are 
        !! priority_domain cells. Also store regular integer grids with the priority domain index
        !! and image.
        class(nc_grid_output_type), intent(in) :: nc_grid_output 
            !! Initialised nc_grid_output type 
        integer(ip), intent(in) :: priority_domain_index(:,:), priority_domain_image(:,:)
            !! Arrays denoting the priority domain index and image, i.e., the index/image of
            !! the domain that is consider to contain the 'real' solution on each cell. Generally the priority domain
            !! corresponds to the finest resolution domain covering that cell. Each domain will
            !! hold arrays like this, which cover its own grid points only, and are used to determine
            !! whether a given part of the flow solution should be considered the 'true' solution, or merely
            !! a halo region. It is important to know this, e.g. for mass conservation calculations, for
            !! the creation of the nesting communication data structures, and for making seamless outputs using
            !! only priority domain values.
        integer(ip), intent(in) :: is_priority_domain_not_periodic(:,:)
            !! Integer array with value 1_ip where the cell is in the priority domain and not in periodic regions,
            !! and value 0_ip otherwise
        integer(ip) :: my_index, my_image
            !! The domain index (i.e. index of the domain in md%domains(:)) and image (i.e. always 1 for non-coarray programs, or
            !! equal to this_image() for coarray programs) 

        integer(ip) :: i, iNcid, spatial_start(2), nxy(2), spatial_stride(2)
    
#ifndef NONETCDF 
        iNcid = nc_grid_output%nc_file_id ! Shorthand

        spatial_stride = nc_grid_output%spatial_stride
        spatial_start = nc_grid_output%spatial_start
        nxy(1) = size(priority_domain_index, dim=1)
        nxy(2) = size(priority_domain_index, dim=2)

        ! Loop to avoid making a temporary variable that contains the 0/1 mask
        do i = spatial_start(2), nxy(2), spatial_stride(2)
            ! NOTE: This records cells that are in the priority domain, but it may double-count
            ! cells if there are periodic regions, and the domain receives periodic halos from itself. 
            ! There is another variable (integer(ip) domain%nesting%is_priority_domain_not_periodic) that records
            ! this info - consider making use of that here instead.
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_is_priority_domain_id, &
                merge(1, 0, &
                    mask=(is_priority_domain_not_periodic(spatial_start(1):nxy(1):spatial_stride(1),i) == 1_ip)), &
                start = [1, (i - spatial_start(2))/spatial_stride(2) + 1]), &
                __LINE__)
        end do 

        ! Priority domain index
        do i = spatial_start(2), nxy(2), spatial_stride(2)

            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_priority_ind, &
                priority_domain_index(spatial_start(1):nxy(1):spatial_stride(1),i), &
                start = [1, (i - spatial_start(2))/spatial_stride(2) + 1]), &
                __LINE__)
        end do 

        ! Priority domain image
        do i = spatial_start(2), nxy(2), spatial_stride(2)

            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_priority_img, &
                priority_domain_image(spatial_start(1):nxy(1):spatial_stride(1),i), &
                start = [1, (i - spatial_start(2))/spatial_stride(2) + 1]), &
                __LINE__)
        end do 

#endif
    end subroutine


    subroutine store_max_stage(nc_grid_output, max_stage, elevation, manning_squared)
        !!
        !! Store the 'max-stage' variable. Also store elevation and manning_squared, if provided.
        !!
        class(nc_grid_output_type), intent(in) :: nc_grid_output
        real(dp), intent(in) :: max_stage(:,:)
        real(dp), optional, intent(in) :: elevation(:,:)
        real(dp), optional, intent(in) :: manning_squared(:,:)

        integer :: iNcid, spatial_start(2), nxy(2), spatial_stride(2)

        ! Gracefully compile in case we don't link with netcdf
#ifndef NONETCDF
        spatial_start = nc_grid_output%spatial_start
        spatial_stride = nc_grid_output%spatial_stride
        nxy(1) = size(max_stage, dim=1)
        nxy(2) = size(max_stage, dim=2)
 
        if(nc_grid_output%record_max_U) then  
 
            ! Shorthand
            iNcid = nc_grid_output%nc_file_id

            ! Save max stage
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_max_stage_id, &
                max_stage(spatial_start(1):nxy(1):spatial_stride(1), spatial_start(2):nxy(2):spatial_stride(2))), &
                __LINE__)

            if(present(elevation)) then
                ! Save elevation
                call check(nf90_put_var(&
                    iNcid, &
                    nc_grid_output%var_elev0_id, &
                    elevation(spatial_start(1):nxy(1):spatial_stride(1), spatial_start(2):nxy(2):spatial_stride(2))), &
                    __LINE__)
            end if

            if(present(manning_squared)) then
                ! Save manning coefficient squared
                call check(nf90_put_var(&
                    iNcid, &
                    nc_grid_output%var_manningsq_id, &
                    manning_squared(spatial_start(1):nxy(1):spatial_stride(1), spatial_start(2):nxy(2):spatial_stride(2))), &
                    __LINE__)
            end if

        end if

#endif

    end subroutine

    subroutine finalise(nc_grid_output)
        !! Close the netcdf file. 

        class(nc_grid_output_type), intent(inout) :: nc_grid_output

    ! Compile gracefully in case we don't use netcdf
#ifndef NONETCDF

        ! Close the netcdf file 
        call check(nf90_close(nc_grid_output%nc_file_id))
        
#endif

    end subroutine

end module
