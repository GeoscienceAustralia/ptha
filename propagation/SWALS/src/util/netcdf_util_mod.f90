module netcdf_util
    !!
    !! Class to write domain variables to a netcdf file.
    !!

    use global_mod, only: charlen, ip, dp
    use stop_mod, only: generic_stop
    use logging_mod, only: log_output_unit
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
    public :: nc_grid_output_type, check, gridded_output_variables_time, &
        gridded_output_variables_max_U, gridded_output_variables_other

    ! Indices of stage, uh, vh, elevation in domain%U
    integer(ip), parameter :: STG=1, UH=2, VH=3, ELV=4
    ! Note -- the above definitions must agree with those in domain_mod.f90

    !
    ! Names of gridded variables that SWALS can output
    !
    character(len=charlen), parameter :: gridded_output_variables_time(4) = &
        [character(len=charlen) :: 'stage', 'uh', 'vh', 'elev']
        !! Gridded variables that can be output as time-series. Must match [STG, UH, VH, ELV] in order.
    character(len=charlen), parameter :: gridded_output_variables_max_U(6) = &
        [character(len=charlen) :: 'max_stage', 'max_speed', 'max_flux', 'arrival_time', 'time_of_max_stage', 'min_stage']
        !! Gridded variables that can be output, are not time-series, but require tracking flow maxima.
        !! Keep 'max_stage' as the first entry (for now) to support code that assumes domain%max_U(:,:,1) contains the max-stage.
    character(len=charlen), parameter :: gridded_output_variables_other(3) = &
        [character(len=charlen) :: 'manning_squared', 'elevation0', 'elevation_source_file_index']
        !! Other gridded variables that can be output and are not time-series

    integer, parameter :: max_temporal_var = size(gridded_output_variables_time)
        !! Maximum number of temporal variables we could possibly store (in practice we might store less)
    integer, parameter :: max_static_var = size(gridded_output_variables_max_U) + &
                                           size(gridded_output_variables_other)
        !! Maximum number of static variables we could possibly store (in practice we might store less)

    type :: nc_grid_output_type
        !!
        !! Type for saving domain data to netcdf. Holds unit numbers, etc
        !!

        character(len=charlen) :: filename
        integer :: nc_file_id
            !! Netcdf file id should be 'plain' integer

        integer :: dim_time_id, dim_x_id, dim_y_id
        integer :: var_time_id, var_x_id, var_y_id

#ifdef DEBUG_ARRAY
        integer :: var_debug_array_id
            !! Time-varying array for debugging
#endif
        integer :: var_is_priority_domain_id, var_priority_ind, var_priority_img
            !! Variables describing the priority domain

        integer :: flush_every_n_output_steps = 25
            !! It can be slow to flush every time we write to files. But
            !! this can be useful for interactively looking at model results. So
            !! use this variable to control the frequency of file-flushing.

        integer :: num_output_steps = 0
            ! Track number of output steps

        !
        ! Time-varying variables to store.
        !
        integer :: time_var_id(max_temporal_var)
            !! netcdf variable id's for stage, uh, vh, elev
        character(len=charlen) :: time_var_names(max_temporal_var) = gridded_output_variables_time
            !! Names of the time-varying variables - do not change this. We turn on/off storage by
            !! adjusting the corresponding entry of time_var_store_flag
        logical :: time_var_store_flag(max_temporal_var) = .true. 
            !! Options to store stage/uh/vh/elev over time.
            !! .true. or .false. for STG, UH, VH, ELV
            !! For example, if we didn't want to store ELV every time-step, we
            !! would do "nc_grid_output%time_var_store_flag(ELV) = .false."

        !
        ! Static variables that can be stored. These describe the solution, are real, and are only
        ! stored once (i.e. no time evolution)
        !
        integer :: static_var_id(max_static_var)
        character(len=charlen) :: static_var_names(max_static_var) = &
            [gridded_output_variables_max_U, gridded_output_variables_other]
            !! Names of solution variables that can be stored once. We turn on/off storeage
            !! by adjusting the static_var_store_flag
        logical :: static_var_store_flag(max_static_var) = .true.
            !! True if we want to store the corresponding variable in static_var_names, false otherwise.

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
        procedure, non_overridable :: initialise => initialise
        procedure, non_overridable :: write_grids => write_grids
        procedure, non_overridable :: store_static_variable => store_static_variable
        procedure, non_overridable :: store_priority_domain_cells => store_priority_domain_cells
        procedure, non_overridable :: finalise => finalise

    end type


    contains

    !
    ! Convenience subroutine for netcdf error handling
    !
    subroutine check(error_status, line, stop_on_error)
        integer, intent(in) :: error_status
            !! Integer returned by a netcdf function call indicating the error.
        integer, optional, intent(in) :: line
            !! Line of the file (typically provided using the __LINE__ macro)
        logical, optional, intent(in) :: stop_on_error
            !! If .true. (default) then if(error_status != nf90_noerr) the program will halt.
            !! Useful to set to .false. when writing files, particularly if we store double precision
            !! variables with single precision, since if the double variables are outside
            !! the range of single precision then netcdf will provide an error message
            !! but everything will otherwise work.

        logical :: stop_on_err

        stop_on_err = .true.
        if(present(stop_on_error)) stop_on_err = stop_on_error

    ! Gracefully compile in case we don't link with netcdf
#ifndef NONETCDF

        if(error_status /= nf90_noerr) then
          write(log_output_unit,*) 'netcdf error_status: ', error_status
          write(log_output_unit,*) trim(nf90_strerror(error_status))
          if(present(line)) then
              write(log_output_unit,*) 'line ', line
          end if
          if(stop_on_err) call generic_stop()
        end if

#endif

    end subroutine

    subroutine initialise(nc_grid_output, filename, output_precision, &
        xs, ys, &
        time_grids_to_store, nontemporal_grids_to_store, &
        attribute_names, attribute_values)
        !!
        !! Create nc file to hold gridded outputs
        !!

        class(nc_grid_output_type), intent(inout) :: nc_grid_output 
            !!type with nc_grid_output info
        character(len=charlen), intent(in) :: filename 
            !! Name of the output netcdf file
        integer(ip), intent(in) :: output_precision 
            !! Fortran real kind denoting precision of output (either C_FLOAT or C_DOUBLE)
        real(dp), intent(in) :: xs(:), ys(:) 
            !! array of grid x/y coordinates
        character(charlen), allocatable, intent(inout) :: time_grids_to_store(:) 
            !! Array with names of time-variables to store. Values should match those in nc_grid_output%time_var_names, or ''.
            !! To support some legacy behaviour we allow this to be unallocated on entry. Some older code specified
            !! variables to store by directly manipulating the nc_grid_output variable, and that approach is still supported
            !! if time_grids_to_store is unallocated (although it is not recommended -- in future code please just specify
            !! the variables to store in domain%time_grids_to_store).
        character(charlen), intent(inout) :: nontemporal_grids_to_store(:) 
            !! Array with names of variables without a time dimension. Note that irrespective of the value of this variable, 
            !! we also store other "priority-domain" variables. Values should match those in nc_grid_output%static_var_names, or ''.
        character(charlen), optional, intent(in):: attribute_names(:) 
            !! global attribute names for netcdf file
        character(charlen), optional, intent(in):: attribute_values(:)
            !! values for global attributes of netcdf file (same length as attribute_names)

        integer:: iNcid, output_prec, i, output_byte, output_int4, output_int2, spatial_stride, spatial_start(2)
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

        if(spatial_stride > minval([nx, ny])) then
            write(log_output_unit,*) ''
            write(log_output_unit,*) 'Error in nc_grid_output initialisation for filename: '
            write(log_output_unit,*) '    ', trim(filename)
            write(log_output_unit,*) '  Cannot have spatial_stride > minval(domain%nx)'
            write(log_output_unit,*) '    spatial_stride: ', spatial_stride, &
                '; minval(domain%nx): ', minval([nx, ny])
            flush(log_output_unit)
            call generic_stop()
        end if

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
        output_int2 = NF90_INT2 ! Can get away with this for priority_domain_index/priority_domain_image

        !
        ! Create output file
        !
        call check(nf90_create(filename, netcdf_file_type, nc_grid_output%nc_file_id), __LINE__)

        iNcid = nc_grid_output%nc_file_id ! Shorthand

        ! Define dimensions
        call check(nf90_def_dim(iNcid, "x", len=nc_grid_output%spatial_count(1), dimid=nc_grid_output%dim_x_id), __LINE__)
        call check(nf90_def_dim(iNcid, "y", len=nc_grid_output%spatial_count(2), dimid=nc_grid_output%dim_y_id), __LINE__)
        call check(nf90_def_dim(iNcid, "time", NF90_UNLIMITED, nc_grid_output%dim_time_id), __LINE__)

        ! Define variables of dimensions
        call check( nf90_def_var(iNcid, "x", output_prec_force_double, &
            (/ nc_grid_output%dim_x_id /), nc_grid_output%var_x_id), __LINE__)
        call check( nf90_def_var(iNcid, "y", output_prec_force_double, &
            (/ nc_grid_output%dim_y_id /), nc_grid_output%var_y_id), __LINE__ )
        call check( nf90_def_var(iNcid, "time", output_prec_force_double, &
            (/ nc_grid_output%dim_time_id /), nc_grid_output%var_time_id), __LINE__ )

        if(allocated(time_grids_to_store)) then
            ! If time_grids_to_store has been specified, use it to define which variables to store
            ! This will over-ride any specification of nc_grid_output%time_var_store_flag (which some
            ! legacy programs do). So make sure the latter has not been set already
            if(.not. all(nc_grid_output%time_var_store_flag)) then
                write(log_output_unit, *)&
                    "Cannot specify BOTH time_grids_to_store AND time_var_store_flag. Use one or the other."
                call generic_stop
            end if

            ! Check for typos in time_grids_to_store. It can contain either an empty string '', or
            ! a value in nc_grid_output%time_var_names
            do i = 1, size(time_grids_to_store)
                if(.not. (time_grids_to_store(i) == '' .or. &
                      any(time_grids_to_store(i) == nc_grid_output%time_var_names))) then
                    write(log_output_unit, *)&
                        "Unknown value in time_grids_to_store:", trim(time_grids_to_store(i))
                    call generic_stop
                end if
            end do

            ! Set the time_var_store_flag based on the presence/absence of the variable in time_grids_to_store
            do i = 1, size(nc_grid_output%time_var_names)
                nc_grid_output%time_var_store_flag(i) = &
                    any(nc_grid_output%time_var_names(i) == time_grids_to_store)
            end do
        else
            ! If unallocated, make time_grids_to_store consistent with time_var_store_flag
            ! This gives backward compatibility with older application codes that directly 
            ! set nc_grid_output%time_var_store_flag, which was previously the only way to 
            ! control these things.
            time_grids_to_store = pack(nc_grid_output%time_var_names, nc_grid_output%time_var_store_flag)
        end if

        ! Store the time-varying variables
        do i = STG, ELV

            if( .not. nc_grid_output%time_var_store_flag(i)) cycle

            ! Store the variable
            call check( &
                nf90_def_var(iNcid, &
                    trim(nc_grid_output%time_var_names(i)), & ! Name of the variable
                    output_prec, & ! Precision
                    [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id, nc_grid_output%dim_time_id], & ! Dimensions
                    nc_grid_output%time_var_id(i) & ! File ID of the variable
                    ), &
                __LINE__ )

            ! Compress if we are using netcdf4
            USING_NETCDF4_GUARD( call check(nf90_def_var_deflate(iNcid, nc_grid_output%time_var_id(i), 1, 1, 3), __LINE__) )
        end do

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
        call check(nf90_def_var(iNcid, "priority_domain_index", output_int2, &
            [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
            nc_grid_output%var_priority_ind), __LINE__ )
        ! priority domain image
        call check(nf90_def_var(iNcid, "priority_domain_image", output_int2, &
            [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
            nc_grid_output%var_priority_img), __LINE__ )


        ! Use nontemporal_grids_to_store to define which variables to store.
        !
        ! Check for typos. It can contain either an empty string '', or
        ! a value in nc_grid_output%static_var_names
        do i = 1, size(nontemporal_grids_to_store)
            if(.not. (nontemporal_grids_to_store(i) == '' .or. &
                  any(nontemporal_grids_to_store(i) == nc_grid_output%static_var_names))) then
                write(log_output_unit, *)&
                    "Unknown value in nontemporal_grids_to_store:", trim(nontemporal_grids_to_store(i))
                call generic_stop
            end if
        end do
        ! Set the value of static_var_store_flag so that the variables that we want to store are stored.
        do i = 1, size(nc_grid_output%static_var_names)
            nc_grid_output%static_var_store_flag(i) = &
                any(nc_grid_output%static_var_names(i) == nontemporal_grids_to_store)
        end do

        ! Storage for the static variables
        do i = 1, size(nc_grid_output%static_var_names)

            if(.not. nc_grid_output%static_var_store_flag(i)) cycle

            call check( &
                nf90_def_var(iNcid, &
                    trim(nc_grid_output%static_var_names(i)), &
                    output_prec, &
                    [nc_grid_output%dim_x_id, nc_grid_output%dim_y_id], &
                    nc_grid_output%static_var_id(i)), &
                __LINE__ )
        end do

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

        integer:: iNcid, nxy(2), spatial_start(2), spatial_stride(2), i
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
            start=[nc_grid_output%num_output_steps]), __LINE__, stop_on_error=.false.)


        ! Save the flow variables
        do i = STG, ELV
            if( .not. nc_grid_output%time_var_store_flag(i)) cycle

            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%time_var_id(i), &
                U(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2), i), &
                start=[1,1,nc_grid_output%num_output_steps]),&
                __LINE__, stop_on_error=.false.)
        end do

#ifdef DEBUG_ARRAY
        if(present(debug_array)) then
            call check(nf90_put_var(&
                iNcid, &
                nc_grid_output%var_debug_array_id, &
                debug_array(spatial_start(1):nxy(1):spatial_stride(1),spatial_start(2):nxy(2):spatial_stride(2)), &
                start=[1,1,nc_grid_output%num_output_steps]), &
                __LINE__, stop_on_error=.false.)
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


    subroutine store_static_variable(nc_grid_output, variable_name, values)
        !!
        !! Store a static flow variable (e.g. 'max-stage'), which is a variable desribing the flow that we only store once.
        !!
        class(nc_grid_output_type), intent(in) :: nc_grid_output
        character(*), intent(in) :: variable_name
        real(dp), intent(in) :: values(:,:)

        integer :: iNcid, spatial_start(2), nxy(2), spatial_stride(2), j, my_j
        logical :: found_var

        ! Gracefully compile in case we don't link with netcdf
#ifndef NONETCDF
        spatial_start = nc_grid_output%spatial_start
        spatial_stride = nc_grid_output%spatial_stride
        nxy(1) = size(values, dim=1)
        nxy(2) = size(values, dim=2)

        ! Shorthand
        iNcid = nc_grid_output%nc_file_id

        found_var = .false.
        do j = 1, size(nc_grid_output%static_var_names)
            ! Only one 'j' will correspond to 'variable_name'
            if((trim(variable_name) == nc_grid_output%static_var_names(j)) .and. &
               nc_grid_output%static_var_store_flag(j)) then

                ! Save the variable
                call check(nf90_put_var(&
                    iNcid, &
                    nc_grid_output%static_var_id(j), &
                    values(spatial_start(1):nxy(1):spatial_stride(1), spatial_start(2):nxy(2):spatial_stride(2))), &
                    __LINE__, stop_on_error=.false.)

                found_var = .true. ! If this is never set, we know that no variable was saved.
                exit
            end if
        end do

        if(.not. found_var) then
            write(log_output_unit, *) "Cannot store static variable: ", trim(variable_name)
            call generic_stop
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
