module read_raster_mod
    !! Read GDAL supported rasters. Both a type-based interface and a basic procedural interface are provided.
    !!
    !! The type-based interface (using the multi_raster_type) is more flexible and should
    !! be preferred. Given a set of raster filenames (in order of preference) it can
    !! extract raster data at a given set of coordinates. At each coordinates, it will use
    !! the highest-preference raster that has non-missing data. Constant and bilinear interpolation
    !! are supported. The code does not fully read the files into memory.
    !!
    !! The procedural interface is not recommended anymore, because it does not separate "opening and closing
    !! the files" from "reading the data". But see 'read_gdal_raster' and 'get_raster_dimensions'.
    !!
    use global_mod, only: charlen, ip, dp
    use stop_mod, only: generic_stop
    use logging_mod, only: log_output_unit
    use iso_c_binding

    implicit none

    private

    ! object based interface
    public:: gdal_raster_dataset_type, multi_raster_type
    ! procedural interface
    public:: read_gdal_raster, get_raster_dimensions
    ! tests
    public:: test_read_raster1

    ! Flag to use bilinear interpolation (1) or not (0) by default
    integer(C_INT), parameter :: use_bilinear_default = 1_C_INT

    type :: gdal_raster_dataset_type
        !!
        !! Type to work with a single raster dataset.
        !!

        type(c_ptr) :: hDataset !! Pointer to GDAL raster data structure
        character(charlen):: inputFile !! Filename
        logical:: isOpen = .FALSE. !! Record if the file is open
        integer(C_INT) :: xydim(2) !! x-y dimensions of the raster
        real(C_DOUBLE) :: lowerleft(2) !! Lower left (corner) of raster
        real(C_DOUBLE) :: upperright(2) !! Upper right (corner) of raster
        real(C_DOUBLE) :: adfGeoTransform(6) !! GDAL's geotransform information
        real(C_DOUBLE) :: dx(2) !! Raster cell size
        real(C_DOUBLE) :: nodata_value !! The raster's nodata value

        contains
        procedure:: initialise => initialise_gdal_raster_dataset !! Initialise the object
        procedure:: get_xy => get_xy_values !! Interpolate at coordinates
        procedure:: finalise => close_gdal_raster !! Close the object
        procedure:: print => print_summary_gdal_raster_dataset !! Print summary information

    end type

    type :: multi_raster_type
        !!
        !! Type to interpolate data from multiple rasters with different resoltuions and spatial extents.
        !!
        !! Often in applications we have multiple raster datasets for interpolation, and have an order in which
        !! we would like to use them. This type will manage that.
        !!
        !! **The raster_files should all have the same projection, but they can have different resolutions,
        !! extents, and they can overlap. Where the raster overlap, their preference order (for interpolation)
        !! corresponds to the order of the input file names.**
        !!

        character(len=charlen), allocatable :: raster_files(:) !! Raster filenames in order of preference
        type(gdal_raster_dataset_type), allocatable :: raster_datasets(:) !! Objects corresponding to each file

        real(C_DOUBLE) :: lowerleft(2) !! minimum lower left (corner) of all rasters combined
        real(C_DOUBLE) :: upperright(2) !! maximum upper right (corner) of all rasters combined

        contains
        procedure:: initialise => initialise_multi_raster !! Initialise the object
        procedure:: get_xy => get_xy_values_multi_raster !! Interpolate at coordinates
        procedure:: finalise => close_multi_raster !! Close the object

    end type


    interface

    !
    ! Here we interface many C routines
    !

    subroutine open_gdal_raster_Cfun(inputFile, hDataset) bind(C, name='open_gdal_raster')
        !! Open raster via GDAL
        use iso_c_binding
        implicit none
        character(kind=C_CHAR), intent(in) :: inputFile(*) !! Input filename
        type(c_ptr), intent(out) :: hDataset !! GDAL raster object
    end subroutine

    subroutine close_gdal_raster_Cfun(hDataset) bind(C, name='close_gdal_raster')
        !! Close the raster via GDAL
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: hDataset !! GDAL raster object
    end subroutine

    subroutine get_gdal_raster_dimensions_Cfun(hDataset, xydim, lowerleft, upperright, &
        adfGeoTransform, dx, nodata_value) bind(C, name='get_gdal_raster_dimensions')
        !! Store the raster dimensions, lower-left and upper-right coordinates, cellsize, nodata value, and some GDAL geo-transform
        !! information
        use iso_c_binding
        implicit none
        type(C_PTR), intent(in) :: hDataset !! GDAL raster object
        integer(C_INT), intent(inout):: xydim(2) !! x/y dimensions
        real(C_DOUBLE), intent(inout):: lowerleft(2), upperright(2), adfGeoTransform(6), dx(2), nodata_value
    end subroutine

    subroutine get_values_at_xy_Cfun(hDataset, adfGeoTransform, x, y, z, N, verbose, bilinear, &
        band) bind(C, name='get_values_at_xy')
        !! Get the raster values at the given x/y coordinates
        use iso_c_binding
        implicit none
        type(C_PTR), intent(in) :: hDataset !! GDAL raster object
        integer(C_INT), value, intent(in):: N !! Number of x/y/z points
        real(C_DOUBLE) :: adfGeoTransform(6) !! GDAL geo-transform data
        real(C_DOUBLE) :: x(N), y(N), z(N)
        integer(C_INT), value, intent(in):: verbose, bilinear, band
    end subroutine

    subroutine read_gdal_raster_Cfun(inputFile, x, y, z, N, verbose, bilinear) bind(C, name='read_gdal_raster')
        !! Procedural interface for reading the raster. This routine opens the inputFile, reads the data and closes the file before
        !! returning.
        use iso_c_binding
        implicit none
        character(kind = C_CHAR), intent(in) :: inputFile(*) !! Raster file name
        integer(C_INT), value, intent(in) :: N, verbose, bilinear
        real(C_DOUBLE), intent(in) :: x(N), y(N)
        real(C_DOUBLE), intent(out) :: z(N)

    end subroutine

    subroutine get_raster_dimensions_Cfun(inputFile, xydim, lowerleft, upperright) bind(C, name='get_raster_dimensions')
        use iso_c_binding
        !! Procedural interface for reading raster dimensions, lower-left and upper-right coordinates.
        implicit none
        character(kind = C_CHAR), intent(in) :: inputFile(*)
        real(C_DOUBLE), intent(out):: lowerleft(2), upperright(2)
        integer(C_INT), intent(out):: xydim(2)

    end subroutine

    end interface

    interface read_gdal_raster
        !! read_gdal_raster is a procedural interface. Each call opens the file, extracts the data, and closes the file,
        !! whereas the object based interface (multi_raster_type & gdal_raster_dataset_type) can separate opening & closing from reading.
        module procedure read_gdal_raster_C_DOUBLE, read_gdal_raster_C_FLOAT, &
            read_gdal_raster_C_DOUBLE_2D, read_gdal_raster_C_FLOAT_2D
    end interface

    contains

    subroutine open_gdal_raster(inputFile, gdal_raster_dataset)
        character(len=charlen), intent(in) :: inputFile
        type(gdal_raster_dataset_type), intent(inout):: gdal_raster_dataset

        character(kind = C_CHAR, len=len_trim(inputFile) + 1) :: inputFile_c

        inputFile_c = trim(inputFile) // C_NULL_CHAR

        gdal_raster_dataset%inputFile = inputFile

        call open_gdal_raster_Cfun(inputFile_c, gdal_raster_dataset%hDataset)

        gdal_raster_dataset%isOpen = .TRUE.

    end subroutine

    subroutine close_gdal_raster(gdal_raster_dataset)
        !! Finalise routine for gdal_raster_dataset_type
        class(gdal_raster_dataset_type), intent(inout):: gdal_raster_dataset

        call close_gdal_raster_Cfun(gdal_raster_dataset%hDataset)
        gdal_raster_dataset%isOpen = .FALSE.

    end subroutine

    subroutine get_gdal_raster_dimensions(gdal_raster_dataset)
        type(gdal_raster_dataset_type), intent(inout):: gdal_raster_dataset

        call get_gdal_raster_dimensions_Cfun(&
            gdal_raster_dataset%hDataset,&
            gdal_raster_dataset%xydim, &
            gdal_raster_dataset%lowerleft, gdal_raster_dataset%upperright, &
            gdal_raster_dataset%adfGeoTransform, gdal_raster_dataset%dx, &
            gdal_raster_dataset%nodata_value)

    end subroutine

    subroutine initialise_gdal_raster_dataset(gdal_raster_dataset, inputFile)
        !! This is the main initialiser of gdal_raster_dataset_type
        class(gdal_raster_dataset_type), intent(inout):: gdal_raster_dataset
        character(len=charlen), intent(in) :: inputFile

        call open_gdal_raster(inputFile, gdal_raster_dataset)
        call get_gdal_raster_dimensions(gdal_raster_dataset)

    end subroutine

    subroutine get_xy_values(gdal_raster_dataset, x, y, z, N, verbose, bilinear, band)
        !! This is the key useful routine for gdal_raster_dataset_type
        class(gdal_raster_dataset_type), intent(in):: gdal_raster_dataset
        integer(ip), intent(in):: N
        real(dp), intent(in):: x(n), y(n)
        real(dp), intent(out):: z(n)
        integer(ip), intent(in), optional:: verbose
        integer(ip), intent(in), optional:: bilinear
        integer(ip), intent(in), optional:: band

        real(C_DOUBLE):: x_c(N), y_c(N), z_c(N)
        integer(C_INT):: verbose_c, bilinear_c, band_c, N_c

        ! Convert to C_DOUBLE
        x_c = real(x, C_DOUBLE)
        y_c = real(y, C_DOUBLE)
        z_c = x*0.0_C_DOUBLE
        N_c = N

        if(present(verbose)) then
            verbose_c = verbose
        else
            verbose_c = 0_C_INT
        end if

        if(present(bilinear)) then
            bilinear_c = bilinear
        else
            bilinear_c = use_bilinear_default
        end if

        if(present(band)) then
            band_c = band
        else
            band_c = 1_C_INT
        end if


        call get_values_at_xy_Cfun(gdal_raster_dataset%hDataset, &
            gdal_raster_dataset%adfGeoTransform, x_c, y_c, z_c, N_c, &
            verbose_c, bilinear_c, band_c)

        z = real(z_c, dp)

    end subroutine

    subroutine get_raster_dimensions(inputFile, xydim, lowerleft, upperright)
        !! Non-object-based interface
        !! Avoid this
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(out):: xydim(2)
        real(dp), intent(out):: lowerleft(2), upperright(2)

        real(C_DOUBLE):: lowerleft_c(2), upperright_c(2)
        integer(C_INT):: xydim_c(2)
        character(kind = C_CHAR, len = len_trim(inputFile) + 1) :: inputFile_c

        inputFile_c = trim(inputFile) // C_NULL_CHAR

        call get_raster_dimensions_Cfun(inputFile_c, xydim_c, lowerleft_c, upperright_c)

        xydim = xydim_c
        lowerleft = real(lowerleft_c, dp)
        upperright = real(upperright_c, dp)

    end subroutine

    ! Read gdal raster, assuming x,y,z are rank 1 arrays with N entries of type C_DOUBLE
    subroutine read_gdal_raster_C_DOUBLE(inputFile, x, y, z, N, verbose, bilinear)
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(in) :: N, verbose
        integer(ip), optional, intent(in) :: bilinear
        real(C_DOUBLE), intent(in) :: x(N), y(N)
        real(C_DOUBLE), intent(out) :: z(N)

        character(kind = C_CHAR, len = len_trim(inputFile) + 1) :: inputFile_c
        integer(C_INT) :: N_c, verbose_c, use_bilinear_local


        if(present(bilinear)) then
            use_bilinear_local = bilinear
        else
            use_bilinear_local = use_bilinear_default
        end if

        ! Use C types
        inputFile_c = trim(inputFile) // C_NULL_CHAR
        N_c = size(x, kind=C_INT)
        verbose_c = 1_C_INT * verbose

        call read_gdal_raster_Cfun(inputFile_c, x, y, z, N_c, verbose_c, use_bilinear_local)

    end subroutine

    ! Read gdal raster, assuming x,y,z are rank 2 arrays with N entries of type C_DOUBLE
    subroutine read_gdal_raster_C_DOUBLE_2D(inputFile, x, y, z, N, verbose, bilinear)
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(in) :: N, verbose
        integer(ip), optional, intent(in) :: bilinear
        real(C_DOUBLE), intent(in) :: x(:,:), y(:,:)
        real(C_DOUBLE), intent(out) :: z(:,:)

        integer(C_INT) :: use_bilinear_local

        if(present(bilinear)) then
            use_bilinear_local = bilinear
        else
            use_bilinear_local = use_bilinear_default
        end if

        call read_gdal_raster_C_DOUBLE(inputFile, x, y, z, N, verbose, use_bilinear_local)

    end subroutine

    ! Read gdal raster, assuming x,y,z are rank 1 arrays with N entries of type C_FLOAT
    subroutine read_gdal_raster_C_FLOAT(inputFile, x, y, z, N, verbose, bilinear)
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(in) :: N, verbose
        integer(ip), optional, intent(in) :: bilinear
        real(C_FLOAT), intent(in) :: x(N), y(N)
        real(C_FLOAT), intent(out) :: z(N)

        character(kind = C_CHAR, len = len_trim(inputFile) + 1) :: inputFile_c
        integer(C_INT) :: N_c, verbose_c, use_bilinear_local
        ! Copy x,y,z into these C_DOUBLE arrays to allow the routine to work
        ! with either single or double precision
        real(C_DOUBLE), allocatable :: x0(:), y0(:), z0(:)

        if(present(bilinear)) then
            use_bilinear_local = bilinear
        else
            use_bilinear_local = use_bilinear_default
        end if


        ! Use C types
        inputFile_c = trim(inputFile) // C_NULL_CHAR
        N_c = size(x, kind=C_INT)
        verbose_c = 1_C_INT * verbose

        ! If x,y,z are not C_DOUBLE, we copy C_DOUBLE versions
        allocate(x0(N_c), y0(N_c), z0(N_c))
        x0 = real(x, C_DOUBLE)
        y0 = real(y, C_DOUBLE)
        z0 = real(z, C_DOUBLE)
        call read_gdal_raster_Cfun(inputFile_c, x0, y0, z0, N_c, verbose_c, use_bilinear_local)
        z = real(z0, C_FLOAT)
        deallocate(x0, y0, z0)

    end subroutine

    ! Read gdal raster, assuming x,y,z are rank 2 arrays with N entries of type C_FLOAT
    subroutine read_gdal_raster_C_FLOAT_2D(inputFile, x, y, z, N, verbose, bilinear)
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(in) :: N, verbose
        integer(ip), optional, intent(in) :: bilinear
        real(C_FLOAT), intent(in) :: x(:,:), y(:,:)
        real(C_FLOAT), intent(out) :: z(:,:)
        integer(C_INT) :: use_bilinear_local

        if(present(bilinear)) then
            use_bilinear_local = bilinear
        else
            use_bilinear_local = use_bilinear_default
        end if

        call read_gdal_raster_C_FLOAT(inputFile, x, y, z, N, verbose, use_bilinear_local)

    end subroutine

    subroutine initialise_multi_raster(multi_raster, raster_files)
        !!
        !! Read a number of rasters into a single multi_raster object. The
        !! rasters should all have the same projection, but can have different resolutions
        !! and extents.
        !!
        class(multi_raster_type), intent(inout) :: multi_raster
        character(len=charlen) :: raster_files(:)

        integer(ip) :: n, i

        if(allocated(multi_raster%raster_datasets)) then
            write(log_output_unit, *) 'multi_raster already allocated'
            call generic_stop
        end if

        n = size(raster_files, kind=ip)
        allocate(multi_raster%raster_files(n))
        multi_raster%raster_files = raster_files

        allocate(multi_raster%raster_datasets(n))

        do i = 1, n
            call multi_raster%raster_datasets(i)%initialise(raster_files(i))
        end do

        !! Make the lowerleft/upperright of the multiraster correspond to the extremes of the individual rasters
        multi_raster%lowerleft = [HUGE(1.0_C_DOUBLE) , HUGE(1.0_C_DOUBLE)]
        multi_raster%upperright = [-HUGE(1.0_C_DOUBLE), -HUGE(1.0_C_DOUBLE)]
        do i = 1, n
            multi_raster%lowerleft  = min(multi_raster%lowerleft , multi_raster%raster_datasets(i)%lowerleft )
            multi_raster%upperright = max(multi_raster%upperright, multi_raster%raster_datasets(i)%upperright)
        end do

    end subroutine

    subroutine close_multi_raster(multi_raster)
        !!
        !! Close/cleanup a multi raster object
        !!
        class(multi_raster_type), intent(inout) :: multi_raster

        integer(ip) :: i, n

        do i = 1, size(multi_raster%raster_datasets, kind=ip)
            call multi_raster%raster_datasets(i)%finalise()
        end do

        deallocate(multi_raster%raster_datasets)
        deallocate(multi_raster%raster_files)

    end subroutine

    subroutine get_xy_values_multi_raster(multi_raster, x, y, z, N, verbose, bilinear, band, na_below_limit)
        !!
        !! Interpolation from multi_raster
        !!
        class(multi_raster_type), intent(in) :: multi_raster !! a multi_raster object
        integer(ip), intent(in) :: N !! Number of x/y/z coordinates
        real(dp), intent(in) :: x(N), y(N) !! Input coordinates
        real(dp), intent(inout):: z(N) !! Output raster values at x/y
        integer(ip), optional, intent(in) :: verbose !! Print information
        integer(ip), optional, intent(in) :: bilinear !! Use bilinear interpolation (1) or nearest-cell interpolation (0)
        integer(ip), optional, intent(in) :: band !! Integer giving the raster band to interpolate from
        real(dp), optional, intent(in) :: na_below_limit
        !! Treat raster values below this number as NA. This can be useful if nodata values are not preserved exactly (e.g. due to
        !! changes in precision). Often nodata values are large negative numbers, in which case we can be sure that all numbers below
        !! some threshold (e.g. `na_below_limit = -1.0e+10`) should be treated as NA.

        real(dp) :: empty_value, ll(2), ur(2), lower_limit_l
        integer(ip) :: i, j, verbose_l, bilinear_l, band_l

        if(present(verbose)) then
            verbose_l = verbose
        else
            verbose_l = 0
        end if

        if(present(bilinear)) then
            bilinear_l = bilinear
        else
            bilinear_l = use_bilinear_default
        end if

        if(present(band)) then
            band_l = band
        else
            band_l = 1
        end if

        if(present(na_below_limit)) then
            lower_limit_l = na_below_limit
        else
            lower_limit_l = -HUGE(1.0_dp)
        end if

        ! Flag for unset 'z'
        empty_value = -huge(1.0_dp)
        z = empty_value

        ! Read rasters until no values are missing
        do j = 1, size(multi_raster%raster_datasets, kind=ip)

            ll = real(multi_raster%raster_datasets(j)%lowerleft, dp)
            ur = real(multi_raster%raster_datasets(j)%upperright, dp)

            do i = 1, N

                ! If z is populated with a value, we are done
                if(z(i) /= empty_value) cycle

                ! Read values inside the raster extent
                if(x(i) >= ll(1) .and. x(i) <= ur(1) .and. y(i) >= ll(2) .and. y(i) <= ur(2)) then
                    call multi_raster%raster_datasets(j)%get_xy(x(i), y(i), z(i), 1, verbose_l, bilinear_l, band_l)
                end if

                ! Set 'nodata' values back to empty values
                if(z(i) == real(multi_raster%raster_datasets(j)%nodata_value, dp)) z(i) = empty_value
                if(z(i) < lower_limit_l) z(i) = empty_value

            end do

        end do

    end subroutine

    subroutine print_summary_gdal_raster_dataset(raster_data)
        !! Print summary information about the raster data
        class(gdal_raster_dataset_type) :: raster_data

        write(log_output_unit, *) 'inputFile   : ', trim(raster_data%inputFile)
        write(log_output_unit, *) '  isOpen     : ', raster_data%isOpen
        write(log_output_unit, *) '  xydim      : ', raster_data%xydim
        write(log_output_unit, *) '  lowerleft  : ', raster_data%lowerleft
        write(log_output_unit, *) '  upperright : ', raster_data%upperright
        write(log_output_unit, *) '  dx         : ', raster_data%dx
        write(log_output_unit, *) '  nodata     : ', raster_data%nodata_value

    end subroutine

    subroutine test_read_raster1()
        !! Unit tests for read_raster_mod
        implicit none

        character(len=charlen) :: inputFile, inputFiles(2)
        integer(ip), parameter:: N = 10, verbose = 0
        integer(ip) :: i
        integer, parameter:: cdp = C_DOUBLE, csp = C_FLOAT
        real(cdp), allocatable :: x(:), y(:), z(:), real_z(:), real_z_bl(:)
        real(csp), allocatable :: xf(:), yf(:), zf(:), real_zf(:), real_zf_bl(:)

        integer(ip):: xydim(2)
        real(dp):: lowerleft(2), upperright(2)
        type(gdal_raster_dataset_type):: test_file
        type(multi_raster_type) :: test_multi_raster

        real(dp), allocatable:: x_dp(:), y_dp(:), z_dp(:)

        ! This requires that the file exists in the subdirectory 'data' where
        ! the program is run
        inputFile = "./data/test_rast.tif"

        allocate(x(N), y(N), z(N), real_z(N), real_z_bl(N), x_dp(N), y_dp(N), z_dp(N))
        allocate(xf(N), yf(N), zf(N), real_zf(N), real_zf_bl(N))

        ! Points in the raster we will look up
        x = [502523.122040_cdp, 506073.391951_cdp, 504937.667827_cdp, 501782.285213_cdp, 502901.704270_cdp, &
             505457.970585_cdp, 500256.751295_cdp, 502305.114871_cdp, 501575.432335_cdp, 502512.685946_cdp ]

        y = [1617967.96301_cdp, 1614589.30853_cdp, 1618488.11768_cdp, 1609423.47632_cdp, 1618325.86200_cdp, &
            1613200.33784_cdp, 1616185.67459_cdp, 1613888.10362_cdp, 1618757.31551_cdp, 1617615.75308_cdp]

        ! Values with 'nearest cell' interpolation
        real_z = [3725100, 3735700, 3729300, 3731900, 3725500, 3735700, 3722500, 3728700, 3722300, 3725300]

        ! Values with bilinear interpolation
        real_z_bl = [3725078.28107_cdp, 3735557.47537_cdp, 3729387.21797_cdp, 3732141.09411_cdp, 3725477.54654_cdp, &
                     3735715.60333_cdp, 3722327.82800_cdp, 3728722.12613_cdp, 3722393.54915_cdp, 3725409.61881_cdp]

        x_dp = real(x, dp)
        y_dp = real(y, dp)
        z_dp = real(x, dp)*0

        ! Single precision variants ('float')
        xf = real(x, csp)
        yf = real(y, csp)
        real_zf = real(real_z, csp)
        real_zf_bl = real(real_z_bl, csp)


        !
        !
        ! Test gdal raster dataset class
        !
        !

        call test_file%initialise(inputFile)

        if((test_file%inputFile == inputFile) .AND. &
            c_associated(test_file%hDataset) .AND. &
            test_file%isOpen) then
            print*, 'PASS'
        else
            print*, 'FAIL'

        end if

        if( (test_file%dx(1) == 200.0_C_DOUBLE) .AND. (test_file%dx(2) == -200.0_C_DOUBLE) .AND. &
            (test_file%xydim(1) == 50_C_INT) .AND. (test_file%xydim(2) == 50_C_INT) .AND. &
            (test_file%lowerleft(1) == 499000.0_C_DOUBLE) .AND. (test_file%lowerleft(2) == 1609000.0_C_DOUBLE) .AND. &
            (test_file%upperright(1) == 509000.0_C_DOUBLE) .AND. (test_file%upperright(2) == 1619000.0_C_DOUBLE)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            print*, test_file%dx
            print*, test_file%xydim
            print*, test_file%lowerleft
            print*, test_file%upperright
        end if

        call test_file%get_xy(x_dp, y_dp, z_dp, size(z_dp), verbose=0, bilinear=0, band=1)
        if(all(abs(z_dp - real_z) < (1.0e-6_dp * real_z))) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        call test_file%get_xy(x_dp, y_dp, z_dp, size(z_dp), verbose=0, bilinear=1, band=1)
        if(all(abs(z_dp - real_z_bl) < (1.0e-6_dp*real_z_bl))) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            print*, z_dp
            print*, real_z_bl
            print*, abs(z_dp - real_z_bl)
        end if

        call test_file%finalise()

        if( (test_file%isOpen .eqv. .FALSE.) .AND. &
            (c_associated(test_file%hDataset) .eqv. .FALSE.)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        !
        ! Test reading of raster summary data
        !
        !
        call get_raster_dimensions(inputFile, xydim, lowerleft, upperright)
        if( (xydim(1) == 50) .and. (xydim(2) == 50) ) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        if( (lowerleft(1) == 499000.0_dp) .and. (lowerleft(2) == 1609000.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        if( (upperright(1) == 509000.0_dp) .and. (upperright(2) == 1619000.0_dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        !
        !
        ! Test of raster IO and bilinear interpolation
        !
        !


        ! Test double precision version
        call read_gdal_raster(inputFile, x, y, z, N, verbose, bilinear=0_ip)

        ! Got these values from R's extract
        if(all(abs(z - real_z) < 1.0e-10_cdp*real_z)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            do i = 1, N
                print*, real_z(i), z(i)
            end do
        end if

        ! Test bilinear version -- got these values from R's extract with
        ! method='bilinear'.
        call read_gdal_raster(inputFile, x, y, z, N, verbose, bilinear=1_ip)
        if(all(abs(z - real_z_bl) < 1.0e-10_cdp*real_z_bl)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            do i = 1, N
                print*, real_z_bl(i), z(i)
            end do
        end if


        ! Now try single precision version
        call read_gdal_raster(inputFile, xf, yf, zf, N, verbose, bilinear=0_ip)

        if(all(abs(zf - real_zf) < 1.0e-6_csp*real_zf)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            do i = 1, N
                print*, real_zf(i), zf(i)
            end do
        end if

        ! Test bilinaer single precision version
        CALL read_gdal_raster(inputFile, xf, yf, zf, N, verbose, bilinear=1_ip)

        if(all(abs(zf - real_zf_bl) < 1.0e-6_csp*real_zf_bl)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            do i = 1, N
                print*, real_zf_bl(i), zf(i)
            end do
        end if


        !
        ! Test reading a raster with NaN values
        !

        inputFile = 'data/test_rast_nans.tif'
        call test_file%initialise(inputFile)
        !print*, test_file%nodata_value
        if( abs(test_file%nodata_value + 1.7e+308_C_DOUBLE) < 1.0e-06*(1.7e+308_C_DOUBLE)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Point right in a nodata cell. Check it is identified as nodata
        x_dp(1) = 508000.0_dp
        y_dp(1) = 1614000.0_dp
        call test_file%get_xy(x_dp(1:1), y_dp(1:1), z_dp(1:1), size(z_dp(1:1)), verbose=0, bilinear=1, band=1)
        if(z_dp(1) == real(test_file%nodata_value, dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        ! Point not exactly in a nodata cell. Check that it is still
        ! correctly set to 'nodata', rather than interpolated.
        x_dp(1) = 500901.0_dp
        y_dp(1) = 1617111.0_dp
        call test_file%get_xy(x_dp(1:1), y_dp(1:1), z_dp(1:1), size(z_dp(1:1)), verbose=0, bilinear=1, band=1)
        if(z_dp(1) == real(test_file%nodata_value, dp)) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if

        call test_file%finalise()

        !
        ! Test of multi-raster. Make the first raster contain NaN values, which
        ! 'pass through' to the second raster
        !
        inputFiles(1) = 'data/test_rast_nans.tif'
        inputFiles(2) = 'data/test_rast.tif'

        ! Test values
        x_dp = real(x, dp)
        y_dp = real(y, dp)
        z_dp = real(x, dp)*0
        ! ... with one point in an nan region
        x_dp(1) = 500901.0_dp
        y_dp(1) = 1617111.0_dp
        real_z_bl(1) = 3722691.0_dp

        call test_multi_raster%initialise(inputFiles)
        call test_multi_raster%get_xy(x_dp, y_dp, z_dp, size(z_dp, kind=ip), verbose=0, bilinear=1)

        if(all(abs(z_dp - real_z_bl) < (1.0e-6_dp*real_z_bl))) then
            print*, 'PASS'
        else
            print*, 'FAIL'
            print*, z_dp
            print*, real_z_bl
            print*, abs(z_dp - real_z_bl)
        end if

        call test_multi_raster%finalise()

        if(allocated(test_multi_raster%raster_files) .or. &
            allocated(test_multi_raster%raster_datasets)) then
            print*, 'FAIL'
        else
            print*, 'PASS'
        end if

    end subroutine
end module read_raster_mod


