module read_raster_mod 
    ! Fortran 2008 code to call C and read from a gdal raster
    !
    ! Both a type-based interface and a basic procedural interface are provided.
    !
    ! The type-based interface is a bit more flexible in that you can open
    ! the file once, then read many times, and then close it, and also make various
    ! enquiries about the file. See 'gdal_raster_dataset_type'.
    !
    ! The procedural interface involves opening/closing the input file each
    ! time a subroutine is called. This is not necessarily a problem, but could
    ! cause issues. See 'read_gdal_raster' and 'get_raster_dimensions'.
    !
    use global_mod, only: charlen, ip, dp
    use iso_c_binding

    implicit none

    private

    ! object based interface
    public:: gdal_raster_dataset_type
    ! procedural interface
    public:: read_gdal_raster, get_raster_dimensions
    ! tests
    public:: test_read_raster1

    ! Flag to use bilinear interpolation (1) or not (0) by default
    integer(C_INT), parameter :: use_bilinear_default = 1

    ! Type to hold c pointers and related info for the raster
    type :: gdal_raster_dataset_type
        type(c_ptr) :: hDataset
        character(charlen):: inputFile
        logical:: isOpen = .FALSE.
        integer(C_INT) :: xydim(2)
        real(C_DOUBLE):: lowerleft(2), upperright(2), adfGeoTransform(6), dx(2)
    
        contains
        procedure:: initialise => initialise_gdal_raster_dataset
        procedure:: get_xy => get_xy_values
        procedure:: finalise => close_gdal_raster
    end type

    interface

    ! Here we interface many C routines

    ! Object based raster access
    subroutine open_gdal_raster_Cfun(inputFile, hDataset) bind(C, name='open_gdal_raster')
        use iso_c_binding
        implicit none
        character(kind=C_CHAR), intent(in) :: inputFile(*)
        type(c_ptr), intent(out) :: hDataset
    end subroutine
    
    ! Object based raster access
    subroutine close_gdal_raster_Cfun(hDataset) bind(C, name='close_gdal_raster')
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: hDataset
    end subroutine

    ! Object based raster access
    subroutine get_gdal_raster_dimensions_Cfun(hDataset, xydim, lowerleft, upperright, &
        adfGeoTransform, dx) bind(C, name='get_gdal_raster_dimensions')
        use iso_c_binding
        implicit none
        type(C_PTR), intent(in) :: hDataset
        integer(C_INT), intent(out):: xydim(2)
        real(C_DOUBLE), intent(out):: lowerleft(2), upperright(2), adfGeoTransform(6), dx(2)
    end subroutine

    ! Object based raster access
    subroutine get_values_at_xy_Cfun(hDataset, adfGeoTransform, x, y, z, N, verbose, bilinear, &
        band) bind(C, name='get_values_at_xy')
        use iso_c_binding
        implicit none
        type(C_PTR), intent(in) :: hDataset
        integer(C_INT), value, intent(in):: N
        real(C_DOUBLE) :: adfGeoTransform(6)
        real(C_DOUBLE) :: x(N), y(N), z(N)
        integer(C_INT), value, intent(in):: verbose, bilinear, band
    end subroutine

    ! Non-object based interface
    subroutine read_gdal_raster_Cfun(inputFile, x, y, z, N, verbose, bilinear) bind(C, name='read_gdal_raster')
        use iso_c_binding
        implicit none
        character(kind = C_CHAR), intent(in) :: inputFile(*)
        integer(C_INT), value, intent(in) :: N, verbose, bilinear
        real(C_DOUBLE), intent(in) :: x(N), y(N)
        real(C_DOUBLE), intent(out) :: z(N)
        
    end subroutine

    ! Non-object based interface
    subroutine get_raster_dimensions_Cfun(inputFile, xydim, lowerleft, upperright) bind(C, name='get_raster_dimensions')
        use iso_c_binding
        implicit none
        character(kind = C_CHAR), intent(in) :: inputFile(*)
        real(C_DOUBLE), intent(out):: lowerleft(2), upperright(2)
        integer(C_INT), intent(out):: xydim(2)

    end subroutine
    
    end interface

    !! read_gdal_raster is a non-object based interface
    !! It can handle C_DOUBLE and C_FLOAT as input
    !! In the latter case we need to copy input data, which increases memory
    !! so it's good to use a different routine. 
    !! Also, it is based around passing x,y (at which z is desired) in one go,
    !! which can force us to use lots of memory (temporaries for x,y)
    interface read_gdal_raster
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

    !! Finalise routine for gdal_raster_dataset_type
    subroutine close_gdal_raster(gdal_raster_dataset)
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
            gdal_raster_dataset%adfGeoTransform, gdal_raster_dataset%dx) 

    end subroutine

    !! This is the main initialiser of gdal_raster_dataset_type
    subroutine initialise_gdal_raster_dataset(gdal_raster_dataset, inputFile)
        class(gdal_raster_dataset_type), intent(inout):: gdal_raster_dataset
        character(len=charlen), intent(in) :: inputFile

        call open_gdal_raster(inputFile, gdal_raster_dataset)
        call get_gdal_raster_dimensions(gdal_raster_dataset)

    end subroutine

    !! This is the key useful routine for gdal_raster_dataset_type
    subroutine get_xy_values(gdal_raster_dataset, x, y, z, N, verbose, bilinear, band)
        class(gdal_raster_dataset_type), intent(in):: gdal_raster_dataset
        integer(ip), intent(in):: n
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
            verbose_c = 0
        end if

        if(present(bilinear)) then
            bilinear_c = bilinear
        else
            bilinear_c = 1 
        end if

        if(present(band)) then
            band_c = band
        else
            band_c = 1
        end if


        call get_values_at_xy_Cfun(gdal_raster_dataset%hDataset, &
            gdal_raster_dataset%adfGeoTransform, x_c, y_c, z_c, N_c, &
            verbose_c, bilinear_c, band_c)

        z = real(z_c, dp)

    end subroutine

    !! Non-object-based interface
    subroutine get_raster_dimensions(inputFile, xydim, lowerleft, upperright)
        character(len=charlen), intent(in) :: inputFile
        integer(ip), intent(out):: xydim(2)
        real(dp), intent(out):: lowerleft(2), upperright(2)

        real(C_DOUBLE):: lowerleft_c(2), upperright_c(2)
        integer(C_INT):: xydim_c(2)
        character(kind = C_CHAR, len = len_trim(inputFile) + 1) :: inputFile_c

        inputFile_c = trim(inputFile) // C_NULL_CHAR

        call get_raster_dimensions_Cfun(inputFile_c, xydim_c, lowerleft_c, upperright_c)

        xydim = xydim_c
        lowerleft = lowerleft_c
        upperright = upperright_c

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
        N_c = size(x)
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
        N_c = size(x)
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


    subroutine test_read_raster1()
        ! Test code for read_raster_mod
        implicit none

        character(len=charlen) :: inputFile
        integer(ip), parameter:: N = 10, verbose = 0
        integer(ip) :: i
        integer, parameter:: cdp = C_DOUBLE, csp = C_FLOAT
        real(cdp), allocatable :: x(:), y(:), z(:), real_z(:), real_z_bl(:)
        real(csp), allocatable :: xf(:), yf(:), zf(:), real_zf(:), real_zf_bl(:)

        integer(ip):: xydim(2)
        real(dp):: lowerleft(2), upperright(2)
        type(gdal_raster_dataset_type):: test_file

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
       

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Test gdal raster dataset class
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        call close_gdal_raster(test_file)

        if( (test_file%isOpen .eqv. .FALSE.) .AND. &
            (c_associated(test_file%hDataset) .eqv. .FALSE.)) then
            print*, 'PASS' 
        else
            print*, 'FAIL' 
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Test reading of raster summary data
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Test of raster IO and bilinear interpolation
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

    end subroutine
end module read_raster_mod 


