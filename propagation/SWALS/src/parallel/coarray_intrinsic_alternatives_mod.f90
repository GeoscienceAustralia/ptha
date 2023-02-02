module coarray_intrinsic_alternatives
    !! Module to allow use of Fortran 2018 intrinsics co_sum and co_broadcast etc in limited situations, even
    !! if the compiler doesn't support them (e.g. true of ifort 2019).
    !! This code will only do something if compiled with -DCOARRAY_PROVIDE_CO_ROUTINES

    use iso_fortran_env, only: int32, real32, real64, real128
    use iso_c_binding, only: c_float, c_double
#if defined(COARRAY_PROVIDE_CO_ROUTINES) || defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
    use mpi
#endif

    implicit none

    private
    public sync_all_generic, this_image2, num_images2, swals_mpi_init, swals_mpi_finalize

#ifdef COARRAY_PROVIDE_CO_ROUTINES
    public :: co_broadcast, co_sum, co_max, co_min

    ! Can add more kinds/ranks as required
    interface co_broadcast
        module procedure co_broadcast_int32, co_broadcast_int32_rank1, co_broadcast_int32_rank2, &
            co_broadcast_c_float, co_broadcast_c_float_rank1, co_broadcast_c_float_rank2, &
            co_broadcast_c_double, co_broadcast_c_double_rank1, co_broadcast_c_double_rank2
    end interface

    interface co_sum
        module procedure co_sum_int32, co_sum_real32, co_sum_real64, co_sum_real128
    end interface

    interface co_max
        module procedure co_max_int32, co_max_real32, co_max_real64, co_max_real128
    end interface

    interface co_min
        module procedure co_min_int32, co_min_real32, co_min_real64, co_min_real128
    end interface
#endif

    contains

    subroutine sync_all_generic
        !! Alternative to "sync all" which uses mpi directly if COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined
        integer :: ierr
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#elif defined(COARRAY)
        sync all
#endif
    end subroutine

    integer function this_image2()
        !! Alternative to this_image() which uses mpi directly if COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined
        integer :: ierr
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call mpi_comm_rank(MPI_COMM_WORLD, this_image2, ierr)
        this_image2 = this_image2 + 1
#elif defined(COARRAY)
        this_image2 = this_image()
#else
        this_image2 = 1
#endif
    end function


    integer function num_images2()
        !! Alternative to num_images() which uses mpi directly if COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined
        integer :: ierr
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call mpi_comm_size(MPI_COMM_WORLD, num_images2, ierr)
#elif defined(COARRAY)
        num_images2 = num_images()
#else
        num_images2 = 1
#endif
    end function

    subroutine swals_mpi_init()
        !! Calls mpi_init(ierr) if COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined, and mpi is not initialised. Otherwise do nothing
        integer :: ierr
        logical :: flag
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call mpi_initialized(flag, ierr)
        if(.not. flag) call mpi_init(ierr)
#endif
    end subroutine

    subroutine swals_mpi_finalize()
        !! Calls mpi_finalize(ierr) if COARRAY_USE_MPI_FOR_INTENSIVE_COMMS is defined. Otherwise do nothing
        integer :: ierr
        logical :: flag
#if defined(COARRAY_USE_MPI_FOR_INTENSIVE_COMMS)
        call mpi_finalized(flag, ierr)
        if (.not. flag) call mpi_finalize(ierr)
#endif
    end subroutine



#ifdef COARRAY_PROVIDE_CO_ROUTINES

    !
    ! Co-broadcast from a source image
    !
    subroutine co_broadcast_int32(var, source_image)
        integer(int32), intent(inout) :: var
        integer(int32), intent(in) :: source_image

        integer :: ierr

        call MPI_Bcast(var, 1, MPI_INTEGER4, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_int32_rank1(var, source_image)
        integer(int32), intent(inout) :: var(:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_INTEGER4, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_int32_rank2(var, source_image)
        integer(int32), intent(inout) :: var(:,:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_INTEGER4, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_float(var, source_image)
        real(c_float), intent(inout) :: var
        integer(int32), intent(in) :: source_image

        integer :: ierr

        call MPI_Bcast(var, 1, MPI_REAL, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_float_rank1(var, source_image)
        real(c_float), intent(inout) :: var(:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_REAL, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_float_rank2(var, source_image)
        real(c_float), intent(inout) :: var(:,:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_REAL, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_double(var, source_image)
        real(c_double), intent(inout) :: var
        integer(int32), intent(in) :: source_image

        integer :: ierr

        call MPI_Bcast(var, 1, MPI_DOUBLE_PRECISION, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_double_rank1(var, source_image)
        real(c_double), intent(inout) :: var(:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_DOUBLE_PRECISION, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    subroutine co_broadcast_c_double_rank2(var, source_image)
        real(c_double), intent(inout) :: var(:,:)
        integer(int32), intent(in) :: source_image

        integer :: ierr, var_size

        var_size = size(var)
        call MPI_Bcast(var, var_size, MPI_DOUBLE_PRECISION, source_image-1, MPI_COMM_WORLD, ierr)

    end subroutine

    !
    ! co-sum
    !

    subroutine co_sum_int32(var)
        integer(int32), intent(inout) :: var

        integer :: ierr, size_var

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_sum_real32(var)
        real(real32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_sum_real64(var)
        real(real64), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_sum_real128(var)
        real(real128), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL16, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    !
    ! co-max
    !

    subroutine co_max_int32(var)
        integer(int32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_max_real32(var)
        real(real32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_max_real64(var)
        real(real64), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_max_real128(var)
        real(real128), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL16, MPI_MAX, MPI_COMM_WORLD, ierr)
    end subroutine

    !
    ! co-min
    !

    subroutine co_min_int32(var)
        integer(int32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_INTEGER4, MPI_MIN, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_min_real32(var)
        real(real32), intent(inout) :: var

        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_min_real64(var)
        real(real64), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_min_real128(var)
        real(real128), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL16, MPI_MIN, MPI_COMM_WORLD, ierr)
    end subroutine

#endif

end module
