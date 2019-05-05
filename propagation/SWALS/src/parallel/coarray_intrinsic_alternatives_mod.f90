!
! Module to allow use of co_sum and co_broadcast, in limited situations, even
! if the compiler doesn't have them (e.g. currently true of ifort).
!
! The point is to work around ifort's lack of these routines
!
! This will only do something if compiled with -DCOARRAY_PROVIDE_CO_ROUTINES
!
module coarray_intrinsic_alternatives

    use iso_fortran_env, only: int32, real32, real64
#ifdef COARRAY_PROVIDE_CO_ROUTINES
    use mpi

    implicit none

    private
    public :: co_broadcast, co_sum, co_max, co_min

    ! Can add more kinds/ranks as required
    interface co_broadcast
        module procedure co_broadcast_int32, co_broadcast_int32_rank1, co_broadcast_int32_rank2
    end interface

    interface co_sum
        module procedure co_sum_int32, co_sum_real32, co_sum_real64
    end interface

    interface co_max
        module procedure co_max_int32, co_max_real32, co_max_real64
    end interface

    interface co_min
        module procedure co_min_int32, co_min_real32, co_min_real64
    end interface

    contains

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

    !
    ! co-sum
    !

    subroutine co_sum_int32(var)
        integer(int32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_sum_real32(var)
        real(real32), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine

    subroutine co_sum_real64(var)
        real(real64), intent(inout) :: var

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
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
#endif
end module
