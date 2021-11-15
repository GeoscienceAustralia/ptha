module ragged_array_mod
    !! Make a type to hold 'ragged' arrays
    !! (i.e. an allocatable array where each element is an allocatable
    !! array)

    use global_mod, only: ip, dp, force_double
    implicit none

    public

    ! integer case

    type allocatable_array_ip_type
        integer(ip), allocatable:: i1(:)
    end type

    type ragged_array_2d_ip_type
        type(allocatable_array_ip_type), allocatable:: i2(:)
    end type

    type ragged_array_3d_ip_type
        type(ragged_array_2d_ip_type), allocatable:: i3(:)
    end type

    ! real case

    type allocatable_array_dp_type
        real(dp), allocatable:: i1(:)
    end type

    type ragged_array_2d_dp_type
        type(allocatable_array_dp_type), allocatable:: i2(:)
    end type

    type ragged_array_3d_dp_type
        type(ragged_array_2d_dp_type), allocatable:: i3(:)
    end type

    !
    ! real case that is always double precision
    !
    type allocatable_array_force_double_type
        real(force_double), allocatable:: i1(:)
    end type

    type ragged_array_2d_force_double_type
        type(allocatable_array_force_double_type), allocatable:: i2(:)
    end type

    type ragged_array_3d_force_double_type
        type(ragged_array_2d_force_double_type), allocatable:: i3(:)
    end type

end module
