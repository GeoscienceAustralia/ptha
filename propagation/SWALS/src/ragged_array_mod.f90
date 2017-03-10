MODULE ragged_array_mod
    !
    ! Make a type to hold 'ragged' arrays
    ! (i.e. an allocatable array where each element is an allocatable
    ! array)
    ! 

    USE global_mod, only: ip, dp
    IMPLICIT NONE

    PUBLIC

    ! Integer case

    TYPE allocatable_array_ip_type
        INTEGER(ip), ALLOCATABLE:: i1(:)
    END TYPE

    TYPE ragged_array_2d_ip_type
        TYPE(allocatable_array_ip_type), ALLOCATABLE:: i2(:)
    END TYPE
    
    TYPE ragged_array_3d_ip_type
        TYPE(ragged_array_2d_ip_type), ALLOCATABLE:: i3(:)
    END TYPE

    ! Real case

    TYPE allocatable_array_dp_type
        REAL(dp), ALLOCATABLE:: i1(:)
    END TYPE

    TYPE ragged_array_2d_dp_type
        TYPE(allocatable_array_dp_type), ALLOCATABLE:: i2(:)
    END TYPE
    
    TYPE ragged_array_3d_dp_type
        TYPE(ragged_array_2d_dp_type), ALLOCATABLE:: i3(:)
    END TYPE

END MODULE
