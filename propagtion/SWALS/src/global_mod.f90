MODULE global_mod

!USE iso_fortran_env, only: REAL32, INT16, INT32
USE iso_c_binding, only: C_FLOAT, C_INT, C_DOUBLE, C_LONG

IMPLICIT NONE

! Default character length, real / integer precision
INTEGER(8), PARAMETER:: charlen = 1024, ip = C_INT !dp = C_FLOAT, ip = C_INT

! If -DREALFLOAT is passed to the compiler, then reals are single precision, otherwise
! we use double
#ifdef REALFLOAT
INTEGER(8), PARAMETER:: dp = C_FLOAT
#else
INTEGER(8), PARAMETER:: dp = C_DOUBLE
#endif
INTEGER(8), PARAMETER:: output_precision = C_FLOAT

! I sometimes need work arrays, which I by default declare to this length,
! rather than having magic numbers everywhere
! INTEGER(ip), PARAMETER, PUBLIC:: worklen = 100

! Physical constants
REAL(dp), PARAMETER:: gravity = 9.8_dp ! m/s**2
REAL(dp), PARAMETER:: advection_beta = 1.0_dp  ! Used to rescale advective terms
REAL(dp), PARAMETER:: radius_earth = 6371000.0_dp ! 6378137.0_dp ! 6371000.0_dp ! Radius of the earth

! Wetting and drying
REAL(dp), PARAMETER:: minimum_allowed_depth = 1.0e-05_dp
! Walls (e.g. reflective boundary) are assigned this elevation
REAL(dp), PARAMETER:: wall_elevation = 1.0e+6_dp

! Timestepping
REAL(dp), PARAMETER:: cfl = 0.9_dp
REAL(dp), PARAMETER:: maximum_timestep = 1.0e+20_dp
CHARACTER(len=charlen), PARAMETER:: default_timestepping_method = 'euler'

! Spatial extrapolation
REAL(dp), PARAMETER:: extrapolation_theta = 1.0_dp ! 0 for first order, 1 for second order

! Output folder
CHARACTER(len=charlen), PARAMETER:: default_output_folder = 'OUTPUTS'

! pi
REAL(dp), PARAMETER:: pi = acos(-1.0_dp)

END MODULE global_mod
