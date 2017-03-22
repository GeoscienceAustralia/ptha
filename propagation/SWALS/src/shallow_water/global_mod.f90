module global_mod

!use iso_fortran_env, only: REAL32, INT16, INT32
use iso_c_binding, only: C_FLOAT, C_INT, C_DOUBLE, C_LONG

implicit none

! Default character length, real / integer precision
integer(8), parameter:: charlen = 1024, ip = C_INT !dp = C_FLOAT, ip = C_INT

! If -DREALFLOAT is passed to the compiler, then reals are single precision, otherwise
! we use double
#ifdef REALFLOAT
integer(8), parameter:: dp = C_FLOAT
#else
integer(8), parameter:: dp = C_DOUBLE
#endif
integer(8), parameter:: output_precision = C_FLOAT

! Physical constants
real(dp), parameter:: gravity = 9.8_dp ! m/s**2
real(dp), parameter:: advection_beta = 1.0_dp  ! Used to rescale advective terms
real(dp), parameter:: radius_earth = 6371000.0_dp ! 6378137.0_dp ! 6371000.0_dp ! Radius of the earth

! Wetting and drying
real(dp), parameter:: minimum_allowed_depth = 1.0e-05_dp
! Walls (e.g. reflective boundary) are assigned this elevation
real(dp), parameter:: wall_elevation = 1.0e+6_dp

! Timestepping
real(dp), parameter:: cfl = 0.9_dp
real(dp), parameter:: maximum_timestep = 1.0e+20_dp
character(len=charlen), parameter:: default_timestepping_method = 'euler'

! Spatial extrapolation
real(dp), parameter:: extrapolation_theta = 1.0_dp ! 0 for first order, 1 for second order

! Output folder
character(len=charlen), parameter:: default_output_folder = 'OUTPUTS'

! pi
real(dp), parameter:: pi = acos(-1.0_dp)

end module global_mod
