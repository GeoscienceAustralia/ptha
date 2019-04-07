module global_mod

use iso_c_binding, only: C_FLOAT, C_INT, C_DOUBLE, C_LONG, C_SIZEOF
use iso_fortran_env, only: REAL128, REAL32 !, INT32, INT64

implicit none

! Default character length, real / integer precision
integer(C_INT), parameter:: charlen = 1024, ip = C_INT !dp = C_FLOAT, ip = C_INT

! If -DREALFLOAT is passed to the compiler, then reals are single precision, otherwise
! we use double
#ifdef REALFLOAT
integer(ip), parameter:: dp = C_FLOAT
#else
integer(ip), parameter:: dp = C_DOUBLE
#endif
integer(ip), parameter:: output_precision = C_FLOAT
! In some parts of the code we need to be sure to use double precision, even if
! 'dp' is single precision. For example, this is needed to get reasonable mass
! conservation tracking in some models (e.g. where integrating the volume 
! involves subtracing the stage from the elevation, where these differ by several
! km). The following constant is used for that purpose
integer(ip), parameter:: force_double = C_DOUBLE

! Physical constants
real(dp), parameter:: gravity = 9.8_dp ! m/s**2
real(dp), parameter:: advection_beta = 1.0_dp  ! Used to rescale advective terms, e.g. to account for different assumed vertical profiles of horizontal velocity when integrating the shallow water equations
real(dp), parameter:: radius_earth = 6371000.0_dp ! 6378137.0_dp ! 6371000.0_dp ! Radius of the earth

! Wetting and drying
real(dp), parameter:: minimum_allowed_depth = 1.0e-05_dp
! Walls (e.g. reflective boundary) are assigned this elevation
real(dp), parameter:: wall_elevation = 1.0e+6_dp

! Timestepping
real(dp), parameter:: cfl = 0.9_dp
real(dp), parameter:: maximum_timestep = 1.0e+20_dp
character(len=charlen), parameter:: default_timestepping_method = 'euler'

! Turn on/of sending of boundary flux data and flux correction (for nesting)
logical, parameter :: send_boundary_flux_data = .true.

! Output folder
character(len=charlen), parameter:: default_output_folder = 'OUTPUTS'

! pi
real(dp), parameter:: pi = acos(-1.0_dp)

integer(ip), parameter :: real_bytes = c_sizeof(1.0_dp) 
integer(ip), parameter :: integer_bytes = c_sizeof(1_ip) 
integer(ip), parameter :: force_double_bytes = c_sizeof(1.0_force_double) 

end module global_mod
