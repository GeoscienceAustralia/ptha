module global_mod
!! Defines globally useful parameters (e.g. real and integer precisions, physical constants, defaults)

use iso_c_binding, only: C_FLOAT, C_INT, C_DOUBLE, C_LONG, C_SIZEOF, C_LONG_DOUBLE, C_LONG_LONG
use iso_fortran_env, only: REAL128 

implicit none

integer(C_INT), parameter:: charlen = 1024 !! Default character length
integer(C_INT), parameter:: ip = C_INT  !! Default integer precision

integer(C_INT), parameter :: long_long_ip = C_LONG_LONG !! Very long integer precision

#ifdef REALFLOAT
integer(ip), parameter:: dp = C_FLOAT !! Default real precision if -DREALFLOAT is passed to the compiler. Sometimes OK, but not generally recommended.
#else
integer(ip), parameter:: dp = C_DOUBLE !! Default real precision unless -DREALFLOAT is passed to the compiler. Appropriate in most cases. 
#endif
!! If -DREALFLOAT is passed to the compiler, then most reals are single precision, otherwise
!! use double.

integer(ip), parameter:: output_precision = C_FLOAT !! Store real output at this precision in netcdf files.

! In some parts of the code we need to be sure to use double precision, even if
! 'dp' is single precision. For example, this is needed to get reasonable mass
! conservation tracking in some models (e.g. where integrating the volume 
! involves subtracting the stage from the elevation, where these differ by several
! km). The following constants are used for that purpose
integer(ip), parameter:: force_double = C_DOUBLE !! Double precision irrespective of whether -DREALFLOAT was passed to compiler.
#ifdef PGI_COMPILER
integer(ip), parameter:: force_long_double = C_DOUBLE !C_LONG_DOUBLE !! Long double precision
#else
integer(ip), parameter:: force_long_double = REAL128 !C_LONG_DOUBLE !! Long double precision
#endif

! Physical constants
real(dp), parameter:: gravity = 9.8_dp !! Gravitational acceleration in m/s**2
real(dp), parameter:: advection_beta = 1.0_dp  !! Used to rescale momentum advection terms, e.g. to account for different assumed vertical profiles of horizontal velocity when integrating the shallow water equations
real(dp), parameter:: radius_earth = 6371000.0_dp !! Radius of the earth ! 6378137.0_dp ! 6371000.0_dp 

! Wetting and drying
real(dp), parameter:: minimum_allowed_depth = 1.0e-05_dp !! Depth at which velocities are zeroed. 
!! For the 'cliffs' solver, the minimum_allowed_depth needs problem-specific tuning so is set with
!! domain%cliffs_minimum_allowed_depth

real(dp), parameter:: wall_elevation = 1.0e+6_dp !! Walls (e.g. reflective boundary) may be given this elevation

! Timestepping
!real(dp), parameter:: cfl = 0.9_dp ! FIXME: Is this still used? Defaults now set in domain%allocate_quantities
real(dp), parameter:: maximum_timestep = 1.0e+20_dp !! Default upper limit on the time-step. 
character(len=charlen), parameter:: default_nonlinear_timestepping_method = 'rk2' !! Used in test suite
character(len=charlen), parameter:: default_linear_timestepping_method = 'linear' !! Used in test suite

logical, parameter :: send_boundary_flux_data = .true. !! Turn on/of sending of boundary flux data and flux correction (for nesting)

character(len=charlen), parameter:: default_output_folder = 'OUTPUTS' !! Default base location for outputs

real(dp), parameter:: pi = acos(-1.0_dp) !! pi

integer(ip), parameter :: real_bytes = (storage_size(1.0_dp)/8_ip) !! Number of bytes in a real(dp)
integer(ip), parameter :: integer_bytes = (storage_size(1_ip)/8_ip) !! Number of bytes in integer(ip)
integer(ip), parameter :: force_double_bytes = (storage_size(1.0_force_double)/8_ip) !! Number of bytes in real(force_double)

end module global_mod
