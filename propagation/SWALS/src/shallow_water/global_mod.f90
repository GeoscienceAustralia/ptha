module global_mod
!! Defines globally useful parameters (e.g. real and integer precisions, physical constants, defaults)

use iso_c_binding, only: C_FLOAT, C_INT, C_DOUBLE, C_LONG, C_SIZEOF, C_LONG_DOUBLE, C_LONG_LONG
use iso_fortran_env, only: REAL128

implicit none


! Default character length and precisions for reals, integers, and real outputs
integer(C_INT), parameter:: charlen = 1024 !! Default character length
integer(C_INT), parameter:: ip = C_INT  !! Default integer precision
integer(C_INT), parameter :: long_long_ip = C_LONG_LONG !! Very long integer precision

! If -DREALFLOAT is passed to the compiler then most reals are single precision, otherwise double
#ifdef REALFLOAT
integer(ip), parameter:: dp = C_FLOAT !! Default real precision if -DREALFLOAT is passed to the compiler. Sometimes OK, but not generally recommended.
#else
integer(ip), parameter:: dp = C_DOUBLE !! Default real precision unless -DREALFLOAT is passed to the compiler. Appropriate in most cases.
#endif

integer(ip), parameter:: output_precision = C_FLOAT !! Store real output at this precision in netcdf files.

! In some parts of the code we need to be sure to use double precision, even if
! 'dp' is single precision. For example, this is needed to get reasonable mass
! conservation tracking in some models (e.g. where integrating the volume
! involves subtracting the stage from the elevation, where these differ by several
! km). The following constants are used for that purpose
integer(ip), parameter:: force_double = C_DOUBLE !! Double precision irrespective of whether -DREALFLOAT was passed to compiler.
#ifdef PGI_COMPILER
integer(ip), parameter:: force_long_double = C_DOUBLE !! If compiled with -DPGI_COMPILER, deal with lack of support for long double precision
#else
integer(ip), parameter:: force_long_double = REAL128 !! Long double precision (unless compiled with -DPGI_COMPILER)
#endif

! Physical constants
real(dp), parameter:: gravity = 9.8_dp !! Gravitational acceleration in m/s**2
real(dp), parameter:: advection_beta = 1.0_dp  !! Used to rescale momentum advection terms, e.g. to account for different assumed vertical profiles of horizontal velocity when integrating the shallow water equations
real(dp), parameter:: radius_earth = 6371000.0_dp !! Radius of the earth

! Wetting and drying
real(dp), parameter:: minimum_allowed_depth = 1.0e-05_dp !! Depth at which velocities are zeroed. For the 'cliffs' solver, the minimum_allowed_depth needs problem-specific tuning so is set with domain%cliffs_minimum_allowed_depth

real(dp), parameter:: wall_elevation = 1.0e+6_dp !! Walls (e.g. reflective boundary) may be given this elevation

! Timestepping
real(dp), parameter:: maximum_timestep = 1.0e+20_dp !! Default upper limit on the time-step.
character(len=charlen), parameter:: default_nonlinear_timestepping_method = 'rk2' !! Used in test suite
character(len=charlen), parameter:: default_linear_timestepping_method = 'linear' !! Used in test suite

! Flux correction
logical, parameter :: send_boundary_flux_data = .true.  !! Turn on/off sending boundary flux data. Details of applied flux correction vary depending on the algorithm, see  timestepping_metadata_mod.f90

! Output dir default
character(len=charlen), parameter:: default_output_folder = 'OUTPUTS' !! Default base directory for outputs

! pi
real(dp), parameter:: pi = acos(-1.0_dp) !! pi

! Various constants for computing memory
integer(ip), parameter :: real_bytes = (storage_size(1.0_dp)/8_ip) !! Number of bytes in a real(dp)
integer(ip), parameter :: integer_bytes = (storage_size(1_ip)/8_ip) !! Number of bytes in integer(ip)
integer(ip), parameter :: force_double_bytes = (storage_size(1.0_force_double)/8_ip) !! Number of bytes in real(force_double)

end module global_mod
