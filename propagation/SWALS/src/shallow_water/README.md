# Main codes for solving variants of the shallow water equations on a multidomain.
----------------------------------------------------------------------------------

The key codes in this folder are:

* [global_mod.f90](global_mod.f90) defines various "globally relevant" constants (e.g. default real precision, gravity, pi, etc).

* [multidomain_mod.f90](multidomain_mod.f90) defines the multidomain_type, containing multiple nested domains

* [domain_mod.f90](domain_mod.f90) defines the domain_type (a single grid on which we solve the shallow water equations in various ways) and implements the various numerical solvers. The domain_mod.f90 also uses:
  * [extrapolation_limiting_mod.f90](extrapolation_limiting_mod.f90) for extrapolation to edges in the finite volume solvers.
  * [cliffs_tolkova_mod.f90](cliffs_tolkova_mod.f90) for the CLIFFS solver
  * All the files starting with `domain_mod_*.f90` are included in some way, to facilitate generation of efficient code and prevent those details from overwhelming [domain_mod.f90](domain_mod.f90). For instance, at the time of writing the domain%compute_fluxes routine include various logical variables that lead to "if" statements in the inner loop - and the code runs much faster if we make those variables constants (parameters in Fortran). This can be achieved by repeating code using "#include" statements. 

* [boundary_mod.f90](boundary_mod.f90) provides some subroutines that implement boundary conditions. Users can provide their own subroutine (`domain%boundary_subroutine => my_boundary_subroutine`). 

* [timestepping_metadata_mod.f90](timestepping_metadata_mod.f90) defines several key constants for each of the timestepping methods -- for instance, their allowed CFL number, and the "halo thickness" that they require to evolve one step. 

* [linear_dispersive_mod.f90](linear_dispersive_mod.f90) to include dispersive terms in some solvers.

* [nested_grid_comms_mod.f90](nested_grid_comms_mod.f90) defines a couple of classes that are useful for communicating between nested domains.

* [point_gauge_mod.f90](point_gauge_mod.f90) defines a class to take care of writing gauges (i.e. point time-series) to a netcdf file.

* [spherical_mod.f90](spherical_mod.f90) has routines that help with spherical coordinates.

* [forcing_mod.f90](forcing_mod.f90) defines a type that is helpful in applying forcings over subsets of a domain (such as a time-varying free-surface perturbation; rainfall; etc). In general users can provide forcings by setting `domain%forcing_subroutine` (which can refer to data in `domain%forcing_context`). To use multiple forcing routines, the latter variables can be set for the first forcing term followed by `call domain%store_forcing`, then set for the second forcing term followed by `call domain%store_forcing`, and so on.
