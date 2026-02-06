project: SWALS
output_dir: ./doc/
project_github: https://github.com/GeoscienceAustralia/ptha/tree/master/propagation/SWALS/
summary: Shallow WAter Like Solvers (SWALS) computes solutions to several variants of the 2D shallow-water equations (linear/nonlinear) in cartesian and spherical coordinates, on domains represented as a connected set of structured grids. 
author: G Davies
graph: true
src_dir: ./src/raster/
src_dir: ./src/util/
src_dir: ./src/parallel/
exclude_dir: ./src/parallel/tmp/
src_dir: ./src/shallow_water/
exclude_dir: ./src/shallow_water/tmp/
src_dir: ./tests/unit_tests/
exclude_dir: ./tests/unit_tests/OUTPUTS/
src_dir: ./tests/parallel_tests/
exclude_dir: ./tests/parallel_tests/OUTPUTS/
exclude_dir: ./examples
exclude: point2point_include_recv_p2p.f90
exclude: point2point_include_send_p2p.f90
exclude: domain_mod_compute_fluxes_alternatives_include.f90
exclude: domain_mod_compute_fluxes_DE1_inner_include.f90  
exclude: domain_mod_routines_vectorized_include.f90
exclude: domain_mod_compute_fluxes_EEC_include.f90 
exclude: domain_mod_compute_fluxes_EEC_inner_include.f90 
exclude: domain_mod_leapfrog_nonlinear_include.f90 
exclude: domain_mod_linear_solver_include.f90
exclude: domain_mod_leapfrog_solver_friction_include.f90
exclude: domain_mod_update_U_DE1_alternatives_include.f90
exclude: domain_mod_update_U_DE1_inner_include.f90
exclude: domain_mod_update_U_DE1_inner_include_openacc.f90
exclude: domain_mod_timestepping_alternatives_include.f90
exclude: domain_mod_with_tasks.f90
exclude: sort_index_template.f90
exclude: sort_index_template2.f90

Shallow WAter Like Solvers (SWALS) computes solutions to several variants of the 2D shallow-water equations (linear/nonlinear) in cartesian and spherical coordinates, on domains represented as a connected set of structured grids. A number of different numerical methods are implemented, suitable for a range of flow regimes. This includes classical leapfrog schemes, shock-capturing finite volume schemes, and the [CLIFFS](https://github.com/Delta-function/cliffs-src) solver developed by Elena Tolkova (which is similar to MOST but uses a different wetting and drying scheme). 

Two-way nesting allows for the use of higher-resolution domains in specified areas. Where domains overlap, the highest-resolution domain is the "priority domain" meaning it contains the solution. Domains communicate with each other via halo exchange, and only receive data from "priority-domain" regions. Different domains may use different numerical solvers, and take different timesteps. For some solvers, flux correction is used to enforce the conservation of mass and advected momentum through nested domain interfaces. 

Parallel computation (shared and distributed memory CPU) is supported with a mixture of MPI (or Fortran coarrays) and/or openmp. Static load balancing can be used to improve the efficiency of large parallel jobs. The code includes a unit test suite, a parallel unit test suite, and a validation test suite (the latter focussing on tsunami type problems).
