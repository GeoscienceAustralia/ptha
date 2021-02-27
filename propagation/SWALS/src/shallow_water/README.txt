Here the main codes are:

    - global_mod.f90 defines various "globally relevant" constants (e.g. default real precision, gravity, pi, etc).

    - multidomain_mod.f90 defines the multidomain_type, containing multiple nested domains

    - domain_mod.f90 defines the domain_type (a single grid on which we solve the shallow water equations in various ways) and implements the various numerical solvers. The domain_mod.f90 also uses extrapolation_limiting_mod.f90, cliffs_tolkova_mod.f90, and boundary_mod.f90. In addition, all the files starting with domain_mod_*.f90 are included in some way, to facilitate generation of efficient code and prevent those details from overwhelming domain_mod.f90. For instance, at the time of writing the domain%compute_fluxes routine include various logical variables that lead to "if" statements in the inner loop - and the code runs much faster if we make those variables constants (parameters in Fortran). This can be achieved by repeating code using "#include" statements. In future if Fortran develops generic programming facilities, our use of "#include" here could be simplified.

    - timestepping_metadata_mod.f90 defines several key constants for the timestepping methods -- for instance, their allowed CFL number, and the "halo thickness" that they require to evolve one step. It is easiest to have them all defined in one place.

    - nested_grid_comms_mod.f90 defines a couple of classes that are useful for communicating between nested domains.

    - point_gauge_mod.f90 defines a class to take care of writing gauges (i.e. point time-series) to a netcdf file.

    - spherical_mod.f90 defines helps with working in spherical coordinates.

    - forcing_mod.f90 defines a type to help with applying forcings (e.g. time-varying free-surface perturbation; rainfall; etc) over subsets of a domain.
