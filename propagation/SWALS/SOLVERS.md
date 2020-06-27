SWALS solvers
--------------

SWALS contains a range of solvers for variants of the 2D shallow water equations. For models that involve multiple domains (i.e. using a `multidomain` object, typically denoted `md`), one can use different solvers on different domains, and they may take different time-steps. 

In general the different solvers have different strengths and limitations, which are discussed below. For instance some of our solvers are good for global-scale tsunami propagation, but totally inappropriate for nearshore inundation simulation (or the converse). Furthermore not all solvers are equally well supported by the boundary conditions and the nesting formulation. For instance: 

* Some solvers are more prone to long-time nesting instabilities if they receive data from a coarser domain; others much less so. While this is partly related to the numerical schemes, it also reflect the degree-of-effort to make nesting well-behaved in different cases.
* Some solvers are not well supported in complex boundary conditions (e.g. a radiation condition with a stage forcing). Again this just reflects the developer's priorities.

Generally speaking the SWALS code development has focussed most on:

* A range of shock-capturing finite volume schemes. These are mainly used for flows with significant nonlinearity, medium to high Froude-numbers, etc. They are quite robust for this purpose, but tend to be more diffusive than the leapfrog schemes (and thus less well suited to global scale tsunami propagation). Some variants of these schemes have better performance at low Froude-numbers than others.

* A range of leapfrog schemes. These solve the linear shallow water equations; variants on the latter with some additional nonlinear terms (e.g. nonlinear pressure gradient terms; manning friction); and there is a variant with all nonlinear terms. The leapfrog schemes have mostly been used for flows domainated by low Froude-numbers (especially global-scale tsunami propagation simulation). The fully nonlinear leapfrog scheme can perform well in some more strongly nonlinear problems (e.g. inundation), but we haven't focussed on making it extremely robust for such purposes (and generally prefer the finite-volume schemes such as `rk2`).

We also include the CLIFFS solver developed by Elena Tolkova, but haven't worked extensively with that.

# Setting the solver type.

A SWALS `multidomain` object (denoted `md`) contains an allocatable array of domains (`md%domains(:)`). The is allocated at the beginning of an application code (e.g. `allocate(md%domains(10))` will allow 10 domains). The solver for each domain is set by specifying its timestepping method, e.g. `md%domains(3)%timestepping_method = "linear"` would make the 3rd domain using the linear leapfrog solver. Alternatively this information can be provided to the subroutine `md%domains(3)%match_geometry_to_parent` which would be used to nest the 3rd domain inside another domain; this commonly done in the example programs (such as [here](./examples/nthmp/BP09/BP09.f90)). 

The allowed `timestepping_method` values are specified in the file [timestepping\_metadata\_mod.f90](./src/shallow_water/timestepping_metadata_mod.f90), which sets solver specific variables such as the halo thickness required to advance one time-step. For some solvers, additional aspects of their behavior are controlled by other variables inside the domain. 

# Details on the solvers

Although many variants of numerical methods are supported in SWALS (detailed below), the validation test-suite most often uses `rk2` for nonlinear problems, and `linear` for linear problems. These defaults can be changed by setting `default_nonlinear_timestepping_method` and `default_linear_timestepping_method` in [global\_mod.f90](./src/shallow_water/global_mod.f90). Many test problems use the latter variables to specify the solver type; but beware not all problems do this (because we deliberately exercise non-default solvers for some problems). 

This mean you should not assume that all schemes are equally tested -- they definitely are not. However it is generally straightforward to adjust validation test codes to use other schemes, either by changing the `default_XX` variables mentioned previously, or by directly changing the test-suite source-code in case other solvers are used.

## The finite-volume solvers

SWALS has a number of classical shock-capturing finite-volume schemes with accuracy up to second order in space and time. These work well for flows with moderate or high Froude-numbers, but may be too dissipative at very low Froude-numbers (e.g. they are not well suited to global-scale tsunami propagation, but can work well the nearshore and inundation simulation). 

In general the finite-volume schemes in SWALS have good stability when used in conjunction with nesting. Occasionally these solvers can be subject to artifical vortices at coarse-to-fine nesting regions, especially where the elevation has rapid variation. This is not particularly common though, and can usually be solved by moving the domain boundary to an area with weaker elevation variation.

The finite-volume schemes in SWALS use a range of different timestepping schemes, flux functions, and spatial extrapolations to derive the flow state at edges. Each scheme has a default limiter coefficient, flux function, and friction model, but these can be changed by the user prior to calling `md%setup`:

The finite-volume solvers are supported by setting `md%domains(j)%timestepping_method` to:

* `rk2`. This is a good default choice, and is heavily exercised in the validation test suite. It has second order accuracy in space and time. It uses 2nd-order Runge-kutta timestepping (i.e. advancing two first-order timesteps and taking the average). By default it uses second order spatial extrapolation with a default limiter coefficient `md%domains(j)%theta = 1.6`.

* `midpoint`. This is like `rk2` but uses a midpoint timestepping method. Results are generally very similar to `rk2`.

* `euler`. This is a variant on the above schemes with first order (explicit Euler) timestepping and lower spatial accuracy. It is less accurate and more dissipative than the aforementioned schemes, but can be computationally cheaper. For stability we find it requires a more dissipative spatial extrapolation ( default limiter coefficient `md%domains(j)%theta = 0.9`, which is less than second order accurate in space). It is not recommended for general usage, although may be fine for some problems. Regardless the underlying code is an important building block for the more advanced schemes.

* `rk2n`. This is an unusual timestepping scheme with second order accuracy in space and time. Internally it uses 5 first-order timesteps to advance 4 timesteps, while retaining second order accuracy in time. In principle it can be faster than the `rk2` and `midpoint` schemes, which internally use 2 first-order timesteps to advance 1 timestep. However the scheme is less stable, and not recommended for general usage. 

For all schemes above the spatial extrapolation of the flow-state to cell edges employs the water-surface, depth, and velocity (x and y components), using the now classical approach of Audusse (*Audusse, E.; Bouchut, F.; Bristeau, M.-O.; Klein, R. & Perthame, B. A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows SIAM Journal of Scientific Computing, 2004, 25, 2050-2065*). The extrapolation employs a slope-limiter which controls the spatial order of accuracy. It is more-or-less dissipative depending on a domain-specific coefficient `theta`; the value for the `j`th domain is stored in `md%domains(j)%theta` and can be any non-negative real number: 

* `theta=0.0` corresponds to purely first-order extrapolation in space, and is quite dissipative. Higher values will make use of gradients when extrapolation edge values, and can have up-to second-order accuracy for smooth flows.

* `theta=1.0` is the smallest value that is second order accurate for sufficiently smooth flows. In practice it is still relatively dissipative (Minmod limiter). 

* `theta=2.0` is the upper threshold for which the extrapolation to edges is "total-variation-diminishing (TVD)" (here we use the term loosely to mean that `left-cell-value <= edge-value <= right-cell-value`, or `left-cell-value >= edge-value >= right-cell-value`). For some other PDEs and numerical methods, this will prevent the creation of new extrema. This is not true for the shallow water equations because they are not TVD, but the term is often used in this way to discuss limiters and the technique remains useful. 

* `theta>2.0` is even less dissipative and is not TVD (i.e. edge values can exceed neighbouring cell values). In practice this remains stable for many problems, and can be important in cases where the limiter otherwise leads to excessive dissipation (e.g. low Froude-number flows with many minima and maxima).

In all cases `theta` is tapered to zero in flows that are near a wet/dry boundary or extremely shallow -- for details on how this is done see [domain_mod.f90](./src/shallow_water/domain_mod.f90) including comments near the definition of `limiter_coef`1-3.

The friction model can be controlled by setting the variable `md%domains(j)%friction_type`. Values are:

* `"manning"`. This is Manning friction, so the array `md%domains(j)%manning_squared(:,:)` is interpreted as the (manning friction)**2

* `"chezy"`. This is Chezy friction. In this case the array `md%domains(j)%manning_squared(:,:)` is interpreted as (1/chezy_friction)**2

In addition one can set a linear friction coefficient `md%domains(j)%linear_drag_coef` which is zero by default. If non-zero this implements a non-standard linear friction model following *Fine, I. V.; Kulikov, E. A. & Cherniawsky, J. Y. Japans 2011 Tsunami: Characteristics of Wave Propagation from Observations and Numerical Modelling Pure and Applied Geophysics, Springer Science and Business Media LLC, 2012, 170, 1295-1307*. In practice you probably do not want to use this for the finite-volume schemes, which are already somewhat numerically dissipative; it is more likely to be useful to add slow friction to the leapfrog schemes.

The numerical fluxes are derived from the left and right cell edge values using an approximate Riemann solver. Various options can be controlled by setting the character string `md%domains(j)%compute_fluxes_inner_method`. Values are:

* `"DE1"`. This is similar to the DE1 solver in ANUGA [see *Davies and Roberts, 2015, Open source flood modelling with a 2D discontinuous elevation hydrodynamic model. Proceedings of MODSIM 2015*]. It uses a HLL approximate Riemann solver, with a treatment of the pressure gradient term. Variants in this approach are popular in the literature.

* `"DE1_low_fr_diffusion"`. This is the default approach. It is similar to `DE1` but reduces the diffusive component of the UH/VH fluxes based on the Froude-number, similar to *Xue-song Li and Chun-wei Gu (2013) Mechanism and Improvement of the Harten–Lax–van Leer Scheme for All-Speed Flows. Computers & Fluids Volume 86, 5, Pages 56-70*.

* `"DE1_upwind_transverse"`. This is like `DE1` but uses an upwind approach for the transverse momentum flux. This approach is quite popular in the literature, for instance it is used in the shallow water solver provided by the well-known Basilisk code.

* `"DE1_low_fr_diffusion_upwind_transverse"`. This is like `DE1_low_fr_diffusion` but uses the upwind approach for the transverse momentum flux. 

Broadly speaking the solvers with `_low_fr_diffusion` might perform better in low Froude-number flows, although often results will be similar. The solvers with `_upwind_transverse` may have less numerical dissipation in shear flows, which may be a good or bad thing, but often the results will be similar. 


## The leapfrog schemes

Leapfrog schemes are classically used for deep-ocean tsunami propagation, and are also popular for inundation modelling (although this is has not been a focus in SWALS). They have good energy conservation properties which makes them much better suited to long-distance deep ocean tsunami propagation, as compared with the finite-volume schemes discussed above. 

In general the leapfrog schemes in SWALS do not have good long-time stability when used in conjunction with nesting. They work well as the coarsest grid in a multidomain, but if used for a refined grid (which receives halo data from a coarser domain) then they often develop instability during long-time integration. This is not always a problem, but more common than we'd like. 

The leapfrog solver variants provided by SWALS are:

* `"linear"`. This solves the linear shallow water equations if `md%domains(j)%linear_solver_is_truely_linear=.true.` (default). In that case the pressure gradient term is linearized by setting the depth equal to the depth at mean-sea-level; the mean-sea level should be set by specifying `md%domains(j)%msl_linear` (which is 0.0 by default). With these settings we solve the standard linear shallow water equations. If you set `md%domains(j)%linear_solver_is_truely_linear=.false.` then the pressure gradient term is not linearized, and the equations are nonlinear (but without nonlinear advection terms or nonlinear friction).

* `"leapfrog_linear_plus_nonlinear_friction"`. This is similar to `"linear"` above with the addition of a nonlinear friction term (either Manning or Chezy as discussed below). The same parameters as discussed above should be specified to control nonlinearity in the pressure gradient term.

* `"leapfrog_nonlinear"`. This solves the nonlinear shallow water equations with a leapfrog scheme. 


For all leapfrog solvers except `"linear"`, the friction model can be controlled by setting the variable `md%domains(j)%friction_type`. Values are:

* `"manning"`. This is Manning friction, so the array `md%domains(j)%manning_squared(:,:)` is interpreted as the (manning friction)**2

* `"chezy"`. This is Chezy friction. In this case the array `md%domains(j)%manning_squared(:,:)` is interpreted as (1/chezy_friction)**2

For all solvers (including `"linear"`), one can set the linear friction coefficient `md%domains(j)%linear_drag_coef` which is 0.0 by default. If this is set to a non-zero value then it implements a non-standard linear friction model following *Fine, I. V.; Kulikov, E. A. & Cherniawsky, J. Y. Japans 2011 Tsunami: Characteristics of Wave Propagation from Observations and Numerical Modelling Pure and Applied Geophysics, Springer Science and Business Media LLC, 2012, 170, 1295-1307*. 


## The CLIFFS solver

This solver was developed by Elena Tolkova and is provided [this git repository](https://github.com/Delta-function/cliffs-src). It is similar the well-known MOST solver but uses a different wetting and drying technique. Tolkova has reported considerable benchmarking of CLIFFS; see links in the README of the aforementioned repository, and also *Tolkova, E. Land--Water Boundary Treatment for a Tsunami Model With Dimensional Splitting Pure and Applied Geophysics, 2014, 171, 2289-2314*.

The CLIFFS solver was added to SWALS by extracting the inner computational routine, with only minor modifications to get it to compile in our framework. CLIFFS solves the shallow water equations using characteristic variables; in SWALS we transfer to/from these variables every time-step because our nesting/parallel framework requires the solver to work with the stage, UH, VH and bed elevation. For this reason the implementation in SWALS might be slower than a pure CLIFFS implementation (on the other hand CLIFFS is quite efficient anyway, and one can employ mixed openmp/mpi in SWALS which may permit larger parallel runs).

In general CLIFFS needs a wetting and drying threshold to be tuned to the problem at hand. This can be done by setting the variable `md%domains(j)%cliffs_minimum_allowed_depth`; typical values might be 0.001 (1 mm) for lab experiments, around 0.05 - 0.1 for high-resolution field cases, and 1 m or more for global-propagation problems. Tuning may be required for stability on any particular problem.

The CLIFFS solver requires sufficiently smooth bathymetry for stability; see the analysis by Tolkova in the CLIFFS user manual. A bathymetry-smoothing routine is provided by Tolkova for this purpose and has been integrated into SWALS; one must generally use this to avoid instabilities. To apply the bathymetry smoothing, after the elevation has been set in `md%domains(j)%U(:,:,ELV)` you would call the bathymetry smoothing routine like `call md%domains(j)%smooth_elevation(smooth_method='cliffs')`.

Although we have solved various know problems accurately using the CLIFFS solver in SWALS, in general limited work has been done. The SWALS boundary conditions and two-way nesting could likely be better tailored to the CLIFFS solver, which has different numerical properties to our other solvers. If using SWALS with CLIFFS solvers in conjunction with two-way nesting or non-trivial boundary-conditions (e.g. radiation with stage forcing), one should double check that it is performing OK. 

By default our CLIFFS solver uses `md%domains(j)%friction_type = "manning"`. This is Manning friction, and can be spatially variable. It is set using the array `md%domains(j)%manning_squared(:,:)` which interpreted as the (manning friction)**2. Thus far we have not implemented Chezy friction in the CLIFFS solver (although that would not be difficult). We have implemented the non-standard linear drag model associated with `md%domains(j)%linear_drag_coef`, which is zero by default. In practice this model is unlikely to be desirable for use with CLIFFS, which is numerically dissipative in any case.

The CLIFFS solver is not mass conservative; see Tolkova papers for discussion of this especially in relation to wetting and drying. For this reason in SWALS nesting framework, no flux correction is applied between neighbouring domains if one uses CLIFFS, and mass conservation errors are expected.
