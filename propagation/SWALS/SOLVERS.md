# SWALS solvers
----------------

SWALS contains a range of solvers for variants of the 2D shallow water equations (optionally with dispersion). For models that involve multiple domains (using a `multidomain` object), one can use different solvers on different domains, and they may take different time-steps. Both Cartesian and Spherical coordinates are supported for all solver types and the choice is made when compiling the code, so all domains must use the same coordinate system. [See here for a study that uses a range of SWALS solvers](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full).

In general the different solvers have different strengths and limitations, which are discussed below. For instance, some solvers work well for global-scale tsunami propagation (at least for long-waves that satisfy the shallow water equations) but are inappropriate for nearshore inundation simulation. The converse holds for some other solvers.

Generally speaking the SWALS code development has focussed most on:

* A range of shock-capturing finite volume schemes. These are mainly used for flows with significant nonlinearity, medium to high Froude-numbers, etc. They are quite robust for this purpose, but tend to be more diffusive than the leapfrog schemes (and thus less well suited to global scale tsunami propagation). Some variants of these schemes have better performance at low Froude-numbers than others. Dispersion is supported (but more accurately with `midpoint` than other schemes).

* A range of leapfrog schemes. These solve the linear shallow water equations; variants on the latter with some additional nonlinear terms (e.g. nonlinear pressure gradient terms; Manning friction); and there is a variant with all nonlinear terms. Dispersion is also supported. Our applications of the leapfrog schemes have focussed on low Froude-number flows (especially global-scale tsunami propagation simulation). The fully nonlinear leapfrog scheme can perform well in some more strongly nonlinear problems (e.g. inundation), but it is only first order accurate in the nonlinear advection term, and so for flows with higher Froude numbers we generally prefer second order finite-volume schemes such as `rk2` (or `midpoint` if also using dispersion).

We also include the CLIFFS solver developed by Elena Tolkova (with the option of a dispersive term).

## Setting the solver type.

A SWALS `multidomain` object (denoted `md`) contains one or more structured grids (`domains`) on which we solve the shallow water equations. These domains are stored in an allocatable array (`md%domains(:)`) of type `domain_type`. This is allocated at the beginning of an application code (e.g. `allocate(md%domains(10))` will allow 10 domains). 

The solver for each domain is set by specifying its timestepping method. For example by setting `md%domains(3)%timestepping_method = "linear"`, we make the 3rd domain use the linear leapfrog solver. Many examples are included in the [validation test suite](./examples), such as [here](./examples/nthmp/BP09/BP09.f90). 

The solver types supported by SWALS are specified in the file [timestepping\_metadata\_mod.f90](./src/shallow_water/timestepping_metadata_mod.f90). This also sets solver specific constants, such as the halo thickness required to advance one time-step. Additional aspects of the solver behavior are controlled by variables in the `domain_type`, for example, `md%domains(j)%lower_left` defines the lower-left coordinate for the `j`th domain. The `domain_type` variables are well documented in the source code, see definitions under `type :: domain_type` in [domain_mod.f90](./src/shallow_water/domain_mod.f90).

The solvers support use of either Cartesian or Spherical coordinates by specifying the appropriate [compile time option](README.md#advanced-compilation). By default Cartesian coordinates are used. To use spherical coordinates we set the preprocessor varable `-DSPHERICAL` when compiling. [Many examples](./examples) are provided in the test-problems; for instance [this makefile](examples/nthmp/Hilo_Tohoku_tsunami/make_hilo) uses spherical coordinates (and unrelated code timing) by defining the variable `SWALS_PREPROCESSOR_FLAGS := -DTIMER -DSPHERICAL`. If the variable `-DSPHERICAL` were not set then Cartesian coordinates would be used.


# The shallow water equations

SWALS solves variants of the 2D shallow water equations:

$$ \frac{\partial \eta}{\partial t} + \nabla \cdot \mathbf{q} = 0 $$

$$ \frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) + g h \nabla \eta + g h \mathbf{S_f} + \mathbf{\Psi} + \mathbf{\Omega_s} + \mathbf{M_s} + \mathbf{\Upsilon} = 0 $$

Here $\eta = h + z$ is the free surface elevation (m), $t$ is time (s), $h$ is the depth (m), $z$ is the bed elevation (m), $\mathbf{q}=h\mathbf{u}$ is the 2D flux vector $(m^2/s)$, $\mathbf{u} = (u,v)$ is the 2D velocity vector (m/s), $\otimes$ is the [outer product](https://en.wikipedia.org/wiki/Outer_product), $g$ is gravity $(m/s^2)$, $\mathbf{S_f}$ is the friction slope vector, $\mathbf{\Psi}$ is a turbulent diffusion term, $\mathbf{\Omega_s}$ is the Coriolis force for spherical coordinates, $\mathbf{M_s}$ gives some additional spherical coordinate terms related to the momentum flux divergence, and $\mathbf{\Upsilon}$ is a dispersive term. 

In Cartesian coordinates $(x,y)$ are in meters. The "grad" operator $\nabla f$ applied to a scalar function $f$ gives $( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y})$. The divergence operator $\nabla \cdot \mathbf{u}$ applied to a vector $\mathbf{u} = (u,v)$ gives $( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y})$, and applied to a matrix (representing $\mathbf{u}\otimes\mathbf{q}$) it gives the vector resulting from application to each column separately. 

In Spherical coordinates $(\theta, \phi)$ are (longitude, latitude) in degrees. The "grad" operator $\nabla f$ applied to a scalar function $f$ gives $( \frac{1}{R \cos(\phi_r)} \frac{\partial f}{\partial \theta_r}, \frac{1}{R \phi_r} \frac{\partial f}{\partial \phi_r})$ where the subscript $r$ indicates angles in radians rather than degrees, and $R$ is the earth radius in meters. The divergence operator $\nabla \cdot \mathbf{u}$ applied to a vector $\mathbf{u}$ gives $( \frac{1}{R \cos(\phi_r)} (\frac{\partial u}{\partial \theta_r} + \frac{\partial \cos(\phi_r) v}{\partial \phi_r}))$, and applied to a matrix (representing $\mathbf{u}\otimes\mathbf{q}$) it gives the vector resulting from application to each column, plus an additional term written here as $\mathbf{M_s}$ to avoid tensor calculus notation. We stress that SWALS uses units of degrees (not radians) for longitude and latitude, although the latter operators are easier to write in radians.

The friction slope term by default uses Manning friction $\mathbf{S_f} = n^2\mathbf{u} |\mathbf{u}| h^{-4/3}$, with Manning coefficient $n$ that is zero by default. But it can also use Chezy friction, which can be written in a similar form using a different power of the depth $(\mathbf{S_f} = n^2\mathbf{u} |\mathbf{u}| h^{-1})$, where $n^{-2}$ is an unusual way of writing the squared Chezy coefficient. The friction model can be controlled by setting the variable `md%domains(j)%friction_type`. Values are:
* `"manning"`. This is Manning friction, so the array `md%domains(j)%manning_squared(:,:)` is interpreted as the (Manning friction)^2
* `"chezy"`. This is Chezy friction. In this case the array `md%domains(j)%manning_squared(:,:)` is interpreted as (1/Chezy friction)^2
* These can be combined with a non-standard linear friction coefficient `md%domains(j)%linear_drag_coef` which is zero by default. If non-zero this implements a linear-friction model following *Fine, I. V.; Kulikov, E. A. & Cherniawsky, J. Y. Japans 2011 Tsunami: Characteristics of Wave Propagation from Observations and Numerical Modelling Pure and Applied Geophysics, Springer Science and Business Media LLC, 2012, 170, 1295-1307*. In practice you probably do not want to use this for the finite-volume schemes, which are already somewhat numerically dissipative; it is more likely to be useful to add slow friction to the leapfrog schemes. See [here](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) for a paper that uses this friction model in SWALS.

The turbulent diffusion term $\mathbf{\Psi}$ is by default set to zero. But an rudimentary eddy viscosity formulation can also be applied *for the finite volume schemes only* by specifying `md%domains(j)%use_eddy_viscosity=.true.`. Then $\mathbf{\Psi} = -\nabla \cdot (\epsilon h (\nabla \otimes \mathbf{u}))$. Here $\epsilon$ is the eddy viscosity model and $\nabla \otimes \mathbf{u}$ is a 2x2 matrix with columns $\nabla u$, $\nabla v$ respectively. The eddy viscosity can be controlled and adjusting the values of `md%domains(j)%eddy_visc_constants(1:2)`. 
* The eddy viscosity $\epsilon$ is equal to `md%domains(j)%eddy_visc_constants(1) + md%domains(j)%eddy_visc_constants(2) * depth * shear_velocity`. 
    * The first term allows a simple constant eddy viscosity
    * The second term corresponds to a very common eddy viscosity model (e.g. *Shiono, K. & Knight, D. W. 1991. Turbulent open-channel flows with variable depth across the channel. Journal of Fluid Mechanics, 222, 617-646*) 
* In Spherical coordinates, SWALS neglects some terms that are second order (except perhaps near the poles) related to the basis vectors changing direction from site to site.
* A simple second-order-in-space explicit scheme is used to discretise the equations. The eddy viscosity $\epsilon$ is clipped if needed to prevent violation of the time-step constraint for explicit diffusion (`eddy_viscosity <= (min_grid_cell_side_length**2 / (2 * timestep))`). This approach can work well in practical cases where turbulent diffusion is not very strong. For models with high resolution and/or strong turbulent diffusion, you may need to reduce the `timestep` to prevent clipping from affecting the solution. The same approach to limiting the eddy-viscosity is used in the model SWASH (*Zijlema, M.; Stelling, G. & Smit, P. SWASH: An operational public domain code for simulating wave fields and rapidly varied flows in coastal waters Coastal Engineering, 2011, 58, 992 - 1012*).
* See [here](examples/overbank_flow/) for a test of the turbulent diffusion model in SWALS, which compares with a separate solution for steady-state overbank flow.

In Cartesian coordinates there is no Coriolis term $(\mathbf{\Omega_s} = 0)$. In Spherical coordinates the Coriolis force is $\mathbf{\Omega_s} = 2 \sin(\phi_r) \gamma_r (-vh, uh)$, where $\gamma_r = 7.292115 \times 10^{-5}$ gives the earth angular frequency in radians per second.

In Cartesian coordinates there are no additional spherical coordinate terms associated with the divergence of the momentum flux $(\mathbf{M_s} = 0)$. In Spherical coordinates the momentum flux divergence leads to $\mathbf{M_s} = \frac{\tan(\phi_r)}{R} u (-vh, uh)$. This term is associated with the changing direction of basis vectors in Spherical coordinates. Using tensor calculus notation we may write the flux divergence term $\nabla_i(u^i q^j)\mathbf{Z_j}$ where $\nabla_i$ is the covariant derivative, which includes a sum of partial derivatives and other terms due to the non-zero Christoffel symbols, while $u^i$ and $q^j$ are the contravariant components of velocity and flux, and $\mathbf{Z_j}$ is the covariant basis. Expanding and converting to velocity components in m/s leads to the extra terms represented by $\mathbf{M_s}$ here. The leads to a very small correction away from the poles and is often ignored in the tsunami literature. For examples of papers which discuss it, see [Williamson et al. 1992](https://doi.org/10.1016/S0021-9991(05)80016-6), [Titov et al. 2016](http://dx.doi.org/10.1061/(ASCE)WW.1943-5460.0000357), and [Popinet 2011](https://doi.org/10.1007/s10236-011-0438-z). 

By default the dispersive term $\mathbf{\Upsilon}=0$. However if `md%domains(j)%use_dispersion = .true.` then dispersive terms are included on domain `j`. This causes SWALS to use a different subroutine to advance the solution on all domains (whether or not they use dispersion), which needs substantially more communication and computation, places more restrictions on the allowed timestepping, and tends to reduce the stability of nesting (although factor of 3 refinement seems to work OK). The dispersive terms have the form $\mathbf{\Upsilon} = -\frac{h_0^2}{3} \nabla (\nabla \cdot \frac{\partial \mathbf{q}}{\partial t})$ which is widely used in the tsunami literature. Here $h_0$ is the depth below `md%domains(j)%msl_linear`, so when using dispersion it is critical to set the latter to a sensible value. In practice we taper the dispersive terms to zero linearly between two user-specified values of $h_0$ (below). There is a preprocessor option (`-DDISPERSIVE_PEREGRINE_IN_FLUX_FORM`) that modifies the dispersive terms to $\mathbf{\Upsilon}=-\frac{h_0^2}{2} \nabla (\nabla \cdot \frac{\partial \mathbf{q}}{\partial t}) + \frac{h_0^3}{6} \nabla (\nabla \cdot \frac{\partial \mathbf{q/h_0}}{\partial t})$ . In the linear regime the latter model is equivalent to the dispersive model of [Peregrine (1967)](https://doi.org/10.1017/S0022112067002605) although in flux form rather than the velocity form presented by Peregrine ( the difference matters in the nonlinear regime, e.g. [Filippini et al. 2015](http://dx.doi.org/10.1016/j.coastaleng.2015.02.003)). 
* The dispersive model is mainly valid in the linear regime, and should be turned off when nonlinearities become too great. SWALS allows the user to specify two depths (`domain%ds%td1` $>$ `domain%ds%td2`) so that dispersion is linearly tapered to zero as $h_0$ varies between `td1` and `td2`. If dispersion is really needed in the nearshore then you should carefully consider whether the model here is suitable.
* Energy conservation is not expected using dispersion (it is no longer a property of the underlying equations). 
* If using finite volume methods with dispersion, it is suggested to use `midpoint` which is more accurate than the other finite volume solvers because the dispersive terms are included in both the predictor and the main timestep.
* Solution method (for all numerical schemes)
    * The time-derivatives are represented as the difference between terms involving purely spatial derivatives at the current time level and the next time level. There is one discretisation for the finite-volume and cliffs solvers, and another for the leapfrog solvers.
    * To solve the resulting implicit equations, imagine splitting the multidomain global timestep $dt$ into a $N_{t}$ substeps of size $\delta t$, where $\delta t$ is the smallest timestep required by any *dispersive* domain. Thus $dt = N_{t}\delta t$ and the substeps can be indexed as $1, 2, 3, \ldots, N_{t}$. Domain $j$ has its own timestep (`alpha_j`$\delta t$), usually determined by its CFL condition. SWALS requires dispersive domains to have `alpha_j` being an integer that EXACTLY divides $N_{t}$ (e.g. if $N_{t} = 8$ then dispersive domains could have `alpha_j = 1, 2, 4, 8`). Non-dispersive domains have the same constraint, unless they take a timestep smaller than $\delta t$ (in which case we could have `alpha_j = 1./2, 1./3, 1./4, ...`).  
    * For each substep $s = 1, 2, \ldots, N_{t}$, we loop over domains
        * Skip domains where `alpha_j > 1` and `mod(s-1, alpha_j) != 0`. (These have a timestep large enough that they do not need to evolve every substep).
        * For domains that remain
            * Backup the solution for dispersive domains. Then take one shallow water step using the domain specific timestep (`alpha_j`$\delta t$), unless `alpha_j < 1` (non-dispersive domains only) in which case we take enough steps to advance time by $\delta t$. 
            * Communicate halo data, but do not immediately update the solution.
            * Implicit solve
                * First *outer iteration*
                    * Dispersive domains receive halos from domains that have advanced to the same time. (The constraints on the timesteps mean that this will include data from all domains with the same timestep. On some substeps it will also include data from domains that take a larger timestep).
                    * Dispersive domains compute the right-hand-side required for the implicit solve.
                    * Dispersive domains guess the solution the end of the substep using quadratic-in-time extrapolation from the solution each global timestep. 
                    * Dispersive domains iteratively update the solution with implicit solves applied sequentially in the x-and-y directions (tridiagonal solves). These *inner iterations* are repeated once by default (corresponding to `md%domains(j)%tridiagonal_inner_iter=2`).
                    * Dispersive domains communicate updated halo data, but do not immediately update the solution
                * Repeat the *outer iteration* once by default (`md%dispersive_outer_iterations_count = 2`), except after the first time these *outer iterations* only receive halos from domains that have advanced to the same time AND have the same grid size. Unequal grid sizes would always disagree, due to interpolation, leading to worse numerical behaviour if they are included beyond the first iteration. 
                * The default number of inner and outer iterations are adequate for problems tried so far, but can be adjusted if needed.
            * Communicate halo data between domains (both dispersive and non-dispersive) that have advanced to the same time.
        * Advance to the next substep
    * After the final substep, apply flux-correction and communicate once more. At this point all domains will have advanced one global timestep $dt = N_{t}\delta t$.
* Models using dispersion thus require much more communication and computation (about 3x compute for `midpoint`). In nested grid models, the implicit method may need halos to be thicker than SWALS provides by default, particularly on domains that only take 1 timestep for every global timestep (which are less likely to have thick halos by default). To address this one can increase the halo thickness by setting the integer `md%domains(j)%minimum_nesting_layer_thickness` (e.g. a value of 30 is used on the coarsest grids for some test problems).

Elevation smoothing is sometimes useful to manage numerical instabilities (especially for the CLIFFS scheme, or for models with complex nesting). It can be applied a specific location with the subroutine `md%domains(j)%smooth_elevation_near_point`, or along all coarse-to-fine nesting boundaries with `md%domains(j)%smooth_elevation_near_nesting_fine2coarse_boundaries`.


# Details on the solvers

Multiple numerical methods are supported in SWALS (detailed below). The validation test-suite most often uses `default_nonlinear_timestepping_method` (currently equal to `'rk2'`) for nonlinear problems, and `default_linear_timestepping_method` (currently equal to `'linear'`) for linear problems. For problems involving dispersion, `midpoint` seems to be the most accurate finite volume method and is emphasised in the dispersive test problems, although the leapfrog schemes can also work well (but tend to be less accurate if nonlinearities are strong).

With these defaults, other schemes are less heavily exercised by the validation test-suite. If you wish to better understand the performance of other schemes, the default nonlinear and linear solvers used in the test suite can be changed by setting `default_nonlinear_timestepping_method` and `default_linear_timestepping_method` in [global\_mod.f90](./src/shallow_water/global_mod.f90). This will affect many but not all test problems. A number of tests deliberately exercise non-default solvers by directly specifying `md%domains(j)%timestepping_method`, rather than using the defaults. In such cases one should directly change `md%domains(j)%timestepping_method` in the code used to setup the problem.

## The Finite-Volume solvers

SWALS has a number of classical shock-capturing finite-volume schemes with accuracy up to second order in space and time. The SWALS finite-volume solvers approximate the shallow water equations in flux conservative form (but use finite differences to treat the dispersive terms). They are shock capturing, and well suited to flows with moderate or high Froude-numbers and wetting/drying, but may be too numerically dissipative to efficiently model flow at very low Froude-numbers. For example they can work well for nearshore and inundation simulation, but are not as well suited to global-scale tsunami propagation as the leapfrog schemes (discussed below). In practice we often develop global-to-local nested grid models by using a leapfrog solver on the coarsest global grid, and finite-volume solvers on the nested grids (such as in [this paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full)). 
In general the finite-volume schemes in SWALS have good stability when used in conjunction with nesting. Occasionally these solvers can be subject to artifical vortices at coarse-to-fine nesting regions, especially where the elevation has rapid variation. This is not particularly common though, and can usually be solved by moving the domain boundary to an area with weaker elevation variation, or by slightly smoothing the elevation data in the neighbourhood of the spurious flow. 

The finite-volume schemes in SWALS use a range of different timestepping schemes, flux functions, and spatial extrapolations to derive the flow state at edges. They are essentially variants of the scheme described [here](https://www.researchgate.net/publication/289882157_Open_source_flood_simulation_with_a_2D_discontinuous-elevation_hydrodynamic_model), but on a structured grid rather than an unstructured mesh. Each scheme has a default limiter coefficient, flux function, and friction model, but these can be changed by the user prior to calling `md%setup`:

The finite-volume solvers are supported by setting `md%domains(j)%timestepping_method` to:

* `"rk2"`. This is a good default choice, and is heavily exercised in the validation test suite. It has second order accuracy in space and time. It uses 2nd-order Runge-kutta timestepping (i.e. advancing two first-order timesteps and taking the average). By default it uses second order spatial extrapolation with a default limiter coefficient `md%domains(j)%theta = 1.6`.

* `"midpoint"`. This is like `"rk2"` but uses a midpoint timestepping method. Results are generally very similar to `"rk2"`, unless using dispersion where `midpoint` is more accurate.

* `"euler"`. This is a variant on the above schemes with first order (explicit Euler) timestepping and lower spatial accuracy. It is less accurate and more dissipative than the aforementioned schemes, but can be computationally cheaper. For stability we find it requires a more dissipative spatial extrapolation ( default limiter coefficient `md%domains(j)%theta = 0.9`, which is less than second order accurate in space). It is not recommended for general usage, although may be fine for some problems. Regardless the underlying code is an important building block for the more advanced schemes.

* `"rk2n"`. This is an unusual timestepping scheme with second order accuracy in space and time. Internally it uses 5 first-order timesteps to advance 4 timesteps, while retaining second order accuracy in time. In principle it can be faster than the `"rk2"` and `"midpoint"` schemes, which internally use 2 first-order timesteps to advance 1 timestep. However the scheme is less stable, does not yet support dispersion, and is not recommended for general usage. 

For all schemes above the spatial extrapolation of the flow-state to cell edges employs the water-surface, depth, and velocity (x and y components), using the now classical approach of Audusse (*Audusse, E.; Bouchut, F.; Bristeau, M.-O.; Klein, R. & Perthame, B. A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows SIAM Journal of Scientific Computing, 2004, 25, 2050-2065*). The extrapolation employs a slope-limiter which controls the spatial order of accuracy. It is more-or-less dissipative depending on a domain-specific coefficient `theta`; the value for the `j`th domain is stored in `md%domains(j)%theta` and can be any non-negative real number: 

* `theta=0.0` corresponds to purely first-order extrapolation in space, and is quite dissipative. Higher values will make use of gradients when extrapolating edge values, and can have up-to second-order accuracy for smooth flows.

* `theta=1.0` is the smallest value that is second order accurate in space for sufficiently smooth flows. In practice it is still relatively dissipative (MinMod limiter). 

* `theta=2.0` is the upper threshold for which the extrapolation to edges is bounded by the neighbouring cell average values. This means that either `left-cell-value <= edge-value <= right-cell-value`, or `left-cell-value >= edge-value >= right-cell-value`.   

* `theta>2.0` is even less dissipative, and extrapolated edge values are not necessarily bound by neighbouring cell values. However in practice higher values (e.g. `theta=4.0`) retains good stability for many problems, even with shocks. Adjusting `theta` like this can be very useful in cases where the limiter otherwise leads to excessive dissipation, for instance some low Froude-number flows with many minima and maxima which are excessively clipped by limiters.

In all cases `theta` is tapered to zero in flows that are near a wet/dry boundary or extremely shallow, reducing the scheme to first order accuracy there. For details on how this is done see [domain_mod.f90](./src/shallow_water/domain_mod.f90) including the comments near the definition of `limiter_coef1 - limited_coef4`.

The numerical fluxes are derived from the left and right cell edge values using an approximate Riemann solver. Various options can be controlled by setting the character string `md%domains(j)%compute_fluxes_inner_method`. Values are:

* `"DE1"`. This is similar to the DE1 solver in ANUGA [see *Davies and Roberts, 2015, Open source flood modelling with a 2D discontinuous elevation hydrodynamic model. Proceedings of MODSIM 2015*]. It uses a HLL approximate Riemann solver, with a treatment of the pressure gradient term. Variants in this approach are popular in the literature.

* `"DE1_low_fr_diffusion"`. This is the default approach. It is similar to `DE1` but reduces the diffusive component of the UH/VH fluxes based on the Froude-number, similar to *Xue-song Li and Chun-wei Gu (2013) Mechanism and Improvement of the Harten–Lax–van Leer Scheme for All-Speed Flows. Computers & Fluids Volume 86, 5, Pages 56-70*.

* `"DE1_upwind_transverse"`. This is like `DE1` but uses an upwind approach for the transverse momentum flux. This approach is quite popular in the literature, for instance it is used in the shallow water solver provided by the well-known Basilisk code.

* `"DE1_low_fr_diffusion_upwind_transverse"`. This is like `DE1_low_fr_diffusion` but uses the upwind approach for the transverse momentum flux. 

Broadly speaking the solvers with `_low_fr_diffusion` might perform better in low Froude-number flows, although often results will be similar. The solvers with `_upwind_transverse` may have less numerical dissipation in shear flows, which may be a good or bad thing, but more often will be similar. 


## The Leapfrog schemes

The SWALS leapfrog schemes solve the nonlinear shallow water equations, with a particular focus on linear or quasi-linear variants that are efficient for low-Froude-number flow. The most complete scheme is `leapfrog_nonlinear`, which solves the nonlinear shallow water equations (although without any eddy-viscosity term).

Leapfrog schemes are classically used for deep-ocean tsunami propagation, and are also popular for inundation modelling (although this has not been a focus for leapfrog schemes in SWALS). They have good energy conservation properties for the nonlinear shallow water equations, which makes them much better suited to long-distance deep ocean tsunami propagation, as compared with the finite-volume schemes discussed above. All the SWALS leapfrog schemes also support dispersion (similar caveats as discussed above). Some references that describe classical leapfrog schemes for the linear and nonlinear shallow water equations are:

* *IOC Numerical method of tsunami simulation with the leap-frog scheme IUGG/IOC Time Project, IUGG/IOC Time Project, 1997*
* *Liu, P. L. F.; Cho, Y.-S.; Briggs, M. J.; Kanoglu, U. & Synolakis, C. E. Runup of solitary waves on a circular Island Journal of Fluid Mechanics, Cambridge University Press, 1995, 302, 259–285*

In general the leapfrog schemes in SWALS do not have good long-time stability when used in conjunction with nesting. They work well as the coarsest grid in a multidomain, but if used in a refined grid (which receives halo data from a coarser domain) then they often develop instability during long-time integration (even more so if using dispersion). This is not always a problem, especially for short model runs. If it is a problem then consider using a finite-volume scheme such as `midpoint` on nested domains (excepting the coarsest grids). 

The leapfrog solver variants provided by SWALS are:

* `"linear"`. This solves the linear shallow water equations (with Coriolis in spherical coordinates). The nonlinear shallow water equations are simplified by ignoring the nonlinear advection term $( \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) )$, the associated $\mathbf{M_s}$ term, and the nonlinear friction term $(g h \mathbf{S_f})$. The pressure gradient term $( g h \nabla \eta )$ is linearized to $g h_0 \nabla \eta$ by replacing the depth $h$ with the time-invariant depth at mean-sea-level $h_0$. In this case the mean-sea level must be specified in the variable `md%domains(j)%msl_linear` (=0.0 by default) so that SWALS can compute $h_0$. Alternatively one can retain nonlinearity in the pressure gradient term by setting `md%domains(j)%linear_solver_is_truely_linear=.false.`.

* `"leapfrog_linear_plus_nonlinear_friction"`. This is similar to `"linear"` above but also retains a nonlinear friction term (either Manning or Chezy as discussed below). The same parameters as discussed above should be specified to control nonlinearity in the pressure gradient term. This is useful as a cheap way to include friction in global scale tsunami models.

* `"leapfrog_nonlinear"`. This solves the nonlinear shallow water equations with a leapfrog scheme, including wetting and drying.

One should not try to specify nonlinear friction terms with the linear solvers. However, the non-starndard linear friction model can still be used. 


## The CLIFFS solver

This solver was developed by Elena Tolkova and is provided [this git repository](https://github.com/Delta-function/cliffs-src). It is similar the well-known MOST solver but uses a different wetting and drying technique. CLIFFS solves the nonlinear shallow water equations, but without any eddy-viscosity or Coriolis terms. SWALS adds support for dispersion using the same finite difference methods as used for the finite-volume schemes (although when using dispersion, the `midpoint` finite volume method appears to be typically more accurate than CLIFFS).

Elena Tolkova has reported considerable benchmarking of CLIFFS; see links in the README of the aforementioned repository, and also *Tolkova, E. Land--Water Boundary Treatment for a Tsunami Model With Dimensional Splitting Pure and Applied Geophysics, 2014, 171, 2289-2314*.

The CLIFFS solver was added to SWALS by extracting the inner computational routine, with only minor modifications to get it to compile in our framework. CLIFFS solves the shallow water equations using characteristic variables; in SWALS we transfer to/from these variables every time-step because our nesting/parallel framework requires the solver to work with the stage, UH, VH and bed elevation. For this reason the implementation in SWALS might be slower than a pure CLIFFS implementation. On the other hand CLIFFS is quite efficient anyway, and one can employ mixed openmp/mpi in SWALS which may permit larger parallel runs.

In general CLIFFS needs a wetting and drying threshold to be tuned to the problem at hand. This can be done by setting the variable `md%domains(j)%cliffs_minimum_allowed_depth`; typical values might be 0.001 (1 mm) for lab experiments, around 0.05 - 0.1 for high-resolution field cases, and 1 m or more for global-propagation problems. Tuning may be required for stability on any particular problem.

The CLIFFS solver requires sufficiently smooth bathymetry for stability; see the analysis by Tolkova in the CLIFFS user manual. A bathymetry-smoothing routine is provided by Tolkova for this purpose and has been integrated into SWALS; one must generally use this to avoid instabilities. To apply the bathymetry smoothing, after the elevation has been set in `md%domains(j)%U(:,:,ELV)` you would call the bathymetry smoothing routine, e.g. `call md%domains(j)%smooth_elevation(smooth_method='cliffs')`.

By default our CLIFFS solver uses `md%domains(j)%friction_type = "manning"`. Thus far we have not implemented Chezy friction in the CLIFFS solver (although that would not be difficult). We have implemented the linear-friction model associated with `md%domains(j)%linear_drag_coef` (discussed above), which is zero by default. In practice this model is unlikely to be desirable for use with CLIFFS, which is numerically dissipative in any case.

The CLIFFS solver is not mass conservative; see Tolkova papers for discussion of this, especially in relation to wetting and drying. For this reason in SWALS's nesting framework, no flux correction is applied between neighbouring domains if one uses CLIFFS, and mass conservation errors are expected.
