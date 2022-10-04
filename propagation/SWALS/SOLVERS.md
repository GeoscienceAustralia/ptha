# SWALS solvers
----------------

SWALS contains a range of solvers for variants of the 2D shallow water equations. For models that involve multiple domains (i.e. using a `multidomain` object, typically denoted `md`), one can use different solvers on different domains, and they may take different time-steps. [See here for a study that uses a range of SWALS solvers](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full).

In general the different solvers have different strengths and limitations, which are discussed below. For instance some of our solvers work well for global-scale tsunami propagation (at least for long-waves that satisfy the shallow water equations), but are inappropriate for nearshore inundation simulation. The converse holds for some other solvers. Furthermore, not all solvers in SWALS are equally well supported by the boundary conditions and the nesting formulation. 

Generally speaking the SWALS code development has focussed most on:

* A range of shock-capturing finite volume schemes. These are mainly used for flows with significant nonlinearity, medium to high Froude-numbers, etc. They are quite robust for this purpose, but tend to be more diffusive than the leap-frog schemes (and thus less well suited to global scale tsunami propagation). Some variants of these schemes have better performance at low Froude-numbers than others.

* A range of leap-frog schemes. These solve the linear shallow water equations; variants on the latter with some additional nonlinear terms (e.g. nonlinear pressure gradient terms; manning friction); and there is a variant with all nonlinear terms. The leap-frog schemes have mostly been used for flows domainated by low Froude-numbers (especially global-scale tsunami propagation simulation). The fully nonlinear leap-frog scheme can perform well in some more strongly nonlinear problems (e.g. inundation), but we haven't focussed on making it extremely robust for such purposes (and generally prefer the finite-volume schemes such as `rk2`).

We also include the CLIFFS solver developed by Elena Tolkova, but have not worked extensively with that.


# Setting the solver type.

A SWALS `multidomain` object (denoted `md`) contains one or more structured grids (`domains`) on which we solve the shallow water equations. These domains are stored in an allocatable array (`md%domains(:)`) of type `domain_type`. This is allocated at the beginning of an application code (e.g. `allocate(md%domains(10))` will allow 10 domains). 

The solver for each domain is set by specifying its timestepping method. For example by setting `md%domains(3)%timestepping_method = "linear"`, we make the 3rd domain use the linear leap-frog solver. Many examples are included in the validation test suite, such as [here](./examples/nthmp/BP09/BP09.f90). 

The available solver types are specified in the file [timestepping\_metadata\_mod.f90](./src/shallow_water/timestepping_metadata_mod.f90). This also sets solver specific variables such as the halo thickness required to advance one time-step. Additional aspects of the solver behavior are controlled by other variables in the domain. 

# Details on the solvers

Multiple numerical methods are supported in SWALS (detailed below). The validation test-suite most often uses `default_nonlinear_timestepping_method` (currently equal to `rk2`) for nonlinear problems, and `default_linear_timestepping_method` (currently equal to `linear`) for linear problems. 

With these defaults, other schemes are less heavily exercised by the validation test-suite. If you wish to better understand the performance of other schemes, the default nonlinear and linear solvers used in the test suite can be changed by setting `default_nonlinear_timestepping_method` and `default_linear_timestepping_method` in [global\_mod.f90](./src/shallow_water/global_mod.f90). This will affect many but not all test problems. A number of tests deliberately exercise non-default solvers by directly specifying `md%domains(j)%timestepping_method`, rather than using the defaults. In such cases one should directly change `md%domains(j)%timestepping_method` in the code used to setup the problem.

## The Finite-Volume solvers

SWALS has a number of classical shock-capturing finite-volume schemes with accuracy up to second order in space and time. These solve the nonlinear shallow water equations with Manning (or Chezy) friction. 

$$ \frac{\partial \eta}{\partial t} + \nabla \cdot \mathbf{q} = 0 $$

$$ \frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) + g h \nabla \eta + g h \mathbf{S_f} + \mathbf{\Psi} + \mathbf{\Omega_s} + \mathbf{M_s} = 0 $$

Here $\eta = h + z$ is the free surface elevation (m), $t$ is time (s), $h$ is the depth (m), $z$ is the bed elevation (m), $\mathbf{q}=h\mathbf{u}$ is the 2D flux vector $(m^2/s)$, $\mathbf{u} = (u,v)$ is the 2D velocity vector (m/s), $\otimes$ is the [outer product](https://en.wikipedia.org/wiki/Outer_product), $g$ is gravity $(m/s^2)$, $\mathbf{S_f}$ is the friction slope vector, $\mathbf{\Psi}$ is a turbulent diffusion term, $\mathbf{\Omega_s}$ is the Coriolis force for spherical coordinates, and $\mathbf{M_s}$ gives some additional metric terms for spherical coordinates. 

In Cartesian coordinates $(x,y)$ are in meters. The "grad" operator $\nabla f$ applied to a scalar function $f$ gives $( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y})$. The divergence operator $\nabla \cdot \mathbf{u}$ applied to a vector $\mathbf{u} = (u,v)$ gives $( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y})$, and applied to a matrix it gives the vector resulting from application to each column separately. 

In Spherical coordinates $(\theta, \phi)$ are (longitude, latitude) in degrees. The "grad" operator $\nabla f$ applied to a scalar function $f$ gives $( \frac{1}{R \cos(\phi_r)} \frac{\partial f}{\partial \theta_r}, \frac{1}{R \phi_r} \frac{\partial f}{\partial \phi_r})$ where the subscript $r$ indicates angles in radians rather than degrees, and $R$ is the earth radius in meters. The divergence operator $\nabla \cdot \mathbf{u}$ applied to a vector $\mathbf{u}$ gives $( \frac{1}{R \cos(\phi_r)} (\frac{\partial u}{\partial \theta_r} + \frac{\partial \cos(\phi_r) v}{\partial \phi_r}))$, and applied to a matrix it gives the vector resulting from application to each column. We stress that SWALS uses units of degrees (not radians) for longitude and latitude.

The friction slope term by default uses Manning friction $\mathbf{S_f} = n^2\mathbf{u} |\mathbf{u}| h^{-4/3}$, with Manning coefficient $n$ that is zero by default. But it can also use Chezy friction, which can be written in a similar form using a different power of the depth $(\mathbf{S_f} = n^2\mathbf{u} |\mathbf{u}| h^{-1})$, where $n^{-2}$ is an unusual way of writing the Chezy coefficient. A non-standard linear friction term is also available (as discussed below). 

The turbulent diffusion term $\mathbf{\Psi}$ is by default set to zero. But an eddy viscosity formulation can also be applied, $\mathbf{\Psi} = -\nabla \cdot (\epsilon h (\nabla \otimes \mathbf{u}))$. Here $\epsilon$ is the eddy viscosity model and $\nabla \otimes \mathbf{u}$ is a 2x2 matrix with columns $\nabla u$, $\nabla v$ respectively. More detail is provided on the turbulent diffusion term below. 

In Cartesian coordinates there is no Coriolis term $(\mathbf{\Omega_s} = 0)$. In Spherical coordinates the Coriolis force is $\mathbf{\Omega_s} = 2 \sin(\phi_r) \gamma_r (-vh, uh)$, where $\gamma_r = 7.292115 \times 10^{-5}$ gives the earth angular frequency in radians per second.

In Cartesian coordinates there are no spherical coordinate metric terms $(\mathbf{M_s} = 0)$. In Spherical coordinates $\mathbf{M_s} = \frac{\tan(\phi_r)}{R} u (-vh, uh)$; this leads to a very small correction away from the poles and is often ignored in the tsunami literature. For examples of papers which discuss it, see [Williamson et al. 1992](https://doi.org/10.1016/S0021-9991(05)80016-6), [Titov et al. 2016](http://dx.doi.org/10.1061/(ASCE)WW.1943-5460.0000357), and [Popinet 2011](https://doi.org/10.1007/s10236-011-0438-z). 

The SWALS finite-volume solvers approximate the above equations in flux conservative form. They are shock capturing, and well suited to flows with moderate or high Froude-numbers and wetting/drying, but may be too numerically dissipative to efficiently model flow at very low Froude-numbers. For example they can work well for nearshore and inundation simulation, but are not as well suited to global-scale tsunami propagation as the leapfrog schemes (discussed below). In practice we often develop global-to-local nested grid models by using a leapfrog solver on the coarsest global grid, and finite-volume solvers on the nested grids (such as in [this paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full)). 

In general the finite-volume schemes in SWALS have good stability when used in conjunction with nesting. Occasionally these solvers can be subject to artifical vortices at coarse-to-fine nesting regions, especially where the elevation has rapid variation. This is not particularly common though, and can usually be solved by moving the domain boundary to an area with weaker elevation variation, or by slightly smoothing the elevation data in the neighbourhood of the spurious flow (the subroutine `md%domains(j)%smooth_elevation_near_point` can do this).

The finite-volume schemes in SWALS use a range of different timestepping schemes, flux functions, and spatial extrapolations to derive the flow state at edges. They are essentially variants of the scheme described [here](https://www.researchgate.net/publication/289882157_Open_source_flood_simulation_with_a_2D_discontinuous-elevation_hydrodynamic_model), but on a structured grid rather than an unstructured mesh. Each scheme has a default limiter coefficient, flux function, and friction model, but these can be changed by the user prior to calling `md%setup`:

The finite-volume solvers are supported by setting `md%domains(j)%timestepping_method` to:

* `"rk2"`. This is a good default choice, and is heavily exercised in the validation test suite. It has second order accuracy in space and time. It uses 2nd-order Runge-kutta timestepping (i.e. advancing two first-order timesteps and taking the average). By default it uses second order spatial extrapolation with a default limiter coefficient `md%domains(j)%theta = 1.6`.

* `"midpoint"`. This is like `"rk2"` but uses a midpoint timestepping method. Results are generally very similar to `"rk2"`.

* `"euler"`. This is a variant on the above schemes with first order (explicit Euler) timestepping and lower spatial accuracy. It is less accurate and more dissipative than the aforementioned schemes, but can be computationally cheaper. For stability we find it requires a more dissipative spatial extrapolation ( default limiter coefficient `md%domains(j)%theta = 0.9`, which is less than second order accurate in space). It is not recommended for general usage, although may be fine for some problems. Regardless the underlying code is an important building block for the more advanced schemes.

* `"rk2n"`. This is an unusual timestepping scheme with second order accuracy in space and time. Internally it uses 5 first-order timesteps to advance 4 timesteps, while retaining second order accuracy in time. In principle it can be faster than the `"rk2"` and `"midpoint"` schemes, which internally use 2 first-order timesteps to advance 1 timestep. However the scheme is less stable, and not recommended for general usage. 

For all schemes above the spatial extrapolation of the flow-state to cell edges employs the water-surface, depth, and velocity (x and y components), using the now classical approach of Audusse (*Audusse, E.; Bouchut, F.; Bristeau, M.-O.; Klein, R. & Perthame, B. A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows SIAM Journal of Scientific Computing, 2004, 25, 2050-2065*). The extrapolation employs a slope-limiter which controls the spatial order of accuracy. It is more-or-less dissipative depending on a domain-specific coefficient `theta`; the value for the `j`th domain is stored in `md%domains(j)%theta` and can be any non-negative real number: 

* `theta=0.0` corresponds to purely first-order extrapolation in space, and is quite dissipative. Higher values will make use of gradients when extrapolating edge values, and can have up-to second-order accuracy for smooth flows.

* `theta=1.0` is the smallest value that is second order accurate in space for sufficiently smooth flows. In practice it is still relatively dissipative (MinMod limiter). 

* `theta=2.0` is the upper threshold for which the extrapolation to edges is bounded by the neighbouring cell average values. This means that either `left-cell-value <= edge-value <= right-cell-value`, or `left-cell-value >= edge-value >= right-cell-value`.   

* `theta>2.0` is even less dissipative, and extrapolated edge values are not necessarily bound by neighbouring cell values. However in practice higher values (e.g. `theta=4.0`) retains good stability for many problems, even with shocks. Adjusting `theta` like this can be very useful in cases where the limiter otherwise leads to excessive dissipation, for instance some low Froude-number flows with many minima and maxima which are excessively clipped by limiters.

In all cases `theta` is tapered to zero in flows that are near a wet/dry boundary or extremely shallow, reducing the scheme to first order accuracy there. For details on how this is done see [domain_mod.f90](./src/shallow_water/domain_mod.f90) including the comments near the definition of `limiter_coef`1-3.

The numerical fluxes are derived from the left and right cell edge values using an approximate Riemann solver. Various options can be controlled by setting the character string `md%domains(j)%compute_fluxes_inner_method`. Values are:

* `"DE1"`. This is similar to the DE1 solver in ANUGA [see *Davies and Roberts, 2015, Open source flood modelling with a 2D discontinuous elevation hydrodynamic model. Proceedings of MODSIM 2015*]. It uses a HLL approximate Riemann solver, with a treatment of the pressure gradient term. Variants in this approach are popular in the literature.

* `"DE1_low_fr_diffusion"`. This is the default approach. It is similar to `DE1` but reduces the diffusive component of the UH/VH fluxes based on the Froude-number, similar to *Xue-song Li and Chun-wei Gu (2013) Mechanism and Improvement of the Harten–Lax–van Leer Scheme for All-Speed Flows. Computers & Fluids Volume 86, 5, Pages 56-70*.

* `"DE1_upwind_transverse"`. This is like `DE1` but uses an upwind approach for the transverse momentum flux. This approach is quite popular in the literature, for instance it is used in the shallow water solver provided by the well-known Basilisk code.

* `"DE1_low_fr_diffusion_upwind_transverse"`. This is like `DE1_low_fr_diffusion` but uses the upwind approach for the transverse momentum flux. 

Broadly speaking the solvers with `_low_fr_diffusion` might perform better in low Froude-number flows, although often results will be similar. The solvers with `_upwind_transverse` may have less numerical dissipation in shear flows, which may be a good or bad thing, but more often will be similar. 

The friction model can be controlled by setting the variable `md%domains(j)%friction_type`. Values are:

* `"manning"`. This is Manning friction, so the array `md%domains(j)%manning_squared(:,:)` is interpreted as the (manning friction)^2

* `"chezy"`. This is Chezy friction. In this case the array `md%domains(j)%manning_squared(:,:)` is interpreted as (1/chezy_friction)^2

In addition one can set a linear friction coefficient `md%domains(j)%linear_drag_coef` which is zero by default. If non-zero this implements a linear-friction model following *Fine, I. V.; Kulikov, E. A. & Cherniawsky, J. Y. Japans 2011 Tsunami: Characteristics of Wave Propagation from Observations and Numerical Modelling Pure and Applied Geophysics, Springer Science and Business Media LLC, 2012, 170, 1295-1307*. In practice you probably do not want to use this for the finite-volume schemes, which are already somewhat numerically dissipative; it is more likely to be useful to add slow friction to the leap-frog schemes.

By default there is no additional turbulence model. However a simple eddy viscosity model can be turned on by specifying `md%domains(j)%use_eddy_viscosity=.true.` and adjusting the values of `md%domains(j)%eddy_visc_constants(1:2)`. 
* With this option, turbulence-terms of the form $\nabla \cdot (\epsilon h (\nabla \otimes \mathbf{u}))$ are used on the right-hand-side of the depth integrated velocity equations (for the evolution of $\mathbf{q} = h\mathbf{u}$). 
* The eddy viscosity $\epsilon$ is set equal to `md%domains(j)%eddy_visc_constants(1) + md%domains(j)%eddy_visc_constants(2) * depth * shear_velocity`. 
    * The first term allows a simple constant eddy viscosity
    * The second term corresponds to a very common eddy viscosity model (e.g. *Shiono, K. & Knight, D. W. 1991. Turbulent open-channel flows with variable depth across the channel. Journal of Fluid Mechanics, 222, 617-646*) 
* A simple second-order-in-space explicit scheme is used to discretise the equations. The eddy viscosity $\epsilon$ is clipped if needed to prevent violation of the time-step constraint for explicit diffusion (`eddy_viscosity <= (min_grid_cell_side_length**2 / (2 * timestep))`). This approach can work well in practical cases where turbulent diffusion is not very strong. For models with high resolution and/or strong turbulent diffusion, you may need to reduce the `timestep` to prevent clipping from affecting the solution. The same approach to limiting the eddy-viscosity is used in the model SWASH (*Zijlema, M.; Stelling, G. & Smit, P. SWASH: An operational public domain code for simulating wave fields and rapidly varied flows in coastal waters Coastal Engineering, 2011, 58, 992 - 1012*).

## The Leap-frog schemes

The SWALS leap-frog schemes solve the nonlinear shallow water equations, with a particular focus on linear or quasi-linear variants that are efficient for low-Froude-number flow. The most complete scheme is `leapfrog_nonlinear`, which solves the nonlinear shallow water equations without any eddy-viscosity term:

$$ \frac{\partial \eta}{\partial t} + \nabla \cdot \mathbf{q} = 0 $$

$$ \frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) + g h \nabla \eta + g h \mathbf{S_f} + \mathbf{\Omega_s} + \mathbf{M_s} = 0 $$

For definitions of these variables see the previous section on finite-volume schemes. 

Leap-frog schemes are classically used for deep-ocean tsunami propagation, and are also popular for inundation modelling (although this has not been a focus for leap-frog schemes in SWALS). They have good energy conservation properties which makes them much better suited to long-distance deep ocean tsunami propagation, as compared with the finite-volume schemes discussed above. Some references that describe classical leap-frog schemes for the linear and nonlinear shallow water equations are:

* *IOC Numerical method of tsunami simulation with the leap-frog scheme IUGG/IOC Time Project, IUGG/IOC Time Project, 1997*
* *Liu, P. L. F.; Cho, Y.-S.; Briggs, M. J.; Kanoglu, U. & Synolakis, C. E. Runup of solitary waves on a circular Island Journal of Fluid Mechanics, Cambridge University Press, 1995, 302, 259–285*

In general the leap-frog schemes in SWALS do not have good long-time stability when used in conjunction with nesting. They work well as the coarsest grid in a multidomain, but if used in a refined grid (which receives halo data from a coarser domain) then they often develop instability during long-time integration. This is not always a problem, especially for short model runs. If it is a problem then consider using a finite-volume scheme such as `rk2` on problematic domains. 

The leap-frog solver variants provided by SWALS are:

* `"linear"`. This solves the linear shallow water equations in the default case that `md%domains(j)%linear_solver_is_truely_linear=.true.`. In that case the pressure gradient term $( g h \nabla \eta )$ is linearized by replacing the depth $h$ with the depth at mean-sea-level $h_0$, so the term becomes $(g h_0 \nabla \eta )$. In this case the mean-sea level should be set by specifying `md%domains(j)%msl_linear` (=0.0 by default) so that SWALS knows how to compute $h_0$. In addition we ignore other nonlinear terms, including the nonlinear advection term $( \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) )$, the nonlinear friction term $(g h \mathbf{S_f})$, and the nonlinear spherical coordinate metric term $(\mathbf{M_s})$. These settings correspond to the standard linear shallow water equations, including the Coriolis force in spherical coordinates. Note that one can retain nonlinearity in the pressure gradient term by setting `md%domains(j)%linear_solver_is_truely_linear=.false.`, while still ignoring the other nonlinear terms.

* `"leapfrog_linear_plus_nonlinear_friction"`. This is similar to `"linear"` above but also retains a nonlinear friction term (either Manning or Chezy as discussed below). The same parameters as discussed above should be specified to control nonlinearity in the pressure gradient term. This is useful as a cheap way to include friction in global scale tsunami models.

* `"leapfrog_nonlinear"`. This solves the nonlinear shallow water equations with a leap-frog scheme, including wetting and drying.


For all leap-frog solvers except `"linear"`, the friction model can be controlled by setting the variable `md%domains(j)%friction_type`. Values are:

* `"manning"`. This is Manning friction, so the array `md%domains(j)%manning_squared(:,:)` is interpreted as the (manning friction)**2

* `"chezy"`. This is Chezy friction. In this case the array `md%domains(j)%manning_squared(:,:)` is interpreted as (1/chezy_friction)**2

For all solvers (including `"linear"`), one can set the linear friction coefficient `md%domains(j)%linear_drag_coef` which is 0.0 by default. If this is set to a non-zero value then it implements the linear friction model of *Fine, I. V.; Kulikov, E. A. & Cherniawsky, J. Y. Japans 2011 Tsunami: Characteristics of Wave Propagation from Observations and Numerical Modelling Pure and Applied Geophysics, Springer Science and Business Media LLC, 2012, 170, 1295-1307*. 


## The CLIFFS solver

This solver was developed by Elena Tolkova and is provided [this git repository](https://github.com/Delta-function/cliffs-src). It is similar the well-known MOST solver but uses a different wetting and drying technique. CLIFFS solves the nonlinear shallow water equations, but without any eddy-viscosity or Coriolis terms.

$$ \frac{\partial \eta}{\partial t} + \nabla \cdot \mathbf{q} = 0 $$

$$ \frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot (\mathbf{u}\otimes\mathbf{q}) + g h \nabla \eta + g h \mathbf{S_f} + \mathbf{M_s} = 0 $$

For definitions of these variables see the above section on Finite Volume Schemes. Tolkova has reported considerable benchmarking of CLIFFS; see links in the README of the aforementioned repository, and also *Tolkova, E. Land--Water Boundary Treatment for a Tsunami Model With Dimensional Splitting Pure and Applied Geophysics, 2014, 171, 2289-2314*.

The CLIFFS solver was added to SWALS by extracting the inner computational routine, with only minor modifications to get it to compile in our framework. CLIFFS solves the shallow water equations using characteristic variables; in SWALS we transfer to/from these variables every time-step because our nesting/parallel framework requires the solver to work with the stage, UH, VH and bed elevation. For this reason the implementation in SWALS might be slower than a pure CLIFFS implementation. On the other hand CLIFFS is quite efficient anyway, and one can employ mixed openmp/mpi in SWALS which may permit larger parallel runs.

In general CLIFFS needs a wetting and drying threshold to be tuned to the problem at hand. This can be done by setting the variable `md%domains(j)%cliffs_minimum_allowed_depth`; typical values might be 0.001 (1 mm) for lab experiments, around 0.05 - 0.1 for high-resolution field cases, and 1 m or more for global-propagation problems. Tuning may be required for stability on any particular problem.

The CLIFFS solver requires sufficiently smooth bathymetry for stability; see the analysis by Tolkova in the CLIFFS user manual. A bathymetry-smoothing routine is provided by Tolkova for this purpose and has been integrated into SWALS; one must generally use this to avoid instabilities. To apply the bathymetry smoothing, after the elevation has been set in `md%domains(j)%U(:,:,ELV)` you would call the bathymetry smoothing routine, e.g. `call md%domains(j)%smooth_elevation(smooth_method='cliffs')`.

By default our CLIFFS solver uses `md%domains(j)%friction_type = "manning"`. This is Manning friction, and can be spatially variable. It is set using the array `md%domains(j)%manning_squared(:,:)` which interpreted as the (manning friction)^2. Thus far we have not implemented Chezy friction in the CLIFFS solver (although that would not be difficult). We have implemented the linear-friction model associated with `md%domains(j)%linear_drag_coef` (discussed above), which is zero by default. In practice this model is unlikely to be desirable for use with CLIFFS, which is numerically dissipative in any case.

The CLIFFS solver is not mass conservative; see Tolkova papers for discussion of this, especially in relation to wetting and drying. For this reason in SWALS's nesting framework, no flux correction is applied between neighbouring domains if one uses CLIFFS, and mass conservation errors are expected.
