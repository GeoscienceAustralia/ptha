# Non-symmetric potential flow: analytical solution vs dispersive solver in spherical coordinates

Compares the dispersive solver with a **non-symmetric** potential flow
solution using spherical coordinates. This is a variation of the radially symmetric case in
[../potential_solution_spherical](../potential_solution_spherical). 

The initial displacement has different scales in the x and y directions
(elliptical Gaussian plus asymmetric oscillations). The domain extent is small
enough that the sphere is well approximated by the Cartesian solution.

One of the tests uses a single grid, and the other uses a nested grid with
non-symmetric placement.

## Single Grid Case (OpenMP)

### Transect comparison at final time
Comparison of numerical solution (SWALS with dispersion) versus analytical potential wave theory along transect at x=0.

![Numerical vs potential solution transects](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_nonsymmetric/numerical_vs_potential_solution_single_grid_OMP.png)

### Final water surface comparison
Spatial distribution of the water surface at the final time, showing the numerical solution (left) and analytical solution (right).

![Final water surface](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_nonsymmetric/numerical_free_surface_single_grid_OMP.png)

## Nested Grid Case (MPI)

### Transect comparison at final time
Comparison of numerical solution (SWALS with dispersion using nested grids) versus analytical potential wave theory along transect at x=0. The nested grid has non-symmetric placement.

![Numerical vs potential solution transects with nested grids](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_nonsymmetric/numerical_vs_potential_solution_nested_with_MPI.png)

### Final water surface comparison
Spatial distribution of the water surface at the final time with nested grids, showing the numerical solution (left) and analytical solution (right).

![Final water surface with nested grids](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_nonsymmetric/numerical_free_surface_nested_with_MPI.png)
