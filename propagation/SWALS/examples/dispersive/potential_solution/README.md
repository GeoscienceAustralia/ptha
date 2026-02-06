# Radially symmetric potential flow: analytical solution vs dispersive solver

Compares the linear dispersive solver with a radially symmetric potential flow
solution. One of the tests uses a single grid, and the other uses a nested grid
with non-symmetric placement.

NB: Extending the domain size and running this for longer will show up differences
between potential flow and the dispersive model in SWALS (which we expect).

## Single Grid Case (OpenMP)

### Transect comparison at final time
Comparison of numerical solution (SWALS with dispersion) versus analytical potential wave theory along transects at y=0 and x=0.

![Numerical vs potential solution transects](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution/numerical_vs_potential_solution_single_grid_OMP.png)

### Final water surface comparison
Spatial distribution of the water surface at the final time, showing the numerical solution (left) and analytical solution (right).

![Final water surface](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution/numerical_free_surface_single_grid_OMP.png)

## Nested Grid Case (MPI)

### Transect comparison at final time
Comparison of numerical solution (SWALS with dispersion using nested grids) versus analytical potential wave theory along transects at y=0 and x=0. The nested grid has non-symmetric placement.

![Numerical vs potential solution transects with nested grids](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution/numerical_vs_potential_solution_nested_with_MPI.png)

### Final water surface comparison
Spatial distribution of the water surface at the final time with nested grids, showing the numerical solution (left) and analytical solution (right).

![Final water surface with nested grids](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution/numerical_free_surface_nested_with_MPI.png)
