# Radially symmetric potential flow in spherical coordinates: analytical solution vs dispersive solver

This problem is just like [../potential_solution](../potential_solution) but
uses spherical coordinates with a domain centred at 45 degrees north. The
domain is very small, so a local Cartesian approximation is accurate.

This provides an opportunity to check the dispersive solver in spherical
coordinates.

## Validation approach

The numerical model results are compared against an analytical potential flow solution. Two test cases are run:
1. **Single grid with OpenMP** - Tests the dispersive solver on a uniform grid
2. **Nested grid with MPI** - Tests the dispersive solver with grid nesting and domain decomposition

## Single grid validation (OpenMP)

### Transect comparison at final time
![Numerical vs analytical solution transects](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_spherical/numerical_vs_potential_solution_EW_single_grid_OMP.png)

Comparison of numerical and analytical solutions along east-west and north-south transects at the final time. The SWALS dispersive solver closely matches the potential wave theory solution.

### Free surface comparison
![Numerical vs analytical free surface](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_spherical/numerical_free_surface_single_grid_OMP.png)

Comparison of the full water surface elevation between the numerical model (left) and analytical solution (right) at the final time.

## Nested grid validation (MPI)

### Transect comparison at final time
![Numerical vs analytical solution transects (nested)](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_spherical/numerical_vs_potential_solution_EW_nested_with_MPI.png)

Comparison of numerical and analytical solutions along east-west and north-south transects for the nested grid case. Results from both the coarse and fine grids are shown, demonstrating accurate solution across the nested domain boundary.

### Free surface comparison
![Numerical vs analytical free surface (nested)](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/dispersive/potential_solution_spherical/numerical_free_surface_nested_with_MPI.png)

Comparison of the full water surface elevation for the nested grid configuration, showing the numerical model (left) and the analytical solution (right). A dotted line shows the high-resolution priority domain region.
