# Radially symmetric potential flow: analytical solution vs dispersive solver

Compares the linear dispersive solver with a radially symmetric potential flow
solution. One of the tests uses a single grid, and the other uses a nested grid
with non-symmetric placement.

NB: Extending the domain size and running this for longer will show up differences
between potential flow and the dispersive model in SWALS (which we expect).
