# Multidomain design namelists
Folder containing the namelist files and the load balance file for the multidomain simulations. When invoking SWALS, input the relevant namelist.

## Load balance files

The first row has the MPI ranks (modular). So this would work for 48 processes, or any divisor of 48 (2, 4, 8, 16, 24, 48).
1, 2, ..., 48
1
2
3
3

Each row has the domain index. This is done by create_boxes.R in the multidomain folder.
