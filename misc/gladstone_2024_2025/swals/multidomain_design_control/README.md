# Multidomain design namelists
Folder containing the namelist files and the load balance file for the multidomain simulations. When invoking SWALS, input the relevant namelist.

## Load balance files

The first row has the MPI ranks (modular). So this would work for 48 processes, or any divisor of 48 (2, 4, 8, 16, 24, 48).
1, 2, ..., 48
1
2
3
3

Each row corresponds to a domain index (e.g. row 1 = domain 1, row 2 = domain 2, ...). Default load balance files are created by the script create_boxes.R in the multidomain_design folder. These give a starting point to run the model but won't usually perform well. Load balance files tuned to the actual model can be constructed using [../post_process/load_balance_script.R](../post_process/load_balance_script.R).
