Test against the analytical solution of "Zhang and Zhu (1994) New solutions for
the propagation of long waves over variable depth. Journal of fluid mechanics
278: 391-406"

To run and make plots comparing the analytical/numerical solution, do

    make -B -f make_circular_island
    # Choose any number of threads >= 1
    export OMP_NUM_THREADS=6 
    ./circular_island
    Rscript plot.R

To make the simulation reasonably fast, it is run on a 2km grid, and small but
detectible differences with the analytical solution remain. These reduce with
grid refinement.
