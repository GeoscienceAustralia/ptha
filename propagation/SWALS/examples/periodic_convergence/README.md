# Test order-of-accuracy with smooth solutions

This problem tests the numerical convergence of flow over a periodic domain with grid refinement. It was previously studied by [de la Asuncion et al (2013)](https://doi.org/10.1016/j.compfluid.2012.01.012) and [Xing et al. (2005)](https://doi.org/10.1016/j.jcp.2005.02.006). The idea is that for problems with smooth flow, we expect the numerical error of the model to decrease as the grid is refined in a predictable way that reflects the order of accuracy of the numerical scheme. 

To estimate the order of accuracy, the [SWALS model](model.f90) is first run on a coarse grid, and then on a sequence of finer grids with cell side-lengths a factor of 2, 4, and 8 times smaller. The [test code](test_convergence.R) uses the result on the finest grid is used as a reference to estimate the error on the other grids. The error is defined based on the maximum absolute value of the difference with the reference result anywhere on the domain. This is done separately for each of the flow quantities Stage, UH, and VH. 

For each flow quantity, the order of accuracy is estimated from the base-2 logarithm of the error ratios on successively finer grids (mesh refinement from 1-to-2, and then from 2-to-4). For the SWALS `rk2` finite-volume scheme, we expect convergence close to order 2 for this problem, and the test code checks that this is the case.
