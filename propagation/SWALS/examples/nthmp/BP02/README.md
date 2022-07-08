# Solitary wave on a composite beach (analytical solution and experiment)

This problem covers benchmarks 2 and 5 from the NTHMP test suite.

Three different solitary waves (cases A, B, C) are propagated over a 1D
composite beach. The analytical solutions are known for the linear shallow water
equations. 

![Figure 1: Beach profile and gauge locations](solution_geometry_caseA_linear.png)

For this problem we expect solutions based on the nonlinear shallow water
equations to deviate from the linear solution. Physically, nonlinearity is
significant for this problem because the wave amplitudes are not negligible compared
to the depth (especially for cases B and C). 

Dispersion also matters physically for this problem because the wave length is
sufficiently short compared to the water depth. It causes the reflected wave to
be delayed in the experimental data, as compared to the linear solution. This
also causes differences between experiments and non-dispersive shallow water
models.

The [SWALS model](BP2_testcases.f90) is setup to take the numerical method as a commandline
argument. This is used to test the linear shallow water solver (`linear`) as well as two
nonlinear shallow water solvers (`rk2` and `leapfrog_nonlinear`).

## Linear solution

Here we test the `linear` scheme, which should give excellent agreement with the analytical solution
derived from the linear shallow water equations.

![Figure 2: Comparison of analytical, experimental, and SWALS linear solutions for case A](solutionA_linear.png)

![Figure 3: Comparison of analytical, experimental, and SWALS linear solutions for case B](solutionB_linear.png)

![Figure 4: Comparison of analytical, experimental, and SWALS linear solutions for case C](solutionC_linear.png)

## Nonlinear solutions - `rk2`

Here we test the `rk2` finite volume scheme, which solves the nonlinear shallow water equations. 

In this problem nonlinear terms are not negligable, so we expect the numerical
solution to differ from the analytical linear solution.  The model also differs
significantly from the experiments due to dispersion, which tends to suppress
shock formation.

![Figure 5: Comparison of analytical, experimental, and SWALS nonlinear solutions for case A](solutionA_rk2.png)

![Figure 6: Comparison of analytical, experimental, and SWALS nonlinear solutions for case B](solutionB_rk2.png)

![Figure 7: Comparison of analytical, experimental, and SWALS nonlinear solutions for case C](solutionC_rk2.png)


## Nonlinear solutions - `leapfrog_nonlinear`

Here we test the `leapfrog_nonliear` finite difference scheme, which solves the
nonlinear shallow water equations. 

The numerical solution is almost identical to the `rk2` solution above, which
is unsurprising as both are solving the nonlinear shallow water equations.

In this problem nonlinear terms are not negligable, so we expect the numerical
solution to differ from the analytical linear solution.  The model also differs
significantly from the experiments due to dispersion, which tends to suppress
shock formation.

![Figure 8: Comparison of analytical, experimental, and SWALS nonlinear solutions for case A](solutionA_leapfrog_nonlinear.png)

![Figure 9: Comparison of analytical, experimental, and SWALS nonlinear solutions for case B](solutionB_leapfrog_nonlinear.png)

![Figure 10: Comparison of analytical, experimental, and SWALS nonlinear solutions for case C](solutionC_leapfrog_nonlinear.png)
