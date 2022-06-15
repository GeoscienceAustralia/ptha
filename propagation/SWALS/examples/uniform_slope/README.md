# Steady uniform flow on a uniform planar slope.

We simulate steady uniform flow down a uniform planar slope. Given a constant upstream discharge and ignoring the effects of the boundaries, we expect the solution to evolve to have a constant velocity and depth. Values of the latter can be computed from Manning's equation. 

The [SWALS model](uniform_slope.f90) simulates this problem in a closed domain, with a fixed discharge input upstream. It is configured to take the numerical solution method as a command line argument, and the run script tests two solution methods: `rk2`, and `leapfrog_linear_plus_nonlinear_friction`. The flow algorithm `leapfrog_linear_plus_nonlinear_friction` is not especially well suited to this problem, so the code gives it a special-case treatment of initial conditions and time-step to ensure stability as this solver evolves to steady state. There are no such issues with the `rk2` solver. 

In both cases we expect the numerical solution in the interior of the domain (where boundary effects are negligable) to accurately reproduce the Manning solution once it reaches steady state. The main program checks that the Manning solution is satisfied to high accuracy, and that mass is conserved in the closed domain. The outputs look like this for `rk2` (and are slightly less accurate for `leapfrog_linear_plus_nonlinear_friction`):
```  
## Testing ##
 
    Theoretical volume:    28800.000000003256     
    Model volume:    28800.000000000451     
      error:    2.8048816602677107E-009
PASS
 
   Theoretical vd:    1.0000000000000000     
   Model vd (near centre):   0.99999999979405507     
PASS
 
   Theoretical d:   0.48559337483020382     
   Model d:   0.48559337472148462     
PASS

```
