# Very shallow steady uniform flow on a uniform planar slope.

We simulate very shallow steady uniform flow (< 2mm deep) down a uniform planar slope. Given a constant upstream discharge and ignoring the effects of the boundaries, we expect the solution to evolve to have a constant velocity and depth. Values of the latter can be computed from Manning's equation. The main numerical challenge of this problem relates to the shallowness of the flow. 

The [SWALS model](uniform_slope.f90) simulates this problem in a closed domain, with a fixed discharge input upstream. The `rk2` solver is used. We expect the numerical solution in the interior of the domain (where boundary effects are negligable) to accurately reproduce the Manning solution once it reaches steady state. The main program checks that the Manning solution is satisfied to high accuracy, and that mass is conserved in the closed domain. The outputs look like this for `rk2`:
```  
## Testing ##
     Froude:   0.37581822142801274     
  
     Theoretical volume:    103.42080000000001     
     Model volume:    103.42079999808000     
       error:    1.9200143697162275E-009
 PASS
  
    Theoretical d:    1.9331820449317630E-003
    Model d:    1.9331813419861260E-003
 PASS
  
    Theoretical vd:    1.0000000000000000E-004
    Model vd (near centre):    9.9999939343947084E-005
PASS
```

If we reduce the order of accuracy of the solver then the error is significant; see comments in the test code that show how to do this.
