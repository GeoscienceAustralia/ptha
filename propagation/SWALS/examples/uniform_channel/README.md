# Steady uniform flow in a grid-aligned rectangular channel.

We simulate steady uniform flow in a channel with rectangular cross-section. The channel has uniform across-channel width, and a constant downstream slope. Given a constant discharge and ignoring the effects of the boundaries, we expect the solution to evolve to have a constant within channel velocity and depth. Values of the latter can be computed from Manning's equation. 

The [SWALS model](uniform_channel.f90) simulates this problem in a closed domain, with a fixed discharge input in the upstream part of the channel. We expect the numerical solution in the interior of the domain (where boundary effects are negligable) to accurately reproduce the Manning solution once it reaches steady state. The main program checks that the Manning solution is satisfied to high accuracy, and that mass is conserved in the closed domain. The outputs look something like this:
```  
## Testing ##

    Theoretical volume:    18000.000000003303     
    Model volume:    18000.000000000255     
      error:    3.0486262403428555E-009
PASS
 
   Theoretical vd:    1.0000000000000000     
   Model vd (near centre):   0.99999999999925315     
PASS
 
   Theoretical d:   0.48559337483020382     
   Model d:   0.48559337482999831     
PASS
```
