This folder can store codes for tsunami propagation modelling. 

Currently a simple linear shallow water equations solver is provided (SWALS),
which is suitable for offshore tsunami propagation (so long as disperision can
be neglected), but not suitable for inundation modelling or nearshore tsunami
modelling (for which the non-linear equations are required).

Users are encouraged to use other propagation codes instead, if desired,
depending on whatever suits their needs. Well known open source examples
include GEOCLAW, JAGURS, COMCOT, ANUGA, and easyWave. Some of the latter
include a range of other solvers, e.g. non-linear, dispersive, etc. They
vary widely in their features, and the required computational effort. 

As much as possible, we try to keep the code logic in other parts of the PTHA
package independent of any particular solver. This is especially enforced in
rptha. Most of our template scripts for PTHA can be adapted to use any solver
which is able to output stage time-series at tide gauges. In general this will
require the user to change function calls which read the solver output data, to
refer to new (user provided) functions which read from the new solver. 
