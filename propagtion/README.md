This folder can store codes for tsunami propagation modelling. 

Currently a simple linear shallow water equations solver is provided (SWALS). 

Users are encouraged to use other propagation codes instead, if desired,
depending on whatever suits their needs. Well known open source examples
include JAGURS and easyWave (with the former providing a range of other
solvers, e.g. non-linear, dispersive, etc). 

As much as possible, we try to keep the code logic in other parts of the PTHA
package independent of any particular solver. This is especially enforced in
rptha. Most of our template scripts for PTHA can be adapted to use any solver
which is able to output stage time-series at tide gauges. In general this will
will require the user to change function calls which read the solver output
data, to refer to new (user provided) functions which read from the new solver. 
