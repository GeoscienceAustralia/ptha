This folder can store codes for tsunami propagation modelling. 

Currently the code SWALS is included (SWALS = Shallow WAter Like Solvers).

Other well known open source shallow water solvers include GEOCLAW, JAGURS,
COMCOT, ANUGA, Basilisk, and easyWave. They vary widely in their features, and
the required computational effort.

We aim to keep the other parts of the PTHA package independent of any
particular solver. This is especially enforced in rptha, but less so for
application-specific template scripts (for instance, assumptions about the file
structure are embedded in the 2018 Australian PTHA template scripts). However
most of our template scripts could be adapted to use any solver which is able
to output stage time-series at tide gauges. In general this will require the
user to change function calls which read the solver output data, to refer to
new (user provided) functions which read from the new solver. 
