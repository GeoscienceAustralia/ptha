# Tests using dispersion

All the test problems here use dispersive solvers. 
* Folders beginning with `nthmp_` contain dispersive versions of some NTHMP test problems, which are adapted from the nonlinear shallow water variants [here](../../nthmp/). 
* [potential_solution](potential_solution) and [potential_solution_spherical](potential_solution_spherical) compare the model with potential flow solution, using both Cartesian and spherical coordinates.
* [shoaling_variable_seabed](shoaling_variable_seabed) models a highly dispersive/nonlinear wave.
* [solitary_shoaling_grilli](solitary_shoaling_grilli) models the shoaling of a solitary wave on a linear slope compares with experimental results from [Gilli et al (1994)](http://dx.doi.org/10.1061/(ASCE)0733-950X(1994)120:6(609)).
* [undular_bore](undular_bore) models the formation of an undular bore and compares the results with a more complex dispersive model.
