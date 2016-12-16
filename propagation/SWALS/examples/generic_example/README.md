Here we attempt to give a fairly generic interface to a spherical coordinates
linear shallow water model.

To run it:
* It is **strongly** suggested that you ensure the unit tests all pass before
attempting to run the current code. See
[../../test/unit_tests](../../test/unit_tests). These will help detect
problems with your computing environment which might be hard to diagnose by
just running a model.
* Assuming all the unit tests pass, please continue .....!
* Copy 'model_global_4min.in' to a new file (e.g. 'model_test.in'), and edit it
to ensure that the input_elevation_raster and input_stage_raster exist on your
filesystem, and that the model extents are as desired. The input rasters should
give the initial elevation and stage in lon-lat coordinates. 
* The elevation raster should cover the desired model domain, but does
not have to have exactly the same north-east-south-west extent -- since the
model will get the data it needs using bilinear interpolation. So for example, 
you can run a small area model using a global DEM, without editing the DEM.
Note however that the model does not account for the periodicity of longitude
(i.e. if the input raster has east-west extent [-180, 180], then the tsunami
model cannot have east-west extent exceeding 180, or less than -180). 
* The stage raster can be either smaller or larger than the desired model domain. The
model will extract stages from this raster where possible, using bilinear
interpolation. It will use a value of 0 elsewhere. 
* If the model extent covers 360 degrees of latitude, then east-west periodic
boundary conditions are used, with reflective north-south boundaries. Land is
always treated as reflective.  In this case, the EW model boundaries should
agree exactly with the EW boundaries of the input elevation data. 
* If the model extent does not cover 360 degrees of longitudde, a transmissive
boundary is used (but land is treated as reflective)

Then compile and run the model with (example here using 6 openmp threads):

    export OMP_NUM_THREADS=6
    make -B -f make_generic_model
    ./generic_model model_global_4min.in

If you have are problems compiling or linking, note that gdal and netcdf have
to be installed, and in particular the latter must be compiled with the same
fortran compiler as used to compile this code. If you need to use custom
gdal/netcdf installs, then you can redefine variables mentioned in
[../../src/src_standard_compiler_var](../../src/src_standard_compiler_var)
(don't edit the latter, just redefine them in the makefile
[make_generic_model](make_generic_model)). For more information on doing this,
see the NCI relevant example in [../../test/unit_tests](../../test/unit_tests).
