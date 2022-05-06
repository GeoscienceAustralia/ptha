# Single grid shallow water model in spherical coordinates.

This code is often convenient for simulating global-scale tsunami propagation.

It enables a variety of shallow water solvers to be used on a single grid in spherical coordinates.

The boundary conditions are transmissive, unless the longitude-extent is 360
degrees. In the latter case, the model uses a combination of periodic east-west
boundaries and reflective north-south boundaries. 

It you wish to create more complex models (e.g. with nesting, or different
boundary conditions) it is better to start from [../periodic_multidomain](../periodic_multidomain).

## To compile and run the code

* Copy one of the '\*.in' files to a new file (e.g. 'my_model.in'), and edit it
to ensure that the input_elevation_raster and input_stage_raster exist on your
filesystem, and that the model extents are as desired. The input rasters should
give the initial elevation and stage in lon-lat coordinates. 
* The elevation raster should cover the desired model domain, but it can
have a larger north-east-south-west extent than the model. The
model will get the data it needs using bilinear interpolation. So for example, 
you can run a small area model using a global DEM, without editing the DEM.
Furthermore it is not essential for the model resolution to match the data
resolution.  However the model does not correct for the periodicity of
longitude when extracting raster data. For instance, if the input raster has
east-west extent [-180, 180], then the tsunami model cannot have east-west
extent like [0, 360], but it could have [-180, 180], or [-180, 12], or [-50,
100], etc. 
* The stage raster can be either smaller or larger than the desired model domain. The
model will extract stages from this raster where possible, using bilinear
interpolation. It will use a value of 0 elsewhere. 
* If the model extent covers 360 degrees of longitude (e.g. [0, 360] or [-180, 180] or [-40, 320]), 
then east-west periodic boundary conditions are used, with reflective
north-south boundaries.  In this case, the EW model boundaries should agree
exactly with the EW boundaries of the input elevation data. 
* If the model extent does not cover 360 degrees of longitudde, a transmissive
  boundary is used.

Then compile and run the model with (e.g.):

    make -B -f make_generic_model
    OMP_NUM_THREADS=6 ./generic_model my_model.in

The above code uses 6 openmp threads -- you can use any number relevant to your system.
It can also be advantageous to control the process affinity (e.g. include
"OMP_PROC_BIND=true" in the model run command), although this depends on your
hardware and what other jobs are being run.

## The tests

The test-suite runs this code in the source-region of the 2011 Tohoku tsunami,
using a number of different solver options. The earthquake source was selected
from among randomly generated sources, so is not a precise representation of
the real forcing of Tohoku.

We perform a regression test by checking that the difference between the
modelled and observed wave is very close to previously computed values (which
are distinct for each solver). The idea is to detect any changes in the solver
predictions; the model is not supposed to give precise agreement to the data.

We also check that for models using a rise-time (i.e. applying the earthquake
forcing over a finite time), the solutions are equivalent to temporally
smoothed solutions from models using an instantaneous forcing. 
* This is true analytically for the linear shallow water equations. 
* In the specific situation tested, it is also a good approximation for the nonlinear shallow water equations.
