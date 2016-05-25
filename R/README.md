rptha
-----
-----

rptha is the main workhorse R package for doing probabilistic tsunami hazard
assessment (PTHA). 

The main steps in the PTHA methodology that is supported by rptha are:

* Define tsunami sources. Currently rptha only supports earthquake tsunami
sources. The user needs to provide contours defining the interface on
which earthquake slip occurs. 

* Discretize each earthquake source-zone into a grid of unit-sources, and
for each compute the tsunami initial condition using the Okada solution
combined with Kajiura filtering. This produces a (static) tsunami initial
condition for each source-zone.

* Create a synthetic catalogue of earthquake events from linear
combinations of the unit sources. Currently rptha supports uniform slip
earthquakes with dimensions determined to (approximately) agree with the
scaling relations of Strasser et al (2010). 

* Compute the tsunami associated with each earthquake event. rptha does not
provide a tsunami solver, so other software is required for this step. In
realistic PTHA applications this is the most computationally demanding part of
the process. 
  * If a linear propagation code is used then the user can optionally compute
the tsunami for each unit-source separately and then combine then later. This
is relatively efficient, but only theoretically valid for tsunami with a
sufficiently small amplitude-to-depth ratio (or equivalently, with sufficiently
small velocities).
  * If nonlinear solvers are required, then the used must first create the 
initial conditions for each event (by linearly combining the unit source
initial conditions), and then run each through the nonlinear propagation code. 

* Assign an mean annual rate to each event in the earthquake catalogue. This
is based on seismic moment conservation principles, and requires the user to
specify the source-zone convergence rate and various other parameters
controlling seismicity. Uncertainties can be accounted for using a logic tree.

* Compute the tsunami associated with each event in the catalogue by
linearly combining the unit-source tsunami.

In a typical application you would write scripts to call routines in rptha.
Examples of scripts for common tasks are provided here.


Installation
------------

To build rptha, go inside 'rptha', start R, and do:

    source('build_package.R')

This will make an R package file in the directory above the package, which can be installed on the command line with:

    R CMD INSTALL rptha_XXXXX.tar.gz

where the XXXX are adapted to match the file name (and on Ubuntu, you would add 'sudo' to the start).

If the above fails because you are missing packages, then try running this (from inside R) prior to the install to get the required packages:

    install.packages(c('sp', 'rgdal', 'rgeos', 'FNN', 'raster', 'minpack.lm', 'geometry', 'geosphere', 'rgl', 'testthat', 'devtools'))

If you are separately reading ncdf rasters, you will also need the 'ncdf4' package.

source_contours_2_unit_sources
------------------------------

This contains
[example_code](source_contours_2_unit_sources/produce_unit_sources.R) to make
tsunami unit sources from source contours and a
[tutorial](source_contours_2_unit_sources/tutorial.md) on its usage.


make_hazard_points
------------------

Example scripts to make hazard points (i.e. offshore points where the tsunami tide-gauges are recorded)


combine_tsunami_sources
-----------------------
Example script to convert from tsunami unit sources to tsunami event offshore
wave heights (assuming the URSGA solver was used for tsunami propagation)
