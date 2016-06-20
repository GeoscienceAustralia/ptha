rptha
=====

rptha is the main workhorse R package for doing probabilistic tsunami hazard
assessment (PTHA). 

PTHA methodologies can be quite diverse depending on the level of detail
required. The main steps in PTHA methodology supported by rptha are outlined
below. (Note that the code can optionally be used just to [generate tsunami
initial conditions](source_contours_2_unit_sources) as raster output files,
which may then form part of some other workflow):

* Define tsunami sources. Currently rptha only supports earthquake tsunami
sources. The user needs to provide contours defining the interface on
which earthquake slip occurs. 

* Discretize each earthquake source-zone into a grid of unit-sources, and
for each compute the tsunami initial condition using the Okada solution
combined with Kajiura filtering. This produces a (static) tsunami initial
condition for each unit-source. See [source_contours_2_unit_sources](source_contours_2_unit_sources).

* Create a synthetic catalogue of earthquake events from linear
combinations of the unit sources. Currently rptha supports uniform slip
earthquakes with dimensions determined to (approximately) agree with the
scaling relations of Strasser et al. (2010). See the function
`get_all_earthquake_events` in the rptha package. Modifications are required to
treat non-uniform slip earthquakes, but for an example of making tsunami
initial conditions for complex earthquake scenarios, see
[here](combine_tsunami_sources/combine_tsunami_sources.R).

* Compute the tsunami associated with each earthquake event. rptha does not
provide a tsunami solver, so other software is required for this step. In
realistic PTHA applications this is the most computationally demanding part of
the process. 
  * If a linear propagation code is used then the user can optionally compute
the tsunami for each unit-source separately and then [combine them
later](combine_tsunami_sources/combine_tsunami_gauges.R). This is relatively
efficient, but only theoretically valid for tsunami with a sufficiently small
amplitude-to-depth ratio (or equivalently, with sufficiently small velocities).
The tsunami propagation results are generally stored at offshore points, and
scripts to make these are [here](make_hazard_points).
  * If nonlinear solvers are required, then the user must first create the 
initial conditions for each event (by linearly combining the unit source
initial conditions), and then run each through the nonlinear propagation code. 
See [here](combine_tsunami_sources/combine_tsunami_sources.R) for an example.

* Assign an mean annual rate to each event in the earthquake catalogue. This
is based on seismic moment conservation principles, and requires the user to
specify the source-zone convergence rate and various other parameters
controlling seismicity. Uncertainties can be accounted for using a logic tree.
See the function `rate_of_earthquakes_greater_than_Mw` in the rptha package.

In a typical application you would write scripts to call routines in rptha.
Examples of scripts for common tasks are provided here.


Installation from source
------------------------

To build rptha, you firstly need to install the R packages that it depends on.
The usual way to get these is to start R, and do:

    install.packages(c('sp', 'rgdal', 'rgeos', 'FNN', 'raster', 'minpack.lm', 'geometry', 'geosphere', 'rgl', 'ncdf4', 'testthat', 'devtools', 'roxygen2'))

This will ask you to choose a mirror to download from. Just choose something that
is close to your location -- for example in Canberra, Australia, you can first select
'http mirrors' and then select the Canberra-Australia mirror. 

As the packages install, check the printed text to ensure it worked. If any packages
failed to install, use google to troubleshoot. The packages rely on various non-R
libraries already being installed on your system (e.g. gdal, netcdf, geos).
Experience suggests that the installation is straightforward on a standard
ubuntu desktop, but can be tricky on non standard environments (e.g. raijin on
the NCI). Note that the `rgl` install is not essential, as it is only used for
3d interactive graphics, so errors associated with this can optionally be
ignored.

Once the packages have installed, you should cd into the ptha/R/rptha directory, start
R, and then do:

    source('build_package.R')

This will make an R package file in the same directory as this README.md
(ptha/R/). That package can be installed on the command line with:

    R CMD INSTALL rptha_XXXXX.tar.gz

where the XXXX are adapted to match the file name (and on Ubuntu, you would add
'sudo' to the start). 

You can also give other people copies of the rptha_XXXXX.tar.gz file to install
themselves (provided they are running a similar environment to you, and have
the package dependencies installed).

source_contours_2_unit_sources
------------------------------

This contains
[example_code](source_contours_2_unit_sources/produce_unit_sources.R) to make
tsunami unit sources from source contours and a
[tutorial](source_contours_2_unit_sources/tutorial.md) on its usage.


make_hazard_points
------------------

Example scripts to [make hazard points](make_hazard_points/make_hazard_pts.R) (i.e. offshore points where the tsunami tide-gauges are recorded)


combine_tsunami_sources
-----------------------
[Example script](combine_tsunami_sources/combine_tsunami_gauges.R) to convert
from tsunami unit sources to tsunami event offshore wave heights (assuming the
URSGA solver was used for tsunami propagation)
