rptha
=====

rptha is the main workhorse R package for doing probabilistic tsunami hazard
assessment (PTHA). 

Background
----------

PTHA methodologies can be quite diverse depending on the level of detail
required. The main steps in PTHA methodology supported by rptha are outlined
below. In a typical application you would write scripts to call routines in rptha.
Examples of scripts for common tasks are provided here. (Note that the code can
optionally be used just to [generate tsunami initial
conditions](examples/source_contours_2_unit_sources) as raster output files, which may
then form part of some other workflow):

* Define tsunami sources. Currently rptha only supports earthquake tsunami
sources. The user needs to provide contours defining the interface on
which earthquake slip occurs. 

* Discretize each earthquake source-zone into a grid of unit-sources, and
for each compute the tsunami initial condition using the Okada solution
combined with Kajiura filtering. This produces a (static) tsunami initial
condition for each unit-source. See [source_contours_2_unit_sources](examples/source_contours_2_unit_sources).

* Create a synthetic catalogue of earthquake events from linear
combinations of the unit sources. Currently rptha supports uniform slip
earthquakes with dimensions determined to (approximately) agree with the
scaling relations of Strasser et al. (2010). See the function
`get_all_earthquake_events` in the rptha package. Modifications are required to
treat non-uniform slip earthquakes. For an example of making tsunami
initial conditions for complex earthquake scenarios, see
[here](examples/combine_tsunami_sources/combine_tsunami_sources.R). In the
[2018 Australian PTHA](http://dx.doi.org/10.11636/Record.2018.041) we treated non-uniform slip
earthquakes by associating a set of them with a `parent' uniform-slip
earthquake; see codes associated with that project [here](examples/austptha_template) and
papers on the approach [here](https://link.springer.com/article/10.1007/s00024-019-02299-w) 
and [here](https://doi.org/10.1093/gji/ggz260)

* Compute the tsunami associated with each earthquake event. In
realistic PTHA applications this is the most computationally demanding part of
the process. In general, the user is encouraged to provide their own
tsunami propagation code, although one option is provided in [propagation](../propagation).
  * If a linear propagation code is used then the user can optionally compute
the tsunami for each unit-source separately and then [combine them
later](examples/austptha_template/SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS). This is relatively
efficient, but only theoretically valid for tsunami with a sufficiently small
amplitude-to-depth ratio (or equivalently, with sufficiently small velocities).
The tsunami propagation results are generally stored at offshore points, and
scripts to make these are [here](examples/make_hazard_points).
  * If nonlinear solvers are required, then the user must first create the 
initial conditions for each event (by linearly combining the unit source
initial conditions), and then run each through the nonlinear propagation code. 
See [here](examples/combine_tsunami_sources/combine_tsunami_sources.R) for an example
of linearly combining tsunami water surface deformations.

* Assign an mean annual rate to each event in the earthquake catalogue. This
is based on seismic moment conservation principles, and requires the user to
specify the source-zone convergence rate and various other parameters
controlling seismicity. Uncertainties can be accounted for using a logic tree.
See the the example script [here](examples/event_rates/single_source_rate_computation.R),
and a relatively complex global-scale application for the 2018 Australian PTHA
[here](examples/austptha_template/EVENT_RATES).


Installation from source
------------------------

### Standard install

These notes apply to linux. 

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

### Windows 

On Windows, building R packages is more challenging, and the developer does not 
generally build the package on Windows machines. Users attempting to do this
should first seriously consider installing a linux virtual machine instead. If 
you still wish to attempt a windows build, then familiarise yourself with the
general issues in building R packages for windows
(http://cran.us.r-project.org/bin/windows/Rtools/), and proceed along the above
lines once appropriate dependencies are installed.

### Install on NCI

* This may differ from the above instructions -- the approach has varied over
  the years. See some detailed notes in the [install](install) directory


Usage
-----
To learn to use rptha, you can look at the [examples](examples) directory, as
well as the package documentation



