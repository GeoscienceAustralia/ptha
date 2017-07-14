Create hazard points
--------------------

Here we create points at which the tsunami model time-series will be stored.

The model uses the GA_250m DEM where available, and GEBCO_2014 elsewhere. However,
we 'patch' the GA250 DEM in some locations in the SW pacific + Macquarie Island, where
GEBCO is clearly better.

We would like:
- Global points at 100m water depth, 50km spacing.
- Australian points at 100m water depth, 20km spacing
- Australian points at 1000m water depth, 20km spacing (to see if we get better convergence)
- DART buoy locations
- (We also need to manually add some points in Australia, wherever the algorithm 'misses' places)
- (And we also add gridded points, to make it easier to use them for structured grid modelling)

All the above categories of points should have a short character label, to
reduce the risk that we mix up their locations later.

Running
-------

Assuming that you have already made the merged DEM's (see [../ELEV][../ELEV]),
and that R is installed with required packages, then do:

    source make_all.sh
