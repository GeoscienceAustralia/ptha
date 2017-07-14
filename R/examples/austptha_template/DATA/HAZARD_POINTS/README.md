Create hazard points
--------------------

Here we create points at which the tsunami model time-series will be stored.

The model uses the GA_250m DEM where available, and GEBCO_2014 elsewhere. However,
we 'patch' the GA250 DEM in some locations in the SW pacific + Macquarie Island, where
GEBCO is clearly better.

We would like:
- Australian points at 20m water depth, 20km spacing
- Australian points at 100m water depth, 20km spacing
- Australian points at 1000m water depth, 20km spacing (to see if we get better convergence)
- Gridded points around Australia, to make it easier to use them for structured grid modelling
- DART buoy locations
- Global points at 100m water depth, 50km spacing.
- We also need to manually add some points in Australia, wherever the algorithm 'misses' places.

Running
-------

Assuming that you have already made the merged DEM's (see [../ELEV][../ELEV]),
and that R is installed with required packages, then do:

    source make_all.sh
