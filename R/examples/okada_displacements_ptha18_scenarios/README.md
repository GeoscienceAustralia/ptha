These scripts can be used to compute and store the 3D Okada displacement field
for unit-sources in the 2018 Australian Probabilistic Tsunami Hazard Assessment (PTHA18).

To run it, edit [config.R](config.R) for your problem. Principally this
involves specifing your sourcezone and desired region where the outputs are
stored. You also have to provide some input data from the PTHA18. See comments
in that script for further explanation.

Once that is done, the code can be executed with (e.g.):

    Rscript produce_okada3d_unit_sources.R

Currently it is setup to work on the kermadectonga2 source-zone, and extract deformation
in a region in the vicinity of Tonga.
