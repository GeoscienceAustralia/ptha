This folder contains the main modelling runs for each source-zone.
-------------------------------------------------------------------

To make files for a new source-zone, you need to copy the TEMPLATE folder to a
new folder with the name **source_name** (e.g. for the Sunda Arc, the folder is
called 'sunda').

Then the code is ready to run. It assumes that ../DATA/SOURCE_CONTOURS has a
file named **source_name**.shp, and also ../DATA/SOURCE_DOWNDIP_LINES has one
called **source_name**_downdip.shp. (For example, for the 'sunda' case, they
would be 'sunda.shp' and 'sunda_downdip.shp' respectively)

Since some source-zones have few unit sources (or might have a number of unit
sources that is not close to a multiple of 16), we also have a 'run_16.sh'
script in the current directory, which can run jobs from any source-zone. This
means we are less likely to set off a job on 1 node which cannot find 16 models
to run.

