This folder contains the main modelling runs for each source-zone.
-------------------------------------------------------------------

To make files for a new source-zone, you need to copy the TEMPLATE folder to a
new folder with the name **source_name** (e.g. for the Sunda Arc, the folder is
called 'sunda').

Then the code is ready to run. It assumes that ../DATA/SOURCE_CONTOURS has a
file named **source_name**.shp, and also ../DATA/SOURCE_DOWNDIP_LINES has one
called **source_name**_downdip.shp. (For example, for the 'sunda' case, they
would be 'sunda.shp' and 'sunda_downdip.shp' respectively)

## Extras

### run_16.sh

The current directory has a [run_16.sh](run_16.sh) script which can run tsunami
propagation models (after they have been set-up by modifying the TEMPLATE
folder and following the instructions in the TSUNAMI_UNIT_SOURCE sub-folder).

Run it with

    qsub run_16.sh

The [run_16.sh](run_16.sh) script is similar to the script
[./TEMPLATE/TSUNAMI_UNIT_SOURCE/run_16.sh](./TEMPLATE/TSUNAMI_UNIT_SOURCE/run_16.sh),
but a key difference is that **it will run jobs from any source-zone**.

This is nice, because consider a situation in which we have run all but a few
(say 3) unit-sources for a particular source-zone. If we ran the run_16.sh
script from within TSUNAMI_UNIT_SOURCE, then we would only end up running 3
jobs on the compute node (instead of 16). However, if we used the run_16.sh
script in the current directory, the code would look for more tsunami
propogation runs to put on the node (i.e. in other source-zones).  This is more
likely to be computationally efficient. 

### checkruns.R

This code can be used to run the completion checks for the propagation models,
over multiple source-zones. In [checkruns.R](checkruns.R), modify the variable
`source_zones` which points to the folders that will be checked, and then
run it from within R as

```r
    source('checkruns.R')
```


### move_nc_files_and_replace_with_symbolic_links.R

Code to move netcdf files in */TSUNAMI_EVENTS/*.nc onto NCI's gdata filesystem,
and replace the local files with symbolic links. To reduce the risk of
accidental deletions, the script renames the existing nc files (by appending
BACKUP to their name). These files can be manually deleted later if required.

