Merge the GA250 and GEBCO2014 dems
---------------------------------

This folder contains code to 'merge' the GA250 and GEBCO2014 DEMs 
onto a 1-arc-minute DEM for the PTHA tsunami model run.

Before running this script, you need to have run the scripts in 
[../GEBCO_2014_1m](../GEBCO_2014_1m) to downsample the GEBCO data and
adjust its longitudinal extent.

The GA250 dem is used in preference to GEBCO where available, except 
in a few obviously erronious places, where the GA250 DEM is 'patched'
by the code [patch_dem.R](patch_dem.R).


# How to run

## Step 1

Run the codes in [../GEBCO_2014_1m](../GEBCO_2014_1m)

## Step 2

Run the code to make an initial merger of the DEMs (note: the scripts may
require modification to point to the correct files).

    source make_merged_dem.sh

The resulting file uses GA250 where available, and GEBCO2014 elsewhere. 

## Step 3

Patch over a few obvious defects in the GA250 dem using the GEBCO DEM. 

    Rscript patch_dem.R

## [Step 4 -- this was done well after the PTHA18 publication]

It is useful to have a DEM that is '1' in areas below MSL, and '0' elsewhere.

    Rscript make_wet_or_dry_dem.R

This has a very small file size so is easy to distribute -- see [here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/DATA/wet_or_dry_gebco_ga250_dem_patched.tif)
