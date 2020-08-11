# Elevation data used in our modelling 
--------------------------------------

## How to get the data

The data that resides in this folder has been compressed and posted to the NCI THREDDS Server for download here: http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/elevation.tar.bz2

You need to download that file and extract it so that the current folder contains sub-folders like this:

    derived_for_model  orig

as well as the current README file.

## What it is

Following the download above, this folder will contain two sub-folders with elevation data that was created for some nearshore tsunami model testing in Australia.

See the sub-folders for information on the underlying data sources. For the most part they were derived from data available from [ELVIS](https://elevation.fsdf.org.au/) and [AODN](https://portal.aodn.org.au/search) and [Ausseabed](http://www.ausseabed.gov.au/).

The data herein is provided for scientific transparency. It was created to run the model code in [../swals](../swals), using various products discussed in the sub-folders. We make no guarentees about the accuracy or fitness for any purpose, aside from repeating the computations in our study.

## How it is used in our study 

SWALS sets the tsunami model's elevation using a heirachy of rasters with a given preference order. To get an elevation value at any point, the code scans through the rasters until it finds the highest-preference raster with non-missing data at the point of interest, and then does bilinear interpolation to get the desired value. 

The raster preference order is defined in the model application code [../swals/model_local_routines.f90](../swals/model_local_routines.f90) in a  subroutine named `setup_elevation`.

This means that there is no fixed relationship between the raster extents/resolution and the model. In contrast, many tsunami models use the elevation grid to define the model grid. 

Further, by construction we don't use parts of the raster that are covered by higher-preference raster data.

