# Elevation data used in our modelling 
--------------------------------------

## How to get the data

The data that resides in this folder has been compressed and posted to the NCI THREDDS Server for download here: http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/elevation.tar.bz2

You need to download that file and extract it (`tar -zxf elevation.tar.bz2`), and perhaps move folders around to ensure the current folder contains sub-folders like this:

    derived_for_model  orig

as well as the current README file.

## What it is

Following the download above, this folder will contain two sub-folders with elevation data. The elevation data was created for our nearshore tsunami model testing in Australia.

See the sub-folders for information on the underlying data sources and associated licencing information. For the most part they were derived from data available (during mid 2019-early 2020) from [ELVIS](https://elevation.fsdf.org.au/) and [AODN](https://portal.aodn.org.au/search) and [Ausseabed](http://www.ausseabed.gov.au/) and the [WA Department of Transport](https://catalogue.data.wa.gov.au/dataset/composite-surfaces-multibeam-lidar-laser), with much of the NSW data being originally collected by the [NSW Department of Planning, Industry and Environment](https://www.dpie.nsw.gov.au/) which includes the (now renamed) NSW Office of Environment and Heritage. The coarsest-scale grid is based on the 2014 version of [GEBCO](https://www.gebco.net/) and the [2009 Australian Bathymetry and Topography Grid](https://data.gov.au/data/dataset/australian-bathymetry-and-topography-grid-june-2009), and we use the same merger of these datasets as was derived for the [2018 Australian Probabilistic Tsunami Hazard Assessment](http://dx.doi.org/10.11636/Record.2018.041). 

The underlying datasets are almost all distributed under the [Creative-Commons-BY-4.0 licence](https://creativecommons.org/licenses/by/4.0/), except for derivatives of the [VCDEM](https://vmdp.deakin.edu.au/geonetwork/srv/eng/metadata.show?uuid=8d3ccf63-ee85-41cd-917e-933624a50b2e) which is released under the [Creative Commons BY NC 4.0 licence](https://creativecommons.org/licenses/by-nc/4.0/), and so similar licences apply to our derived grids as discussed in the README files.

The data herein is provided for scientific transparency. It was created to run the model code in [../swals](../swals), using various products discussed in the sub-folders. We make no guarentees about the accuracy or fitness for any purpose, aside from repeating the computations in our study, which was the sole focus of our elevation data processing.

## How it is used in our study 

SWALS sets the tsunami model's elevation using a heirachy of rasters with a given preference order. To get an elevation value at any point, the code scans through the rasters until it finds the highest-preference raster with non-missing data at the point of interest, and then does bilinear interpolation to get the desired value. 

The raster preference order is defined in the model application code [../swals/model_local_routines.f90](../swals/model_local_routines.f90) in a  subroutine named `setup_elevation`.

This means that there is no fixed relationship between the raster extents/resolution and the model. In contrast, many tsunami models use the elevation grid to define the model grid. 

Further, by construction we don't use parts of the raster that are covered by higher-preference raster data.

