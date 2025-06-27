This folder contains elevation data used in the tsunami model. As the full set of files is large, it is provided as a 15 Gb bzip2 compressed tar archive for separate download here: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/elevation.tar.bz2

When the latter folder is extracted (e.g. using `tar -jxf elevation.tar.bz2`), the current directory should contain the following subfolders and README.md (which provides more information on the elevation data):
    derived_for_model
    derived_for_model_NWWA_update
    derived_for_model_SWWA_update
    orig
    README.md

Although SWALS internally uses the above datasets to create the elevation seen by the model, you may be more interested to directly see the elevation in SWALS. For convenience it can be downloaded here: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/full_resolution_model_elevation.tar.bz2
