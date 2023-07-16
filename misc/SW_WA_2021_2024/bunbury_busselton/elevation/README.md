# Information to set the model elevation

* [make_swals_elevation_files_preference_list.R](make_swals_elevation_files_preference_list.R) is used to create a list of files (in preference order) that SWALS uses to fix the elevation

* [convert_polygon_subdirs_to_csv.R](convert_polygon_subdirs_to_csv.R) is used to convert shapefiles to csvs. In practice I use:
    * Data in the folder `initial_stage_40cmAHD` to define a polygon inside the Vasse estuary where the initial stage is set to 40cm. This region is behind floodgates, and various reports indicate that the water level is controlled to be 40cm AHD (e.g. Shane Martin's 2014 study on storm surge in Busselton).
    * Data in the folder `bridges_to_remove` to define a set of polygons where the elevation will be clipped to a maximum of zero (but can be less). Initially I used this to remove some bridges that were still in the DEM. Later I added other polygons that don't related to bridges, but other elevation artefacts that the coast where zeroing the elevation is sensible (e.g. eroding coastlines where some of the older datasets suggest high topography that I want to remove).
