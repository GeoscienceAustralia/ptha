# Make walls to burn into the SWALS elevation
---------------------------------------------

Quite a few parts of Tongatapu are built around wetlands, with the roads acting as flow barriers. It is important to not miss these features in the model DEM, and this is most reliably achieved by burning the features into the grid.

We burn walls into the DEM to represent:
* Coastal seawalls
* Breakwaters that might otherwise not be properly captured
* Roads that act as barriers to the flow
* A broad ridge [not super well defined] in the main town.

The code and data here is used to create XYZ files with the max-elevation along these linear features. The max-elevation is defined by searching the elevation around each point, in a specified buffer distance (see the code for details).

Run the code with `Rscript make_breakwalls.R`


