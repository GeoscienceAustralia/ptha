Here are some scripts for making hazard points. They are not in 'package' form,
and they rely on input datasets and tools which are not provided
here. The requirements include a linux-like GDAL install. 

Other methods can be used to make hazard points -- this is just one example

Key code is:
make_hazard_pts.R -- for automatically making hazard points
contour_util.R -- useful GIS functions
point_util.R -- useful functions
country_area_check.R -- check on area of countries (used when preparing
   ISLAND_CLIP_LAYER so we don't cut out small island states)
