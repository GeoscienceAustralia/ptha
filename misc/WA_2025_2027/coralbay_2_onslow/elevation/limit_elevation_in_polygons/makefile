# Create polygons/values defining upper limits on the elevation.
#
# The SWALS model will ensure the elevation inside these polygons does not
# exceed the specified value.
#
# First make a polygon shapefile(s) with field "limit" defining the upper limit 
# of the elevation in that polygon. Then Point hard coded paths in the R file to
# it. Then run as 
#
# make

.PHONY: limit-elevation

limit-elevation:
	. ../../modules_R_431.sh && \
	Rscript create_elevation_limit_csv.R
