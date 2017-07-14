For every source-zone (with name *source_name*), this folder should contain
a shapefile named *source_name*_downdip.shp (and corresponding .prj, .shx, .dbf files). 

The shapefile should define down-dip breaks between unit-sources. It can be
made from the sourcezone contours, using the script
[../make_initial_downdip_lines.R](../make_initial_downdip_lines.R).  Outputs
from the latter script will typically be manually edited somewhat, to remove
any kinks in the downdip lines, and ensure they are reasonably spaced, etc.
