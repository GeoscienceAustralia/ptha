# DEM blending

This script makes a raster that smoothly transitions between 2 other input rasters: one which is high-resolution but spatially limited, and the other which is coarse-resolution but spatially extensive. 

The idea is that where the high-resolution data exists, we use that. Away from the high-resolution data, we transition smoothly to the coarse-resolution data. 

How should the transition be determined? This script does it by creating contour lines on the coarse-resolution dataset, with line-depths defined manually to be denser in shallow water, and coarser in deep water. We then calculate the "error" of the coarse-resolution data (i.e. difference with the fine-resolution data) at points spaced along these contour lines. Where the fine-resolution data does not exist, we interpolate the error along the contour lines, but assume the error transitions to zero at sufficiently long distances from known error points. Finally we interpolate the "error surface" using the points defined above, with triangulation, and use this to correct the coarse-resolution data. 

In general the scripts parameters will need case-specific adjustment for new applications. Also the data is not provided. While it may be useful the script is purely for demonstration purposes, and there is no guarentee it will work well in your application!

The script assumes you can run `gdalwarp` via a call to `system` (which is true when gdal is installed on a linux system), and also that you've installed rptha as [described here](../../README.md). 
