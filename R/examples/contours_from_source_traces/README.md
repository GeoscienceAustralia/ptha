# Geometric method to make source contours

This code creates 'geometrically defined' source contours, based on an input
shapefile defining the source-trace (i.e. the trench).

It is useful when you know the trench location, and are able to estimate the dip and
other key fault parameters, but don't have higher quality information (such as SLAB1.0).

*Note the contours will generally require manual editing to remove kinks and other artefacts (e.g with GIS software like QGIS). This is especially true for the trench contour, which tends to show 'jagged' artefacts due to rasterization. We suggest this code be used as a starting point in the contour definition only.*

## Shapefile format
Examples of the input shapefile are in the directories [timor](timor), [tanimbar](tanimbar), and [tolo_thrust](tolo_thrust). These are examples only, and we do not claim that these inputs lead to a realistic representation of these source-zones. 

Notice how the traces are *ordered along-strike, with end-points of each trace segment matching*. This is a requirement of the input shapefiles. 

The shapefile attribute tables contain some important parameters that are used to define the source:
1. `Max_depth`: The maximum depth of the source-contours in km. 
    * This must be identical on each individual source-zone segment. 
2. `Concavity`: Either 'Linear' or 'Up' or 'Down'. 
    * Depending on the `Concavity`, you need to provide different dip values.
    * In the case of 'Linear', we assume the fault dips linearly. 
        1. In this case, there needs to be an attribute named `Dip` which gives the dip of the fault in degrees 
    * In the case of 'Up' or 'Down', we assume the fault dips parabolically, either concave-up or concave-down. 
        1. In this case you need to provide an attribute named `Dip_0` giving the dip (degrees) at the trench, and another attribute with a name like `Dip_10` giving the dip (degrees) at the depth specified after the underscore (i.e. 10km the latter example -- but you could replace the `10` with any other depth that was less than `Max_depth` and greater than zero). 
        2. These dip values need to respect the concavity. For example if the concavity is 'Up', then the dip at the trench must be greater than the deeper dip; and conversely for concave down contours. Further, the parabola needs to intersect `Max_depth` (not necessarily true for 'Up' contours, see discussion below).

## Running the code
Once the source-trace shapefiles are created, you should open the script [run_convert_traces_to_contours.R](run_convert_traces_to_contours.R), and edit the variable `source_traces` to point to the traces you'd like to convert to contours (this can contain one or more files). One that has been done, you run the code from the command line with:

    Rscript run_convert_traces_to_contours.R

The above command will make contour shapefiles in the CONTOURS directory, which should be further checked and edited to remove jagged artefacts.

## Troubleshooting
The code tries to provide useful error messages. One problem I have seen repeatedly is that, in the case of a `Concavity`='Up' parabolic fault, the user may accidently specify parameters that do not allow the source to reach `Max_depth` (i.e. the parabola turning point occurs at a depth less than `Max_depth`). If that happens, you'll get an error message something like this, which tells you to edit the input data:

    OGR data source with driver: ESRI Shapefile
    Source: "tanimbar", layer: "tanimbar"
    with 3 features
    It has 10 fields
    Integer64 fields read as strings:  id Dip_0 Dip_20 Max_depth
    Error in extend_trace_to_depth_contours(mytrace, maxdepth = maxdepth,  :
      Based on your inputs, the parabolic profile has a turning 
     point near to depth=20.470049629778 which is < max_depth=25
     This suggests a data input error, and is not allowed. You should either
     reduce the dip at the trench, or increase the deeper dip value, or
     decrease the maxdepth. The parameters were:
     mindepth = 0
     dip@mindepth = 30
     dip_depth = 20
     dip@dip_depth = 5
     Leading to a parabola like " depth = mindepth + beta * distance + alpha * distance^2 "
     with parameters alpha=-0.00407098833859726 beta=0.577350269189626
     where distance=0 at depth=mindepth
