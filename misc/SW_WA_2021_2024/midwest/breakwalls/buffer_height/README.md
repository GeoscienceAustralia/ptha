This directory contains a script to create x/y/z data representing the tops of breakwalls.
The elevation is set using the maximum elevation within some radius (buffer) of points along the wall.

These can be 'burned' into the elevation grid in SWALS to ensure that barriers to flow are represented, even if sometimes the model is too coarse to reliably capture the breakwall maxima.

Run with

    Rscript make_breakwalls.R


__NOTE__: In this model we create a breakwall to close the Bunbury floodgate, but deliberately exclude it from the default list of breakwall files (generated by [../make_breakwalls_list.R](../make_breakwalls_list.R)).
Another model separately enforced it depending on a parameter in the code.
