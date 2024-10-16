This directory contains a script to create x/y/z data representing the tops of
breakwalls.

These can be 'burned' into the elevation grid in SWALS to ensure that barriers
to flow are represented, even if sometimes the model is too coarse to reliably
capture the breakwall maxima.

Run with

    Rscript make_breakwalls.R
    Rscript make_breakwalls_list.R


__NOTE__: In this model we create a breakwall to close the Bunbury floodgate, but deliberately exclude it from the default list of breakwall files (generated by [make_breakwalls_list.R](make_breakwalls_list.R)). This is so we can control whether the breakwall is closed from inside the [SWALS model setup code](../swals/model_local_routines.f90). 

__NOTE__: In April 2024 I added in patches for recent work on the Vasse Diversion drain. To avoid overwriting the older model files, instead `make_breakwalls_list.R` was edited to produce an updated output file that is used for that specific model.
