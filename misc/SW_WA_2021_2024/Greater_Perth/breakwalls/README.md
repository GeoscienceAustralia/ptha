This directory contains a script to create x/y/z data representing the tops of
breakwalls.

These can be 'burned' into the elevation grid in SWALS to ensure that barriers
to flow are represented, even if sometimes the model is too coarse to reliably
capture the breakwall maxima.

Run with

    Rscript make_breakwalls.R
    Rscript make_breakwalls_list.R
