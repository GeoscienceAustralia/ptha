This directory contains a script to create x/y/z data representing the bottoms of channels.

These can be 'burned' into the elevation grid in SWALS to ensure that flow paths
are represented, even if sometimes the model is too coarse to reliably
capture the the minima.

The code is very similar to that in ../breakwalls

Run with

    Rscript make_inverts_perth.R
    Rscript make_inverts_list_have_fullpath.R
