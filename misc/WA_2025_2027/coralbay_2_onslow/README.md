# Inundation PTHA for Coral Bay to Onslow

Part 2 of phase 4 of the Western Australia Tsunami Inundation Modelling Project with DFES.

The project sets up the inputs (inverts, breakwalls, gauges, elevation and multidomains) for SWALS.
Once the inputs are ready, compile SWALS with it's [makefile](./swals/makefile) and run from there. 
The project uses [one module file](./modules_R_431.sh) to source for running R scripts,
and another [one](./modules_ifx_2025_llvm.sh) to compile SWALS.
