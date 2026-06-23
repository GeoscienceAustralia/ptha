# Convert max-stage outputs back to AHD

The SWALS model uses an adjusted elevation such that, at the coast, 0.0 m is equivalent to a local high tide.

That means the 'raw' max-stage outputs are all reported relative to a 'local high tide' datum.

To convert these outputs to AHD, we add in the tidal adjustment using the script here.

Run with 
```
Rscript convert_max_stage_outputs_to_AHD.R
```
