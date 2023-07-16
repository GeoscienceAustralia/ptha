# R interface to the TPXO program `predict_tide`

This code provides a basic R interface to the [TPXO tidal model](https://www.tpxo.net/global). It relies on a separate install of the [OPTS software](https://www.tpxo.net/otps) (note we use the version that takes binary input data - not the netcdf version). To run this you'll also need to register and download model files for TPXO9v5a [as described here](https://www.tpxo.net/tpxo-products-and-registration). 

The code is slightly modified from Gareth Davies' [stormwavecluster](https://github.com/GeoscienceAustralia/stormwavecluster) github repository, in particular [this interface to TPXO7.2](https://github.com/GeoscienceAustralia/stormwavecluster/tree/master/R/tpxo7.2). 

## Information on the codes

* [tidal_computations.R](tidal_computations.R) provides an example of usage.
* [OPTS_directory_name_TPXO9.R](OTPS_directory_name_TPXO9.R) defines paths to key TPXO files on your machine. **Users will need to edit this after installing TPXO on their machine**.
    * [OPTS_directory_name.R](OTPS_directory_name.R) shows an earlier example using TPXO7.2 (but is not otherwise used herein).
* [predict_tide.R](predict_tide.R) provides the main computational routines.

## Testing

Once you've got everything installed and updated [OPTS_directory_name_TPXO9.R](OTPS_directory_name_TPXO9.R), you can start `R` from inside the [test](test) directory and run:
```
    Rscript test_tides.R
```
This should print `PASS` and make a plot comparing predicted and observed tides at a site near Coffs Harbour in NSW, Australia (with reasonable agreement).
