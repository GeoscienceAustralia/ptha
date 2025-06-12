This folder contains analysis of the models and data. To run the code, it is necessary for all the code in `../swals` to have been run as described therein.

The main analysis is in the `tsunami_size` directory.

## Details

Once all the code in `../swals/` has been run, this folder will contain several folders with names beginning with `UNTAR_PROCESSING*`. They were created by running `../swals/extract_key_outputs_from_tarred_multidomain_3.R`, as explained in the README in that directory. 

The script `make_directory_structure.R` is then used to copy the contents of all these folders into a single folder `REDUCED_OUTPUTS` for further analysis.

At that point one can go inside the `tsunami_size` directory and run some calculations.
