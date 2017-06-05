create_files_for_all_simulations.R -- used to make directory structures for model runs

run_16.sh -- qsub this to run 16 of the models which are not already slated for running

run_permute_netcdf_gauges.sh -- qsub this to permute the dimension of the tide-gauge netcdf output files, to ensure fast access to single-station data [which is extremely convenient for later analysis]
