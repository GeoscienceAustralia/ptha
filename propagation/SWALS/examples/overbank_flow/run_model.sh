export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./model ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_model > build_outfile.log
# Run the job with a north/south aligned channel.
eval "$OMP_RUN_COMMAND ./model NS_aligned > outfile.log"
eval "$OMP_RUN_COMMAND ./model EW_aligned > outfile.log"
# Plot it and report tests
Rscript plot_results.R
