export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./model outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_model > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./model > outfile.log"
# Plot it and report tests
Rscript plot_results.R testcase
