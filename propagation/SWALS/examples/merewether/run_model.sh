export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./merewether_example outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_merewether_example > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./merewether_example > outfile.log"
# Plot it and report tests
Rscript plot.R
