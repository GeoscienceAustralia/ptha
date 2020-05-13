export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./circular_island_testcase
rm -r ./OUTPUTS
# Build the code
make -B -f make_circular_island > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./circular_island_testcase > outfile.log"

# Plot it and report tests
Rscript plot.R
