export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./benchmark_problem1
rm -r ./OUTPUTS
# Build the code
make -B -f make_benchmark_problem1 > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./benchmark_problem1 > outfile.log"
# Plot it and report tests
Rscript plot_results.R
