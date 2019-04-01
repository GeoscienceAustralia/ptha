# Clean existing binary
rm ./benchmark_problem1
rm -r ./OUTPUTS
# Build the code
make -B -f make_benchmark_problem1 > build_outfile.log
# Run the job
./benchmark_problem1 > outfile.log
# Plot it and report tests
Rscript plot_results.R
