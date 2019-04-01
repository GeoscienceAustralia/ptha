# Clean existing binary
rm ./uniform_slope ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_uniform_slope > build_outfile.log
# Run the job
./uniform_slope > outfile.log
# Plot it and report tests
Rscript plot_results.R
