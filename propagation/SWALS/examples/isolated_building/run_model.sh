# Clean existing binary
rm ./model outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_model > build_outfile.log
# Run the job
./model > outfile.log
# Plot it and report tests
Rscript plot_results.R testcase
