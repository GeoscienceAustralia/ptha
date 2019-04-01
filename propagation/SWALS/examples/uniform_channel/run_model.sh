# Clean existing binary
rm ./uniform_channel ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_uniform_channel > build_outfile.log
# Run the job
./uniform_channel > outfile.log
# Plot it and report tests
Rscript plot_results.R
