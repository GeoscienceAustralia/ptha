# Clean existing binary
rm ./radial_dam_break ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_radial_dambreak > build_outfile.log
# Run the job
./radial_dam_break > outfile.log
# Plot it and report tests
Rscript plot_results.R
