# Clean existing binary
rm ./dam_break
rm -r ./OUTPUTS
# Build the code
make -B -f make_dam_break > build_outfile.log
# Run the job
./dam_break > outfile.log
# Plot it and report tests
Rscript plot.R
