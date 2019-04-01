# Clean existing binary
rm ./monai
rm -r ./OUTPUTS
# Build the code
make -B -f make_monai > build_outfile.log
# Run the code
./monai > outfile.log
# Plot and report tests
Rscript plot.R
