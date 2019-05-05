# Clean existing binary
rm ./nesting_reflection
rm -r ./OUTPUTS
# Build the code
make -B -f make_nesting_reflection > build_outfile.log
# This script runs the job, plots, and reports tests
Rscript plot.R
