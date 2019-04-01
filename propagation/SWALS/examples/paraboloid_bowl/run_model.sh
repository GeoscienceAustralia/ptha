# Clean existing binary
rm ./paraboloid_bowl
rm -r ./OUTPUTS
# Build the code
make -B -f make_paraboloid_bowl > build_outfile.log
# Run the job
./paraboloid_bowl > outfile.log
# Plot it and report tests
Rscript plot.R
