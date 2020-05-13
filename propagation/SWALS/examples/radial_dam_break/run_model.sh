export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./radial_dam_break ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_radial_dambreak > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./radial_dam_break > outfile.log"
# Plot it and report tests
Rscript plot_results.R
