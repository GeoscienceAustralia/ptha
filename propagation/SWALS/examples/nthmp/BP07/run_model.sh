export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./monai
rm -r ./OUTPUTS
# Build the code
make -B -f make_monai > build_outfile.log
# Run the code
eval "$OMP_RUN_COMMAND ./monai > outfile.log"
# Plot and report tests
Rscript plot.R
