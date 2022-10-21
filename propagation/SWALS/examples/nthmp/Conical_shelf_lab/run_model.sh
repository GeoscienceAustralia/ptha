
export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./model
rm -r ./OUTPUTS
# Build the code
make -B -f make_model > build_outfile.log
# Run the code
eval "$OMP_RUN_COMMAND ./model '' > outfile_omp.log"
# Plot and report tests
Rscript plot.R
