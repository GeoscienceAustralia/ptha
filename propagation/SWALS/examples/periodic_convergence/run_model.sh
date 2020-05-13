export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Run the convergence test, so the "test_convergence.R" script can be run
rm -rf ./OUTPUTS model
make -B -f make_model > build_log.log
eval "$OMP_RUN_COMMAND ./model 1"
eval "$OMP_RUN_COMMAND ./model 2"
eval "$OMP_RUN_COMMAND ./model 4"
eval "$OMP_RUN_COMMAND ./model 8"
#./run_model 16
Rscript test_convergence.R
