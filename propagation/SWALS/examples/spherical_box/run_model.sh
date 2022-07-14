export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Run the same model with 3 different spherical solvers
rm -rf ./OUTPUTS model
make -B -f make_model > build_log.log
eval "$OMP_RUN_COMMAND ./model leapfrog_nonlinear"
eval "$OMP_RUN_COMMAND ./model cliffs"
eval "$OMP_RUN_COMMAND ./model rk2"

# Run the test and make figures
Rscript compare_models.R
