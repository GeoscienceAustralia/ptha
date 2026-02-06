export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

#
# Test without MPI
#

rm ./nesting_stability
rm -r ./OUTPUTS
make -B -f make_nesting_stability > build_outfile.log
# Run it
eval "$OMP_RUN_COMMAND ./nesting_stability 'linear'"
Rscript check_speed.R "linear_OMP"
eval "$OMP_RUN_COMMAND ./nesting_stability 'midpoint'"
Rscript check_speed.R "midpoint_OMP"

#
# Test with MPI
#
rm ./nesting_stability
rm -r ./OUTPUTS
make -B -f make_nesting_stability_coarray > build_outfile.log

# Make sure the initial conditions are up to date
eval "$CAF_RUN_COMMAND ./nesting_stability 'linear'"
Rscript check_speed.R "linear_MPI"
eval "$CAF_RUN_COMMAND ./nesting_stability 'midpoint'"
Rscript check_speed.R "midpoint_MPI"

