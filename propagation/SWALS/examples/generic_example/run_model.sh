export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

#
# CASES WITH CORIOLIS
#

# Clean existing binary
rm ./generic_model
rm -r ./OUTPUTS
# Build the code
make -B -f make_generic_model > build_outfile.log

# Run linear
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_linear.in > outfile.log"
Rscript plot.R linear

# Run linear_with_nonlinear_friction
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_linear_with_nonlinear_friction.in > outfile.log"
Rscript plot.R linear_with_nonlinear_friction

# Run almost_linear
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_almost_linear.in > outfile.log"
Rscript plot.R almost_linear

# Run almost_linear_with_nonlinear_friction
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_almost_linear_with_nonlinear_friction.in > outfile.log"
Rscript plot.R almost_linear_with_nonlinear_friction

# Run leapfrog_nonlinear
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_leapfrog_nonlinear.in > outfile.log"
Rscript plot.R leapfrog_nonlinear

#
# CASES WITHOUT CORIOLIS
#

rm ./generic_model
# Build the code
make -B -f make_generic_model_nocoriolis > build_outfile.log

# Run linear
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_linear.in > outfile.log"
Rscript plot.R linear_no_coriolis

# Run linear_with_nonlinear_friction
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_linear_with_nonlinear_friction.in > outfile.log"
Rscript plot.R linear_with_nonlinear_friction_no_coriolis

# Run almost_linear
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_almost_linear.in > outfile.log"
Rscript plot.R almost_linear_no_coriolis

# Run almost_linear_with_nonlinear_friction
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_almost_linear_with_nonlinear_friction.in > outfile.log"
Rscript plot.R almost_linear_with_nonlinear_friction_nocoriolis

# Run leapfrog_nonlinear
rm -r ./OUTPUTS
eval "$OMP_RUN_COMMAND ./generic_model test_model_japan_leapfrog_nonlinear.in > outfile.log"
Rscript plot.R leapfrog_nonlinear_nocoriolis
