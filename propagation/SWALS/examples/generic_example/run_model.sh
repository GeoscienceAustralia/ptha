# Clean existing binary
rm ./generic_model
rm -r ./OUTPUTS
# Build the code
make -B -f make_generic_model > build_outfile.log

# Run linear
rm -r ./OUTPUTS
./generic_model test_model_japan_linear.in > outfile.log
Rscript plot.R linear

# Run linear_with_nonlinear_friction
rm -r ./OUTPUTS
./generic_model test_model_japan_linear_with_nonlinear_friction.in > outfile.log
Rscript plot.R linear_with_nonlinear_friction

# Run almost_linear
rm -r ./OUTPUTS
./generic_model test_model_japan_almost_linear.in > outfile.log
Rscript plot.R almost_linear

# Run almost_linear_with_nonlinear_friction
rm -r ./OUTPUTS
./generic_model test_model_japan_almost_linear_with_nonlinear_friction.in > outfile.log
Rscript plot.R almost_linear_with_nonlinear_friction
