# Clean existing binary
rm ./uniform_slope ./outfile.log
rm -r ./OUTPUTS
# Build the code
make -B -f make_uniform_slope > build_outfile.log
# Run the job
./uniform_slope 'rk2' > outfile.log
# Report tests
Rscript plot_results.R

#
# As above, for leapfrog_linear_plus_nonlinear_friction
# This numerical method is not 'good' for this problem, so
# the test code tries to make it easier using "nice" initial conditions.
#

# Clean existing binary
rm ./outfile.log
rm -r ./OUTPUTS
# Run the job
./uniform_slope 'leapfrog_linear_plus_nonlinear_friction' > outfile.log
# Report tests
Rscript plot_results.R
