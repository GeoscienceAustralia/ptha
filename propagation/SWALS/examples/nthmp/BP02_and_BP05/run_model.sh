export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./BP2_testcases
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP2_testcases > build_outfile.log
# In this case the R script does plotting and reports tests
Rscript plot.R linear
Rscript plot.R leapfrog_nonlinear
Rscript plot.R rk2

# Cliffs result is interesting. We can see it is not shock-capturing (delayed
# shocks, like in dam-break test case).
# Rscript plot.R cliffs 
