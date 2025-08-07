export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./BP2_testcases
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP2_testcases > build_outfile.log
# In this case the R script does plotting and reports tests
#Rscript plot.R linear
Rscript plot.R leapfrog_nonlinear
Rscript plot.R rk2
Rscript plot.R midpoint
#Rscript plot.R cliffs
