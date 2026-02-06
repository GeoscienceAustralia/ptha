export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./BP5_testcases_single_grid
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP5_testcases_single_grid > build_outfile.log
# In this case the R script does plotting and reports tests
#Rscript plot.R linear
Rscript plot_single_grid.R leapfrog_nonlinear
Rscript plot_single_grid.R rk2
Rscript plot_single_grid.R midpoint
#Rscript plot_single_grid.R cliffs
