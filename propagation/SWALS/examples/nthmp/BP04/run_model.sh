# Clean existing binary
rm ./BP4_testcases
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP4_testcases > build_outfile.log
# In this case the R script does plotting and reports tests
Rscript plot.R
