# Clean existing binary
rm ./BP2_testcases
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP2_testcases > build_outfile.log
# In this case the R script does plotting and reports tests
Rscript plot.R
