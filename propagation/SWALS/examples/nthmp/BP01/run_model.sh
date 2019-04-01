# Clean existing binary
rm ./BP1_testcases
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP1_testcases > build_outfile.log
# In this case the R script does plotting and reports tests
Rscript plot.R
