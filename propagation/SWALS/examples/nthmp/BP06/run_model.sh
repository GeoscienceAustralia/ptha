# Clean existing binary
rm ./BP06
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP06 > build_outfile.log
# In this case the R script does plotting and reports tests
Rscript plot.R
