#
# NON-COARRAY
#
# Clean existing binary
rm ./BP09
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP09 > build_outfile.log
# Run the code
./BP09 > outfile.log
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log.log ./multidomain_log_lowresolution_openmp.log

#
# SAME MODEL WITH COARRAY
#
# Clean existing binary
rm ./BP09
rm -r ./OUTPUTS 
# Build the code
make -B -f make_BP09_coarray > build_outfile.log
# Run the code
OMP_NUM_THREADS=2 cafrun -np 6 ./BP09 > outfile.log
# Plot and report tests
echo '# Testing coarray version '
Rscript plot_results.R lowresolution_coarray
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log*01.log ./multidomain_log_lowresolution_coarray.log

# Run the comparison script
echo '# Comparing coarray/openmp versions '
Rscript compare_logs_coarray_openmp.R
