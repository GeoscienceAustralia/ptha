#
# NON-COARRAY
#
export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

echo 'Will run openmp version with: ' $OMP_RUN_COMMAND
echo 'Will run coarray version with: ' $CAF_RUN_COMMAND

# Clean existing binary
rm ./model
rm -r ./OUTPUTS
# Build the code
make -B -f make_hilo > build_outfile.log
# Run the code
eval "$OMP_RUN_COMMAND ./model > outfile.log"
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log.log ./multidomain_log_lowresolution_openmp.log

#
# SAME MODEL WITH COARRAY
#
# Clean existing binary
rm ./model
#rm -r ./OUTPUTS 
# Build the code
make -B -f make_hilo_coarray > build_outfile.log
# Run the code
eval "$CAF_RUN_COMMAND ./model > outfile.log"
# Plot and report tests
echo '# Testing coarray version '
Rscript plot_results.R lowresolution_coarray
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log*01.log ./multidomain_log_lowresolution_coarray.log

# Run the comparison script
echo '# Comparing coarray/openmp versions '
Rscript compare_logs_coarray_openmp.R
