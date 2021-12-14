#
# Here we compare openmp/coarray versions of the models that use the
# same domain partition. This eliminates the small differences between
# the models when they use a different domain partition. 
#
# Before running this, consider storing the output in full precision
# (i.e. in src/global_mod.f90, set "output_precision = C_DOUBLE").
#

export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

echo 'Will run openmp version with: ' $OMP_RUN_COMMAND
echo 'Will run coarray version with: ' $CAF_RUN_COMMAND

#
# Regular openmp, with enforced domain partitioning that will match the coarray
# default partition
#

# Clean existing binary
rm ./BP09
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP09 > build_outfile.log
# Run the code
eval "$OMP_RUN_COMMAND ./BP09 load_balance_6_trivial.txt > outfile_omp.log"
#eval "OMP_NUM_THREADS=1 ./BP09 load_balance_6_trivial.txt > outfile_omp.log"
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp
# Move the openmp results to a new folder
mv $( dirname OUTPUTS/RUN*/multidomain_log.log ) OUTPUTS/openmp_results

#
# Same model with coarray, same partition
#

# Clean existing binary
rm ./BP09
#rm -r ./OUTPUTS 
# Build the code
make -B -f make_BP09_coarray > build_outfile.log
# Run the code
eval "$CAF_RUN_COMMAND ./BP09 load_balance_6_trivial.txt > outfile_ca.log"
#eval "OMP_NUM_THREADS=1 mpiexec -np 6 ./BP09 load_balance_6_trivial.txt > outfile_ca.log"
# Plot and report tests
echo '# Testing coarray version '
Rscript plot_results.R lowresolution_coarray
# Move the coarray results to a new folder
mv $( dirname OUTPUTS/RUN*/mu*001.log ) OUTPUTS/coarray_results

# Run the comparison script
echo '# Comparing coarray/openmp versions '
Rscript compare_logs_coarray_openmp.R
