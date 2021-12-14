#
# Below we run a few versions of the model
# - openmp, without local timestepping
# - coarray, without local timestepping
# - openmp, with local timestepping
#

export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

echo 'Will run openmp version with: ' $OMP_RUN_COMMAND
echo 'Will run coarray version with: ' $CAF_RUN_COMMAND

#
# Run regular openmp model without local timestepping
#

# Clean existing binary
rm ./BP09
rm -r ./OUTPUTS
# Build the code
make -B -f make_BP09 > build_outfile.log
# Run the code
#eval "$OMP_RUN_COMMAND ./BP09 load_balance_6_trivial.txt > outfile_omp.log"
eval "$OMP_RUN_COMMAND ./BP09 '' > outfile_omp.log"
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp
# Move the openmp results to a new folder
mv $( dirname OUTPUTS/RUN*/multidomain_log.log ) OUTPUTS/openmp_results

#
# Same model with coarray -- results not identical as above due to domain
# partitioning (but close)
#

# Clean existing binary
rm ./BP09
#rm -r ./OUTPUTS 
# Build the code
make -B -f make_BP09_coarray > build_outfile.log
# Run the code
#eval "$CAF_RUN_COMMAND ./BP09 load_balance_6_trivial.txt > outfile_ca.log"
eval "$CAF_RUN_COMMAND ./BP09 '' > outfile_ca.log"
# Plot and report tests
echo '# Testing coarray version '
Rscript plot_results.R lowresolution_coarray
# Move the coarray results to a new folder
mv $( dirname OUTPUTS/RUN*/mu*001.log ) OUTPUTS/coarray_results

# Run the comparison script
echo '# Comparing coarray/openmp versions '
Rscript compare_logs_coarray_openmp.R

#
# Openmp model with local timestepping and partitioning -- results not
# identical to previous due to local timestepping (but close)
#
rm ./BP09
#rm -r ./OUTPUTS
# Build the code
make -B -f make_BP09_localtimestep > build_outfile.log
# Run the code
eval "$OMP_RUN_COMMAND ./BP09 'load_balance_6_trivial.txt' > outfile_omp.log"
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp_localtimestep
# Move the results to a new folder
mv $( dirname OUTPUTS/RUN*/mu*001.log ) OUTPUTS/openmp_results_localtimestep

# Run the comparison script
echo '# Comparing openmp/openmp-localtimestep versions '
Rscript compare_logs_openmp_localtimestep.R

