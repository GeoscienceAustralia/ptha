#
# NON-COARRAY
#
# Clean existing binary
rm ./tauranga
rm -r ./OUTPUTS
# Build the code
make -B -f make_tauranga > build_outfile.log
# Run the code
./tauranga > outfile.log
# Plot and report tests
echo '# Testing openmp version '
Rscript plot_results.R lowresolution_omp
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log.log ./multidomain_log_lowresolution_openmp.log

#
# SAME MODEL WITH COARRAY
#
# Clean existing binary
rm ./tauranga
rm -r ./OUTPUTS 
# Build the code
make -B -f make_tauranga_coarray > build_outfile.log
# Run the code
OMP_NUM_THREADS=2 cafrun -np 6 ./tauranga > outfile.log
#OMP_NUM_THREADS=2 OMP_PROC_BIND=true mpiexec -n 6 --map-by core ./tauranga > outfile.log
# Plot and report tests
echo '# Testing coarray version '
Rscript plot_results.R lowresolution_coarray
# Store a log file for later tests
cp ./OUTPUTS/RUN*/multidomain_log*01.log ./multidomain_log_lowresolution_coarray.log

# Run the comparison script
echo '# Comparing coarray/openmp versions '
Rscript compare_logs_coarray_openmp.R
