# Compile the code and run in distributed memory parallel
# We run with many cores since it can be easier to provoke failures

# Can be handy to try other compilers [e.g. for PGI need to point to their variant of mpiexec]
export MPIEXEC='mpiexec --oversubscribe'
rm ./parallel_unit_tests
#source setup_coarray.sh
make -B -f make_test
#export OMP_NUM_THREADS=1
time OMP_NUM_THREADS=1 $MPIEXEC -n 1 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 2 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 3 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 4 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 6 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 8 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 12 ./parallel_unit_tests
time OMP_NUM_THREADS=1 $MPIEXEC -n 16 --oversubscribe ./parallel_unit_tests
#time mpiexec -n 32 ./parallel_unit_tests
#time mpiexec -n 64 ./parallel_unit_tests
