# get the CAF compiler, compile the code and run the tests
# We run with many cores since it can be easier to provoke failures
rm ./parallel_unit_tests
#source setup_coarray.sh
make -B -f make_test
#export OMP_NUM_THREADS=1
time OMP_NUM_THREADS=1 mpiexec -n 1 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 2 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 3 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 4 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 6 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 8 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 12 ./parallel_unit_tests
time OMP_NUM_THREADS=1 mpiexec -n 16 --oversubscribe ./parallel_unit_tests
#time mpiexec -n 32 ./parallel_unit_tests
#time mpiexec -n 64 ./parallel_unit_tests
