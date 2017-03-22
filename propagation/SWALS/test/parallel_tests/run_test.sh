# get the CAF compiler, compile the code and run the tests
# We run with many cores since it can be easier to provoke failures
source ~/caf_setup.sh
make -B -f make_test
export OMP_NUM_THREADS=1
time cafrun -np 1 ./parallel_unit_tests
time cafrun -np 2 ./parallel_unit_tests
time cafrun -np 3 ./parallel_unit_tests
time cafrun -np 4 ./parallel_unit_tests
time cafrun -np 6 ./parallel_unit_tests
time cafrun -np 8 ./parallel_unit_tests
time cafrun -np 12 ./parallel_unit_tests
time cafrun -np 16 ./parallel_unit_tests
time cafrun -np 32 ./parallel_unit_tests
time cafrun -np 64 ./parallel_unit_tests
