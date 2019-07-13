rm -r OUTPUTS
rm paraboloid_bowl
make -B -f make_paraboloid_bowl_coarray > build_outfile.log
# Pure openmp
OMP_NUM_THREADS=4 mpiexec -n 1 ./paraboloid_bowl 'load_balance_partition.txt' > outfile.log
# Pure coarray
OMP_NUM_THREADS=1 mpiexec -n 4 ./paraboloid_bowl 'load_balance_partition.txt' > outfile.log
