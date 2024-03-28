OMP_NUM_THREADS=2 OMP_PROC_BIND=true mpirun -np 8 --map-by ppr:1:core:PE=2 ./model_gfortran '' test03 test 0.0

