#
# Code to define commands to execute test programs
#
# This gives a central place to define e.g. number of openmp threads, or method
# to invoke mpi. Makes it easier to test efficiently on different machines.
#

# FIXME: Currently some of the validation tests ASSUME 6 mpi processes (or
# coarray images). This includes the tests with a prescribed
# load_balance_partition file. In future the tests might be adapted to generate
# those files from the information below. For now keep 6 mpi processes, but
# optionally change the number of openmp threads

export OMP_RUN_COMMAND='OMP_NUM_THREADS=48 OMP_PROC_BIND=true '

export CAF_RUN_COMMAND='OMP_NUM_THREADS=8 mpiexec -n 6 -ppr:3:socket:PE=8 '

# These calls are used for a parallel reproducibility check that splits the domain in 4.
export CAF_RUN_COMMAND_ONE_IMAGE='OMP_NUM_THREADS=48 mpiexec -n 1 -ppr:1:node:PE=48 '
export CAF_RUN_COMMAND_FOUR_IMAGES='OMP_NUM_THREADS=12 mpiexec -n 4 -ppr:2:socket:PE=12 '

