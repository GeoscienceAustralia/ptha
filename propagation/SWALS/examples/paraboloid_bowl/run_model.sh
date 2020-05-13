export SWALS_SRC='../../src'
source ${SWALS_SRC}/test_run_commands

# Clean existing binary
rm ./paraboloid_bowl
rm -r ./OUTPUTS
# Build the code
make -B -f make_paraboloid_bowl > build_outfile.log
# Run the job
eval "$OMP_RUN_COMMAND ./paraboloid_bowl '' > outfile.log"
# Plot it and report tests
Rscript plot.R

#
# New test -- check we get the same result with openmp and mpi
#
source parallel_reproduce.sh
Rscript parallel_reproduce_check.R
