export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

rm -r ./OUTPUTS
rm ./solitary_shoaling_grilli ./libSWE.a
make -B -f make_solitary_shoaling_grilli > build_outfile.log

eval "$OMP_RUN_COMMAND ./solitary_shoaling_grilli midpoint 0.44 0.2 > outfile.log"
Rscript plot.R 'midpoint'

eval "$OMP_RUN_COMMAND ./solitary_shoaling_grilli 'rk2' 0.44 0.2 > outfile.log"
Rscript plot.R 'rk2'

eval "$OMP_RUN_COMMAND ./solitary_shoaling_grilli 'leapfrog_nonlinear' 0.44 0.2 > outfile.log"
Rscript plot.R 'leapfrog_nonlinear'

eval "$OMP_RUN_COMMAND ./solitary_shoaling_grilli 'cliffs' 0.44 0.2 > outfile.log"
Rscript plot.R 'cliffs'
