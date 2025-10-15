export SWALS_SRC='../../../src'
source ${SWALS_SRC}/test_run_commands

rm -r ./OUTPUTS
rm ./solitary_shoaling_grilli ./libSWE.a
make -B -f make_solitary_shoaling_grilli > build_outfile.log

eval "$OMP_RUN_COMMAND ./solitary_shoaling_grilli midpoint 0.44 0.2 > outfile.log"

Rscript plot.R
