# Run in an interactive job on one node
#
# Make sure you have first done
#    source SWALS_ifx_modules_2025_llvm.sh
#
model_binary=/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/examples/generic_example/generic_model

OMP_NUM_THREADS=48 $model_binary scenario_1.in > outfile.log
