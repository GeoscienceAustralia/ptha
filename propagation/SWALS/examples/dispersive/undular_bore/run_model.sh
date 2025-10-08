rm -r ./OUTPUTS
rm ./undular_bore ./libSWE.a
make -B -f make_undular_bore
Rscript run_and_plot.R
