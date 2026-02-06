rm -r ./OUTPUTS
rm ./shoaling_variable_seabed ./libSWE.a
make -B -f make_shoaling_variable_seabed
Rscript run_and_plot.R
