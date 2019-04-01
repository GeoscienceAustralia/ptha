# Run the convergence test, so the "test_convergence.R" script can be run
rm -rf ./OUTPUTS model
make -B -f make_model > build_log.log
./model 1
./model 2
./model 4
./model 8
#./run_model 16
Rscript test_convergence.R
