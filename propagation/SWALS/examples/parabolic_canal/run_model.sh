# Clean existing binary
rm ./parabolic_canal
rm -r ./OUTPUTS
# Build the code
make -B -f make_parabolic_canal > build_outfile.log

## Run with rk2
./parabolic_canal 'rk2' > outfile.log
# Make plot and report tests
Rscript plot.R 'rk2'

## Run with leapfrog_nonlinear
./parabolic_canal 'leapfrog_nonlinear' > outfile.log
# Make plot and report tests
Rscript plot.R 'leapfrog_nonlinear'

### Run with CLIFFS -- it needs more work
#./parabolic_canal 'cliffs' > outfile.log
## Make plot and report tests
#Rscript plot.R 'cliffs'
