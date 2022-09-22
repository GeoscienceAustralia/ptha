rm ./radial_potential
rm -r ./OUTPUTS
make -B -f make_radial_potential

# Make sure the initial conditions are up to date
Rscript potential_solution.R
OMP_NUM_THREADS=12 ./radial_potential

# Compare the numerical/analytical solutions
Rscript plot.R
