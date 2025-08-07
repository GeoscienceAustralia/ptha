rm ./radial_potential_nested
rm -r ./OUTPUTS
make -B -f make_radial_potential

# Make sure the initial conditions are up to date
Rscript potential_solution.R
eval "$OMP_RUN_COMMAND ./radial_potential_nested '' 1"

# Compare the numerical/analytical solutions
Rscript plot.R "single_grid_OMP"

#
# Nested grid version + MPI
#
rm ./radial_potential_nested
rm -r ./OUTPUTS
make -B -f make_radial_potential_coarray > build_outfile.log

# Make sure the initial conditions are up to date
Rscript potential_solution.R
eval "$CAF_RUN_COMMAND ./radial_potential_nested '' 2"

# Compare the numerical/analytical solutions
Rscript plot.R "nested_with_MPI"
