Make the random sources with:
```
Rscript sample_random_scenarios_simple.R
Rscript select_random_scenarios.R
```
Then use the following, on a machine with 48 cores, to make the initial conditions.
```
qsub run_generate_initial_conditions.PBS
```
