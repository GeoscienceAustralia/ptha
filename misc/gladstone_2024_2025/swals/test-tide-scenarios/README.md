# Testing validity of tidal adjustment

Run a mini PTHA at a static tide and a varying tide to see if the results match.

1. Create PBS scripts to run the (n_sites + 1) x (n_scenarios=50) simulations using [create_tidal_tests.sh](create_tidal_tests.sh).
2. Submit the jobs using [submit_jobs.sh](submit_jobs.sh).
3. Inspect the match using [tide_gof.R](tide_gof.R).
