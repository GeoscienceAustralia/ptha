# 1 Exceedance rate calculations
Modify the test code to suit your case, then run it and check that all cases PASS
```
Rscript test_exceedance_rate_raster_calculations_offshore.R
Rscript test_exceedance_rate_raster_calculations_rosslyn_bay.R
```
Should both print a few "PASS"

# 2 exceedance rates with uncertainty
Modify the test code to suit your case then run and check that all checks PASS 
This uses multiple cores so probably needs an interactive job on NCI.
```
Rscript test_compute_exceedance_rates_at_epistemic_uncertainty.R
```
Should print "PASS"

# 3. Thresholds
Modify the test code below to suit your case, then run it and check that it prints PASS.
This will require a node, with as many cores needed in application_specific_inputs.R::MC_CORES.
The first couple of tests take about 10-15 mins to run each, but the 3rd test works on a whole domain and takes about 2 hours to run. 
```
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R > test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.log
```
It will output some logs but each test should print PASS at the end.
If the test works, proceed with calculations of interest.
