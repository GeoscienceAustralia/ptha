# 1 Exceedance rate calculations
Modify the test code to suit your case, then run it and check that all cases PASS.
These depend on having run the max/min_stage_at_a_point/ptha calculations. 
```
Rscript test_max_stage_exceedance_rate_raster_calculations_offshore.R
Rscript test__max_stage_exceedance_rate_raster_calculations_rosslyn_bay.R
Rscript test_min_stage_exceedance_rate_raster_calculations_offshore.R
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
```
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R 2>&1 | tee log/test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.log
```
If the test works, proceed with calculations of interest.
