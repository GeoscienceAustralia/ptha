# 1 Exceedance rate calculations
Modify the test code to suit your case, then run it and check that all cases PASS.
These depend on having run the max_stage_at_a_point/ptha calculations. 
```bash
Rscript test_max_stage_exceedance_rate_raster_calculations_offshore.R
Rscript test__max_stage_exceedance_rate_raster_calculations_rosslyn_bay.R
```
Should both print a few "PASS"

# 2 exceedance rates with uncertainty
Modify the test code to suit your case then run and check that all checks PASS 
This uses multiple cores so probably needs an interactive job on NCI.
```bash
Rscript test_compute_exceedance_rates_at_epistemic_uncertainty.R
```
Should print "PASS"

# 3. Thresholds
Modify the test code below to suit your case, then run it and check that it prints PASS.
This will require a node, with as many cores needed in application_specific_inputs.R::MC_CORES.
```bash
Rscript test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.R 2>&1 | tee log/test_compute_threshold_at_exceedance_rate_of_epistemic_uncertainty.log
```
If the test works, proceed with calculations of interest.

# 4. Negative max-stage thresholds
Computing values of the min-stage for given exceedance rates requires working with negative numbers. To ensure that the [../threshold_given_exrate/](../threshold_given_exrate/) calculations work with negative numbers, we negated the max-stages and computed the thresholds as normal. These should be the negative of the max_stage rasters. First both the `max_stage` and `neg_max_stage` threshold_given_exrate rasters need to be created. To do this, run `max_stage` and `neg_max_stage` results generated in [../threshold_given_exrate/make_threshold_epistemic_uncertainty_jobs.R](../threshold_given_exrate/make_threshold_epistemic_uncertainty_jobs.R), which requires compute. Unlike the previous tests mentioned, this test must be run after computing the `threshold_given_exrate` results. Then run the test as:
``bash
Rscript test_neg_max_stage_vs_max_stage_threshold_rasters.R
```
