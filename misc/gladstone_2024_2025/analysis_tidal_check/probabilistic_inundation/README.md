# Probabilistic Inundation calculations

This folder contains scripts for probabilistic inundation calculation, which can be run after all the models in [../../../swals](../../../swals) have been simulated.

Before running anything you'll need to modify the [application_specific_file_metadata.R](application_specific_file_metadata.R) for your case.

Once set up run every script in the [test](test) directory (tests with `compute` in the name run on multiple cores). They will need file paths adjusted and depends on the [../max_stage_at_a_point](../max_stage_at_a_point/) calculations. Then you may:

- compute the [exrate_given_threshold](exrate_given_threshold),
- calculate a [threshold_given_exrate](threshold_given_exrate) corresponding to an exceedance-rate and epistemic uncertainty percentile,
- compute the min and mean [arrival_time](arrival_time)
