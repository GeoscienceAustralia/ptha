Compare the PTHA18 maximum-stage at a given offshore site with the result of the nonlinear SWALS models. This is a useful way to check the high-resolution Monte Carlo results. 
* Basically the Monte Carlo results should agree well with the offshore PTHA in deep water far from the coast (i.e. at sites where limitations in the offshore PTHA's hydrodynamic model are not important) for tsunami wave heights that are well covered by the Monte Carlo sampling regime. 
  * At sites closer to the coast or in shallower water, we do not expect the offshore PTHA's hydrodynamic model to work so well -- and we expect differences.
* Our Monte Carlo sampling methods deliberately avoid putting much effort on sampling very small scenarios, or very large scenarios. So we also expect differences for these.

It is of interest to separate the effects of A) insufficient sampling, and B) differences in the hydrodynamic models. 
* To help focus on A),  we provide figures (with `_Sampled_` in the title) where the max-stage in the sampled scenarios is replaced with the PTHA18 max-stage.
* To help focus on B), we provide figures that compare the max-stage in the PTHA18 and the nonlinear model, without regard to exceedance rates. 

The key code is [extract_max_stage_at_a_point.R](extract_max_stage_at_a_point.R). 
* This has many hard-coded input parameters that need modification for each case.

When running this on NCI, you need access to a range of filesystems. I tend to run with an interactive job, requested like: `qsub -I -X -P w85 -q express -l ncpus=48 -l mem=192gb -l wd -l jobfs=10gb -l walltime=02:00:00 -lstorage=scratch/w85+gdata/w85+gdata/fj6`
