# Estimating Monte-Carlo uncertainties for stratified and stratified/importance sampling.
-----------------------------------------------------------------------------------------

This tutorial follows on from the
[random_scenario_sampling.md](random_scenario_sampling.md) tutorial. It shows
how to estimate the statistical properties of the Monte-Carlo exceedance-rate errors. 

Make sure you understand the content in 
[random_scenario_sampling.md](random_scenario_sampling.md) before proceeding.


# Background
------------

Monte-Carlo sampling techniques have some approximation error, associated with
use of a finite random sample of scenarios. In theory as more and more
scenarios are sampled, the solution will converge to the "correct" solution
(i.e. that which would be obtained by using every PTHA scenario). 

The theory of Monte-Carlo sampling is well understood, and estimates of the
Monte-Carlo errors can be obtained in 2 situations:
1. Prior to sampling, the variance of the Monte-Carlo exceedance-rate (for a given
   site and threshold) can be analytically computed at PTHA18 hazard points
   (i.e. offshore sites).
2. After sampling, when high-resolution models have been run for every random
   scenario, the variance of the Monte-Carlo exceedance-rate can be estimated
   using only the random scenarios (at any site covered by the high-resolution
   model).

Above the "variance of the Monte-Carlo exceedance-rate" is the number that we would
get by repeating Monte-Carlo sampling many many times, storing the exceedance-rate
each time, and then computing their variance. It is a measure of the "typical accuracy"
of the estimated Monte-Carlo exceedance-rate. With a smaller variance, any given Monte-Carlo
estimate is more likely to be accurate.

In practical applications there are 2 situations where this is useful:
* When designing the Monte-Carlo sampling scheme. Then approach 1 is helpful
  to cheaply estimate the errors associated with a particular scheme and
  sampling effort. By checking the performance of a Monte-Carlo sampling scheme
  at offshore sites, prior to expensive inundation computation, it is easier to
  make good decisions on which scheme to use. For example if a scheme leads to high
  errors offshore of the site-of-interest, we might assume it will also be inaccurate
  at the onshore site of interest.
* After high-resolution computations have been conducted for a random sample of
  scenarios. In this situation we often want an estimate of the Monte-Carlo
  uncertainties at coastal and onshore sites (covered by the high-resolution
  model, but not PTHA18). Approach 2 enables the Monte-Carlo uncertainty to be
  estimated at such sites, using just the random Monte-Carlo sample.

## Get the source-zone event data, and some maximum-stage data.
---------------------------------------------------------------

Here we get the source-zone events data and tsunami maxima at a site of
interest, following the approach used in the
[random_scenario_sampling.md](random_scenario_sampling.md) tutorial.


```r
# Get the scripts to access the PTHA18
ptha18 = new.env()
source('../../get_PTHA_results.R', local=ptha18, chdir=TRUE)

# Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in PTHA18)
source_zone = 'kermadectonga2'
source_zone_scenarios = ptha18$get_source_zone_events_data(source_zone,  slip_type='stochastic')

# Get the tsunami maxima at a point of interest
event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    target_point = c(185.1239, -21.0888), # Known location of PTHA18 hazard point
    slip_type='stochastic', # Matches the slip-type used to create scenarios
    all_source_names=source_zone)
# Convenient shorthand
event_peak_stage_ref = event_peak_stage_at_refpoint[[source_zone]]$max_stage

# Convenient shorthand for the magnitudes and rates in the event table
event_Mw = source_zone_scenarios$events$Mw 
event_rates = source_zone_scenarios$events$rate_annual # Logic-tree mean model
```

The above variables will be used in the code below.

## 1. Confidence interval for Monte-Carlo exceedance-rates at offshore sites, prior to Monte-Carlo sampling
-----------------------------------------------------------------------------------------------------------

Suppose we want to understand the likely errors in the exceedance-rate at an offshore point of interest in the PTHA18. 

Below we do the calculation with a stage-threshold of 3m.

### Stratified-sampling

We firstly illustrate the calculations for stratified-sampling. In this example
we assume 12 samples per magnitude bin are used up to magnitude 9.6, although
any non-uniform sampling effort could also be specified.

```r
samples_per_Mw_stratified<-function(Mw){ 12*(Mw < 9.65) } # Sample size function used below
```

Below we compute the exact exceedance-rate (using all scenarios), and the
variance of Monte-Carlo exceedance-rates that would result from the proposed
stratified-sampling scheme.

```r
stage_threshold_example = 3 # 3m at the site where event_peak_stage_ref was defined.

# Get the "exact" exceedance-rate and the analytical Monte-Carlo variance, stratifed-sampling
exrate_and_var_stratified = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
    event_Mw, event_rates, 
    event_peak_stage_ref, 
    stage_threshold=stage_threshold_example, 
    samples_per_Mw=samples_per_Mw_stratified)

# Print the exceedance-rate and its variance
exrate_and_var_stratified
```

```
## [1] 6.573855e-04 4.912610e-08
```

An approximate 95% confidence interval for the Monte-Carlo exceedance-rate can
be obtained, assuming a normal distribution

```r
# Quantiles of a normal distribution
stratified_confint = qnorm(c(0.025, 0.5, 0.975), 
    mean=exrate_and_var_stratified[1], 
    sd=sqrt(exrate_and_var_stratified[2]))
stratified_confint
```

```
## [1] 0.0002229711 0.0006573855 0.0010918000
```
This means that about 95\% of the time, we expect the exceedance-rate with a stage-threshold= 
3
to be within 
2.23 &times; 10<sup>-4</sup> 
and 
0.00109. 


### Stratified/importance sampling

Below we show the same calculations using stratified/importance-sampling. 

We assume 12 samples per magnitude bin are used, although any non-uniform sampling effort
can be specified.

```r
samples_per_Mw_stratified_importance<-function(Mw){ 12*(Mw < 9.65) } # Sample size function used below
```
Because we are using stratified/importance-sampling, we also need to know the `event_importance`.

```r
event_importance = event_peak_stage_ref # Better sample scenarios where this is high
```

Below we compute the exact exceedance-rate (using all scenarios), and the
variance of Monte-Carlo exceedance-rates using the proposed sampling scheme: 

```r
stage_threshold_example = 3 # 3m at the site where event_peak_stage_ref was defined.

# Get the "exact" exceedance-rate and the analytical Monte-Carlo variance
exrate_and_var_stratified_importance = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
    event_Mw, event_rates, 
    event_peak_stage_ref, 
    stage_threshold=stage_threshold_example, 
    samples_per_Mw=samples_per_Mw_stratified,
    event_importance_weighted_sampling_probs = (event_rates * event_importance) # Importance sampling
    )

# Print the exceedance-rate and its variance
exrate_and_var_stratified_importance
```

```
## [1] 6.573855e-04 7.946050e-09
```
Note the exceedance-rate is the same as for stratified-sampling (because it
is exact). However the variance has reduced.

An approximate 95% confidence interval for the result obtained by Monte-Carlo
sampling can be obtained if we assume a normal distribution

```r
# Quantiles of a normal distribution
stratified_importance_confint = qnorm(c(0.025, 0.5, 0.975), 
    mean=exrate_and_var_stratified_importance[1], 
    sd=sqrt(exrate_and_var_stratified_importance[2]))
stratified_importance_confint
```

```
## [1] 0.0004826731 0.0006573855 0.0008320980
```
This means that about 95\% of the time, we expect the exceedance-rate with a stage-threshold= 
3
to be within 
4.83 &times; 10<sup>-4</sup> 
and 
8.32 &times; 10<sup>-4</sup>.

Notice this interval is substantially narrower than the previous interval.

### Summary

Above we computed the variance of Monte-Carlo exceedance-rates that would be
expected at a given site, for a given stage-threshold. This was also used to
compute a 95% confidence interval. If the Monte-Carlo sampling were repeated
many times, then approximately 95% of them would be contained in the confidence
interval; it is approximate only because we assumed a normal distribution. 

Notice that no Monte-Carlo sampling was involved in these calculations. This is
possible at PTHA18 hazard points because the statistical properties of
Monte-Carlo sampling are well understood.

The computations were performed for both stratified and
stratified/importance-sampling, both using 12 scenarios in each magnitude-bin.
One way to quantify the improvement gained via stratified/importance-sampling
is to look at the ratios of the variances in each case.

```r
# Ratio of variances of Monte-Carlo exceedance-rates at the stage-threshold
variance_reduction_factor = exrate_and_var_stratified[2]/exrate_and_var_stratified_importance[2]
variance_reduction_factor
```

```
## [1] 6.182456
```
For these schemes, the variance reduces in proportion to the sampling effort.
The result indicates we would need approximately 6.18 
times as many samples with stratified-sampling to get the same accuracy as
stratified/importance-sampling. This is specific to the chosen site,
scenario-frequency model, and stage threshold.

In both cases, the errors can also be reduced by using a larger Monte-Carlo
sample. The average error scales with the square-root of the sample size, so
e.g., to halve the errors, use 4x as many scenarios.

## 2. Confidence interval for Monte-Carlo exceedance-rates at coastal sites after high-resolution simulation. 
-------------------------------------------------------------------------------------------------------------

Suppose that high-resolution simulations have been performed for every scenario
in a Monte-Carlo sample. That means that we can compute onshore tsunami
statistics for those scenarios (which we cannot do for other PTHA18 scenarios). 

Below we show how to estimate the Monte-Carlo errors using only the random
sample. This is typically useful to estimate the Monte-Carlo error for onshore
quantities of interest.

### Stratified-sampling

The implementation requires a Monte-Carlo sample:

```r
# Make a reproducible random seed to make the code reproducible (this is optional)
set.seed(12345)

# Make the random scenarios -- stratified-sampling
random_scenarios_stratified = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    samples_per_Mw=samples_per_Mw_stratified
    )
```

We also need to know the quantity-of-interest at some onshore site of interest, for
only the random tsunami scenarios. Here, for demonstration purposes, we assume that the
quantity of interest is equal to twice the offshore `event_peak_stage_ref` value. But in general
this value should be extracted from the tsunami model.

```r
# The code takes a quantity of interest for EVERY PTHA SCENARIO, but only uses the values
# for the randomly sampled scenarios. 
onshore_quantity_of_interest = rep(NA, length(event_peak_stage_ref)) # One for every PTHA scenario
# Now set values for the randomly sampled scenarios (would normally come from a tsunami model),
# using their indices
random_inds = na.omit(random_scenarios_stratified$inds)
onshore_quantity_of_interest[random_inds] = 
    2* event_peak_stage_ref[random_inds]
```
Please note the `onshore_quantity_of_interest` contains ONE VALUE FOR EVERY PTHA SCENARIO, even though
they are missing for all except the random scenarios.

Below we estimate the exceedance-rate and its variance at a threshold of 6 m

```r
exrate_estimate_stratified = ptha18$estimate_exrate_uncertainty(
    random_scenarios_stratified,
    onshore_quantity_of_interest,
    threshold_stage = 6
    )
# Estimate of the exceedance-rate mean, and the Monte-Carlo variance.
exrate_estimate_stratified
```

```
## [1] 1.026517e-03 1.118763e-07
```

From this we can estimate a 95% confidence interval for the true exceedance-rate, using the
estimated mean and variance, assuming a normal distribution.

```r
stratified_onshore_confint = qnorm(c(0.025, 0.5, 0.975), 
    mean=exrate_estimate_stratified[1], 
    sd=sqrt(exrate_estimate_stratified[2]))
# Print the confidence interval -- [lower, mean, upper]
stratified_onshore_confint
```

```
## [1] 0.0003709495 0.0010265167 0.0016820839
```

In this particular case, because we used a made-up definition of the
quantity-of-interest, we can compute the exact answer: 
6.57 &times; 10<sup>-4</sup> 

The exact answer is contained in the confidence interval in this case, 
which is expected approximately 95% of the time.

### Stratified/importance sampling

Here we repeat the above calculations, using stratified/importance sampling.

The implementation requires a Monte-Carlo sample:

```r
event_importance = event_peak_stage_ref # Suggest this is located "near site of interest"

# Make the random scenarios -- stratified/importance-sampling
random_scenarios_stratified_importance = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    samples_per_Mw=samples_per_Mw_stratified_importance,
    event_importance_weighted_sampling_probs = (event_rates * event_importance), # Importance sampling
    )
```

We also need to know the quantity-of-interest at some onshore site of interest, for
only the random tsunami scenarios. Here, for demonstration purposes, we assume that the
quantity of interest is equal to twice the offshore `event_peak_stage_ref` value. But in general
this value should be extracted from the tsunami model.

```r
# The code takes a quantity of interest for EVERY PTHA SCENARIO, but only uses the values
# for the randomly sampled scenarios. 
onshore_quantity_of_interest = rep(NA, length(event_peak_stage_ref)) # One for every PTHA scenario
# Now set values for the randomly sampled scenarios (would normally come from a tsunami model),
# using their indices
random_inds = na.omit(random_scenarios_stratified_importance$inds)
onshore_quantity_of_interest[random_inds] = 
    2* event_peak_stage_ref[random_inds]
```
Please note the `onshore_quantity_of_interest` contains ONE VALUE FOR EVERY PTHA SCENARIO, even though
they are missing for all except the random scenarios.

Below we estimate the exceedance-rate and its variance at a threshold of 6 m.

```r
exrate_estimate_stratified_importance = ptha18$estimate_exrate_uncertainty(
    random_scenarios_stratified_importance,
    onshore_quantity_of_interest,     
    threshold_stage = 6
    )
# Estimate of the exceedance-rate mean, and the Monte-Carlo variance.
exrate_estimate_stratified_importance
```

```
## [1] 7.008044e-04 6.806318e-09
```

From this we can estimate a 95% confidence interval for the true exceedance-rate, using the
estimated mean and variance, assuming a normal distribution.

```r
stratified_importance_onshore_confint = qnorm(c(0.025, 0.5, 0.975), 
    mean=exrate_estimate_stratified_importance[1], 
    sd=sqrt(exrate_estimate_stratified_importance[2]))
# Print the confidence interval -- [lower, mean, upper]
stratified_importance_onshore_confint
```

```
## [1] 0.0005391065 0.0007008044 0.0008625022
```
Notice this confidence interval is substantially narrower than the one obtained
with stratified-sampling. 

In this particular case, because we used a made-up definition of the
quantity-of-interest, we can compute the exact answer: 
6.57 &times; 10<sup>-4</sup> 

The exact answer is contained in the confidence interval in this case, 
which is expected approximately 95% of the time.

