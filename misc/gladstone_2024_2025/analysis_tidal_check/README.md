Analysis scripts for the [mini-PTHA](../swals/run_ptha_tidal_check/README.md) with the spatially varying tidal adjustment.
This is to compare the results with and without the adjustment technique.

# 1. Check the max-stage exceedance-rate curves vs PTHA18
If these aren't good, the mini-PTHA won't estimate the true hazard well.
However, as long as the scenarios are kept the same between tidal adjustment methods (static vs spatially varying), then it's a fair comparison.  

# 2. Compute the probabilistic inundation
Depends on max-stage-at-a-point for testing.
Required for testing in step 4.

# 3. JATWC to inundation
Using JATWC warning rules, compute the maximum expected inundation.
Optional comparison.

# 4. Check that PTHA analysis is equivalent
Test the spatially varying tide results are similar to static tide results in [ptha_tidal_adjustment_tests/](ptha_tidal_adjustment_tests).

Results like these are published in a conference paper in Coasts and Ports 2025:
*Accommodating Spatially Varying Tidal Planes in Tsunami Hazard Assessments*, available at https://pid.geoscience.gov.au/dataset/ga/150394.
