# Timeline of updates to the PTHA 

Below we note updates to the PTHA following initial publication, and
differences with subsequent related publications.


## 2018

* The [PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) was released in November 2018. 

## 2019

* While drafting [the GJI publication](https://doi.org/10.1093/gji/ggz260) we identified a minor 
bug in the p-value calculation for the median-coverage-statistic. This led to
small changes in the median-coverage-statisic p-values (with updates reported
in the GJI publication). The changes are too small to affect the
interpretation. For example, the FAUS-with-fixed-rigidity p-value changed from
0.017 (PTHA18 report) to 0.003 (GJI Paper). 

* The [PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) tests a range
  of scenario generation models, and suggests that (p96): *For applications we
recommend use of the HS scenarios, with either constant or variable
shear-modulus depending on the perceived characteristics of the site of
interest*. However following further analysis and review, for [the GJI
publication](https://doi.org/10.1093/gji/ggz260) we suggested it is also worth
considering use of the bias-adjusted VAUS scenarios (p1957): 
*In applications it may be beneficial to use both HS and bias-adjusted VAUS
approaches, to represent epistemic uncertainties in scenario generation.
Furthermore, the VAUS approach may facilitate greater computational efficiency
because it involves fewer degrees of freedom*. Considering the large uncertainties
in scenario generation and the fact that both the (bias-adjusted) VAUS and HS
methods lead to a similar offshore hazard for Australia, it seems reasonable to
try both approaches (notwithstanding that real earthquakes do have heterogeneous slip).

* Both the [PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041) and 
[the PAGEOPH paper](https://link.springer.com/article/10.1007/s00024-019-02299-w) 
compare the modelled magnitude-exceedance-rates with site-specific results from
other studies (Section 3.7.8.2 of the PTHA18 report, and Table 1 of the PAGEOPH
paper). While developing the PAGEOPH paper we adjusted some of the
`region-for-comparison` definitions to better agree with the regions in the
cited studies. We also adopted a slightly better approach to interpolating
modelled magnitude-exceedance-rate curves at the magnitudes of interest, which
is less sensitive to the magnitude binning details. Combined this leads to some
differences in the reported earthquake ARIs (e.g. compare Table 1 in the
PAGEOPH paper to Section 3.8.7.2 in the PTHA18 report). Code for both variants
can be accessed [here](../R/examples/austptha_template/EVENT_RATES/README.md),
see the links to scripts with names like `earthquake_rate_comparisons_XXXX.R`.

* The review of [the PAGEOPH paper](https://link.springer.com/article/10.1007/s00024-019-02299-w)
prompted us to improve the stage-vs-exceedance-rate percentile uncertainty calculation method from the
[PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041). This is discussed in Section 3.5 of the PAGEOPH paper. 
This leads to some (generally small) changes in the uncertainty percentiles for
stage-vs-exceedance-rate. The codes to do the revised calculations can be accessed
[here](../R/examples/austptha_template/EVENT_RATES/README.md), see the Section
`Updated stage-vs-exceedance-rate percentile uncertainty calculations` for the
links and context. The online results have been updated to reflect this,
although the older results can still be accessed (as discussed in the relevant
sections of [ptha_access/README.md](README.md) and
[ptha_access/DETAILED_README.md](DETAILED_README.md)).


