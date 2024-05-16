Updates to PTHA, and suggestions for future updates
---------------------------------------------------

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
try both approaches (notwithstanding that real earthquakes have heterogeneous slip). 

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
can be accessed [here](../R/examples/austptha_template/EVENT_RATES/),
see the links to scripts with names like `earthquake_rate_comparisons_XXXX.R`.

* The review of [the PAGEOPH paper](https://link.springer.com/article/10.1007/s00024-019-02299-w)
prompted us to improve the stage-vs-exceedance-rate percentile uncertainty calculation method from the
[PTHA18 report](http://dx.doi.org/10.11636/Record.2018.041). This is discussed in Section 3.5 of the PAGEOPH paper. 
This leads to some (generally small) changes in the uncertainty percentiles for
stage-vs-exceedance-rate. The codes to do the revised calculations can be accessed
[here](../R/examples/austptha_template/EVENT_RATES/), see the README Section
`Updated stage-vs-exceedance-rate percentile uncertainty calculations` for the
links and context. The online results have been updated to reflect this,
although the older results can still be accessed (as discussed in the relevant
sections of [ptha_access/README.md](README.md) and
[ptha_access/DETAILED_README.md](DETAILED_README.md)).


## 2020

* To speed up remote reads (from the NCI THREDDS Server) we have created netcdf
files containing ONLY the max_stage variable for each tsunami event. The code in ptha_access
now often directly reads max_stage from such files, rather than reading it from 
other files that contain additional variables (which tends to be slower).

* To speed up remote reads (from the NCI THREDDS Server) we have created a single
netcdf file containing the ONLY the gauge variables: lon,lat,elev,gaugeID. Now in a
number of situations the ptha_access codes read from this file directly to get
coordinate variables. Previously we would read from other netcdf files, which also
contain other variables - but this can be slow. 

## 2021

* An update to the PTHA server caused the ptha_access scripts to break. This was easily
fixed (by using https links, rather than http). Users will have to update their codes to
get the database working again.


# Ideas for future updates 

* Advice on tsunami dissipation models, and integration of something like this into the peak-stage exceedance-rate products.

* We missed one of the DART buoys for the 2013 Santa Cruz event -- add it in.

* Consider there were strike-slip components in some historical events (e.g. Sumatra 2004, and recent literature suggests also perhaps Chile 1960). One could either add this in to the model, or re-interpret the data to focus on the thrust component.

* Consider making the Mw-frequency model have a stronger relation with "the capacity of a source-zone to host an earthquake". Currently we have hard-threshold geometric limits on earthquake size, which are not particularly restrictive because we don't want to rule out large compact ruptures. However the threshold is actually uncertain. Perhaps this information should have greater influence on the Mw-frequency model, say via the Mw-max priors (e.g. some tapering rather than a hard-threshold)?

* Draw on JG's revised ideas on large magnitude normal fault sources.

* There is a mistake in the `arutrough` and `tanimbar` scenarios -- both of which are normal faults. The scenarios appear to have been created with a rigidity of 30 GPa. Apparently later the rigidity was changed to 60 GPa, which is how we treat normal faults in general - but it appears we did not update the scenarios. The impact is likely that the sources behave like they have a rigidity of 30 GPa, but we have a more aggresive peak-slip-limit than usual. Although these source-zones are minor, and the rigidity poorly constrained, this should be fixed. The other normal fault sources seem fine. 

* The `Trobriand` source-zone is prescribed a much lower convergence rate than it would have been prescribed if we had used Bird 2003, based on arguments in Griffin and Davies (2018). From the longer-term history of that source-zone, it seems likely that a higher convergence rate might give more realistic results. This should be analysed and updated if necessary.

* We should come up with a better method for integrating paleotsunami and longer-term historical data. For earthquake magnitudes, this should be relatively straightforward (adjust the likelihood in the Bayesian update to allow for different data periods, with different observed uncertainties). Dealing with onshore historical/paleo data seems more challenging, but super important. It might be worth doing extensive inundation modelling for sites with good onshore records, so we can "take the PTHA onshore" (e.g. with Monte-Carlo methods) and then use the onshore model & observations in the Bayesian update.

* We could extend the DART testing using records in [this paper](http://dx.doi.org/10.1093/gji/ggt328).

* PTHA18 treats fault segmentation in only a limited way (as an optional local Mw-frequency model, but no real boundaries on rupture location). Consider adding in methods to treat segment boundaries (partially weighted). For example see a simple discussion [here](https://doi.org/10.1785/0120200219). We could also look at geometric methods for source placement, which imply non-trivial spatial magnitude distributions (e.g. [see here and related papers](10.1785/0120220175) ). 

* Determine a strategy to deal with the following situation: PTHA18 often uses SLAB2 geometry and max-depths based on upper limits from Berryman et al. 2015. In some cases (New Guinea, Makran) this implies a very large fault area. Combined with our modelled lower-bound of 10% coupling, this forces substantial seismic moment accumulation. Then, the model prefers high Mw-max (and also low coupling) to try to enforce seismic moment balance. In my subjective opinion, although these results might be possible, it would also be good to have a method that placed some weight on the possibility of very low coupling and/or a smaller active fault area.

* The outer-rise treatment should be considered more deeply and perhaps revised.
