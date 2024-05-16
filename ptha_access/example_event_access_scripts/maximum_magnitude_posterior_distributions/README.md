# Uncertainties in maximum magnitude on each source zone

## Quickstart

The file [Mw_max_percentiles_from_PTHA18.csv](Mw_max_percentiles_from_PTHA18.csv) provides a summary of the maximum-magnitude distribution on each source-zone in PTHA18. It was generated using the script [compute_Mw_max_uncertainty.R](compute_Mw_max_uncertainty.R). It contains the following columns:
* `source_zone` - the name of the source-zone in PTHA18. This might represent a full unsegmented source-zone (e.g. `puysegur2`, `kurilsjapan`) or one of several segments on a source zone (e.g. `sunda2_java`, `southamerica_peru`). For segments, the leading part of the name corresponds to the source-zone. 
* `is_a_segment` - TRUE/FALSE depending on whether the source representation is a segment. FALSE corresponds to unsegmented source representations.
* Columns with names like `p_0`, `p_0.01`, ... `p_0.99`, `p_1` that collectively represent the cumulative distribution function of the maximum magnitudes. e.g.
  * The column `p_0.1` contains, for each source-zone, the largest magnitude `X` such that `PR(Mw_max <= X) = 0.1` where `PR` denotes the posterior probability according to PTHA18. 
  * The column `p_0.5` contains, for each source-zone, the largest magnitude `X` such that `PR(Mw_max <= X) = 0.5` where `PR` denotes the posterior probability according to PTHA18. 
  * The column `p_0.9` contains, for each source-zone, the largest magnitude `X` such that `PR(Mw_max <= X) = 0.9` where `PR` denotes the posterior probability according to PTHA18. 
  * ... and so on

## Background 

The maximum earthquake magnitude on any given source zone is often of interest. Its value is usually uncertain and often controversial. We know that the smallest plausible value is greater than the largest historically observed earthquake. But this is often a very weak constraint, and history shows that the real maximum magnitude could be much greater. 

For example, prior to the 2004 Sumatra Andaman earthquake (Mw 9.2), the largest historical earthquake in that region was smaller than Mw 8 ( [Satake and Atwater 2007](https://doi.org/10.1146/annurev.earth.35.031306.140302) ). On the Cascadia subduction zone, the largest historical earthquake is currently well below Mw 8, yet sedimentary evidence suggests events nearing Mw 9 may occur every 500-1000 years on average (e.g. [Rong et al., 2014](http://www.bssaonline.org/content/104/5/2359.abstract), [Westby et al., 2022](https://doi.org/10.2138/gselements.18.4.251) ). On any particular source-zone, even if the largest known events are clearly high magnitude, there is no guarantee that they will reflect the maximum magnitude - which is potentially larger still. 

PTHA18 attempts to reflect uncertainties in maximum magnitudes (and other source-zone parameters) by using multiple alternative models with Bayesian weights. Details are in [this paper](https://doi.org/10.1007/s00024-019-02299-w) and [this report](http://dx.doi.org/10.11636/Record.2018.041). Our lower limits on Mw-max are based on historical events (or available paleo interpretations), while our upper limits on Mw-max are inferred from scaling relations using prediction errors to allow for compact ruptures (but capped to a limit of 9.6).

The approach leads to a (posterior) probability distribution for the maximum magnitude on each source zone. On source-zones with optional segmentation there are different distributions for the unsegmented model, and each segment.



