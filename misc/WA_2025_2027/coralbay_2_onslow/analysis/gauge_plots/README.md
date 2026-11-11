# Validation tests

Compare published tsunami sources propagated through this model set up to the tide gauge records for validation.

There are four historical scenarios to run for:
    - Sumatra 2004 (Fuji and Satake 2007, also Fujii 2021),
    - Sumatra 2005 (Fujii et al. 2020),
    - Java 2006 (Fujii and Satake 2006) and
    - Sumatra 2007 (Fujii and Satake 2007).

First read in and process the gauge data for each event using
```{bash}
Rscript process_gauges_sumatra2004.R path/to/model/multidomain
```
This finds the gets the nearest gauges in the model that are in the water and within 1 km to the saved gauge locations for this event.
It stores the relevant gauge stage and times from the model and historic observations in an `rds` file in the multidomain directory.
Run all such `process_gauges_*.R` files.

Then plot the model's gauges against historic gauges with e.g. `Rscript plot_gauges_sumatra2004.R`.
Modify this (and similar) scripts to change which gauges are plotted and adjust the figure size.
Optionally pass 1, 2 or 3 command arguments `YLIM`, `YLIM XMAX` or `YLIM XMIN XMIN` to specify the plot limits, e.g.
```{bash}
Rscript plot_gauges_sumatra2004.R 0.6
```

Previously this directory was in [../../swals](../../swals/).
