## Plot the tide gauges using the following
Process the gauges with the appropriate variant of:
```
    Rscript process_gauges_sumatra2004.R ../OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20231122_120644963/
```
Then plot them using:
```
Rscript plot_gauges_perth_sumatra2004.R ../OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20231122_120644963/gauges_plot_2004-12-26_Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0.RDS
```
Then do the same thing for the low resolution version to compare.
Then repeat the process for the 2005 event.
