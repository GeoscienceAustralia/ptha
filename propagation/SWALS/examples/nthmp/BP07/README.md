# NTHMP benchmark problem 7: Laboratory model of runup at Monai

We simulate the experiment of Matsuyama and Tanaka (2001)

  * Matsuyama and Tanaka (2001), An experimental study of the highest run-up height in the 1993 Hokkaido Nansei-oki earthquake tsunami. ITS 2001 Proceedings, Session 7, Number 7-21.

The experiment represents runup at Monai during the 1993 Okushiri tsunami event, which was higher than on much of the neighbouring coast. For the full field problem see [../BP09](../BP09). The test data and problem description are available in [Randy LeVeque's repository](https://github.com/rjleveque/nthmp-benchmark-problems/tree/master/BP07-DmitryN-Monai_valley_beach). This includes the flume geometry, boundary time-series, measured internal gauge time-series, runup photograph time-series, and measurements of runup maxima at several sites (the latter provided for multiple experimental runs).

The [SWALS model](monai.f90) simulates this problem using a nested grid in the main runup zone, with boundary forcing from the measured time-series.

Figure 1 compares the modelled and observed time series at several gauges. The agreement is good, although there are some minor differences in the early part of the time series (around 5-10s) suggesting the experiment had background waves that do not feature in the boundary forcing.

![Figure 1: Modelled and observed gauge time-series at three sites.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/gauges_plot.png)

The figures below compare the model wet-dry front with experimental snapshots near the runup maxima. They were taken every 0.5s starting from about 15s. The benchmark description notes some uncertainty in the start time, and herein a value of 15.1s was found to give reasonable consistency of the model and experiment. The model contour is defined from a 2mm depth threshold (also suggested in the [GEOCLAW NTHMP tests for this problem](https://depts.washington.edu/clawpack/links/nthmp-benchmarks/geoclaw-results.pdf)).

![Figure 2: Snapshot of model and observations in main runup area.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/snapshot_time_15.1.png)
![Figure 3: Snapshot of model and observations in main runup area.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/snapshot_time_15.6.png)
![Figure 4: Snapshot of model and observations in main runup area.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/snapshot_time_16.1.png)
![Figure 5: Snapshot of model and observations in main runup area.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/snapshot_time_16.6.png)
![Figure 6: Snapshot of model and observations in main runup area.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/BP07/snapshot_time_17.1.png)

The test script also checks that the modelled runup maxima at three sites is within the range of the experimental results. Those results are written to [this csv file](model_vs_experiment_test_result.csv).
