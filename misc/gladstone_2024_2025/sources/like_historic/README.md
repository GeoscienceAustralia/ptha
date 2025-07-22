# Tsunami source generation

Sources a generated from tif files for each tsunami scenario. This folder contains tsunami source models that are similar to historic events. After inspecting many tide gauges, the following three events exhibited tsunami signal near Gladstone: Chile 2010, Solomon 2007 and Tohoku 2011. They all have published source inversions, as described in  [the Global Dissipation paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) with [source code in the ptha repository](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/nearshore_testing_2020/sources). However, those source inversions are not always used here.

A Kajiura filter was applied to these using [apply_kajiura_to_rasters.R](apply_kajiura_to_rasters.R) to make a filtered model, such as
```
Rscript apply_kajiura_to_rasters.R Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM.tif
```
This makes a raster in the input tif directory. The results with kajiura smoothing is in `./Chile2010/FujiSatake2013/Fuji_chile2010_sources_SUM_KAJIURA_SMOOTHED.tif`, with `_SMOOTHED` appended to the filename.

PTHA18 scenarios similar to the Solomon 2007 and Chile 2010 events are included. In the past they've had okay performance. For Solomon 2007, one scaled scenario performs well at Gladstone.
