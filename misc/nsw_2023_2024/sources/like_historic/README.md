This folder contains tsunami source models that are similar to historic events. 

Some are from [the Global Dissipation paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) with [source code in the ptha repository](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/nearshore_testing_2020/sources).

Others are from PTHA18 sources that I expect to perform OK (based on past experience).

For non-PTHA18 sources that are based on fault models and need a Kajiura filter applied, I use the script [apply_kajiura_to_rasters.R](apply_kajiura_to_rasters.R) to make a filtered model, such as
```
# This makes a raster in the input tif directory
Rscript apply_kajiura_to_rasters.R Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM.tif
# Result with kajiura smoothing is in
#   ./Chile1960/FujiSatake2013/Fuji_chile1960_sources_SUM_KAJIURA_SMOOTHED.tif
```


