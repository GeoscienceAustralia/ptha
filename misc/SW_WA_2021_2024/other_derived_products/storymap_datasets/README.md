# Datasets for storymap

This folder contains data to depict the Sumatra 2004 tsunami in the storymap.
* [Sumatra2004_max_stage_greater_perth_model](Sumatra2004_max_stage_greater_perth_model) contains rasters depicting the maximum modelled waterlevel. A vrt file is included so they can be viewed as a single file.
* [sunda2_unit_source_grid](sunda2_unit_source_grid) contains a shapefile depicting the Sunda Arc. This originates from the 2018 Australian PTHA.


## How to get the datasets (assuming working on the w85 project on NCI)

1. Clip the rasters in dry areas with 
```
Rscript get_sumatra2004_rasters_clipped_in_dry_regions.R
```

2. Go inside `Sumatra2004_max_stage_greater_perth_model`, and combine the tifs in a vrt with 
```
gdalbuildvrt -resolution highest All_Sumatra2004_max_stage_dry_clipped.vrt DRY_CLIPPED_max_stage_domain_*.tif
```

3. Copy the Sunda Arc unit source grid shapefile from PTHA18 to `sunda2_unit_source_grid` folder
```
cp -r /g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/sunda2/EQ_SOURCE/unit_source_grid sunda2_unit_source_grid
```
