# Analysis of the tidal adjustment technique

The SWALS model tries to simulate the tsunami atop a static HAT sea level, and
since HAT varies, this is achieved by modifying the elevation as in Macaulay
and Davies (2025). 

Since the elevation adjustment approach is heuristic, we need to test it in each case.

For the kalbarri2onslow model we simulated 4 static tidal levels in addition to the "elevation adjusted" simulation.
- 0.65m (like Kalbarri, and close to Steep Point and Geraldton)
- 0.95m (like Carnarvon and Coral Bay)
- 1.17m (Exmouth)
- 1.42m (Onslow)

A number of the runs failed because the timestep wasn't small enough (this was
dealt with in the hazard runs by re-running with shorter timesteps, but I
forgot to apply that logic to the tidal tests). Scenarios that failed in 1 or more
tides had indices [6, 8, 12, 16, 18, 19], with dirnames below (for a single tide only)
- `ptha18_tidal_testing_sunda2_row_0107985_Mw_94_HS-full-ambient_sea_level_0`
- `ptha18_tidal_testing_sunda2_row_0108100_Mw_94_HS-full-ambient_sea_level_0`
- `ptha18_tidal_testing_sunda2_row_0110660_Mw_96_HS-full-ambient_sea_level_0`
- `ptha18_tidal_testing_sunda2_row_0110734_Mw_96_HS-full-ambient_sea_level_0`
- `ptha18_tidal_testing_sunda2_row_0110907_Mw_96_HS-full-ambient_sea_level_0`
- `ptha18_tidal_testing_sunda2_row_0109369_Mw_95_HS-full-ambient_sea_level_0.65`

The simplest check would be to just use the other 14 cases. 

## How to run
```
# If the multidomain folders are tarred, use this to untar them.
Rscript untar_run_folders.R

# Extract the rasters at chosen sites and sea levels
Rscript extract_rasters_at_sites.R

# Make site-specific PDF's comparing each scenario with static vs varying elevation
# (stored in `./rasters_at_sites/SITE_NAME/*.pdf`) as well as summaries of the
# depth difference (stored in the current folder).
Rscript plot_results.R 
```
