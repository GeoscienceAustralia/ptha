# WA tsunami modelling (Kalbarri to Onslow region) 2025-2027

The codes here were mostly derived by modifying code from previous projects, including:
* [../SW_WA_2021_2024](../SW_WA_2021_2024) provided a starting point for the scenario sampling
* [../gladstone_2024_2025](../gladstone_2024_2025) and [../nsw_2023_2024](nsw_2023_2024) for the tsunami hydrodynamics among other things.

Key subfolders:

* [kalbarri_2_coralbay](kalbarri_2_coralbay) - Model covering the region from (somewhat south of) Kalbarri to Coral Bay, with high resolution in towns.

* [kalbarri_2_onslow](kalbarri_2_onslow) - A *preliminary and experimental version* of the model, covering the full area but at lower resolution and without breakwalls/inverts/etc. This was used to get the code working with variable rigidity PTHA18 scenarios, de-risk the main calculations implemented in other folders, check the sensitivity of the results to including horizontal components in the tsunami sources, and test the tidal adjustment technique.

* OTHERS HERE
