# Tsunami model runs for Midwest region of WA, FY2023-2024
--------------------------------------------------

This folder contains code to run the tsunami models, including:
* Scenarios like historical events [(see here)](../../Greater_Perth/sources/like_historic/)
* Random PTHA18 scenarios [(see here)](../../Greater_Perth/sources/hazard/)

The tsunami model includes the area from Seabird to Geraldton in high resolution.
It is similar to earlier work on the Greater Perth area, but here we were able to use some new elevation datasets, refined breakwalls, measured river depths and include some previously excluded areas.

In addition the code has been substantially updated to be easier to use and modify, the model setup is more efficient, and the scripts are setup to use the new SapphireRapids nodes on Gadi.

## Help with running the model
- An outline of the steps to run on the National Computational Infrastucture (NCI) cluster is provided in [README.HOW_TO_RUN.md](README.HOW_TO_RUN.md).
- A summary of the key codes is provided in [README.CODE_SUMMARY.md](README.CODE_SUMMARY.md).
- See comments within the code for further documentation.
- You might also want to review [../GENERAL_GUIDANCE_ON_MODEL_SETUP.md](../GENERAL_GUIDANCE_ON_MODEL_SETUP.md).

## Folder structure
- [pre_process](pre_process/): Scripts to prepare the model setup, earthquake sources and initial stages.
- [load_balance_files](load_balance_files/): Folder to store the load balance files.
- [run](run/): Scripts to run the model for different scenarios on the NCI cluster.
- [post_process](post_process/): Scripts for simple post-processing and data handling of the model results. Probabilistic analysis is done separately.
- [plots](plots/): Scripts to generate plots of the model results.
