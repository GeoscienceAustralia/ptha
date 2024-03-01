# Tsunami model runs for Midwest region of WA, FY2023-2024
--------------------------------------------------

This folder contains code to run the tsunami models for earthquake sources defined in:
* Scenarios like historical events [(see here)](../../Greater_Perth/sources/like_historic/)
* Random PTHA18 scenarios [(see here)](../../Greater_Perth/sources/hazard/)

The tsunami model includes the area from Guilderton to Geraldton in high resolution. The model setup is similar to the [2023 version fo the Greater Perth model codes](../../greater_perth_revised2023/) but here we were able to use some new elevation datasets, refined breakwalls, better representations of some rivers, and include some previously excluded areas.

The individual codes are mostly similar to the [2023 version of the Greater Perth model codes](../../greater_perth_revised2023/). But the file structure has been substantially updated by using subfolders to make a clearer conceptual separation between different modelling steps. The documentation has also been significantly updated and split into separate files. In combination this should make it easier reuse and modify the codes. The scripts are setup to use the SapphireRapids nodes on Gadi. Following advice from NCI we now use a different invocation of `mpirun` which more explicitly binds cores to NUMA domains, which might lead to more consistent performance (previously we would occasionally see models running slowly, although it isn't 100% clear whether this was due to the `mpirun` invocation).

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
