# The files that were in this folder on NCI (Gadi) have been moved to NCI's MDSS (a tape archive), under project w85

The location in `mdss` is "tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/Greater_Perth/swals/OUTPUTS"

To view the files (top level) we can do: 
```
[gxd547@gadi-login-06 OUTPUTS]$ mdss ls tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/Greater_Perth/swals/OUTPUTS
```
which shows the files previously in this folder.
```
Fuji_andaman2004_24hrs_domain010322_lowres_timevaryingRealistic-full-ambient_sea_level_0.0
Fuji_andaman2004_24hrs_domain010322_timevaryingRealistic-full-ambient_sea_level_0.0
Fuji_sumatra2005_domains010322_newProcessDataToSend-full-ambient_sea_level_0.0
ptha18-GreaterPerth-sealevel60cm-lowres
ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres
```

The folders starting with `Fuji_` contain simulations like tsunamis due to historical earthquakes in Sumatra 2004 (in low resolution, and regular resolution) and Sumatra 2005. The other folders contain inundation models for random ptha scenarios, with the main results being `ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres` (and the other being used for convergence testing).

To re-run the analysis scripts, or do other analysis with the scenario files, you will have to copy the required data from `mdss` to the current directory.

## An important note on the multidomain log files.

The folder `ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres` contains subfolders for the outerrise-sunda (`random_outerrisesunda`) and sunda2 (`random_sunda2`) sources. 
* Each contains hundreds of sub-folders corresponding to unique model scenarios. 
* Inside those subfolders are tarred model output directories (each containing a multidomain). 

An example subfolder is `ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres/random_sunda2/ptha18_random_scenarios_sunda2_row_0111332_Mw_96_HS-full-ambient_sea_level_0.6`, representing a magnitude 9.6 scenario on the Sunda Arc thrust source. It contains the following files:
```
[gxd547@gadi-login-09 ptha18_random_scenarios_sunda2_row_0111332_Mw_96_HS-full-ambient_sea_level_0.6]$ ls
raster_output_files.tar  raster_output_files_with_speed.tar  RUN_20220402_065604331.tar.bz2
```
* The multidomain folder is `RUN_20220402_065604331.tar.bz2`
* The rasters for analysis are `raster_output_files_with_speed.tar`. The `raster_output_files.tar` also contains rasters but is missing speed.

Analogous files exist in every scenario subfolder.

When processing, the scenario folders also contained a file `multidomain_log_image_00000000000000000001.log`. Unfortunately I deleted these log files before archiving, mistakenly believing that NCI would not allow archiving of these relatively small files (as is stated by their documentation, but apparently not enforced). However, the log file is contained inside the `RUN_20220402_065604331.tar.bz2`, and can be extracted with
```
tar -x -f RUN_20220402_065604331.tar.bz2 ./RUN_20220402_065604331/multidomain_log_image_00000000000000000001.log
mv ./RUN_20220402_065604331/multidomain_log_image_00000000000000000001.log .
rm -r ./RUN_20220402_065604331/
```
and similarly for all other scenario folders. 

If you want to run analyses that assume the existence of these log files, you'll have to recreate them after manually copying from `mdss`.
