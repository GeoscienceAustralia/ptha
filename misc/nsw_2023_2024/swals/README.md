Below we show the initial model setup and running various historical scenarios
and stress tests.

```bash
# Make a load balance file in the regular way
qsub run_test_model_sapphirerapids
# Then
#   cd to output directory, 
#   run the load_balance_script.R
#   move it to the load_balance_files directory
#   update the multidomain_design_control_NNL4_1arcminoffshore_elevationsourceindex.nml

# Test runs
qsub run_test_model_with_loadbalance_and_elevationsource_sapphirerapids.sh

#
# Stress tests
#
qsub run_model_small_source.sh
qsub run_model_extreme_source.sh
# Check the energy stats, max-stage/flux/speed & last_timestep UH/VH for the previous runs

# Convergence testing was done with the initial model setup (for interim report)
# Files are in the folder.

#
# Historical tests
# (ignoring a few alternative tests for which run scripts exist, that were used for other things)
#
qsub run_model_chile1960_FujiiSatake.sh
qsub run_model_chile1960_HoEtAl.sh
qsub run_model_Puysegur2009_PTHA18.sh
qsub run_model_Chile2010_Lorito11.sh
qsub run_model_solomon2007.sh
qsub run_model_Chile2015_Williamson17.sh
qsub run_model_tohoku2011_Yamazaki.sh
qsub run_model_Kermadec2021_Romano.sh
qsub run_model_Newhebrides2021_Gusman_Kajiura.sh

#
# Plot the historical runs, when they finish
#
# FIXME show code here

```


Here we show the process of running the first batch of random scenarios. They
are associated with ID710.5, which is the PTHA18 hazard point ID that was used
in the importance sampling.
```bash
#
# Create PTHA scenario run qsub scripts
#
Rscript create_random_ptha_qsub_scripts_sealevel110cm.R

# Submit the scripts, storing the qsub files in a directory.
# In practice we submit these in batches, not all at once
mkdir submitted_qsub_files_ID710.5_B
for i in run_ptha18_NSW2023b_ID710.5_sealevel110cm_*.sh; do echo $i; qsub $i; mv $i submitted_qsub_files_ID710.5_B; done

# Wait for the above to finish. When they do, you'll see files with similar
# names and '.e123456' or '.o123456' extensions where the number is a process
# ID. 
# Scan the 'dot o' files for any runs that failed (nonzero Exit status)
grep "Exit" *.sh.o*
# or took a surprising length of time
grep "Wall" *.sh.o*

# Also check the log statistics using code in ../analysis/check_log_files/

# For any problematic runs, check the associated 'dot e / dot o' files and the logs,
# and try to work out what happened. 
#
# I had a few failures, but all seemed hardware related (based on the messages)
# and rerunning the models made them work OK. This included
#   - Re-submitting run_ptha18_NSW2023b_ID710.5_sealevel110cm_1094.sh
#   - Running a single job from run 1094 (see the script run_ptha18_NSW2023b_ID710.5_sealevel110cm_1093_fixfailed.sh).
# These all had messages that were indicative of network failure.
# In my case, none of the re-runs demanded tweaking of the model itself.

# Edit the script to create the rasters -- ensure the file path matches your runs
# Then run it
qsub run_create_tarred_rasters_from_tarred_multidomains.sh

# At this point you can proceed with all calculations in ../analysis
```

For the NSW model, we have ran multiple sets of scenarios (using different
offshore sites to guide the importance sampling). These runs proceed similarly
to above, but need some changes to the scripts. In particular:
* [create_random_ptha_qsub_scripts_sealevel110cm_ID1315p5.R](create_random_ptha_qsub_scripts_sealevel110cm_ID1315p5.R) is a variant on the random scenario creation script for scenario batch ID1315.5.
  * If any scenarios were already run in batch ID710.5, then this script will create outputs for scenarios ID1315.5 using a symbolic link to the existing runs.

There is also a distinct script to run the raster creation [run_create_tarred_rasters_from_tarred_multidomains_ID1315.5.sh](run_create_tarred_rasters_from_tarred_multidomains_ID1315.5.sh).
