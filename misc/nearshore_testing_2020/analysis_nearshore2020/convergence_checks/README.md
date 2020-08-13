The scripts here are used to do convergence tests for 2 different source models.

For each source model, we compare time-series gauges of 2 tsunami models. The model setups are identical except that one model has cells with side-length="1/2 of the other model". We focus on nearshore gauges that have good quality data for that particular source model (i.e. focus on the sites where we will use the model results).

Run the code with 
```
    Rscript compare_gauges_fujii.R
    Rscript compare_gauges_yamakazi.R
```
and then visually examine the pdfs plot to compare the model results. 

Some statistics are computed inside the above scripts, and to ease later reproducibility I have copied key results in the scripts as comments.
