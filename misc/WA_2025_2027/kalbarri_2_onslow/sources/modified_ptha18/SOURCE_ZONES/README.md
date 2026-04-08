Here subfolders are of the form `__SOURCE_ZONE_NAME__/EQ_SOURCE_bathyslope2025/` and contain code to create unit sources including the effect of horizontal components. 

To run it, go into each of these directories and 
```
qsub run_produce_unit_sources.PBS
```

Once all have been run, they can be moved to the `fj6` gdata (making the results accessible remotely) with [copy_bathyslope_unit_sources_to_fj6.R](copy_bathyslope_unit_sources_to_fj6.R).

