Here we create walls WITHOUT using a DEM to set the elevation. This is useful when we know a wall elevation via other means.

Run with
```
Rscript make_manual_walls.R
```

We set the wall elevation using the last part of the directory name (after `_`).

There shouldn't be any shapefiles in this directory (but they can be inside subdirectories).

