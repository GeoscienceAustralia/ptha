# Define polygons in which the initial stage is set to a special value. 

Typically this is useful to remove water from areas below the model's initial sea level which are not connected to the ocean so should not be wet initially.

Make polygon shapefiles defining each region in subfolders of the `shapes` folder, then create the SWALS input file using `Rscript shape_to_csv.R`, which sets the initial water level in all polygons to -20 m (equivalent to making them dry).
