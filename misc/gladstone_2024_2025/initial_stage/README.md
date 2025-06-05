# Polygons used to set the initial stage in the model

[The model](../swals/model_multidomain_design_mod.f90) points to the [override_initial_stages.csv](override_initial_stages.csv) file which has has two columns.
The first column is for filenames of csv files containing polygons.
The second column is to set the value of the initial fluid stage within the polygon.

1. Place any shape files in the [shapes](shapes) directory.
2. Call `Rscript shape_to_csv.R` to convert them into a format that swals can read.
3. This will automatically update the `override_initial_stages.csv` file to point to all these files with a stage of -20m. Manually change to your desired alternative stage. 
