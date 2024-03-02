# Polygons used to set the initial stage in the model

The SWALS model points to the file [override_initial_stages.csv](override_initial_stages.csv) to set the initial stage in the model. This file is used to set the initial stage in the model to a specific value within a polygon. The file has two columns. The first column is for filenames of csv files containing polygons. The second column is to set the value of the initial fluid stage within the polygon. The file should not have trailing commas.

To set this up:
1. Place any shape files in the `shapes` directory.
2. Call `Rscript shape_to_csv.R` to convert them into a format that swals can read.
3. Manually alter the `override_initial_stages.csv` file to point to all these files and input your desired stage.
