#
# Ensure the input elevation files are downsampled to the model resolution
# BEFORE we make the hazard points. 
#
export GA250_FILE='../../DATA/ELEV/GA250_1m/GA_250_1mx1m.tif'
#export GEBCO_2014_FILE='../../DATA/ELEV/GEBCO_2014_1m/GEBCO_2014_1minx1min.tif'
#
# Instead of pure GEBCO2014, use the merged global dem [as used in the model]
#
export GEBCO_2014_FILE='../../DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'

# Points in 20m depth with 20km spacing around Australia
# Arguments are 'input raster', 'desired hazard point elevation', 'desired hazard point spacing',
#   'maximum number of pixels offshore ', 'output_folder'
Rscript make_hazard_pts.R $GA250_FILE -20.0 20.0 60 3 OUTPUTS/OUTPUT_GA250_20m ''

# Points in 100m depth with 20km spacing around Australia
# Arguments are 'input raster', 'desired hazard point elevation', 'desired hazard point spacing',
#   'maximum number of pixels offshore ', 'output_folder'
Rscript make_hazard_pts.R $GA250_FILE -100.0 20.0 60 10 OUTPUTS/OUTPUT_GA250_100m ''

# Points in 1000m depth with 20km spacing around Australia. Interesting
# to see if these converge better
Rscript make_hazard_pts.R $GA250_FILE -1000.0 20.0 180 10 OUTPUTS/OUTPUT_GA250_1000m ''

# Global points. Can come in handy for our other work. 50km spacing
Rscript make_hazard_pts.R $GEBCO_2014_FILE -100.0 50.0 60 10 OUTPUTS/OUTPUT_GEBCO2014_100m ./INPUTS/EXTRA_HAZ_LINES/EXTRA_HAZ_LINES.shp

# Gridded points around Australia. 
Rscript make_gridded_hazard_points.R

#
# Then run a script to combine them
#
Rscript combine_hazard_points.R

#
# Then order the points for fast access
#
Rscript reorder_mhp.R
