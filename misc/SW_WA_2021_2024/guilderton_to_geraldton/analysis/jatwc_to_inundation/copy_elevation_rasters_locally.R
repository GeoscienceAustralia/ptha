#
# Make a local copy of the elevation rasters
# This is desirable because:
#   1) In practice we want to provide the elevation as well as the hazard results
#   2) We don't have to refer to relatively obscure elevation file paths, except in this script.
#

ptha_elevation_rasts = Sys.glob('../../swals/OUTPUTS/Sumatra2005_FujiiEtAl2020-full-ambient_sea_level_0.0/RUN_20231127_151302826/elevation0*.tif')

# Make a local elevation folder and copy the rasters to it
elev_dir = './elevation_in_model'
dir.create(elev_dir)
copy_worked = file.copy(from=ptha_elevation_rasts, to=paste0(elev_dir, '/', basename(ptha_elevation_rasts)))
if(!all(copy_worked)) stop('Copy of elevation data failed')

# Make a vrt as well
mydir = getwd()
setwd(elev_dir)
system('gdalbuildvrt -resolution highest all_elevation_combined.vrt elevation0*.tif')
setwd(mydir)
