# These modules can be used to compile and run 'rptha'
module load intel-compiler
module load gcc
module load netcdf
# Default proj conflicts with default gdal, which prevents installation of 'terra' and many other things.
module load proj/6.2.1
module load geos
module load R/4.2.1
module load gdal
module load udunits

# I also needed to patch the 's2' package before install, see here: https://github.com/r-spatial/s2/issues/199
# Then I could install stars
