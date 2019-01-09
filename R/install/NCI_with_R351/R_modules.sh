# These modules were used by NCI's build, see
# /apps/R/3.5.1/README.nci

# Load them to set-up R/3.5.1 for use with the 'ptha' package on NCI
# See also R_351_package_installs.txt for package installation

module unload intel-fc
module unload intel-cc
module unload openmpi
module load intel-fc/2018.3.222
module load intel-cc/2018.3.222
module load intel-mkl/2018.3.222
module load java/jdk-10.0.1
module load zlib/1.2.8
module load bzip2/1.0.6
module load xz/5.2.2
module load pcre/8.38 # Needs version less than 10.2
module load curl/7.49.1

# GD ADDITIONS
module load gdal
module load geos
module load proj
module load netcdf

module load R/3.5.1
