export R_MAKEVARS_USER='~/.R/Makeconf'
export LD_LIBRARY_PATH=/short/w85/gxd547/PTHA_Aust/SOURCE/proj4/proj.4.9.1/proj-4.9.1/build/lib:$LD_LIBRARY_PATH
export PATH=/short/w85/gxd547/PTHA_Aust/SOURCE/proj4/proj.4.9.1/proj-4.9.1/build/bin:$PATH
export _R_CHECK_FORCE_SUGGESTS_=FALSE

module load R/3.3.0 

module load intel-cc
#module swap gcc gcc/5.2.0
module load gcc
module load gdal/1.10.1
## 15/06/2016 we needed a local version of proj
## Need to edit PATH and LD_LIBRARY_PATH -- which is done above
#module load proj/4.8.0
module load netcdf/4.3.0
module load hdf4
module load hdf5
module load geos/3.4.2
module load curl/7.49.1
