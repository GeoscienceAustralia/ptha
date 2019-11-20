# These modules can be used to compile and run 'rptha' on NCI's gadi supercomputer
module load intel-compiler
module load netcdf/4.7.1
module load proj
module load geos
module load R/3.6.1
export PATH=/g/data1a/w85/tsunami/CODE/gadi/gcc_system/gdal/build/bin/:$PATH
export LD_LIBRARY_PATH=/g/data1a/w85/tsunami/CODE/gadi/gcc_system/gdal/build/lib/:$LD_LIBRARY_PATH
