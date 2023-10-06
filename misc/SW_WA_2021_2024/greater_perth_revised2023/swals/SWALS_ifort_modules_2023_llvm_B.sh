module load pbs

#module load intel-mkl/2023.2.0
module load intel-compiler-llvm/2023.2.0
module load openmpi
module load gdal/3.5.0 # Later version is missing ZLIB compression
module load netcdf/4.8.0 # Need to use non-default netcdf (conflict with other libraries?)
#module load proj/8.1.1
module load proj/6.2.1 # Given the above GDAL, R's terra package needs this.
