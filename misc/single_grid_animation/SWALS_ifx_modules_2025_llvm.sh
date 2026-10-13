#
# Note for the intel-llvm compiler
#   To use the llvm-ar linker when compiling SWALS, put
#     SWALS_AR :=/apps/intel-oneapi/compiler/2023.2.0/linux/bin-llvm/llvm-ar rcs
#   in the application's makefile. This will prevent errors when making the library.
#   (albeit I don't always use the library when compiling -- often just link directly with object files).
#

#export PATH=/apps/intel-oneapi/compiler/2023.2.0/linux/bin-llvm/:${PATH}

module load pbs

module load intel-compiler-llvm/2025.2.0
module load openmpi/4.1.7 # Specify the version
module load gdal/3.5.0 # Later version is missing ZLIB compression
module load netcdf/4.8.0 # Need to use non-default netcdf (conflict with other libraries?)
#module load proj/8.1.1
module load proj/6.2.1 # Given the above GDAL, R's terra package needs this.
