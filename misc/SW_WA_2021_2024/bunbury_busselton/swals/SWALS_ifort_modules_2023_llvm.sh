module load pbs

module load intel-mkl/2023.0.0

module load intel-compiler-llvm/2023.0.0
#module load intel-compiler/2021.7.0
#module load intel-compiler

#module load openmpi/4.1.4
module load openmpi

module load gdal/3.5.0
#module load gdal/3.0.2
#module load gdal

# Some issue with default netcdf (conflict with other libraries?)
module load netcdf/4.8.0
#module load netcdf

# Specify proj to avoid incompatibility with default gdal
module load proj/8.1.1
#module load proj/6.2.1
#module load proj
