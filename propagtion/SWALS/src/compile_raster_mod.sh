
# Note: It is important that the include / linking flags go at the end
#       Otherwise the compilation fails. 
FORTRAN=gfortran
C=gcc

# Global fortran variables
$FORTRAN -c global_mod.f90

# C function to read rasters with gdal
$C -c read_raster_c.c `gdal-config --cflags` 

# Fortran interface to C function to read rasters
$FORTRAN -c read_raster_mod.f90

# Test program for raster IO
#$FORTRAN -o test_read_raster test_read_raster.f90 read_raster_mod.o read_raster_c.o global_mod.o `gdal-config --libs`
