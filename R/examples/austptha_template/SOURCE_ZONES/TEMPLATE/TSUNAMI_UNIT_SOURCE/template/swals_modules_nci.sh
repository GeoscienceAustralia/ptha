module swap gcc/system gcc/5.2.0

# Since we link with non-standard netcdf/gdal installs, we need to provide the LD_LIBRARY_PATH at runtime
export LD_LIBRARY_PATH=/short/w85/tsunami/CODE/alternative_gcc_lib/gcc_5.2.0/netcdf/install/lib/:/short/w85/tsunami/CODE/alternative_gcc_lib/gcc_5.2.0/gdal/install/lib/:$LD_LIBRARY_PATH

