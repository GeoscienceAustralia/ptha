# Before running, need to change LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/short/w85/tsunami/CODE/alternative_gcc_lib/gcc_5.2.0/gdal/install/lib:/short/w85/tsunami/CODE/alternative_gcc_lib/gcc_5.2.0/netcdf/install/lib:$LD_LIBRARY_PATH

make -B -f make_test_NCI
./unit_tests
