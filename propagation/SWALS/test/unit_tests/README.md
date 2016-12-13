Compile with:

    make -B -f make_unit_tests

and run with:

    ./unit_tests

Since the output is a bit verbose, you might use

    ./unit_tests | grep "FAIL"

to check just for failures.

You might want to adjust the preprocessing flags to check that the
tests still pass with different compilation options.


Running on NCI
--------------

On the NCI, more effort is required for installation. The key steps are:

1) Load gfortran version > 4.8 (I have been using 5.2.0 recently). Since
gcc/system is 4.4, and is loaded by default, we need to update to the newer
version with

    module swap gcc/system gcc/5.2.0

2) Now, you need to link to gdal and netcdf libraries which are compiled with 
the same version of gcc/gfortran. Also, the gdal install should link to the new
netcdf install. I have built these libraries for gcc/5.2.0 on NCI. The script
run_test_NCI.sh points the LD_LIBRARY_PATH to these libraries, and then builds
the code with make_test_NCI. The latter is a modified version of make_unit_tests, 
with the key change being that environmental variables related to GDAL and NETCDF
point to the locally installed library locations. If this is all set up ok, then
you can compile and run the code with:

    source run_test_NCI.sh

3) To compile application code, the same process is required, i.e.
* Load the correct gcc/gfortran
* Point LD_LIBRARY_PATH to them
* Ensure the makefile overrides the GDAL and NETCDF variables as above
