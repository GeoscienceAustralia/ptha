# SWALS unit tests

These tests perform elementary checks on SWALS modules. They are the first tests you should run when making changes to SWALS or installing on a new machine. Unit tests specifically for MPI are provided in the [parallel_tests](../parallel_tests).

The unit tests are very quick and reduce the chance that we break existing code. They also help detect issues with your compiler or machine setup, which need to be resolved before running real models. 

The unit tests do not check the performance of the flow algorithms. That is done by the [validation tests](../validation_tests), which take considerably longer to run. 

Compile with:

    make -B -f make_tests

and run with:

    ./unit_tests > outfile.log

Since the output is a bit verbose, you might use

    grep FAIL outfile.log | wc -l

to count failures (should be 0), and 

    grep PASS outfile.log | wc -l

to count passes (383 at the time of writing).

You might want to adjust the preprocessing flags to check that the
tests still pass with different compilation options.

