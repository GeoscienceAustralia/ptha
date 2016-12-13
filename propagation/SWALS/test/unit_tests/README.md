Compile with:

    make -B -f make_unit_tests

and run with:

    ./unit_tests

Since the output is a bit verbose, you might use

    ./unit_tests | grep "FAIL"

to check just for failures.

You might want to adjust the preprocessing flags to check that the
tests still pass with different compilation options.
