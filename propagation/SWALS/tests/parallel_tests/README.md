# Parallel unit-tests

Before running these tests you should check that the unit tests pass, see [here](../unit_tests). 

The parallel unit tests can be run with

    source run_test.sh > outfile.log

This will run the same tests multiple times, with different numbers of MPI ranks, and write a large number of `PASS` or `FAIL` statements to `outfile.log`. We do not expect any `FAIL`s, but the log will contain a few other statements (e.g. multidomain summary information). 

These tests check for accuracy, not parallel performance. We do not expect speed-ups proportional to the number of ranks, because the etst code is not written to scale well. But it might detect problems with your MPI setup.

To quickly count the `FAIL` statements I suggest

    cat outfile.log | grep FAIL | wc -w

which should return `0` if there are no failures. The same approach can be applied to `PASS` statements

    cat outfile.log | grep FAIL | wc -w

and at the time of writing returns `4215`.
