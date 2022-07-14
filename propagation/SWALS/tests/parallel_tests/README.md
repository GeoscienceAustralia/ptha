# Parallel unit-tests

Run with

    source run_model.sh > outfile.log

This will run the tests multiple times with different numbers of MPI ranks, and write a large number of `PASS` or `FAIL` statements to `outfile.log`. We do not expect any `FAIL`s, but the log will contain a few other statements (e.g. multidomain summary information). 

To quickly count the `FAIL` statements I suggest

    cat outfile.log | grep FAIL | wc -w

which should return `0` if there are no failures. The same approach can be applied to `PASS` statements

    cat outfile.log | grep FAIL | wc -w

and at the time of writing returns `4163`.
