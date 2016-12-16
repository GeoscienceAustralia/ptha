Test case from the NTHMP suite.

This is an anaytical solution, but includes comparison with experimental data
as well. It is only test case in the NTHMP suite that is appropriate for the
linear shallow water equations (to my knowledge the others all require
non-linearity).

To run the code and plot the results, we assume R is installed on your system. Then do:

    make -B -f make_BP2_testcases
    Rscript plot.R # This runs the code and makes pdf plots
