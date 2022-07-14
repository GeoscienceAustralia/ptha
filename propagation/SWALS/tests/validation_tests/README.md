# SWALS validation tests

The script here can be used to run test problems in `../../examples/*` and `../../examples/nthmp/*` via:

    Rscript run_validations.R

It takes about an hour to run on my desktop (as of 15/07/22). It searches directories for a `run_model.sh` script, and runs those it finds from inside their parent directory. 

Most tests use openmp parallelism alone, but a significant minority also use MPI. The commands controlling openmp and MPI (i.e. number of threads, ranks, etc) are determined by the script [../../src/test_run_commands.sh](../../src/test_run_commands).

For each problem the `run_model.sh` scripts print multiple `PASS` (or `FAIL`) statements (and usually make figures in the test problem directories). `FAIL` statements indicate the test criteria were not satisfied.

Aside from genuine bugs, `FAIL` statements might also reflect the ad-hoc nature of the test criteria. For example, a particular problem might specify an error tolerance of 0.01, but depending on the case, somewhat larger errors might also be acceptable. Typically there is no objective error tolerance. For this reason, `FAIL` statements may arise when testing new flow algorithms, or even when porting the code to a different compiler, even if there are not genuine problems. Thus user judgement is required to interpret `FAIL`s -- are they genuine problems, or simply unimportant violations of the test criteria?

