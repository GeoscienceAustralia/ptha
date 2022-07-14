# SWALS validation tests

The script here can be used to run test problems in `../../examples/*` and `../../examples/nthmp/*` via:

    Rscript run_validations.R

This takes about an hour to run on my desktop (as of July 2022). It searches the test problem directories for scripts named `run_model.sh`, and runs them inside their parent directories. 

Most tests use openmp parallelism alone, but a significant minority also use MPI. The commands controlling openmp and MPI (i.e. number of threads, ranks, etc) are determined by the script [../../src/test_run_commands.sh](../../src/test_run_commands). For easier portablility between machines the latter script usually just sources some machine specific script, such as [this example](../../src/test_run_commands_basic).

For each problem the `run_model.sh` scripts print multiple `PASS` (or `FAIL`) statements and, in most cases, make figures in the test problem directories. `FAIL` statements indicate the test criteria were not satisfied.

Aside from genuine bugs, `FAIL` statements may reflect the ad-hoc nature of the test criteria. For example, a test problem might specify an error tolerance of 0.01, but depending on the case, somewhat larger errors might also be acceptable. Typically there is no objective error tolerance. For this reason `FAIL` statements may arise when testing new flow algorithms, or when porting the code to a different compiler, even if there are not genuine problems. Thus user judgement is required to interpret `FAIL`s -- are they genuine problems, or simply unimportant violations of the test criteria? The test problem figures will likely help to judge this.

