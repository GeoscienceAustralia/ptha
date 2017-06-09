#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

// See if we can do y <- a*x+y faster in C than in pure R
// (Possible that R get's confused by the copy?
SEXP axpy_c(SEXP y, SEXP a, SEXP x, SEXP n){

    double* y_c = REAL(y);
    const double* a_c = REAL(a);
    const double* x_c = REAL(x);
    const int* n_c = INTEGER(n);

    const int ionex = 1;
    const int ioney = 1;

    F77_CALL(daxpy)(n_c, a_c, x_c, &ionex, y_c, &ioney); 

    return R_NilValue;
}
