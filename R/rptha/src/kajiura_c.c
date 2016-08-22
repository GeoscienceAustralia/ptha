#include <R.h>
#include <Rinternals.h>

// SUBROUTINE kajiura_convolution(depth_inv, initial_deformation_padded, filterXYr, &
//        kajiuraGmax, lfx, lfy, lnx, lny, output) BIND(C, name='kajiura_convolution')
void kajiura_convolution(double* depth_inv, double* initial_deformation_padded, double* filterXYr,
    double* kajiuraGmax, int* lfx, int* lfy, int* lnx, int* lny, double* output);


// Functions which we want to use with R's `.Call` interface must take all
// parameters as a SEXP, and return a SEXP.
SEXP kajiura_convolution_c(SEXP depth_inv_SEXP, SEXP initial_deformation_padded_SEXP, 
    SEXP filterXYr_SEXP, SEXP kajiuraGmax_SEXP, SEXP lfx_SEXP, SEXP lfy_SEXP, SEXP lnx_SEXP,
    SEXP lny_SEXP, SEXP output_SEXP){

    // Get pointers to vector data
    double* depth_inv = REAL(depth_inv_SEXP); 
    double* initial_deformation_padded = REAL(initial_deformation_padded_SEXP); 
    double* filterXYr = REAL(filterXYr_SEXP); 
    double* kajiuraGmax = REAL(kajiuraGmax_SEXP); 
    int* lfx = INTEGER(lfx_SEXP);
    int* lfy = INTEGER(lfy_SEXP);
    int* lnx = INTEGER(lnx_SEXP);
    int* lny = INTEGER(lny_SEXP);
    double* output = REAL(output_SEXP);

    // Call the fortran function
    kajiura_convolution(depth_inv, initial_deformation_padded, filterXYr,
        kajiuraGmax, lfx, lfy, lnx, lny, output);

    // Return the R object
    return output_SEXP;
    // return;
}
