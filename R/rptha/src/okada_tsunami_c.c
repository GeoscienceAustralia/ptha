#include <R.h>
#include <Rinternals.h>

//    xout = .Fortran('fault_disp', as.double(alp),as.double(elon),as.double(elat),as.double(edep),
//        as.double(strk),as.double(dip),as.double(lnth),as.double(wdt),
//        as.double(disl1),as.double(disl2),as.double(rlon),as.double(rlat),as.double(dstmx),
//        as.double(edsp),as.double(ndsp),as.double(zdsp),as.integer(m),as.integer(n), 
//        DUP=TRUE, # DUP=FALSE is deprecated :(
//        PACKAGE='rptha')

// Interface for the fortran routine fault_disp
void fault_disp_(double* alp, double* elon, double* elat, double* edep, double* strk,
    double* dip, double* lnth, double* wdth, double* disl1, double* disl2, double* rlon,
    double* rlat, double* dstmx, double* edsp, double* ndsp, double* zdsp, int* m, int* n);

// R interface in C
SEXP fault_disp_c(SEXP alp_SEXP, SEXP elon_SEXP, SEXP elat_SEXP, SEXP edep_SEXP, SEXP strk_SEXP,
    SEXP dip_SEXP, SEXP lnth_SEXP, SEXP wdth_SEXP, SEXP disl1_SEXP, SEXP disl2_SEXP, SEXP rlon_SEXP,
    SEXP rlat_SEXP, SEXP dstmx_SEXP, SEXP edsp_SEXP, SEXP ndsp_SEXP, SEXP zdsp_SEXP, SEXP m_SEXP, SEXP n_SEXP){

    // Get pointers to vector data
    double* alp = REAL(alp_SEXP); 
    double* elon = REAL(elon_SEXP); 
    double* elat = REAL(elat_SEXP); 
    double* edep = REAL(edep_SEXP);
    double* strk = REAL(strk_SEXP); 
    double* dip = REAL(dip_SEXP);
    double* lnth = REAL(lnth_SEXP);
    double* wdth = REAL(wdth_SEXP);
    double* disl1 = REAL(disl1_SEXP);
    double* disl2 = REAL(disl2_SEXP);
    double* rlon = REAL(rlon_SEXP);
    double* rlat = REAL(rlat_SEXP);
    double* dstmx = REAL(dstmx_SEXP);
    double* edsp = REAL(edsp_SEXP);
    double* ndsp = REAL(ndsp_SEXP);
    double* zdsp = REAL(zdsp_SEXP);
    int* m = INTEGER(m_SEXP); 
    int* n = INTEGER(n_SEXP);

    // Call the fortran
    fault_disp_(alp, elon, elat, edep, strk, dip, lnth, wdth, disl1, disl2, rlon, rlat, 
        dstmx, edsp, ndsp, zdsp, m, n); 

    // Return something of type SEXP -- note that a range of things were modified, but this 
    // is just required for R
    return edsp_SEXP;
}

