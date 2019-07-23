#include <R.h>
#include <Rinternals.h>


// Interface to the fortran routine
//    subroutine make_stage_vs_rate_curve_site( rate_curve_Mw, rate_curve_exrate, N1, dMw, &
//            unique_Mw, N2, event_Mw_to_unique_Mw_match, event_conditional_prob, N3, &
//            sorted_event_stages, sorted_stage_inds, &
//            output_stages, output_stage_exrates, N4) bind(C, 'make_stage_vs_rate_curve_site')
void make_stage_vs_rate_curve_site(double* rate_curve_Mw, double* rate_curve_exrate, int* N1, 
    double* dMw, double* unique_Mw, int* N2, int* event_Mw_to_unique_Mw_match, double* event_conditional_prob,
    int* N3, double* sorted_event_stages, int* sorted_stage_inds, double* output_stages, double* output_stage_exrates,
    int* N4);

// Interface for R 
SEXP make_stage_vs_rate_curve_site_c(SEXP rate_curve_Mw, SEXP rate_curve_exrate, SEXP N1,
    SEXP dMw, SEXP unique_Mw, SEXP N2, SEXP event_Mw_to_unique_Mw_match, SEXP event_conditional_prob,
    SEXP N3, SEXP sorted_event_stages, SEXP sorted_stage_inds, SEXP output_stages, SEXP output_stage_exrates,
    SEXP N4){

    // Get pointers to R's data
    double* rate_curve_Mw_c = REAL(rate_curve_Mw);
    double* rate_curve_exrate_c = REAL(rate_curve_exrate);
    int* N1_c = INTEGER(N1); 
    double* dMw_c = REAL(dMw); 
    double* unique_Mw_c = REAL(unique_Mw);
    int* N2_c = INTEGER(N2);
    int* event_Mw_to_unique_Mw_match_c = INTEGER(event_Mw_to_unique_Mw_match);
    double* event_conditional_prob_c = REAL(event_conditional_prob);
    int* N3_c = INTEGER(N3);
    double* sorted_event_stages_c = REAL(sorted_event_stages);
    int* sorted_stage_inds_c = INTEGER(sorted_stage_inds);
    double* output_stages_c = REAL(output_stages);
    double* output_stage_exrates_c = REAL(output_stage_exrates);
    int* N4_c = INTEGER(N4);


    // Call Fortran's make_stage_vs_rate_curve_site
    make_stage_vs_rate_curve_site(rate_curve_Mw_c, rate_curve_exrate_c, N1_c, 
        dMw_c, unique_Mw_c, N2_c, event_Mw_to_unique_Mw_match_c, event_conditional_prob_c,
        N3_c, sorted_event_stages_c, sorted_stage_inds_c, output_stages_c, output_stage_exrates_c,
        N4_c);

    // Return something of class SEXP
    return R_NilValue;

}
