#include <R.h>
#include <Rinternals.h>


// Interface to the fortran routine
void gauge_statistics(double* gauge_times, double* flow_var, 
    double* stage_threshold_for_arrival_time, int* ngauges, int* ntimes, 
    int* nflow_var, double* output_var); 

// R interface
SEXP gauge_statistics_c(SEXP gauge_times, SEXP flow_var, 
    SEXP stage_threshold_for_arrival_time, SEXP ngauges, SEXP ntimes, 
    SEXP nflow_var, SEXP output_var){

    // Get pointers to R's data
    double* gauge_times_c = REAL(gauge_times);
    double* flow_var_c = REAL(flow_var);
    double* stage_threshold_for_arrival_time_c = REAL(stage_threshold_for_arrival_time);
    int* ngauges_c = INTEGER(ngauges);
    int* ntimes_c = INTEGER(ntimes);
    int* nflow_var_c = INTEGER(nflow_var);
    double* output_var_c = REAL(output_var);

    // Call Fortran's gauge_statistics
    gauge_statistics(gauge_times_c, flow_var_c, stage_threshold_for_arrival_time_c,
        ngauges_c, ntimes_c, nflow_var_c, output_var_c);

    // Return something of class SEXP
    //return gauge_times; 
    return R_NilValue;

}
