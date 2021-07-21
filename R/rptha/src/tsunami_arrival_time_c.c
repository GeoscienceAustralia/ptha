#include <R.h>
#include <Rinternals.h>

// Interface to the fortran routine
//    subroutine tsunami_arrival_time_and_maxima(n, stage, time, msl, arrival_fraction_of_maxima, &
//        tsunami_maxima, tsunami_arrival_time) bind(C, name='tsunami_arrival_time_and_maxima')
void tsunami_arrival_time_and_maxima_fortran(int* n, double* stage, double* time, double* msl, double* arrival_fraction_of_maxima, 
    double* tsunami_maxima, double* tsunami_arrival_time);


SEXP tsunami_arrival_time_and_maxima_c(SEXP n, SEXP stage, SEXP time, SEXP msl, SEXP arrival_fraction_of_maxima, 
    SEXP tsunami_maxima, SEXP tsunami_arrival_time){
    
    // Get pointers to R's data
    int* n_c = INTEGER(n);
    double* stage_c = REAL(stage); 
    double* time_c = REAL(time);
    double* msl_c = REAL(msl);
    double* arrival_fraction_of_maxima_c = REAL(arrival_fraction_of_maxima);
    double* tsunami_maxima_c = REAL(tsunami_maxima);
    double* tsunami_arrival_time_c = REAL(tsunami_arrival_time);

    // Call Fortran's tsunami_arrival_time_and_maxima
    tsunami_arrival_time_and_maxima_fortran(n_c, stage_c, time_c, msl_c, arrival_fraction_of_maxima_c, 
        tsunami_maxima_c, tsunami_arrival_time_c);

    // Return something of class SEXP
    return R_NilValue;
}

// Interface to the fortran routine
//    subroutine exceedance_rate_given_maxima_and_arrival_time(ns, scenario_maxima, scenario_arrival_times, scenario_rates, &
//            nmax, tsunami_maxima_in_output, narrival, arrival_times_in_output, output_exceedance_rate_by_maxima_time) &
void exceedance_rate_given_maxima_and_arrival_time_fortran(
    int* ns, double* scenario_maxima, double* scenario_arrival_times, double* scenario_rates, 
    int* nmax, double* tsunami_maxima_in_output, int* narrival, double* arrival_times_in_output, 
    double* output_exceedance_rate_by_maxima_time); 

SEXP exceedance_rate_given_maxima_and_arrival_time_c(
    SEXP ns, SEXP scenario_maxima, SEXP scenario_arrival_times, SEXP scenario_rates,
    SEXP nmax, SEXP tsunami_maxima_in_output, SEXP narrival, SEXP arrival_times_in_output, 
    SEXP output_exceedance_rate_by_maxima_time){
    
    // Get pointers to R's data
    int* ns_c = INTEGER(ns);
    double* scenario_maxima_c = REAL(scenario_maxima); 
    double* scenario_arrival_times_c = REAL(scenario_arrival_times); 
    double* scenario_rates_c = REAL(scenario_rates); 
    int* nmax_c = INTEGER(nmax);
    double* tsunami_maxima_in_output_c = REAL(tsunami_maxima_in_output);
    int* narrival_c = INTEGER(narrival);
    double* arrival_times_in_output_c = REAL(arrival_times_in_output);
    double* output_exceedance_rate_by_maxima_time_c = REAL(output_exceedance_rate_by_maxima_time);

    // Call Fortran's exceedance_rate_given_maxima_and_arrival_time
    exceedance_rate_given_maxima_and_arrival_time_fortran(ns_c, scenario_maxima_c, scenario_arrival_times_c, scenario_rates_c,
        nmax_c, tsunami_maxima_in_output_c, narrival_c, arrival_times_in_output_c, output_exceedance_rate_by_maxima_time_c);

    // Return something of class SEXP
    return R_NilValue;
}
