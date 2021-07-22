#' arrival time calculation
#'
#' Interface to the fortran routines
#'
#' @param time vector of times
#' @param stage vector of stage values
#' @param msl mean sea level. length(msl) == 1. Tsunami size is stage-msl
#' @param arrival_fraction_of_maxima We say the tsunami has arrived when (stage-msl) > arrival_fraction_of_maxima*max(stage-msl)
#' @return A vector of length 2 giving max(stage), and the time of arrival.
#' @export
#' @examples
#' time = seq(0, 10, len=100)
#' stage = sin(time)
#' max_and_arrival = tsunami_maxima_and_arrival_time(time, stage)
#' stopifnot(max(stage) == max_and_arrival[1])
#' ind = min(which(time == max_and_arrival[2]))
#' stopifnot(stage[ind] > 0.2 & stage[ind-1] < 0.2)
#'
tsunami_maxima_and_arrival_time<-function(
    time,
    stage,
    msl=0.0,
    arrival_fraction_of_maxima = 0.2){

    # Make sure the following will be treated as vectors
    stopifnot(length(stage) > 1)
    stopifnot(length(time) == length(stage))

    # The following should not be vectors
    stopifnot(length(msl) == 1)
    stopifnot(length(arrival_fraction_of_maxima) == 1)

    # Ensure these are vectors of length=1
    input_msl = as.double(c(msl))
    input_arrival_fraction_of_maxima = as.double(c(arrival_fraction_of_maxima))
    output_tsunami_maxima = as.double(c(0.0))
    output_tsunami_arrival_time = as.double(c(0.0))

    # This is not a vector
    n = as.integer(length(time))

    .Call("tsunami_arrival_time_and_maxima_c",
          n, stage, time, input_msl, input_arrival_fraction_of_maxima,
          output_tsunami_maxima, output_tsunami_arrival_time)

    return(c(output_tsunami_maxima, output_tsunami_arrival_time))

}
   
#' exceedance rate computation
#'
#' Interface to the fortran routines
#' 
#' @param scenario_maxima tsunami maxima for all scenarios
#' @param scenario_arrival_times tsunami arrival times for all scenarios
#' @param scenario_rates individual rates for all scenarios
#' @param tsunami_maxima_in_output vector giving the tsunami maxima values for which we will compute exceedance-rates
#' @param arrival_times_in_output vector giving tsunami arrival times for which we will compute non-exceedance-rates
#' @return a matrix one row for each tsunami_maxima_in_output, and one column for each arrival_times_in_output. The
#' entries give the rate of tsunamis with maxima exceeding tsunami_maxima_in_output, and arrival time before arrival_times_in_output.
#' @export
#' @examples
#' maxima = runif(100)
#' arrival_times = runif(100)
#' scenario_rates = rep(0.1, length=100)
#' tsunami_maxima_in_output = seq(0, 1, len=11)
#' arrival_times_in_output = seq(0, 1, len=9)
#' exrate_matrix = exceedance_rate_given_maxima_and_arrival_time(
#'     maxima,
#'     arrival_times,
#'     scenario_rates,
#'     tsunami_maxima_in_output,
#'     arrival_times_in_output)
#' 
#' # Compute expected solution
#' back_computed_sol = exrate_matrix * 0
#' for(i in 1:ncol(exrate_matrix)){
#'     for(j in 1:nrow(exrate_matrix)){
#'         back_computed_sol[j,i] = sum(scenario_rates[ 
#'             (arrival_times < arrival_times_in_output[i] &
#'              maxima > tsunami_maxima_in_output[j])])
#'     }
#' }
#' # Check it worked
#' stopifnot(all(abs(back_computed_sol - exrate_matrix) <= 1.0e-06*exrate_matrix))
#' 
exceedance_rate_given_maxima_and_arrival_time<-function(
    scenario_maxima,
    scenario_arrival_times,
    scenario_rates,
    tsunami_maxima_in_output,
    arrival_times_in_output){

    # Make sure the following vars will be interpreted as vectors.
    stopifnot(length(scenario_maxima) > 1)
    stopifnot(length(scenario_maxima) == length(scenario_arrival_times))
    stopifnot(length(scenario_maxima) == length(scenario_rates))
    stopifnot(length(tsunami_maxima_in_output) > 1)
    stopifnot(length(arrival_times_in_output) > 1)

    # The following are not vectors
    ns = as.integer(length(scenario_maxima))
    narrival = as.integer(length(arrival_times_in_output))
    nmax = as.integer(length(tsunami_maxima_in_output))

    output_matrix = matrix(0.0, nrow=nmax, ncol=narrival)

    .Call("exceedance_rate_given_maxima_and_arrival_time_c",
        ns, scenario_maxima, scenario_arrival_times, scenario_rates,
        nmax, tsunami_maxima_in_output, narrival, arrival_times_in_output,
        output_matrix)

    return(output_matrix)

}
