% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauge_statistics.R
\name{gauge_statistics_simple}
\alias{gauge_statistics_simple}
\title{Simple gauge summary statistics}
\usage{
gauge_statistics_simple(
  gauge_times,
  flow_data,
  stage_threshold_for_arrival_time,
  use_fortran = TRUE,
  ...
)
}
\arguments{
\item{gauge_times}{vector giving the times at which each flow-time-series
is measured}

\item{flow_data}{3D array with size [number_of_gauges, 
number_of_timesteps, number_of_flow_variables]. It is assumed that stage is
stored where the 'third index is equal to 1'.}

\item{stage_threshold_for_arrival_time}{The arrival time is defined as
the time at which the absolute stage exceeds this value.}

\item{use_fortran}{logical. If FALSE, uses an R implementation. The latter
is slow but useful for testing}

\item{...}{other arguments}
}
\value{
matrix with 5 columns and one row for each gauge. The matrix columns
are: 1. max stage; 2. zero crossing period; 3. stage range; 4. Time that abs(stage)
exceeds stage_threshold_for_arrival_time; 5. The initial stage
}
\description{
This function computes gauge-summary statistics at one or more gauges,
stored as a 3D array of time-series, with dimensions [number_gauges,
number_time-steps, number_flow_variables]. It is implemented in fortran for
speed. The main purpose is efficiently computing summary statistics for all
gauges in the PTHA (because not all the function operations can be vectorized,
looping over the gauges leads to too much overhead
in pure R)
}
\examples{
t = seq(0,1000)

# Store gauge time-series
#

flow_array = array(runif(9*length(t)), dim=c(3, length(t), 3))

# Set the stage values in [,,1]. The other values will not effect the computation
flow_array[1,,1] = sin(2*pi*(t-3)/30)
flow_array[2,,1] = sin(2*pi*(t-3)/10)
flow_array[3,,1] = sin(2*pi*t/3)

library(rptha)

output = gauge_statistics_simple(t, flow_array, stage_threshold_for_arrival_time = 7.0e-01)

#
# Regression test -- we expect this result, although the period
# computation has slight errors if use_fortran=FALSE (but within a few percent)
#
expected_output = structure(c(0.994521895368277, 0.951056516295184, 0.866025403784579,
                              30, 10, 3, 1.98904379073655, 1.90211303259037, 1.73205080756909, 
                              7, 0, 1, -0.587785252292473, -0.951056516295154, 0), .Dim = c(3L, 
                              5L))

outputB = gauge_statistics_simple(t, flow_array, stage_threshold_for_arrival_time = 7.0e-01, use_fortran=FALSE)
# Errors should be very small
stopifnot(all(abs(output-expected_output) <= 1.0e-06*abs(expected_output)))

# Use a weaker test for the R version for period, because it does not do interpolation
# in the zero-crossing-period computation, so errors end up larger (but < 0.2\%)
stopifnot(all(abs(outputB[,-2]-expected_output[,-2]) <= 1.0e-06*abs(expected_output[,-2])))
stopifnot(all(abs(outputB[,2]-expected_output[,2]) <= 2.0e-03*abs(expected_output[,2])))


}
