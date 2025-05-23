% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauge_statistics.R
\name{gauge_range_filtered}
\alias{gauge_range_filtered}
\title{Gauge time-series filtered range}
\usage{
gauge_range_filtered(
  data_t,
  data_s,
  filter_freq = 1/(2 * 60),
  interp_dt = NULL,
  detailed = FALSE
)
}
\arguments{
\item{data_t}{vector of numeric times (usually in seconds)}

\item{data_s}{vector of numeric stages corresponding to data_t (usually in metres)}

\item{filter_freq}{We remove frequencies above filter_freq from the series
before computing the range.}

\item{interp_dt}{numeric time-step (same units as data_t). We interpolate
the data to a fixed time-step given by interp_dt, before performing the
discrete fourier transform. If NULL, the minimum spacing in data_t is used.}

\item{detailed}{logical. If FALSE, only the filtered series range is returned. If TRUE,
then we return a list, containing the filtered series range, but also A) the
time of the max/min of the filtered series, and B) the frequency of the
dominant component of the series in the times between the min and max (or max
and min). The latter is computed by taking the filtered series between the
two extrema, concatenating it with a reversed version of itself, computing
the fft of the latter, and finding the frequency with dominant spectral
peak. The idea is this can be a robust way to get an idea of the peak tsunami
frequency around the time of the largest waves.}
}
\value{
vector giving the min and max of the filtered stage
}
\description{
Compute the range of a gauge stage time-series after spectral filtering of high frequencies.
This is practically useful for avoiding seismic wave signals in some DART buoys.
The results should always be checked visually though -- since it cannot automatically
treat all 'bad-behaviour' in the observations.
}
\examples{
# Make a wave train with 2 spectral components, having frequencies 1 and 1/50
t = seq(1,1000,by=0.1)
stage = 0.2*sin(2*pi*t) + 1.0*sin(2*pi*t/50)

range(stage) # Should be about +-1.2
# Get the stage range after filtering out frequencies > 1/10. By construction it
# should be about +- 1.0
longperiod_range = gauge_range_filtered(t, stage, filter_freq = 1/10)
stopifnot(all(abs(longperiod_range - c(-1,1)) < 1.0e-02))
# Check that if we filter frequencies > 1, then it still works
full_range = gauge_range_filtered(t, stage, filter_freq = 1/0.5)
stopifnot(all(abs(full_range - range(stage)) < 1.0e-02))

# Check on spectral details
longperiod_range2 = gauge_range_filtered(t, stage, filter_freq = 1/10, detailed=TRUE)
# Should show a peak frequency close to 0.02 ( = 1/50)
stopifnot(abs(longperiod_range2$mintomax_peak_frequency - 1/50) < 0.01/50)

# Try it with the full data -- we should still pick out the 1/50 frequency as dominant
full_range2 = gauge_range_filtered(t, stage, filter_freq = 1/0.5, detailed = TRUE)
stopifnot(abs(longperiod_range2$mintomax_peak_frequency - 1/50) < 0.01/50)
}
