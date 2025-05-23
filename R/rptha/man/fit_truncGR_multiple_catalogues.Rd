% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutenberg_richter_densities.R
\name{fit_truncGR_multiple_catalogues}
\alias{fit_truncGR_multiple_catalogues}
\title{Fit a truncated GR distribution to multiple earthquake catalogues with different
completeness magnitudes.}
\usage{
fit_truncGR_multiple_catalogues(
  catalogue_lists,
  start_par,
  mw_max = Inf,
  return_environment = FALSE,
  ...
)
}
\arguments{
\item{catalogue_lists}{list-of-lists, containing a single list for each catalogue. 
For each catalogue, the list must have components 'mw_min', 'duration' (giving
the minimum magnitude and observation duration), and a vector 'mw' giving the
mw_values (may be NA, in which case it is assumed). See the example}

\item{start_par}{vector of length 2 giving starting values for the rate of events
above the minimum mw_min, and the gutenberg-richter b value}

\item{mw_max}{maximum magnitude for the source zone. Can be Inf.}

\item{return_environment}{logical. If TRUE, return the function environment. 
If FALSE, just return the fit (i.e. output from \code{optim})}

\item{...}{further arguments to \code{optim}}
}
\value{
See return_environment above. By default, returns an object where fit$par
  contains c('estimated rate of events above min(mw_min)', 'estimated b'). The
  units of the rate is ('events-per-unit-of_catalogue-duration) (e.g.
  typically events per year). See the examples for a technique to extract
  parameter uncertainties from the fit
}
\description{
This routine estimates the rate of earthquakes and the Gutenberg Richter b
value from a set of catalogues that are assumed to relate to a single
location, but have different completeness magnitudes. For example, at a given
site we might have 3 different 'catalogues': One with magnitude data that is
complete for the last 40 years above Mw=6; older data that is complete for
160 years before that above Mw 7.5; and even earlier observations which are
complete for 100 years above Mw 9.0 (perhaps with no earthquakes in the latter case). 
These can all be combined to estimate the rate of earthquakes (e.g. above Mw 6),
and the 'b' parameter, assuming that the event timings are Possionian, and the
event magnitudes follow a truncated GR distribution with known mw_max.
}
\examples{
#
# A scenario with 3 catalogue completnesses. Assume the data is complete for 
# the last 40 years above Mw=6; complete for 160 years before that above Mw 7.5,
# and complete for 100 years before that above Mw 9.0 (perhaps with no earthquakes in 
# any of the cases)
#
mw_max = Inf # Could be finite, or use Inf to fit pure Gutenberg Richter
mw_mins = c(9.0, 7.5, 6.0) # Compleness magnitudes
durations = c(100, 160, 40) # Observation times

# We will estimate the following 'true' values from simulated data to
# evaluate the statistical method
true_rate_above_min_mwmin = 0.5 
true_b = 0.8 


# Repeat the fit 'Nreps' times on different simulated data, to get an idea of 
# the accuracy of the method.
Nreps = 500
# Vectors to store results
store_b = rep(NA, Nreps)
store_rate = rep(NA, Nreps)
store_b_sd = rep(NA, Nreps)
store_rate_sd = rep(NA, Nreps)

# Use these to store a fit to 'only' the most recent catalogue, for comparison
store_b_r = rep(NA, Nreps)
store_rate_r = rep(NA, Nreps)
store_b_sd_r = rep(NA, Nreps)
store_rate_sd_r = rep(NA, Nreps)

# Use this to store the analytical Kijko & Smit estimator (2012)
store_ks_b = rep(NA, Nreps)


for(jj in 1:Nreps){

    # Simulate random catalogues
    true_rates_above_mwmin = true_rate_above_min_mwmin * (1-
        ptruncGR(mw_mins, b=true_b, mw_min=min(mw_mins), mw_max=mw_max))
    # Number of events
    nevents = rpois(length(mw_mins), lambda=true_rates_above_mwmin * durations)
    #
    # MAKE catalogue_lists HERE
    #
    catalogue_lists = list()
    for(i in 1:length(mw_mins)){
        # Do the i'th catalogue
        catalogue_lists[[i]] = list()

        # Insert duration value
        catalogue_lists[[i]]$duration = durations[i] 
        # Insert mw_min value
        catalogue_lists[[i]]$mw_min = mw_mins[i] 

        # Insert mw values -- special notation if nevents = 0
        if(nevents[i] > 0){
            # Generate random mw values
            catalogue_lists[[i]]$mw = rtruncGR(nevents[i], b=true_b, 
                mw_min=mw_mins[i], mw_max=mw_max)
        }else{
            # Denote 'no events observed' -- still useful information, because it
            # could reduce the chance that 'very high' rates are correct
            catalogue_lists[[i]]$mw = NA
        }
    }
    # Now the catalogue_lists data is created

    # Estimate the 'true_rate_above_min_mwmin' and 'true_b', with starting guesses of
    # 0.5 and 1 respectively
    myfit = fit_truncGR_multiple_catalogues(catalogue_lists=catalogue_lists, 
        start_par=c(0.5, 1), mw_max=mw_max, hessian=TRUE, control=list(reltol=1e-12))

    # Compare with fit to only recent data 
    basic_fit =  fit_truncGR_multiple_catalogues(catalogue_lists=catalogue_lists[3],
        start_par=c(0.5, 1), mw_max=mw_max, hessian=TRUE, control=list(reltol=1e-12))


    # Store the outputs, so we can study the performance of the estimators
    store_b[jj] = myfit$par[2]
    store_rate[jj] = myfit$par[1]
    sds = sqrt(diag(solve(myfit$hessian)))
    store_b_sd[jj] = sds[2]
    store_rate_sd[jj] = sds[1]


    # Store the basic-fit outputs, for comparison
    store_b_r[jj] = basic_fit$par[2]
    store_rate_r[jj] = basic_fit$par[1]
    sds = sqrt(diag(solve(basic_fit$hessian)))
    store_b_sd_r[jj] = sds[2]
    store_rate_sd_r[jj] = sds[1]

    # Compare with 'b' estimator from Kijko and Smit (2012)
    # Note this requires mw_max = Infinity, whereas the other 2 approaches
    # don't need that.
    ns = unlist(lapply(catalogue_lists, function(x) length(x$mw)))
    bs = unlist(lapply(catalogue_lists, function(x) 1/(log(10) * ( mean(x$mw) - x$mw_min)) ))
    # If a catalogue has no data, remove it
    not_missing = unlist(lapply(catalogue_lists, function(x) !all(is.na(x$mw))))
    bs = bs[which(not_missing)]
    ns = ns[which(not_missing)]
    store_ks_b[jj] = 1/( sum(ns/sum(ns) * 1/bs))

}
#
# Report on results of all Nreps fits
#

# store_b should be close to true_b, at least on average
print(summary(store_b))
# store_rate should be close to true_rate_above_min_mwmin, at least on average
print(summary(store_rate))
# Compute the fraction of cases which had the true value inside
# 2-standard-deviation confidence intervals (should be approximately 0.95)
print(mean(abs(store_b - true_b) < 2*store_b_sd))
print(mean(abs(store_rate - true_rate_above_min_mwmin) < 2*store_rate_sd))
# The mean fitted value should generally be close to the true value -- the
# bounds here should very very rarely fail
stopifnot(abs(mean(store_b) - true_b) < 0.07)
stopifnot(abs(mean(store_rate) - true_rate_above_min_mwmin) < 0.05)

# Compare variance of fit using the current multi-catalogue method, the method
# of Kijko and Smit (2012), and the approach of only fitting the 'best'
# catalogue. 
print(summary(store_b - true_b)) # Current approach -- lowest variance
print(summary(store_ks_b - true_b)) # Kijko-Smit 2012 neat analytical formula
print(summary(store_b_r - true_b)) # Pure Aki applied to the best catalogue


}
