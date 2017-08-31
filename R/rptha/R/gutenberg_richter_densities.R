#
# Gutenberg Richter functions
#



#' Cumulative distribution function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR. See ?rGR for an extended
#' example of fitting with maximum likelihood.
#' 
#' @param q vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
#' @examples
#' # Compute exceedance probability for some random Mw values (i.e. conditional
#' # on the fact that an event occurred)
#' random_Mw = rGR(10, b=1, mw_min=5)
#' exceedance_prob = 1 - pGR(random_Mw, b=1, mw_min=5)
pGR<-function(q, b, mw_min){
    
    output = 1-10^(b*(mw_min - q))
    if(any(q < mw_min)){
        output[q<mw_min] = 0
    }
    return(output)
}

#' Inverse cumulative distribution function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR. See ?rGR for an extended example
#' of fitting with maximum likelihood.
#' 
#' @param p vector of probabilities
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
#' @examples
#' # Compute the 90th percentile Mw
#' Mw_90 = qGR(0.9, b=0.8, mw_min=6.0)
qGR<-function(p, b, mw_min){
    if(any(p < 0 | p > 1)) stop('Invalid p value')

    output = mw_min - log10(1-p)/b
    return(output)
}

#' Probability density function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR. See ?rGR for an extended
#' example of fitting with maximum likelihood.
#' 
#' @param x vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability densities at x
#' @export
#' @examples
#' # Compute the pdf of some random Mw values
#' random_Mw = rGR(10, b=0.8, mw_min = 6.0)
#' density_vals = dGR(random_Mw, b=0.8, mw_min=6.0)
dGR<-function(x, b, mw_min){

    output = b*log(10)*10^(b*(mw_min - x))

    if(any(x < mw_min)){
        output[x < mw_min] = 0    
    }
    return(output)
}

#' Random samples from Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR
#' 
#' @param n integer number of samples
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with random Mw values
#' @export
#' @examples
#' # Basic usage
#' random_Mw_samples = rGR(10, b=1, mw_min=5)
#'
#' #
#' # EXAMPLES OF FITTING WITH MAXIMUM LIKELIHOOD
#' #
#' # Suppose we have some data, and fit the GR above Mw = mw_min_eg. What happens, 
#' # in cases where the data is raw, rounded, and/or perturbed randomly
#' #
#' 
#' # Negative log-likelihood function 
#' negloglik_GR<-function(b, x, mw_min){
#'     -sum(log(dGR(x, b, mw_min)))
#' }
#' 
#' # Negative log-likelihood function for data 'rounded' to the nearest 0.1 
#' negloglik_GR_binned<-function(b, x, mw_min, delta=0.1){
#'     -sum(log(pGR(x+delta/2, b, mw_min) - pGR(pmax(x-delta/2, mw_min), b, mw_min)))
#' }
#' 
#' # Reproducible randomness
#' set.seed(1234)
#' 
#' 
#' # Set mw_min for the example (> 5)
#' mw_min_eg = 6.0
#' # Set 'true' b value
#' true_b = 1 
#' 
#' # Generate random Mw events above mw_min << mw_min_eg
#' random_Mw = rGR(1e+06, b=true_b, mw_min=mw_min_eg - 1.0)
#' # Round the data to 1 dp, i.e. nearest 0.1
#' delta = 0.1 # Consistent with use of round(x, 1)
#' rounded_random_Mw = round(random_Mw, 1)
#' # Jitter the data
#' jittered_random_Mw = random_Mw + rnorm(length(random_Mw), 0, delta)
#' # Jittered and rounded
#' jitterRound_random_Mw = round(jittered_random_Mw, 1)
#' 
#' k1 = which(random_Mw >= mw_min_eg)
#' k2 = which(rounded_random_Mw >= mw_min_eg)
#' k3 = which(jittered_random_Mw >= mw_min_eg)
#' k4 = which(jitterRound_random_Mw >= mw_min_eg)
#' 
#' #
#' # UNROUNDED DATA
#' #
#' 
#' # Unrounded data, naive approach -- PRETTY GOOD
#' fit1 = optimize(negloglik_GR, c(0.2, 2), x=random_Mw[k1], mw_min=mw_min_eg)
#' stopifnot(abs(fit1$minimum - true_b) < 0.01)
#' # Unrounded data -- BIASED if we subtract 0.05 from mw_min -- as there has been
#' # no rounding
#' fit1b = optimize(negloglik_GR, c(0.2, 2), x=random_Mw[k1], mw_min=mw_min_eg - delta/2)
#' stopifnot(abs(fit1b$minimum - true_b) > 0.05)
#' 
#' #
#' # ROUNDED DATA
#' # 
#' 
#' # Rounded, naive approach -- BIASED
#' fit1_rounding = optimize(negloglik_GR, c(0.2, 2.0), x=rounded_random_Mw[k2], mw_min=mw_min_eg)
#' stopifnot(abs(fit1_rounding$minimum - true_b) > 0.05)
#' # Rounded, with mw_min shifted -- OK
#' fit1_roundingb = optimize(negloglik_GR, c(0.2, 2.0), x=rounded_random_Mw[k2], mw_min=mw_min_eg-delta/2)
#' stopifnot(abs(fit1_roundingb$minimum - true_b) < 0.02)
#' 
#' # Integrated likelihood type approach, without mw_min adjustment -- BIASED
#' fitR = optimize(negloglik_GR_binned, c(0.2, 2.0), x=rounded_random_Mw[k2], mw_min=mw_min_eg, delta=delta)
#' stopifnot(abs(fitR$minimum - true_b) > 0.05)
#' # Integrated likelihood type approach, with mw_min adjustment -- GOOD
#' fitRb = optimize(negloglik_GR_binned, c(0.2, 2.0), x=rounded_random_Mw[k2], mw_min=mw_min_eg-delta/2, delta=delta)
#' stopifnot(abs(fitRb$minimum - true_b) < 0.01)
#' 
#' #
#' # JITTERED DATA, NO ROUNDING
#' #
#' 
#' # Naive approach -- OK
#' fitJ = optimize(negloglik_GR, c(0.2, 2.0), x=jittered_random_Mw[k3], mw_min=mw_min_eg)
#' stopifnot(abs(fitJ$minimum - true_b) < 0.01)
#' # Naive approach with adjusted Mw-min -- BIASED
#' fitJb = optimize(negloglik_GR, c(0.2, 2.0), x=jittered_random_Mw[k3], mw_min=mw_min_eg-delta/2)
#' stopifnot(abs(fitJb$minimum - true_b) > 0.05)
#' 
#' # Integrated likelihood type approach, without mw_min adjustement -- OK
#' fitJc = optimize(negloglik_GR_binned, c(0.2, 2.0), x=jittered_random_Mw[k3], mw_min=mw_min_eg, delta=delta)
#' stopifnot(abs(fitJc$minimum - true_b) < 0.01)
#' # Integrated likelihood type approach, with mw_min adjustement -- BIASED
#' fitJd = optimize(negloglik_GR_binned, c(0.2, 2.0), x=jittered_random_Mw[k3], mw_min=mw_min_eg-delta/2, delta=delta)
#' stopifnot(abs(fitJd$minimum - true_b) > 0.05)
#' 
#' 
#' #
#' # JITTERED AND ROUNDED DATA
#' #
#' 
#' # Naive approach -- BIASED
#' fitJR = optimize(negloglik_GR, c(0.2, 2.0), x=jitterRound_random_Mw[k4], mw_min=mw_min_eg)
#' stopifnot(abs(fitJR$minimum - true_b) > 0.05)
#' # Naive approach with adjusted Mw-min  -- OK
#' fitJRb = optimize(negloglik_GR, c(0.2, 2.0), x=jitterRound_random_Mw[k4], mw_min=mw_min_eg-delta/2)
#' stopifnot(abs(fitJRb$minimum - true_b) < 0.02)
#' 
#' # Integrated likelihood type approach, without mw_min adjustement -- BIASED
#' fitJRc = optimize(negloglik_GR_binned, c(0.2, 2.0), x=jitterRound_random_Mw[k4], mw_min=mw_min_eg, delta=delta)
#' stopifnot(abs(fitJRc$minimum - true_b) > 0.05)
#' # Integrated likelihood type approach, with mw_min adjustement -- PRETTY GOOD
#' fitJRd = optimize(negloglik_GR_binned, c(0.2, 2.0), x=jitterRound_random_Mw[k4], mw_min=mw_min_eg-delta/2, delta=delta)
#' stopifnot(abs(fitJRd$minimum - true_b) < 0.01)
#' 
#' 
#' #
#' # SUMMARY: If the data is rounded to 1dp, then offset mw_min by -delta/2, irrespective
#' # of whether the data is also measured with random error.
#' #
rGR<-function(n, b, mw_min){

    qGR(runif(n), b=b, mw_min=mw_min)
}


#
# Truncated Gutenberg Richter functions
#


#' Cumulative distribution function for truncated Gutenberg Richter distribution
#'
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR. See '?rtruncGR'
#' for an extended example using the functions to estimate 'b'.
#' 
#' @param q vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
#' @examples
#' # Get the exceedance probability (given that an event occurred) for 10
#' # random Mw values
#' random_Mw = rtruncGR(10, b=0.7, mw_min=6.0, mw_max = 9.5)
#' exceedance_prob = 1 - ptruncGR(random_Mw, b=0.7, mw_min=6.0, mw_max=9.5)
#' 
ptruncGR<-function(q, b, mw_min, mw_max=Inf){
    
    output = (1 - 10^(b*(mw_min - q)))/pGR(mw_max, b, mw_min)

    if(any(q < mw_min)){
        output[q < mw_min] = 0
    }
    if(any(q > mw_max)){
        output[q > mw_max] = 1
    }

    return(output)
}

#' Inverse cumulative distribution function for truncated Gutenberg Richter distribution
#'
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR. See ?rtruncGR
#' for an example of using the functions for b estimation.
#' 
#' @param p vector of probabilities
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with Mw values having non-exceedance probabilities p
#' @export
#' @examples
#' # 90th percentile of truncated GR 
#' Mw_90 = qtruncGR(0.9, b=1.2, mw_min = 5.5, mw_max = 8.0)
#'
qtruncGR<-function(p, b, mw_min, mw_max=Inf){

    if(any(p < 0 | p > 1)) stop('Invalid p value')

    p = p * pGR(mw_max, b, mw_min)
    qGR(p, b, mw_min)
}

#' Probability density function for truncated Gutenberg Richter distribution
#'
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR. See ?rtruncGR
#' for an example of using the functions for b estimation
#' 
#' @param x vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with Mw values having non-exceedance probabilities p
#' @export
#' @examples
#' # Get the probability density for some random Mw values
#' random_Mw = rtruncGR(10, b=0.7, mw_min=6.0, mw_max = 9.5)
#' random_density = dtruncGR(random_Mw, b=0.7, mw_min=6.0, mw_max=9.5)
#' 
dtruncGR<-function(x, b, mw_min, mw_max=Inf){
    
    output = dGR(x, b, mw_min)/pGR(mw_max, b, mw_min)

    if(any(x < mw_min)){
        output[x<mw_min] = 0
    }

    if(any(x > mw_max)){
        output[x>mw_max] = 0 
    }

    return(output)
}

#' Random samples from truncated Gutenberg Richter distribution
#'
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR.
#' 
#' @param n integer giving desired number of samples
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with random Mw values
#' @export
#' @examples
#' #
#' # Make 300 random Mw samples
#' #
#' random_Mws = rtruncGR(300, b=0.7, mw_min=5.0, mw_max=Inf)
#' #
#' # Estimate 'b' with maximum likelihood
#' #
#' negloglik_truncGR<-function(b, data, mw_min, mw_max){
#'     -sum(log(dtruncGR(data, b=b, mw_min=mw_min, mw_max=mw_max)))
#' }
#' fit_inf_mw_max = optimize(negloglik_truncGR, interval=c(0, 2), data=random_Mws, 
#'     mw_min=5.0, mw_max=Inf, tol=1e-12)
#' #
#' # Because mw_max = Inf, this should give Aki's maximum likelihood estimator
#' # = 1/(ln(10) * (mean(data) - mw_min))
#' aki_estimator = 1/(log(10) * ( mean(random_Mws) - 5.0))
#' stopifnot(abs(fit_inf_mw_max$minimum - aki_estimator) < 1.0e-06)
#' 
#' #
#' # If mw_max < Inf, then the maximum likelihood b-value will not be the same
#' # as Aki's estimator, since we are using a truncated GR distribution
#' #
#' random_Mws = rtruncGR(300, b=0.7, mw_min=5.0, mw_max=8.8)
#' fit_finite_mw_max = optimize(negloglik_truncGR, interval=c(0, 2), data=random_Mws, 
#'     mw_min=5.0, mw_max=8.8)
#' aki_estimator = 1/(log(10) * ( mean(random_Mws) - 5.0))
#' 
rtruncGR<-function(n, b, mw_min, mw_max=Inf){

    qtruncGR(runif(n), b=b, mw_min=mw_min, mw_max=mw_max)

}


#
# Multi-catalogue maximum likelihood fitting
#

#' Fit a truncated GR distribution to multiple earthquake catalogues with different
#' completeness magnitudes.
#'
#' This routine estimates the rate of earthquakes and the Gutenberg Richter b
#' value from a set of catalogues that are assumed to relate to a single
#' location, but have different completeness magnitudes. For example, at a given
#' site we might have 3 different 'catalogues': One with magnitude data that is
#' complete for the last 40 years above Mw=6; older data that is complete for
#' 160 years before that above Mw 7.5; and even earlier observations which are
#' complete for 100 years above Mw 9.0 (perhaps with no earthquakes in the latter case). 
#' These can all be combined to estimate the rate of earthquakes (e.g. above Mw 6),
#' and the 'b' parameter, assuming that the event timings are Possionian, and the
#' event magnitudes follow a truncated GR distribution with known mw_max.
#'
#' @param catalogue_lists list-of-lists, containing a single list for each catalogue. 
#'  For each catalogue, the list must have components 'mw_min', 'duration' (giving
#'  the minimum magnitude and observation duration), and a vector 'mw' giving the
#'  mw_values (may be NA, in which case it is assumed). See the example
#' @param start_par vector of length 2 giving starting values for the rate of events
#'  above the minimum mw_min, and the gutenberg-richter b value  
#' @param mw_max maximum magnitude for the source zone. Can be Inf.
#' @param return_environment logical. If TRUE, return the function environment. 
#'   If FALSE, just return the fit (i.e. output from \code{optim})
#' @param ... further arguments to \code{optim}
#' @return See return_environment above. By default, returns an object where fit$par
#'   contains c('estimated rate of events above min(mw_min)', 'estimated b'). The
#'   units of the rate is ('events-per-unit-of_catalogue-duration) (e.g.
#'   typically events per year). See the examples for a technique to extract
#'   parameter uncertainties from the fit
#' @export
#' @examples
#' #
#' # A scenario with 3 catalogue completnesses. Assume the data is complete for 
#' # the last 40 years above Mw=6; complete for 160 years before that above Mw 7.5,
#' # and complete for 100 years before that above Mw 9.0 (perhaps with no earthquakes in 
#' # any of the cases)
#' #
#' mw_max = 9.5 # Could be Inf to fit pure Gutenberg Richter
#' mw_mins = c(9.0, 7.5, 6.0) # Compleness magnitudes
#' durations = c(100, 160, 40) # Observation times
#' true_rate_above_6 = 0.5 # We will estimate this value from simulated data
#' true_b = 0.8 # We will estimate this value from simulated data
#' 
#' 
#' # Repeat the fit 'Nreps' times on different simulated data, to get an idea of 
#' # the accuracy of the method.
#' Nreps = 500
#' # Vectors to store results
#' store_b = rep(NA, Nreps)
#' store_rate = rep(NA, Nreps)
#' store_b_sd = rep(NA, Nreps)
#' store_rate_sd = rep(NA, Nreps)
#' 
#' for(jj in 1:Nreps){
#' 
#'     # Simulate random catalogues
#'     true_rates_above_mwmin = true_rate_above_6 * (1-
#'         ptruncGR(mw_mins, b=true_b, mw_min=min(mw_mins), mw_max=mw_max))
#'     # Number of events
#'     nevents = rpois(3, lambda=true_rates_above_mwmin * durations)
#'     #
#'     # MAKE catalogue_lists HERE
#'     #
#'     catalogue_lists = list()
#'     for(i in 1:3){
#'         # Do the i'th catalogue
#'         catalogue_lists[[i]] = list()
#' 
#'         # Insert duration value
#'         catalogue_lists[[i]]$duration = durations[i] 
#'         # Insert mw_min value
#'         catalogue_lists[[i]]$mw_min = mw_mins[i] 
#' 
#'         # Insert mw values -- special notation if nevents = 0
#'         if(nevents[i] > 0){
#'             # Generate random mw values
#'             catalogue_lists[[i]]$mw = rtruncGR(nevents[i], b=true_b, 
#'                 mw_min=mw_mins[i], mw_max=mw_max)
#'         }else{
#'             # Denote 'no events observed' -- still useful information, because it
#'             # could reduce the chance that 'very high' rates are correct
#'             catalogue_lists[[i]]$mw = NA
#'         }
#'     }
#'     # Now the catalogue_lists data is created
#' 
#'     # Estimate the 'true_rate_above_6' and 'true_b', with starting guesses of
#'     # 0.5 and 1 respectively
#'     myfit = fit_truncGR_multiple_catalogues(catalogue_lists=catalogue_lists, 
#'         start_par=c(0.5, 1), mw_max=9.5, hessian=TRUE)
#' 
#' 
#'     # Store the outputs, so we can study the performance of the estimators
#'     store_b[jj] = myfit$par[2]
#'     store_rate[jj] = myfit$par[1]
#'     sds = sqrt(diag(solve(myfit$hessian)))
#'     store_b_sd[jj] = sds[2]
#'     store_rate_sd[jj] = sds[1]
#' }
#' #
#' # Report on results of all Nreps fits
#' #
#' 
#' # store_b should be close to true_b, at least on average
#' print(summary(store_b))
#' # store_rate should be close to true_rate_above_6, at least on average
#' print(summary(store_rate))
#' # Compute the fraction of cases which had the true value inside
#' # 2-standard-deviation confidence intervals (should be approximately 0.95)
#' print(mean(abs(store_b - true_b) < 2*store_b_sd))
#' print(mean(abs(store_rate - true_rate_above_6) < 2*store_rate_sd))
#'
fit_truncGR_multiple_catalogues<-function(catalogue_lists, start_par, 
    mw_max=Inf, return_environment = FALSE, ...){

    # Parse catalogues
    ncat = length(catalogue_lists)
    mw_mins = unlist(lapply(catalogue_lists, f<-function(x) x$mw_min))
    durations = unlist(lapply(catalogue_lists, f<-function(x) x$duration))

    if(any(durations <= 0)) stop('Error: Duration must be positive')

    # Record catalogues with no observations 
    no_obs = unlist(lapply(catalogue_lists, 
        f<-function(x) (length(x$mw) == 1 & is.na(x$mw[1]))))

    # Check all mw are consistent with mw_min
    for(i in 1:ncat){
        if(!no_obs[i]){
            if(any(catalogue_lists[[i]]$mw < mw_mins[i])){
                stop(paste0('Catalogue ', i, ' has mw < mw_min'))
            }
            if(any(catalogue_lists[[i]]$mw >= mw_max)){
                stop(paste0('Catalogue ', i, ' has mw >= mw_max'))
            }
        }
    }

    min_mw_min = min(mw_mins)
    mw_max = mw_max

    # Function to optimize
    negloglik_fun<-function(par){

        rate = par[1]
        b = par[2]

        nll = 0
        for(i in 1:ncat){

            local_rate = rate * 
                (1-ptruncGR(mw_mins[i], b=b, mw_min=min_mw_min, mw_max=mw_max))

            if(no_obs[i]){
                nll = nll - dpois(0, lambda=local_rate * durations[i], log=TRUE)
            }else{
                local_mw = catalogue_lists[[i]]$mw
                nll = nll - 
                    dpois(length(local_mw), lambda=local_rate*durations[i], log=TRUE) -
                    sum(log(dtruncGR(local_mw, b=b, mw_min = mw_mins[i], mw_max=mw_max)))
            }
        }

        return(nll)
    }

    fit = optim(start_par, negloglik_fun, ...)

    if(return_environment){
        return(environment())
    }else{
        return(fit)
    }
}

