#
# Gutenberg Richter functions
#



#' Cumulative distribution function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR
#' 
#' @param q vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
#' @examples
#' #
#' # Suppose we have some data, and fit the GR above Mw = mw_min_eg. What happens, 
#' # in cases where the data is raw, rounded, and/or perturbed randomly
#' #
#' 
#' negloglik_GR<-function(b, x, mw_min){
#'     -sum(log(dGR(x, b, mw_min)))
#' }
#' 
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
#' # SUMMARY: If the data is rounded to 1dp, then offset mw_min by -0.05, irrespective
#' # of whether the data is also measured with random error.
#' #

pGR<-function(q, b, mw_min){
    
    output = 1-10^(b*(mw_min - q))
    if(any(q < mw_min)){
        output[q<mw_min] = 0
    }
    return(output)
}

#' Inverse cumulative distribution function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR
#' 
#' @param p vector of probabilities
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
#' @examples
#' #
#' # Suppose we have some data, and fit the GR above Mw = mw_min_eg. What happens, 
#' # in cases where the data is raw, rounded, and/or perturbed randomly
#' #
#' 
#' negloglik_GR<-function(b, x, mw_min){
#'     -sum(log(dGR(x, b, mw_min)))
#' }
#' 
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
#' # SUMMARY: If the data is rounded to 1dp, then offset mw_min by -0.05, irrespective
#' # of whether the data is also measured with random error.
#' #

qGR<-function(p, b, mw_min){
    if(any(p < 0 | p > 1)) stop('Invalid p value')

    output = mw_min - log10(1-p)/b
    return(output)
}

#' Probability density function for Gutenberg Richter distribution
#'
#' Note corresponding functions pGR, qGR, dGR, rGR
#' 
#' @param x vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @return vector with probability densities at x
#' @export
#' @examples
#' #
#' # Suppose we have some data, and fit the GR above Mw = mw_min_eg. What happens, 
#' # in cases where the data is raw, rounded, and/or perturbed randomly
#' #
#' 
#' negloglik_GR<-function(b, x, mw_min){
#'     -sum(log(dGR(x, b, mw_min)))
#' }
#' 
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
#' # SUMMARY: If the data is rounded to 1dp, then offset mw_min by -0.05, irrespective
#' # of whether the data is also measured with random error.
#' #
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
#'
#' @examples
#' #
#' # Suppose we have some data, and fit the GR above Mw = mw_min_eg. What happens, 
#' # in cases where the data is raw, rounded, and/or perturbed randomly
#' #
#' 
#' negloglik_GR<-function(b, x, mw_min){
#'     -sum(log(dGR(x, b, mw_min)))
#' }
#' 
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
#' # SUMMARY: If the data is rounded to 1dp, then offset mw_min by -0.05, irrespective
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
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR
#' 
#' @param q vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with probability(mw < q | mw >= mw_min)
#' @export
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
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR
#' 
#' @param p vector of probabilities
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with Mw values having non-exceedance probabilities p
#' @export
qtruncGR<-function(p, b, mw_min, mw_max=Inf){

    if(any(p < 0 | p > 1)) stop('Invalid p value')

    p = p * pGR(mw_max, b, mw_min)
    qGR(p, b, mw_min)
}

#' Probability density function for truncated Gutenberg Richter distribution
#'
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR
#' 
#' @param x vector of quantiles
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with Mw values having non-exceedance probabilities p
#' @export
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
#' Note corresponding functions ptruncGR, qtruncGR, dtruncGR, rtruncGR
#' 
#' @param n integer giving desired number of samples
#' @param b Gutenberg-Richter b value
#' @param mw_min Minimum mw
#' @param mw_max Maximum mw
#' @return vector with random Mw values
#' @export
rtruncGR<-function(n, b, mw_min, mw_max=Inf){

    qtruncGR(runif(n), b=b, mw_min=mw_min, mw_max=mw_max)

}

