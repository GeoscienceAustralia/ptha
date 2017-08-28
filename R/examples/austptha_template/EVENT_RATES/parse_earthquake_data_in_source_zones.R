#
# Code containing experimental/random analysis of earthquakes in our source zones
# Not part of main scripts. (FIXME: Consider removing from version control)
#

# Get GCMT + utility functions
source('gcmt_subsetter.R', local=TRUE)

eq_events = vector(mode='list', length=length(unit_source_grid_poly))
for(i in 1:length(unit_source_grid_poly)){
    eq_events[[i]] = get_gcmt_events_in_poly(names(unit_source_grid_poly)[i])
}
names(eq_events) = names(unit_source_grid_poly)

#
# Get rate of events > threshold on source-zone
#

source_zone_exceedance_rate = sapply(eq_events, f<-function(x) length(x[,1])/cmt_duration_years)

# time between events [ignoring first/last points issues where the time is
# uncertain]
source_zone_time_differences = lapply(eq_events, f<-function(x) diff(x$julianDay1900)/days_in_year )

#
# EXPLORATORY ANALYSES BELOW HERE
#

# Gutenberg Richter cumulative distribution
GR<-function(x, b, mw_min) 1-10^(b*(mw_min - x))

# Gutenberg Richter pdf
gR<-function(x, b, mw_min) b*log(10)*10^(b*(mw_min - x))

negloglik_GR<-function(b, x, mw_min){
    -sum(log(gR(x, b, mw_min)))
}

## Find events within 0.2 degrees of a source
print_b_estimates<-function(poly){

    xx_c = lonlat_in_poly(gcmt[,c('cent_lon', 'cent_lat')], poly, buffer_width=0.2)
    xx_h = lonlat_in_poly(gcmt[,c('hypo_lon', 'hypo_lat')], poly, buffer_width=0.2)

    xx = (xx_c | xx_h)

    # Fit 'b' to the data for different possible Mw-min's
    for(mm in seq(5., 6.5, by=0.1)){

        kk = which(xx & gcmt$Mw > mm & gcmt$depth <= depth_threshold)

        out = optimize(negloglik_GR, lower=0.5, upper=1.5, x=gcmt$Mw[kk], mw_min=mm)
        out_exact = 1/(log(10)*mean(gcmt$Mw[kk] - mm))

        test_b_CI = seq(0.2, 2.5, by=0.01)

        out_multiple = sapply(test_b_CI, f<-function(b) negloglik_GR(b, x=gcmt$Mw[kk], mw_min=mm) )
        out_fun = approxfun(test_b_CI, out_multiple)
        out_fun_offset<-function(x) (out_fun(x) - out$objective - qchisq(0.95, df=1)/2 )
        out_lower = uniroot(out_fun_offset, lower=min(test_b_CI), upper=out$minimum)
        out_upper = uniroot(out_fun_offset, lower=out$minimum, upper=max(test_b_CI))

        # mag, N, solution, solution, lower profile likelihood CI, upper profile likelihood CI, analytical sd
        print(c(mm, length(kk), out_exact, out$minimum, out_lower$root, out_upper$root, out$minimum/sqrt(length(kk))))

    }
}

# Normalise time differences by rate of events [e.g. Geist et al, 2011, and
# much earthquake literature]
rate_normalised_sztd = source_zone_time_differences
for(i in 1:length(rate_normalised_sztd)){
    rate_normalised_sztd[[i]] = rate_normalised_sztd[[i]]*source_zone_exceedance_rate[i]
}

# Try exponential and gamma fit to the entire 'rate-normalised'
# time-between-events dataset
library(fitdistrplus)
mm = unlist(rate_normalised_sztd)
# One zero -- must be removed for fitting -- but check it
mm[which(mm==0)] = min(mm[which(mm>0)])
# exponential (poisson process)
exp_fit_mle = fitdist(data=mm, distr='exp', method='mle')
# gamma
gamma_fit_mme = fitdist(data=mm, distr='gamma', method='mme')
gamma_fit_mle = fitdist(data=mm, distr='gamma', start = as.list(gamma_fit_mme$estimate))
gamma_fit_mge = fitdist(data=mm, distr='gamma', method='mge')
# weibull
weibull_fit_mle = fitdist(data=mm, distr='weibull', method='mle', start=list(shape=1, scale=1))

#
# Goodness of fit tests
#
gofstat(gamma_fit_mle)
gofstat(weibull_fit_mle)
gofstat(exp_fit_mle)

pdf('test.pdf', width=10, height=8)
par(mfrow=c(2,2))
qqcomp(exp_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Exponential model QQ plot')
qqcomp(gamma_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Gamma model QQ plot')
qqcomp(weibull_fit_mle, xlogscale=TRUE, ylogscale=TRUE, a.ppoints=0, main='Weibull model QQ plot')

par(mfrow=c(3,3))
mw_store = c()
time_diff_store = c()
for(i in 1:length(eq_events)){
    try(plot(eq_events[[i]]$Mw, c(rate_normalised_sztd[[i]], NA), main=names(eq_events)[i]) )
    mw_store = c(mw_store, eq_events[[i]]$Mw)
    time_diff_store = c(time_diff_store, rate_normalised_sztd[[i]], NA)
}

#
# Is there a relation between the normalised time between events, and the magnitude of the first event?
#
par(mfrow=c(1,1))
plot(range(mw_store, na.rm=TRUE), range(time_diff_store[time_diff_store>0], na.rm=TRUE), col=0, 
    xlab='Mw', ylab='Rate normalised time between events', log='y')
for(i in 1:length(eq_events)){
    try(points(eq_events[[i]]$Mw, c(rate_normalised_sztd[[i]], NA), pch=i, col=rainbow(9)[i], cex=2))
}
legend('topright', names(eq_events), pch=1:8, col=rainbow(9)[1:8], pt.cex=2)

# Some relationship ?
mw_time_diff_test = cor.test(mw_store, time_diff_store, method='s')


#
# Spatial plots
#

library(RFOC)

for(i in 1:length(eq_events)){

    source_zone_events_plot(names(unit_source_grid_poly)[i], eq_events[[i]])

}

dev.off()

#'
#' Alternative parameterisation of gamma distribution
#'
dgamma2<-function(x, alpha, mean_rate, log=FALSE){
    shape = alpha
    rate = mean_rate * alpha
    dgamma(x, shape=shape, rate=rate, log=log)
}

pgamma2<-function(q, alpha, mean_rate = 1, lower.tail = TRUE, log.p = FALSE){
    shape = alpha
    rate = mean_rate * alpha
    pgamma(q, shape=shape, rate=rate, lower.tail=lower.tail, log.p=log.p)
}

qgamma2<-function(p, alpha, mean_rate = 1, lower.tail = TRUE, log.p = FALSE){
    shape = alpha
    rate = mean_rate * alpha
    qgamma(q, shape=shape, rate=rate, lower.tail=lower.tail, log.p=log.p)
}

rgamma2<-function(n, alpha, mean_rate = 1){
    shape = alpha
    rate = mean_rate * alpha
    rgamma(n, shape=shape, rate=rate)
}

#'
#' Negative log-likelihood for a model where time-between-events have a gamma distribution,
#' with a global shape parameter, and a source-zone specific mean rate
#
nll_regional_gamma2<-function(par, data=eq_events, 
    obs_start = cmt_start_time_julianDay1900, 
    obs_end = cmt_end_time_julianDay1900, 
    dt_limit_days = 0,
    mean_rate_lower_limit = 0.001,
    return_vector=FALSE){

    # par is a vector with the (global) shape parameter, followed by a 'mean_rate' parameter
    # for each source-zone

    alpha = par[1]
    stopifnot(length(par)-1 == length(data))

    if(any(par[-1] < mean_rate_lower_limit)) return(Inf)

    nll_site = rep(NA, length=length(data))

    for(i in 1:length(data)){

        data_times = data[[i]]$julianDay1900
        mean_rate = par[i+1]

        if(length(data_times) == 0){
            # No events -- so we know the time-between-events is at least obs_end-obs_start
            nll_site[i] = -pgamma2(
                (obs_end - obs_start)/days_in_year, 
                alpha, mean_rate, lower.tail=FALSE, log.p=TRUE)
        }

        if(length(data_times) == 1){
            # We only have bounds on the time between the observed event and the ones before/after
            nll_site[i] = 
                -pgamma2((data_times - obs_start)/days_in_year, 
                    alpha, mean_rate, lower.tail=FALSE, log.p=TRUE) -
                pgamma((obs_end - data_times)/days_in_year, 
                    alpha, mean_rate, lower.tail=FALSE, log.p=TRUE)
        }

        if(length(data_times) > 1){
            # Censored start and end time spacings, standard treatment for data inside
            nll_site[i] = -sum(
                dgamma2(pmax(dt_limit_days, diff(data_times))/days_in_year, 
                    alpha, mean_rate, log=TRUE)) -
                pgamma2((data_times[1] - obs_start)/days_in_year, 
                    alpha, mean_rate, lower.tail=FALSE, log.p=TRUE) -
                pgamma2((obs_end - data_times[length(data_times)])/days_in_year, 
                    alpha, mean_rate, lower.tail=FALSE, log.p=TRUE)

        }
    }

    if(return_vector){
        return(nll_site)
    }else{
        return(sum(nll_site))
    }
}


starting_par = c(0.5, pmax(0.003, 
unlist(lapply(eq_events, nrow))/
    ((cmt_end_time_julianDay1900 -cmt_start_time_julianDay1900)/days_in_year)))

# Basic fitting
if(FALSE){
    fit1$par = starting_par
    for(i in 1:50){
        fit1 = optim(fit1$par, nll_regional_gamma2, method='Nelder-Mead')
        print(fit1$value)
    }

    # Drop izumariana -- any impact?
    fit2 = optim(starting_par[-c(2,6)], f<-function(par) nll_regional_gamma2(par, data=eq_events[-c(1,5)]))
    # Just southamerica
    fit3 = optim(starting_par[c(1,8)], f<-function(par) nll_regional_gamma2(par, data=eq_events[7]))

    starting_par2 = starting_par * 2
    fit2$par = starting_par2
    for(i in 1:50){
        fit2 = optim(fit2$par, nll_regional_gamma2, method='Nelder-Mead')
        print(fit2$value)
    }

}
#
# Simulate some data
# 
simulate_fitted_model<-function(par, n=100){

    random_egs = lapply(as.list(par[-1]), 
        f<-function(x) cumsum(rgamma2(n, alpha=par[1], mean_rate=x)))

    return(random_egs)
}

# Plot synthetic data vs model
if(FALSE){
    random_egs = simulate_fitted_model(fit1$par)

    par(mfrow=c(length(random_egs), 1))
    par(mfrow=c(8,1))
    par(mar=c(0,0,0,0)+1)
    for(i in 1:length(random_egs)){
        plot(random_egs[[i]], rep(1, length(random_egs[[i]])), t='h', xlim=c(0, 50))
        title(names(eq_events)[i], line=-5, cex.main=2, col.main='red')
        points((eq_events[[i]]$julianDay1900 - cmt_start_time_julianDay1900)/days_in_year, 
            rep(1.5, length(eq_events[[i]]$julianDay1900)), 
            t='h', col='green')
    }

    nevents_obs = unlist(lapply(eq_events, nrow))
    nevents_model = replicate(100, 
        { sim = simulate_fitted_model(fit1$par); unlist(lapply(sim, f<-function(x) sum(x<(2017-1976)))) })
}

if(FALSE){
    #
    # Check for bias in parameter estimates
    #
    bias_check = replicate(1000,
        {sim = rgamma2(100, alpha=0.5, mean_rate=1/3); 
        sim_keep = which(cumsum(sim) < (2017-1976));
        if(length(sim_keep) > 0){
            sim = sim[sim_keep]
        }else{
            sim = sim[1]
        };
        fit = optim(par=c(0.4, 1/3), 
            f<-function(par) -sum(dgamma2(sim, par[1], par[2], log=TRUE)))$par;
        fit = c(fit, length(sim)) }
        )
}

#
# More complex testing for bias
#
# Simulate data 'like' the GCMT data (i.e. some given start/end time), with
# known parameters, and try to fit the data to it
#
#stop()
true_par = starting_par
max_j = 100 #1000
# If regular_samples = FALSE, do a 'large sample' test -- in that case, our method has little
# bias -- whereas for small samples, the bias is large
regular_sample = TRUE 
counts = matrix(NA, ncol=max_j, nrow=length(true_par)-1)
fitted_par = matrix(NA, ncol=max_j, nrow=length(true_par))
for(j in 1:max_j){

    print(j)
    random_data = simulate_fitted_model(true_par, n=500)

    for(i in 1:length(random_data)){
        # Make hypothetical time-series with duration = gcmt duration

        # Start 50 years before gcmt obs start time
        cm = random_data[[i]] * days_in_year + cmt_start_time_julianDay1900 - 50*days_in_year

        if(max(cm) < cmt_end_time_julianDay1900) stop()

        if(regular_sample){
            keep = which(cm > cmt_start_time_julianDay1900 & cm < cmt_end_time_julianDay1900)
            if(length(keep) > 0){
               random_data[[i]] = data.frame(julianDay1900 = cm[keep])
            }else{
               random_data[[i]] = data.frame(julianDay1900 = c())
            }
        }else{
            # Use this to keep all data, for large sample test purposes
            random_data[[i]] = data.frame(julianDay1900 = cm)
        }
    }

    counts[,j] = unlist(lapply(random_data, nrow))

    myfit = optim(true_par, nll_regional_gamma2, data=random_data)
    for(k in 1:5){
        myfit = optim(myfit$par, nll_regional_gamma2, data=random_data)
    }
        
    fitted_par[,j] = myfit$par
}

summary(t(fitted_par))
true_par

# Result: In small samples, significant bias

##  #
##  # Try with nls.lm
##  #
##  nll_regional_gamma2_b <-function(par){
##  
##      nll_regional_gamma2(par, return_vector=TRUE)
##  
##  }
##  
##  library(minpack.lm)
##  
##  nls.lm(par=fit1$par, fn=nll_regional_gamma2_b)
