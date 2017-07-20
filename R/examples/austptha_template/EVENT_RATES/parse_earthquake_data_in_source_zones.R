#
# Code for casual analysis of earthquakes in our source zones
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

    plot(unit_source_grid_poly[[i]], main=names(unit_source_grid_poly)[i], border='grey', axes=TRUE)
    points(eq_events[[i]]$cent_lon, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')
    # Make sure we don't miss points due to longitude convention
    points(eq_events[[i]]$cent_lon-360, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')
    points(eq_events[[i]]$cent_lon+360, eq_events[[i]]$cent_lat, cex=(eq_events[[i]]$Mw-6.5)**2, col='red')

    #
    # Same plot as above, with beachballs
    #

    plot(unit_source_grid_poly[[i]], main=names(unit_source_grid_poly)[i], border='grey', axes=TRUE)

    pol_centroid = as.numeric(coordinates(gCentroid(unit_source_grid_poly[[i]])))

    for(j in 1:length(eq_events[[i]]$cent_lon)){

        lon_lat = c(eq_events[[i]]$cent_lon[j], eq_events[[i]]$cent_lat[j])
        lon_lat = adjust_longitude_by_360_deg(lon_lat, pol_centroid)
    
        # Beach-ball details
        mec = SDRfoc(eq_events[[i]]$strk1[j], eq_events[[i]]$dip1[j], eq_events[[i]]$rake1[j], PLOT=FALSE)

        fcol =  c( rgb(1,0,0, alpha=1), 'lightblue', 'green', 'orange', 'yellow', 'purple', 'black')[foc.icolor(mec$rake1)]

        par(lwd=0.2) # Try to avoid strong lines on beachballs
        try(
            justfocXY(mec, lon_lat[1], lon_lat[2],
                focsiz=4*((eq_events[[i]]$Mw[j]-6)/10)**2, xpd=TRUE, fcol=fcol,
                fcolback='white'),
            silent=TRUE
            )
        par(lwd=1)

    }
}

dev.off()

