#
# This code takes the data provided by BOM and produces a de-tided series
# corresponding to the Chile 2010 event
#

source('../spectral_highpass_filter.R')

make_chile2010<-function(){
    # Chile 2010
    x = read.delim('20100226_nuku7.csv')
    names(x)[1] = 'time'
    x$time = strptime(x$time, '%d/%m/%Y %H:%M', tz='Etc/UTC')
    x$juliant = julian(x$time)
    # Remove long-periods from the residual
    x_filt = spectral_highpass_filter(as.numeric(x$juliant - x$juliant[1])*3600*24,
                                      x$Residual)
    output = data.frame(time=x$time, juliant=x$juliant, height=x$Observed, 
                        resid=x_filt$highfreq)
    write.csv(output, 'nukualofa_chile2010_detided.csv', row.names=FALSE)

    time_range = strptime(c('2010-02-27 00:00:00', '2010-03-02 00:00:00'), 
                      format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    k = which(output$time > time_range[1] & output$time < time_range[2])
    png('nukualofa_chile2010.png', width=14, height=7.5, units='in', res=300)
    par(mfrow=c(2,1))
    plot(output$time[k], output$resid[k], t='l', xlab='Time', ylab='Stage')
    grid()
    title("Chile 2010 tsunami @ Nuku'Alofa")
    plot(output$time[k], output$height[k], t='l', xlab='Time', ylab='Stage')
    grid()
    title("As above with tides")
    dev.off()

}
make_chile2010()



## Smaller 2009 event [2009/03/19, Mw 7.65]. Use a small time-window to work
## around any artefacts.
#x = read.csv('nhw20092010.csv')
#names(x) = c('time', 'stage', 'pred', 'resid')
#x$time = strptime(x$time, '%d/%m/%Y %H:%M', tz='Etc/UTC')
#x$juliant = julian(x$time)
#time_range = strptime(c('2009-03-19 12:00:00', '2009-03-21 00:00:00'), 
#                      format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
#k = which(x$time > time_range[1] & x$time < time_range[2])

# Event around 2009/09/29 -- big tsunami in Samoa.
make_samoa2009<-function(){

    # 2009-2010, various
    x = read.csv('nhw20092010.csv')
    names(x) = c('time', 'stage', 'pred', 'resid')
    x$time = strptime(x$time, '%d/%m/%Y %H:%M', tz='Etc/UTC')
    x$juliant = julian(x$time)

    # Filter to about 1 month around the 2009 event.
    event_start = julian(strptime('2009/09/20 00:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC'))
    event_end   = julian(strptime('2009/10/07 00:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC'))
    k0 = which.min(abs(x$juliant - event_start))
    k1 = which.min(abs(x$juliant - event_end))
    x = x[k0:k1,]

    x_filt = spectral_highpass_filter(as.numeric(x$juliant - x$juliant[1])*3600*24,
                                      x$resid)
    output = data.frame(time=x$time, juliant=x$juliant, height=x$stage, 
                        resid=x_filt$highfreq)
    write.csv(output, 'nukualofa_samoa2009_detided.csv', row.names=FALSE)

    # Make a plot
    png('nukualofa_samoa2009.png', width=14, height=7.5, units='in', res=300)
    event_start = julian(strptime('2009/09/29 12:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC'))
    event_end   = julian(strptime('2009/09/30 12:00:00', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC'))
    k0 = which.min(abs(x$juliant - event_start))
    k1 = which.min(abs(x$juliant - event_end))
    k = seq(k0, k1)
    par(mfrow=c(2,1))
    plot(output$time[k], output$resid[k], t='l', xlab='Time', ylab='Stage')
    grid()
    title("Samoa 2009 tsunami @ Nuku'Alofa")
    plot(output$time[k], output$height[k], t='l', xlab='Time', ylab='Stage')
    grid()
    title("As above with tides")
    dev.off()
}
make_samoa2009()
