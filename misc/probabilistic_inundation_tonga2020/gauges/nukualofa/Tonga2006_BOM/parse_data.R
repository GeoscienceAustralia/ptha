#
# This code de-tides a record of the 2006 Tonga tsunami from Nuku'alofa
# tide-gauge, then saves the data as a csv and makes a plot.
#

source('../spectral_highpass_filter.R')

parse_data<-function(){
    x = read.csv('nukualofa_tsu.txt')
    times = strptime(x[,1], format='%d/%m/%Y %H:%M', tz='Etc/UTC')
    juliant = julian(times)
    obs = x$Observed
    resid = x$Residual # The includes both tsunami and 'not-predicted-tide' components.

    # Filter the residual
    times_in_seconds = as.numeric((juliant - juliant[1]))*3600*24

    filtered_resid = spectral_highpass_filter(times_in_seconds, resid) 

    output = data.frame(time=times, juliant=juliant, height=obs,
                        resid=filtered_resid$highfreq)

    return(output)
}

#
t06 = parse_data()

png('Tonga2006_residual.png', width=10, height=4, units='in', res=300)
plot(t06$juliant, t06$resid, t='l', 
     xlim=c(t06$juliant[3500], t06$juliant[5000]), 
     xlab='Time (julian day)', ylab='Stage residual (m)',
     main="Tonga 2006 tsunami @ Nuku'alofa")
earthquake_start_time = strptime('2006/05/03 15:26:40', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC')
abline(v=julian(earthquake_start_time), col='red')
abline(v=julian(earthquake_start_time) + seq(0, 2, by=1/24), col='grey', lty='dotted')
abline(h=seq(-0.2, 0.2, by=0.1), col='grey', lty='dotted')
dev.off()

write.csv(t06, file='nukualofa_2006_tsunami_detided.csv', row.names=FALSE)
