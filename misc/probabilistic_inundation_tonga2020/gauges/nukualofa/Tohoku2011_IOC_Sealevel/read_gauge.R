#
# This code de-tides a record of the 2011 Tohoku tsunami at Nuku'alofa, and
# writes an output file and a plot.
#

source('../spectral_highpass_filter.R')

# Nuku'alofa in Tonga
nuku = read.delim('nukualofa_harbour_gauge.txt', skip=3, stringsAsFactors=FALSE)
names(nuku) = c('chartime', 'stage')
nuku$time = strptime(nuku[,1], format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
nuku$time_s = as.numeric(julian(nuku$time)-julian(nuku$time[1]))*24*3600
nuku_filt = spectral_highpass_filter(nuku$time_s, nuku$stage)

trange = strptime(c('2011-03-11 06:00:00', '2011-03-14 00:00:00'), 
                  format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

# Plot
png('Tohoku_tsunami_obs_Nukualofa.png', width=10, height=7, units='in', res=300)
par(mfrow=c(2,1))

inds = which(nuku$time > trange[1] & nuku$time < trange[2])
plot(nuku$time[inds], nuku_filt$highfreq[inds], t='l', 
     main="Nuku'alofa tsunami wave-train following Tohoku tsunami [detided]",
     xlab='Time', ylab='De-tided stage (m)')
grid(col='orange')

plot(nuku$time[inds], nuku$stage[inds], main='As above with tides', t='l',
     xlab='Time', ylab='Stage (m)')
grid(col='orange')

dev.off()

# Export data in a time-range
trange_2 = strptime(c('2011-03-09 00:00:00', '2011-03-19 00:00:00'), 
                    format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
inds = which(nuku$time > trange_2[1] & nuku$time < trange_2[2])
output = data.frame(time=nuku$time[inds], 
                    juliant=julian(nuku$time[inds]), 
                    height=nuku$stage[inds],
                    resid=nuku_filt$highfreq[inds])
write.csv(output, file='nukualofa_tohoku_tsunami_detided.csv', row.names=FALSE)
