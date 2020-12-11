library(rptha)

run_series_name = commandArgs(trailingOnly=TRUE)[1]

#rds_name = 'ptha18_tonga_MSL0_depth_and_stage_exrate_curve_at_parliament.RDS'
#rds_name = 'ptha18_tonga_MSL0_meshrefine2_depth_and_stage_exrate_curve_at_parliament.RDS'
#rds_name = 'ptha18_tonga_MSL0.8_depth_and_stage_exrate_curve_at_parliament.RDS'
rds_name = paste0(run_series_name, '_depth_and_stage_exrate_curve_at_parliament.RDS')

parliament_data = readRDS(rds_name)
x = parliament_data$stage_exrate_curves
coord = parliament_data$results_df[1,c('lon', 'lat')]

# Here we use the self-normalised importance-sampling weights. These conserve the 
# overall exceedance-rate (because they are forced to sum to 1) but they lead to 
# finite-sample bias in integrals (although asymptotically unbiased)
unsegmented_vals = cbind(x$unsegmented_HS$depth, x$unsegmented_HS$exrate_depth)
segmented_vals = cbind(x$unsegmented_HS$depth, 
    (x[["hukurangi_segment_HS"]]$exrate_depth + 
     x[["kermadec_segment_HS" ]]$exrate_depth +
     x[["tonga_segment_HS"    ]]$exrate_depth ))

# Here we use the importance-sampling weights that do not conserve the exceedance-rates,
# but in theory are unbiased irrespective of sample size.
alternate_unsegmented_vals = cbind(x$unsegmented_HS$depth, x$unsegmented_HS$alternate_exrate_depth)
alternate_segmented_vals = cbind(x$unsegmented_HS$depth, 
    (x[["hukurangi_segment_HS"]]$alternate_exrate_depth + 
     x[["kermadec_segment_HS" ]]$alternate_exrate_depth +
     x[["tonga_segment_HS"    ]]$alternate_exrate_depth ))

png(paste0('Tsunami_inundation_depth_vs_frequency_parliament_', rds_name, '.png'),
    width=7.5, height=6, units='in', res=300)
par(mar=c(4.5, 6.5, 3, 2))
par(family='serif')
plot(alternate_unsegmented_vals, t='l', lty='dashed',
     xlab='', ylab='', log='xy', col='green',
     xlim=c(0.1, 10), ylim=c(1.0e-05, 1/100), axes=FALSE)
points(alternate_segmented_vals, t='l', col='blue', lty='dashed')
axis(side=1, cex.axis=1.3)
axis(side=2, at=10**seq(-5, -2), cex.axis=1.3,
     labels=c('1/100000', '1/10000', '1/1000', '1/100'), 
     las=1)
add_log_axis_ticks(side=2)
add_log_axis_ticks(side=1)
abline(v=c(0.1, 0.2, 0.5, 1, 2, 5, 10), col='darkgrey', lty='dotted')
abline(h=10**seq(-5, -2), col='darkgrey', lty='dotted')
mtext('Tsunami inundation depth at Parliament (m)', side=1, line=2.5, cex=1.4)
mtext('Mean exceedances per year', side=2, line=5.5, cex=1.4)

mean_vals = 0.5*(unsegmented_vals + segmented_vals)
#points(mean_vals, t='l', col='purple')

alternate_mean_vals = 0.5*(alternate_unsegmented_vals + alternate_segmented_vals)
points(alternate_mean_vals, t='l', lwd=3, col='black')

title(paste0('Tsunami inundation depth exceedance-rate at Parliament \n (', 
   'lon=', round(coord[1],4), ', lat=', round(coord[2],4), ')'))

legend('bottomleft', 
       c('PTHA18: logic-tree mean exceedance-rate', 
         '(assuming unsegmented source-zone)', 
         '(assuming segmented source-zone)',
         '2% in 50 years'),
       lty=c('solid', 'dashed', 'dashed', 'solid'), 
       lwd=c(3, 1, 1, 1), col=c('black', 'green', 'blue', 'red'),
       cex=1.1)

target_exrate = uniroot(f<-function(x){(1 - exp(-50*x)) - 0.02}, lower=1/3000, upper=1/2000, tol=1e-12)$root
target_depth = approx(alternate_mean_vals[,2], alternate_mean_vals[,1], xout=target_exrate)$y
points(c(1.0e-04, target_depth), c(target_exrate, target_exrate), col='red', t='l')
points(c(target_depth, target_depth), c(1.0e-10, target_exrate), col='red', t='l')
dev.off()

# Some useful stats
depth_at_target_exrate = approx(mean_vals[,2], mean_vals[,1], xout=target_exrate, ties='min')$y
print(c('depth-at-target-exrate (self-normalised importance-sampling): ', depth_at_target_exrate))
depth_at_target_exrate = approx(alternate_mean_vals[,2], alternate_mean_vals[,1], xout=target_exrate, ties='min')$y
print(c('depth-at-target-alternate_exrate (regular importance sampling): ', depth_at_target_exrate))

exrate_10cm = approx(mean_vals[,1], mean_vals[,2], xout=0.1)$y
print(c('exrate-of-10cm-inundation (self-normalised importance-sampling): ', exrate_10cm))
print(c('   chance in 50 years    : ', 1 - exp(-50 * exrate_10cm)))

exrate_10cm = approx(alternate_mean_vals[,1], alternate_mean_vals[,2], xout=0.1)$y
print(c('alternate_exrate-of-10cm-inundation (regular importance-sampling): ', exrate_10cm))
print(c('   chance in 50 years    : ', 1 - exp(-50 * exrate_10cm)))
