#
# Offshore PTHA18 point -- plot max-stage vs exceedance-rate and compare with PTHA18
#
library(rptha)

run_series_name = commandArgs(trailingOnly=TRUE)[1]

#
# Get the PTHA18 access codes, needed for functions below
#
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

#
# Get the PTHA18 file scenaro rates for kermadectonga2 only
#
source_zone = 'kermadectonga2'
nc_file  = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
    'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', 'stochastic',
    '_slip_earthquake_events_', source_zone, '.nc')
fid = nc_open(nc_file, readunlim=FALSE)
rates_full_source = ncvar_get(fid, 'rate_annual')
nc_close(fid)

event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3458.3,
    slip_type='stochastic',
    all_source_names=source_zone)
max_stage_ptha18 = event_peak_stage_at_refpoint$kermadectonga2$max_stage
target_stages = c(0.001, seq(0.01, 20, by=0.01))
ptha18_curve = cbind(target_stages, sapply(target_stages, f<-function(x) sum(rates_full_source*(max_stage_ptha18 > x))))

#
# Get data from this study out of the RDS file
#
ptha18_3458_point_data = readRDS(paste0(run_series_name, '_depth_and_stage_exrate_curve_at_ptha18_point_3458.3.RDS'))
x = ptha18_3458_point_data$stage_exrate_curves

# Note stages will be the same for all entries of the list
stages = x[[4]]$stage
# Here we get the combined 'unsegmented + union-of-segments' exceedance-rate
exrates = 0.5*x[[4]]$exrate_stage + 
          0.5*(x[[1]]$exrate_stage + x[[2]]$exrate_stage + x[[3]]$exrate_stage)
# Here we use the importance-sampling weights, that do not precisely conserve the
# exceedance rates, but in theory have less bias.      
alternate_exrates = 0.5*x[[4]]$alternate_exrate_stage + 
          0.5*(x[[1]]$alternate_exrate_stage + x[[2]]$alternate_exrate_stage + x[[3]]$alternate_exrate_stage)

# Here is the PTHA18 mean-rate curve [all source-zones]
# This was manually copied from a PTHA18 shapefile output.
ptha18_values = matrix(
    c(0.133, 1/10,
      0.276, 1/25,
      0.450, 1/50,
      0.704, 1/100,
      1.190, 1/250,
      1.720, 1/500,
      2.520, 1/1000,
      4.000, 1/2500,
      5.520, 1/5000,
      7.230, 1/10000),
    ncol=2, byrow=TRUE)

# Plot comparing the results

png(paste0(run_series_name, '_stage_vs_rate_validation_at_ptha18_point_3458.png'), width=8, height=6, units='in', res=300)
plot(stages, alternate_exrates, log='y', t='o', 
     xlab='Max-stage (m)', ylab='Exceedance-rate (events/year)',
     main=paste0('Max-stage vs exceedance-rate at PTHA18 point 3458: PTHA18 vs Tonga study\n',
         'Tonga study includes a subset of kermadectonga2 scenarios with a 5 hour model. \n ',
         'The PTHA18 uses a 36 hour model, optionally with additional source-zones.'),
     cex.main=1.0, xlim=c(0, 10), ylim=c(1e-05, 1))
points(ptha18_values, col='red', pch=19, t='o') # Straight out of PTHA18, all source-zones
points(ptha18_curve, t='l', col='green', lwd=2) # Straight out of PTHA18, kermadectonga2 only
#points(stages, exrates, t='l', col='purple')
abline(v=seq(0, 20), col='orange', lty='dotted')
abline(h=10**seq(-6, 0), col='orange', lty='dotted')
add_log_axis_ticks(side=2)
legend('topright', 
       c('Tonga inundation study [only scenarios from "kermadectonga2"]', 
         'PTHA18 [all source-zones]', 'PTHA18 [kermadectonga2 only]'),
       col=c('black', 'red', 'green'), pch=c(1, 19, NA), lty=c(1, 1, 1))
dev.off()

# save.image('quick_plot_image.RDS')
