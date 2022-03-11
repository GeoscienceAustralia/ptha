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
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
ptha18 = new.env()
source(file_here, local=ptha18, chdir=TRUE)

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
event_Mw = round(event_peak_stage_at_refpoint$kermadectonga2$Mw, 1)
target_stages = c(0.001, seq(0.01, 20, by=0.01))
ptha18_curve = cbind(target_stages, sapply(target_stages, function(x) sum(rates_full_source*(max_stage_ptha18 > x))))

event_importance_weighted_sampling_probs = max_stage_ptha18*rates_full_source

#
# Get data from this study out of the RDS file
#
ptha18_3458_point_data = readRDS(paste0(run_series_name, 
   '_depth_and_stage_exrate_curve_at_ptha18_point_3458.3.RDS'))
x = ptha18_3458_point_data$stage_exrate_curves
# 
# Get the analytical variance
#
scenario_mw = round(ptha18_3458_point_data$scenarios_databases$unsegmented_HS$mw,1)
unique_Mw = sort(unique(scenario_mw))
unique_Mw_counts = sapply(unique_Mw, function(x) sum(scenario_mw == x))
samples_per_Mw_FUN = approxfun(unique_Mw, unique_Mw_counts, method='constant')
analytical_variance = sapply(target_stages, function(x){
    ptha18$analytical_Monte_Carlo_exrate_uncertainty(
        event_Mw = event_Mw, 
        event_rates = rates_full_source, 
        event_peak_stage = max_stage_ptha18, 
        stage_threshold=x, 
        samples_per_Mw=samples_per_Mw_FUN,
        event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)[2]
    })
# 95% analytical confidence intervals for the Monte Carlo result, based on the offshore PTHA and theory.
# We will compare the high-res model results with this
lower_CI = ptha18_curve[,2] + qnorm(0.025)*sqrt(analytical_variance)
upper_CI = ptha18_curve[,2] + qnorm(0.975)*sqrt(analytical_variance)

# Note stages will be the same for all entries of the list
stages = x[[4]]$stage
# Here we get the combined 'unsegmented + union-of-segments' exceedance-rate
exrates = 0.5*x[[4]]$exrate_stage + 
          0.5*(x[[1]]$exrate_stage + x[[2]]$exrate_stage + x[[3]]$exrate_stage)

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
for(i in 1:2){
    png(paste0(run_series_name, '_stage_vs_rate_validation_at_ptha18_point_3458_v', i, '.png'), 
        width=8, height=6, units='in', res=300)
    options(scipen=5)
    par(mar=c(5,6,3,2))
    plot(stages, exrates, log=ifelse(i==1, 'xy', 'y'), t='l', lwd=3,
         xlab="", ylab="",
         main='Tsunami-maxima exceedance-rate curve @ Site P',
         #main=paste0('Max-stage vs exceedance-rate at PTHA18 point 3458: PTHA18 vs Tonga study\n',
         #    'Tonga study includes a subset of kermadectonga2 scenarios with a 5 hour model. \n ',
         #    'The PTHA18 uses a 36 hour model, optionally with additional source-zones.'),
         cex.main=1.7, xlim=c(0.05, 10), ylim=c(1e-04, 1e-01), col='red', las=1)
    mtext(side=2, 'Exceedance-rate (events/year)', line=4, cex=1.6)
    mtext(side=1, 'Tsunami-maxima (m)', line=2, cex=1.6)
    #points(ptha18_values, col='red', pch=19, t='o') # Straight out of PTHA18, all source-zones
    points(ptha18_curve, t='l', col='blue', lwd=2) # Straight out of PTHA18, kermadectonga2 only
    points(ptha18_curve[,1], lower_CI, t='l', col='blue', lty='dashed')
    points(ptha18_curve[,1], upper_CI, t='l', col='blue', lty='dashed')
    #abline(v=seq(0, 20), col='orange', lty='dashed')
    abline(v=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50), col='orange', lty='dashed')
    abline(h=10**seq(-6, 0), col='orange', lty='dotted')
    add_log_axis_ticks(side=2)
    add_log_axis_ticks(side=1)
    legend('topright',  
           c('PTHA18 (Kermadec-Tonga only)',
             '95% CI (analytical) for the \n Monte-Carlo exceedance-rate'),
           col=c('blue', 'blue'), pch=c(NA, NA), lty=c('solid', 'dashed'), 
           pt.cex=c(NA, NA), bg=rgb(1,1,1,alpha=0.5), box.col=rgb(1,1,1,alpha=0.5), cex=1.5)
    legend('bottomleft',  
           c('High-resolution model with \n a single Monte-Carlo sample'),
           col=c('red'), lty=c('solid'),  lwd=3, 
           pt.cex=c(0.3), bg=rgb(1,1,1,alpha=0.5), box.col=rgb(1,1,1,alpha=0.5), cex=1.5)
    dev.off()
}




#
# COMPUTE PERCENTILES FOR THE HIGH-RES MODEL OFFSHORE
#




# Make arrays that nominally hold stage/elev for ALL SCENARIOS ON THIS
# SOURCE-ZONE IN THE PTHA. Actually they will only be non-NA for randomly
# sampled scenarios, however, the function we use them in
# (estimate_exrate_uncertainty) requires this structure
all_model_results = ptha18_3458_point_data$results_df

# The following line could alternatively use any segment with no change to results.
random_scenarios  = ptha18_3458_point_data$scenarios_databases$unsegmented_HS
scenarios_to_results_inds = match(random_scenarios$md_dir, all_model_results$md_dir)
random_scenario_stage = rep(NA, length(max(random_scenarios$inds)))
# The random_scenario stages which follow are only "!is.na()" at scenarios we
# actually simulated -- but that is all that is required.
random_scenario_stage[random_scenarios$inds] = all_model_results$max_stage[scenarios_to_results_inds] 
rm(random_scenarios, scenarios_to_results_inds, all_model_results) # Emphasise latter variables were just temporary 

peak_stages = 10**seq(log10(0.01), log10(10), by=0.005)

# Get the R session resulting from "compute_rates_all_sources.R", needed to 
# work with alternative logic-tree branches
ptha18_source_rate_env = new.env()
source(paste0(dirname(file_here), '/get_detailed_PTHA18_source_zone_info.R'),
       local=ptha18_source_rate_env, chdir=TRUE)

# All logic-tree-branches for the kermadectonga 2 unsegmented, and the segments
all_source_names = c('kermadectonga2', 'kermadectonga2_tonga', 'kermadectonga2_kermadec', 'kermadectonga2_hikurangi')

scenario_base = ifelse(grepl('_repeated', run_series_name),
    '../../sources/random_repeated/',
    '../../sources/random/')

all_source_samples = list(
    'kermadectonga2' = paste0(scenario_base, 'random_scenarios_kermadectonga2_unsegmented_HS.csv'),
    'kermadectonga2_tonga' = paste0(scenario_base, 'random_scenarios_kermadectonga2_tonga_segment_HS.csv'),
    'kermadectonga2_kermadec' = paste0(scenario_base, 'random_scenarios_kermadectonga2_kermadec_segment_HS.csv'),
    'kermadectonga2_hikurangi' = paste0(scenario_base, 'random_scenarios_kermadectonga2_hukurangi_segment_HS.csv'))
stopifnot(all(all_source_names == names(all_source_samples)))

# Read the random scenarios (note: scenarios are the same in all cases, but the nominal rates differ)
all_samples = list()
for(nm_i in all_source_names){
    all_samples[[nm_i]] = read.csv(all_source_samples[[nm_i]])
}

# Get the stage exceedance rates over all logic-tree branches, for unsegmented and segmented models
all_stage_exrates_all_logic_trees = list()
# For other functions below it is useful to store some components of these results differently.
# The following lists do that.
stage_exceedances_logic_tree_branch = list()
back_computed_mean_curve = list()
logic_tree_posterior_probs = list()
for(nm_i in all_source_names){
 
     # Make "source_zone" and "segment" to match the ptha18_source_rate_env information
 
     source_zone = 'kermadectonga2'
     if(source_zone == nm_i){
         # Unsegmented branch
         segment = ''
     }else{
         # Segmented branch -- extract the segment name directly
         segment = strsplit(nm_i, '_')[[1]][2]
     }
 
     # Compute exceedance rates for all peak_stages, for all logic-tree branches
     all_stage_exrates_all_logic_trees[[nm_i]] = 
         ptha18_source_rate_env$random_scenario_exceedance_rates_all_logic_tree_branches(
             source_zone = source_zone,
             segment = segment,
             random_scenarios = all_samples[[nm_i]],
             all_scenario_stage = random_scenario_stage, 
             threshold_stages=peak_stages)
 
     #
     # Alternative storage for some results, which is used in code below
     #
 
     # Exceedance rates for each logic tree branch, and each peak stage
     stage_exceedances_logic_tree_branch[[nm_i]] = 
         all_stage_exrates_all_logic_trees[[nm_i]]$logic_tree_branch_exceedance_rates
     # Average of the stage exceedance-rate curves over all logic tree branches
     # (accounting for posterior probs)
     back_computed_mean_curve[[nm_i]] = 
         all_stage_exrates_all_logic_trees[[nm_i]]$logic_tree_mean_exceedance_rates
     logic_tree_posterior_probs[[nm_i]] = 
         all_stage_exrates_all_logic_trees[[nm_i]]$logic_tree_branch_posterior_prob
 
}


segment_inds = 2:4 # List entries corresponding to segments
percentile_probs = c(0.025, 0.16, 0.5, 0.84, 0.975)

# Reproducible randomness
set.seed(123)
# Mixed segmented, unsegmented, comonotonic segments
mean_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_stage_exrates_all_logic_trees$kermadectonga2,
    all_stage_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=0.5,
    union_of_segments_wt=0.5,
    segments_copula_type='comonotonic',
    percentile_probs=percentile_probs)

# Get the PTHA18 stage exceedance-rate curve, for this source zone only, at the target point
# This contains percentile information that we will test against
stage_exrate_curve = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
    target_index=event_peak_stage_at_refpoint[[source_zone]]$target_index,
    source_name=source_zone)

plot_panel<-function(mean_curve_env, plot_title=''){
    # Used to plot the epistemic uncertainties below

    plot(peak_stages, mean_curve_env$mean_exrate, 
        t='l', log='xy', ylim=c(1e-04, 3e-01), xlim=c(0.05, 10), lwd=3, las=1, 
        xlab="", ylab="",  cex.lab=1.4, cex.axis=1.4, axes=FALSE, xaxs='i')
    axis(side=1, cex.axis=1.4)
    axis(side=2, at=10**seq(-5, 0), cex.axis=1.4,
         labels=c('1/100000', '1/10000', '1/1000', '1/100', '1/10', '1'), 
         las=1)
    add_log_axis_ticks(side=2)
    mtext(side=1, "Tsunami maxima (m)", cex=1.7, line=2.5)
    #mtext(side=2, "Mean exceedances per year ", cex=1.7, line=4.5)
    mtext(side=2, 'Exceedance-rate (events/year)', line=5.5, cex=1.6)
    title(main=plot_title, cex.main=2.2)
    points(peak_stages, mean_curve_env$percentile_exrate[3,],  t='l', col='brown', lwd=3)
    points(peak_stages, mean_curve_env$percentile_exrate[2,],  t='l', col='blue', lwd=2)
    points(peak_stages, mean_curve_env$percentile_exrate[4,],  t='l', col='blue', lwd=2)
    points(peak_stages, mean_curve_env$percentile_exrate[1,], t='l', col='green', lwd=2)
    points(peak_stages, mean_curve_env$percentile_exrate[5,], t='l', col='green', lwd=2)
    #abline(v=seq(0, 20), col='grey', lty='dotted')
    #abline(h=1/c(500, 2500, 10000), col='orange')
    abline(v=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50), col='orange', lty='dashed')
    abline(h=10**seq(-6, 0), col='orange', lty='dotted')
    add_log_axis_ticks(side=2)
    add_log_axis_ticks(side=1)
    legend('topright', c('Mean', 'Median', '16/84 %', '2.5/97.5 %'), 
           col=c('black', 'brown', 'blue', 'green'),
           lwd=c(3,3,2,2), pch=c(NA, NA, NA, NA), cex=1.3, bty='n')

    text(0.055, 1/10000 * 2, 'Solid lines: current model \nDashed lines: PTHA18', 
         adj=c(0.,0.5), cex=1.4)

}

#
# 2 panel plot -- one showing the PTHA18 mean curve and analytical
# uncertainties vs the new model, the other showing percentiles for PTHA18 and the new model
#
png(paste0(run_series_name, '_stage_vs_rate_validation_at_ptha18_point_3458_COMBINED.png'), 
    width=12, height=5, units='in', res=300)

# Panel 1
par(mfrow=c(1,2))
options(scipen=5)
par(mar=c(5,7,3,0.5))
plot(stages, exrates, log='xy', t='l', lwd=3,
     xlab="", ylab="",
     main='', #'Tsunami-maxima exceedance-rate curve @ Site P',
     cex.main=1.7, xlim=c(0.05, 10), ylim=c(1e-04, 3e-01), col='red', 
     las=1, axes=FALSE)
axis(side=1, cex.axis=1.4)
axis(side=2, at=10**seq(-5, 0), cex.axis=1.4,
     labels=c('1/100000', '1/10000', '1/1000', '1/100', '1/10', '1'), 
     las=1)
mtext(side=2, 'Exceedance-rate (events/year)', line=5.5, cex=1.6)
mtext(side=1, 'Tsunami-maxima (m)', line=2.5, cex=1.6)
points(ptha18_curve, t='l', col='blue', lwd=2) # Straight out of PTHA18, kermadectonga2 only
points(ptha18_curve[,1], lower_CI, t='l', col='blue', lty='dashed')
points(ptha18_curve[,1], upper_CI, t='l', col='blue', lty='dashed')
abline(v=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50), col='orange', lty='dashed')
abline(h=10**seq(-6, 0), col='orange', lty='dotted')
add_log_axis_ticks(side=2)
add_log_axis_ticks(side=1)
legend('bottomleft',  
       c('High-res Monte-Carlo',
         'PTHA18: Logic-tree mean',
         'PTHA18: 95% CI for Monte-\nCarlo exceedance-rate'
         ),
       col=c('red', 'blue', 'blue'), pch=c(NA, NA, NA), 
       lty=c('solid', 'solid', 'dashed'), lwd = c(3,1,1),
       pt.cex=c(0.3, NA, NA), bg=rgb(1,1,1,alpha=0.0), 
       box.col=rgb(1,1,1,alpha=0.0), cex=1.4, yjust=1)
title(main = 'A) Monte-Carlo uncertainties', cex.main = 1.8)

# Panel 2
plot_panel(mean_curve_env)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate, t='l', lty='dashed', col='black', lwd=1)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_median, t='l', lty='dashed', col='brown', lwd=1)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_16pc, t='l', lty='dashed', col='blue', lwd=1)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_84pc, t='l', lty='dashed', col='blue', lwd=1)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_lower_ci, t='l', lty='dashed', col='green', lwd=1)
points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_upper_ci, t='l', lty='dashed', col='green', lwd=1)

title(main = 'B) Epistemic uncertainties', cex.main = 1.8)

dev.off()

