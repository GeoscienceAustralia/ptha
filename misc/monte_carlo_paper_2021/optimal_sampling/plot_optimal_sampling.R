#
# Explore variance reductions associated with importance sampling, and varying
# the number of scenarios in each magnitude bin.
#

sink(file='log_of_plot_optimal_sampling.log')

# Get the PTHA18 access script in its own environment
ptha18_access_script = '../../../ptha_access/get_PTHA_results.R'
#ptha18_access_script = '../../../../../AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
if(!file.exists(ptha18_access_script)) stop('You need to update the path to the get_PTHA_results.R script')
ptha18 = new.env()
source(ptha18_access_script, local=ptha18, chdir=TRUE)

#
# Get scenario data for the testing
#

# Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in PTHA18)
source_zone = 'kermadectonga2'
kt2_scenarios = ptha18$get_source_zone_events_data(source_zone,  slip_type='stochastic')

# Convenient shorthand for the magnitudes and rates in the event table
event_Mw = kt2_scenarios$events$Mw
event_rates = kt2_scenarios$events$rate_annual

# Useful wrapper for extracting peak-stage values for all scenarios at a given site
get_peak_stage_at_target_point<-function(target_point){
    event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
        target_point = target_point,
        slip_type='stochastic',
        all_source_names=source_zone)
    return(event_peak_stage_at_refpoint[[source_zone]]$max_stage)
}

# Get the peak stage value at a point east of Tongatapu
target_point = c(185.1239, -21.0888)
event_peak_stage_target_point = get_peak_stage_at_target_point(target_point)

TOTAL_SAMPLES = 1200 # How many samples in the Monte-Carlo scheme?

DEFAULT_MC_REPS = 10000 # When we repeat MC-sampling many times, by default this defines 'many'

# If we were to make a synthetic catalogue with the same number of scenarios, how
# many years would it cover?
EQUIVALENT_SYNTHETIC_CATALOGUE_DURATION = (TOTAL_SAMPLES / sum(event_rates))

# How many samples to take if we sample the same amount in all Mw bins (excluding impossible bins)?
unique_mw = seq(7.2, 9.8, by=0.1)
const_samples = (unique_mw < 9.65)
const_samples = const_samples/sum(const_samples)*TOTAL_SAMPLES

#
# Plot the hazard curve, and show convergence of our estimators in the case of
# stratified sampling with constant sampling effort in "possible" Mw bins
#
plot_hazard_curve<-function(
    sampling_type='stratified',
    event_rates = kt2_scenarios$events$rate_annual,
    Number_MC_reps=DEFAULT_MC_REPS,
    hist_xlim=c(6e-04, 2e-03),
    event_importance_weighted_sampling_probs = NULL,
    event_peak_stage_local = event_peak_stage_target_point,
    mw_sampling_fun = NULL,
    fig_title=NULL,
    add_equivalent_synthetic_catalogue=FALSE,
    add_95pc_analytical_CI=TRUE,
    add_hardcoded_normal_distribution_to_second_panel = FALSE){

    if(is.null(fig_title)) stop('Must provide fig title')

    #
    # Setup defaults for stratified and stratified_importance sampling
    #
    if(sampling_type == 'stratified'){

        # No importance sampling
        if(!is.null(event_importance_weighted_sampling_probs)){
            stop('Do not specify event_importance_weighted_sampling_probs with stratified sampling')
        }
        # Weights equivalent to regular stratified sampling
        event_importance_weighted_sampling_probs = event_rates

    }else if(sampling_type == 'stratified_importance'){

        # By default use stage*rate at the target point
        if(is.null(event_importance_weighted_sampling_probs)){
            event_importance_weighted_sampling_probs = event_rates*event_peak_stage_target_point
        }

    }else{
        stop('unknown sampling_type')
    }


    if(is.null(mw_sampling_fun)){
        # How many samples in each Mw bin? By default use a constant for 'possible' scenarios.
        mw_sampling_fun = function(Mw){ const_samples[1]*(Mw < 9.65) }
    }

    png(fig_title, width=9, height=9, units='in', res=300)

    par(mfrow=c(2,1))
    par(mar=c(4,7,2,1))
    options(scipen=5)

    peak_stage_vals = seq(0.01, 15, by=0.01)
    hazard_curve = sapply(peak_stage_vals,
        function(x) sum(event_rates * (event_peak_stage_local > x)))
    plot(peak_stage_vals, hazard_curve, log='xy', t='l', lwd=2, las=1,
        xlab="",
        ylab='', #'Exceedance-rate (events/year)',
        cex.lab=1.5, cex.axis=1.3, 
        ylim=c(1/10000, 1/10), xlim=c(0.02, 10))
    mtext(side=1, "Tsunami maximum-stage (from offshore PTHA)", line=2.5, cex=1.5)
    mtext(side=2, 'Exceedance-rate (events/year)', line=5.5, cex=1.5)
    add_log_axis_ticks(side=1)
    add_log_axis_ticks(side=2)
    #abline(h=c(1/10, 1/25, 1/100, 1/500, 1/2500, 1/10000), col='darkgreen', lty='dashed')
    abline(h=c(1/10, 1/100, 1/1000, 1/10000), col='darkgreen', lty='dashed')

    #
    # Compute multiple Monte-Carlo samples, store their statistics at target_stage, and plot some
    #
    target_stage = 2

    exrate_ts_store = list(mean=rep(NA, Number_MC_reps), var=rep(NA, Number_MC_reps), analytical_mean_var=rep(NA, 2))

    # Store the analytical mean/variance as well, for cross-checking
    exrate_ts_store$analytical_mean_var = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
        event_Mw,
        event_rates,
        event_peak_stage_local,
        stage_threshold=target_stage,
        samples_per_Mw=mw_sampling_fun,
        event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)

    set.seed(1234) # Reproducible randomness
    for(i in 1:Number_MC_reps){

        random_scenarios = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
            event_rates,
            event_Mw,
            samples_per_Mw = mw_sampling_fun,
            event_importance_weighted_sampling_probs = event_importance_weighted_sampling_probs)

        # For both cases this will compute the correct mean/variance
        mean_and_variance = ptha18$estimate_exrate_uncertainty(
            random_scenarios, event_peak_stage_local, target_stage)

        exrate_ts_store$mean[i] = mean_and_variance[1]
        exrate_ts_store$var[i] = mean_and_variance[2]

        # Plot the first 500 samples (not too many to avoid cluttering the plot)
        if(i <= 500){
            random_scenario_stages = event_peak_stage_local[random_scenarios$inds]
            random_hazard_curve = sapply(peak_stage_vals,
                function(x){
                    sum(random_scenarios$importance_sampling_scenario_rates_basic*
                        (random_scenario_stages > x), na.rm=TRUE)
                })

            points(peak_stage_vals, random_hazard_curve, t='l', lwd=0.2, col='grey')
        }
    }
    points(peak_stage_vals, hazard_curve, t='l', lwd=2) # Put the hazard curve back on (it is covered by gray curves)


    if(!add_equivalent_synthetic_catalogue){
        # Different legends depending on 95% analytical CIs
        if(!add_95pc_analytical_CI){
            legend('bottomleft', c('All scenarios in offshore PTHA',
                   paste0('Monte-Carlo estimates (500 only)')), 
                   lty=c(1, 1), col=c('black', 'grey'), lwd=c(2, 2), cex=1.3)
        }else{
            legend('bottomleft', 
                   c('All scenarios in offshore PTHA',
                       paste0('Monte-Carlo estimates (500 only)'),
                       '95% interval (analytical)'), 
                   lty=c('solid', 'solid', 'dotdash'), col=c('black', 'grey', 'darkred'), 
                   lwd=c(2, 2, 2), cex=1.3)
        }

    }else{

        # Different legends depending on 95% analytical CIs
        if(!add_95pc_analytical_CI){

            legend('bottomleft', 
                   c('All scenarios in offshore PTHA', 
                     paste0('Monte-Carlo estimates (500 only)'), 
                     'Equivalent Synthetic Catalogue 95% interval'),
                   lty=c('solid', 'solid', 'dashed'), col=c('black', 'grey', 'darkblue'), 
                   lwd=c(2, 1, 2), cex=1.3, bg=rgb(1,1,1,alpha=0.7), bty='o', box.col=rgb(1,1,1,alpha=0.3))
        }else{

            legend('bottomleft',
                   c('All scenarios in offshore PTHA', 
                     paste0('Monte-Carlo estimates (500 only)'), 
                     '95% interval (analytical)', 
                     'Equivalent Synthetic Catalogue 95% interval'),
                   lty=c('solid', 'solid', 'dotdash', 'dashed'), col=c('black', 'grey', 'darkred', 'darkblue'), 
                   lwd=c(2, 1, 2, 2), cex=1.3, bg=rgb(1,1,1,alpha=0.7), bty='o', box.col=rgb(1,1,1,alpha=0.3))

        }

        # Include 95% intervals for a synthetic catalogue
        equivalent_synthetic_lower = sapply(hazard_curve, function(x){
            qpois(0.025, lambda=x*EQUIVALENT_SYNTHETIC_CATALOGUE_DURATION)/EQUIVALENT_SYNTHETIC_CATALOGUE_DURATION
               })
        equivalent_synthetic_upper = sapply(hazard_curve, function(x){
            qpois(0.975, lambda=x*EQUIVALENT_SYNTHETIC_CATALOGUE_DURATION)/EQUIVALENT_SYNTHETIC_CATALOGUE_DURATION
               })


        # Add lines to the plot -- because it is log-log we should not use zeros -- instead use a very small
        # threshold
        points(peak_stage_vals, pmax(equivalent_synthetic_lower, 1e-100), t='l', col='darkblue', lty='dashed', lwd=2)
        points(peak_stage_vals, pmax(equivalent_synthetic_upper, 1e-100), t='l', col='darkblue', lty='dashed', lwd=2)

    }


    abline(v=target_stage, col='purple', lwd=2)

    mean_stoc = mean(exrate_ts_store$mean)
    var_stoc = var(exrate_ts_store$mean)

    # Get the 'analytical' mean and variance expected for the sampling method
    analytical_mean_variance = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
        event_Mw, event_rates, event_peak_stage_local, stage_threshold=target_stage,
        samples_per_Mw = mw_sampling_fun,
        event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)
    mean_analytical = analytical_mean_variance[1]
    var_analytical = analytical_mean_variance[2]

    if(add_95pc_analytical_CI){
        # Add an analytical confidence interval to the plot, using the normal approximation

        lower_CI = rep(NA, length(peak_stage_vals))
        upper_CI = rep(NA, length(peak_stage_vals))

        for(i in 1:length(peak_stage_vals)){
            tmp = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
                event_Mw, event_rates, event_peak_stage_local, stage_threshold=peak_stage_vals[i],
                samples_per_Mw = mw_sampling_fun,
                event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)
            lower_CI[i] = tmp[1] + qnorm(0.025)*sqrt(tmp[2])
            upper_CI[i] = tmp[1] + qnorm(0.975)*sqrt(tmp[2])
        }
        points(peak_stage_vals, lower_CI, t='l', col='darkred', lty='dotdash', lwd=2)
        points(peak_stage_vals, upper_CI, t='l', col='darkred', lty='dotdash', lwd=2)
    }

    #
    # Histogram of errors at the target_stage
    #

    hist(exrate_ts_store$mean, breaks=100, freq=FALSE, main='',
         xlab="", ylab="", xlim=hist_xlim, cex.axis=1.3, las=1)
    mtext(side=1,
          bquote(paste('Distribution of ', .(Number_MC_reps),
              ' Monte-Carlo exceedance-rates @ ', Q^T, '=',
              .(target_stage), ' m')),
          line=2.8, cex=1.5)
    mtext(side=2, 'Probability Density', line=4, cex=1.5)
    # Add normal distribution
    xs_local = seq(min(exrate_ts_store$mean), max(exrate_ts_store$mean), len=500)
    points(xs_local, dnorm(xs_local, mean=mean_analytical, sd=sqrt(var_analytical)), t='l',
           col='darkred', lty='dotdash', lwd=3)

    if(!add_hardcoded_normal_distribution_to_second_panel){
        legend('topright', 
               paste0("Normal distribution (analytical\n",
                      "mean and variance)"),
               lwd=3, col='darkred', lty='dotdash', pch=NA, cex=1.25, bty='n')
    }else{
        # Useful when we want to compare the spread of Monte-Carlo results with other results
        x_vals = seq(hist_xlim[1], hist_xlim[2], len=201)
        y_vals = dnorm(x_vals, mean=0.00122335093286598, sd=sqrt(0.0000000304813319339911))
        points(x_vals, y_vals, t='l', col='skyblue', lwd=3, lty='dotted')

        legend('topright', c(
               paste0("Normal distribution (analytical\n",
                      "mean and variance)"),
               paste0("Stratified-sampling (Figure 2)")),
               lwd=c(3,3), col=c('darkred', 'skyblue'), lty=c('dotdash', 'dotted'), 
               pch=c(NA,NA), cex=1.25, bty='n')
    }
    dev.off()

    #
    # Summary results
    #

    print(c('Summary results for ', sampling_type))
    print(c('  mean_analytical: ', mean_analytical))
    print(c('  mean_stoc : ', mean_stoc))
    print(c('      ratio : ', mean_analytical/mean_stoc))
    print(c('      %err  : ', (1 - mean_stoc/mean_analytical)*100))
    print(c('  var_analytical: ', var_analytical))
    print(c('  var_stoc : ', var_stoc))
    print(c('     ratio : ', var_analytical/var_stoc))
    print(c('     %err  : ', (1-var_stoc/var_analytical)*100))
    print(c('  sd_analytical: ', sqrt(var_analytical)))
    print(c('  sd_stoc : ', sqrt(var_stoc)))
    print(c('    ratio : ', sqrt(var_analytical/var_stoc)))
    print(c('    %err  : ', (1 - sqrt(var_stoc/var_analytical))*100))

    # Empirical confidence interval 95% true coverage
    coverage_CI = mean( (mean_analytical > exrate_ts_store$mean + qnorm(0.025)*sqrt(exrate_ts_store$var)) &
                        (mean_analytical < exrate_ts_store$mean + qnorm(0.975)*sqrt(exrate_ts_store$var)) )
    print(c('  Empirical confidence interval coverage (ideal 0.95): ', coverage_CI))

}

plot_hazard_curve('stratified', fig_title = 'Exceedance_rate_stratified_target_point.png', 
    add_equivalent_synthetic_catalogue=TRUE)
plot_hazard_curve('stratified_importance', fig_title = 'Exceedance_rate_stratified_importance_target_point.png',
    add_hardcoded_normal_distribution_to_second_panel=TRUE)

#stop()
#
# Optimal sampling
#

# Which threshold-stage values should we do the computation for?
threshold_stages = c(1, 2, 4, 8)

# Determine the optimal number of samples for each threshold
optimal_samples = lapply(threshold_stages, function(x){
    ptha18$get_optimal_number_of_samples_per_Mw(
        event_Mw,
        event_rates,
        event_peak_stage_target_point,
        stage_threshold=x,
        total_samples=TOTAL_SAMPLES) # No importance sampling
    })
# As above, using importance sampling
optimal_samples_IS = lapply(threshold_stages, function(x){
    ptha18$get_optimal_number_of_samples_per_Mw(
        event_Mw,
        event_rates,
        event_peak_stage_target_point,
        stage_threshold=x,
        total_samples=TOTAL_SAMPLES,
        event_importance_weighted_sampling_probs = (event_peak_stage_target_point*event_rates)) # With importance sampling
    })

#
# Develop a 'compromise' optimal number of samples, by averaging the results above
# with the constant-sampling result
#
sum_sampling    = optimal_samples[[1]]$Nsamples*0
sum_sampling_IS = optimal_samples_IS[[1]]$Nsamples*0
for(i in 1:length(optimal_samples)){
    sum_sampling    = sum_sampling    + optimal_samples[[i]]$Nsamples
    sum_sampling_IS = sum_sampling_IS + optimal_samples_IS[[i]]$Nsamples
}
mean_optimal = sum_sampling/length(optimal_samples)
mean_optimal_IS = sum_sampling_IS/length(optimal_samples)

# Proposed approach puts % weight on constant samples, and remaining weight on the mean_optimal samples
WEIGHT_ON_CONSTANT_SAMPLES = 0.25
r1 = (1-WEIGHT_ON_CONSTANT_SAMPLES)
r2 = WEIGHT_ON_CONSTANT_SAMPLES
mean_optimal_including_constant    = (r1*mean_optimal    + r2*const_samples)
mean_optimal_IS_including_constant = (r1*mean_optimal_IS + r2*const_samples)

# Write the case with importance-sampling to a file so we can use it later
write.csv(
    data.frame(unique_mw=unique_mw,
               mean_optimal_samples_IS_including_constant=mean_optimal_IS_including_constant),
    file='Non_uniform_sampling_effort_compromise_stratifiedImportance.csv',
    row.names=FALSE)

# Sanity check -- we always have TOTAL_SAMPLES
stopifnot(isTRUE(all.equal(sum(mean_optimal), TOTAL_SAMPLES)))
stopifnot(isTRUE(all.equal(sum(mean_optimal_IS), TOTAL_SAMPLES)))
stopifnot(isTRUE(all.equal(sum(mean_optimal_including_constant), TOTAL_SAMPLES)))
stopifnot(isTRUE(all.equal(sum(mean_optimal_IS_including_constant), TOTAL_SAMPLES)))

#
# How good are the above choices? Print some statistics
#

# Convenience functions for printing stats below
get_variances<-function(optimal_samples_i, mean_optimal){

    # Get the variance if we optimizes the sampling strategy to minimise it, for this particular
    # site and stage-threshold. This is the 'best-case' in terms of "number of samples' but in
    # general we wouldn't do it (because we need to consider a range of return periods)
    variance_best = sum(optimal_samples_i$variance_numerator/optimal_samples_i$Nsamples, na.rm=TRUE)

    # Get the variance if the comprimise sampling strategy is used
    variance_chosen = sum(optimal_samples_i$variance_numerator/mean_optimal, na.rm=TRUE)

    # Get the variance if we used constant sampling in all bins
    variance_const = sum(optimal_samples_i$variance_numerator/const_samples, na.rm=TRUE)

    return(list(variance_chosen=variance_chosen,
                chosen_on_best=variance_chosen/variance_best,
                constant_on_chosen=variance_const/variance_chosen,
                constant_on_best=variance_const/variance_best))
}
print_variances<-function(optimal_samples, mean_optimal, mean_optimal_including_constant){

    for(i in 1:length(optimal_samples)){
        vars = get_variances(optimal_samples[[i]], mean_optimal)
        print(c('No constant  ', signif(unlist(vars), 4)), width=999)
        vars = get_variances(optimal_samples[[i]], mean_optimal_including_constant)
        print(c('With constant', signif(unlist(vars), 4)), width=999)
    }

}

#
# How does the variance change with the above sampling effort, vs 'constant sampling', and
# vs 'the best possible for the chosen threshold'
#

print('Variances (no importance sampling) with different scenario counts')
print_variances(optimal_samples, mean_optimal, mean_optimal_including_constant)

print('Variances (IMPORTANCE SAMPLING) with different scenario counts')
print_variances(optimal_samples_IS, mean_optimal_IS, mean_optimal_IS_including_constant)

#
# Do some plotting of alternative optimal sampling efforts
#
plot_sampling_effort<-function(
    optimal_samples, mean_optimal, mean_optimal_including_constant,
    optimal_samples_IS, mean_optimal_IS, mean_optimal_IS_including_constant,
    const_samples, threshold_stages){

    # Various plot pars
    PLOT_YLIM = c(0, 350)
    BAR_GROUP = 80
    BAR_LWD = 4
    MIN_BAR_HT = 0*(const_samples > 0) # To show a bar can try using a small non-zero value
    BAR_H_OFFSET = 1.5

    library(cptcity)
    COLZ = cpt(pal='cb_seq_YlOrRd_06', n=6)[3:6]

    png('Optimal_sampling_effort.png', width=9, height=6, units='in', res=300)
    par(mar=c(4,4.5,2,1))
    par(mfrow=c(2,1))
    # Plot without importance sampling
    for(i in 1:length(threshold_stages)){
        # Stagger the x location of the bars, so we fit several bars close to the
        # desired Mw value

        if(i == 1){

            plot(optimal_samples[[i]]$Mw + (i-BAR_H_OFFSET)/BAR_GROUP,
                 pmax(optimal_samples[[i]]$Nsamples, MIN_BAR_HT),
                 t='h', lwd=BAR_LWD, lend=1, col=COLZ[i], ylim=PLOT_YLIM,
                 xlab='', ylab='Optimal # Samples', las=1,
                 main='Stratified sampling', cex.main=1.8, cex.lab=1.5)

            mtext("Mw", side=1, line=2, cex=1.5)
            axis(side=1, at=seq(7.2, 9.6, by=0.1), labels=FALSE)

        }else{

            points(optimal_samples[[i]]$Mw + (i-BAR_H_OFFSET)/BAR_GROUP,
                   pmax(optimal_samples[[i]]$Nsamples, MIN_BAR_HT),
                   t='h', lwd=BAR_LWD, lend=1, col=COLZ[i])
        }
    }
    points(optimal_samples[[1]]$Mw + (0-BAR_H_OFFSET)/BAR_GROUP,
           pmax(const_samples, MIN_BAR_HT),
           t='h', lwd=BAR_LWD, lend=1, col='black')
    grid(col='orange')

    # Report the variance reduction that would be obtained by optimising the
    # number of samples for the specific threshold (independent of the
    # mean_optimal argument)
    legend_VR_local_optimised = unlist(lapply(optimal_samples,
        function(x) get_variances(x, mean_optimal)$constant_on_best))
    legend('topleft', 
           c('Equal in all bins    (VR = 1.00)', 
              paste0('Threshold = ', threshold_stages, 
              ' m   (VR = ', round(legend_VR_local_optimised ,2),')')), 
           fill=c('black', COLZ))

    # Plot with importance sampling
    for(i in 1:length(threshold_stages)){
        # Stagger the x location of the bars, so we fit several bars close to the
        # desired Mw value

        if(i == 1){

            plot(optimal_samples_IS[[i]]$Mw + (i-2.5)/BAR_GROUP,
                 pmax(optimal_samples_IS[[i]]$Nsamples, MIN_BAR_HT),
                 t='h', lwd=BAR_LWD, lend=1, col=COLZ[i], ylim=PLOT_YLIM,
                 xlab='', ylab='Optimal # Samples', las=1,
                 main='Stratified/importance-sampling', cex.main=1.8, cex.lab=1.5)

            mtext("Mw", side=1, line=2, cex=1.5)
            axis(side=1, at=seq(7.2, 9.6, by=0.1), labels=FALSE)

        }else{

            points(optimal_samples_IS[[i]]$Mw + (i-2.5)/BAR_GROUP,
                   pmax(optimal_samples_IS[[i]]$Nsamples, MIN_BAR_HT),
                   t='h', lwd=BAR_LWD, lend=1, col=COLZ[i])
        }
    }
    grid(col='orange')
    points(optimal_samples_IS[[1]]$Mw + (0-2.5)/BAR_GROUP, pmax(const_samples, MIN_BAR_HT),
           t='h', lwd=BAR_LWD, lend=1, col='black')

    # Report the variance reduction that would be obtained by optimising the
    # number of samples for the specific threshold (independent of the
    # mean_optimal_IS)
    legend_VR_local_optimised = unlist(lapply(optimal_samples_IS,
        function(x) get_variances(x, mean_optimal_IS)$constant_on_best))
    legend('topleft', 
           c('Equal in all bins    (VR = 1.00)', 
              paste0('Threshold = ', threshold_stages, 
                     ' m   (VR = ', round(legend_VR_local_optimised,2),')')), 
           fill=c('black', COLZ))
    dev.off()


    #
    # Plot the chosen number of samples, and the variance reductions.
    #
    png('Chosen_sampling_effort.png', width=9, height=6, units='in', res=300)
    BAR_EPS = 0.012
    plot(optimal_samples[[1]]$Mw - BAR_EPS, mean_optimal_including_constant,
         ylim=c(0, 150), xlim=c(7.15, 9.65),
         t='h', lwd=BAR_LWD*1.3, lend=1, col='black',
         xlab='Mw', ylab = '# Samples', cex.lab=1.4, cex.axis=1.3)

    abline(h=const_samples[1], col='purple', lty='dashed')
    points(optimal_samples[[1]]$Mw + BAR_EPS, mean_optimal_IS_including_constant, ylim=c(0, 150),
        t='h', lwd=BAR_LWD*1.3, lend=1, col='skyblue4')

    title(main='Selected non-uniform sampling effort and extra variance-reduction', cex.main=1.5)

    legend_VR_chosen_IS = unlist(lapply(optimal_samples_IS,
        function(x) get_variances(x, mean_optimal_IS_including_constant)$constant_on_chosen))
    legend_VR_chosen = unlist(lapply(optimal_samples,
        function(x) get_variances(x, mean_optimal_including_constant)$constant_on_chosen))

    white_t = rgb(1,1,1,alpha=0.0)
    legend('topleft', 
           paste0('Threshold = ', rep(threshold_stages, 1), 
                  ' (VR = ', format(round(legend_VR_chosen, 2)), ')'),
           title = 'Stratified sampling', 
           box.col=white_t, cex=1.1, bg=white_t)
    legend('top', 
           paste0('Threshold = ', rep(threshold_stages, 1), 
                  ' (VR = ', format(round(legend_VR_chosen_IS, 2)), ')'),
           title = 'Stratified/importance-sampling', 
           text.col='skyblue4',
           title.col='skyblue4',
           box.col=white_t, cex=1.1, bg=white_t)

    text(7.5, 53, 'Uniform sampling', col='purple', cex=1.4)
    dev.off()

    return(invisible(0))
}
# Make the plots
plot_sampling_effort(
    optimal_samples, mean_optimal, mean_optimal_including_constant,
    optimal_samples_IS, mean_optimal_IS, mean_optimal_IS_including_constant,
    const_samples, threshold_stages)

#
# Study how this performs at a bunch of other points
#

test_point_list = list(
    'target_point' = target_point,
    test_point_1 = c(184.8292, -20.9286),
    test_point_2 = c(184.6758, -21.3586),
    test_point_3 = c(184.47748, -21.041364),
    #test_point_4 = c(185.05833, -21.44512), ## This is very shallow, should avoid in offshore PTHA, although the sampling works OK
    test_point_NZ =  c(178.3945, -37.3940)
    )

event_peak_stage_list = list()
optimal_samples_list = list()
optimal_samples_IS_list = list()

for(nm_i in names(test_point_list)){
    #
    # Generalisation to nearby sites?
    event_peak_stage_list[[nm_i]] = get_peak_stage_at_target_point(test_point_list[[nm_i]])

    optimal_samples_list[[nm_i]] = lapply(threshold_stages, function(x){
        ptha18$get_optimal_number_of_samples_per_Mw(
            event_Mw,
            event_rates,
            event_peak_stage_list[[nm_i]],
            stage_threshold=x,
            total_samples=TOTAL_SAMPLES) # No importance sampling
        })

    optimal_samples_IS_list[[nm_i]] = lapply(threshold_stages, function(x){
        ptha18$get_optimal_number_of_samples_per_Mw(
            event_Mw,
            event_rates,
            event_peak_stage_list[[nm_i]],
            stage_threshold=x,
            total_samples=TOTAL_SAMPLES,
            # With importance sampling based on the previous stage
            event_importance_weighted_sampling_probs = (event_peak_stage_target_point*event_rates))
        })

    print(paste0('Variances (no importance sampling) with different scenario counts at ', nm_i))
    print_variances(optimal_samples_list[[nm_i]], mean_optimal, mean_optimal_including_constant)

    print(paste0('Variances (IMPORTANCE SAMPLING) with different scenario counts at ', nm_i))
    print_variances(optimal_samples_IS_list[[nm_i]], mean_optimal_IS, mean_optimal_IS_including_constant)

}

#
# Compare analytical confidence intervals at multiple sites -- 'stratified sampling with uniform-sampling
# effort in each bin, vs stratified/importance sampling with non-uniform sampling in each bin'
#
plot_curve_comparison_multiple_sites<-function(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    event_rates,
    event_importance_weighted_sampling_probs,
    nonuniform_sampling_fun,
    uniform_sampling_fun){

    par(mfrow=c(2,3))
    par(oma=c(0, 1, 0, 0))
    par(mar=c(4,6,2,1))
    options(scipen=5)
    # Target_point
    all_titles = list()
    all_titles$target_point = 'Site P'
    all_titles$test_point_1 = 'Nearby site 1 '
    all_titles$test_point_2 = 'Nearby site 2 '
    all_titles$test_point_3 = 'Nearby site 3 '
    all_titles$test_point_NZ = 'Distant site NZ'

    # Plotting order
    names_plot_order = c('target_point', 'test_point_1', 'test_point_NZ', 'test_point_2', 'test_point_3')
    # Check that the values of names_plot_order are all in event_peak_stage_list
    stopifnot(all(is.finite(match(names_plot_order, names(event_peak_stage_list)))))

    # Geometrically spaced peak-stage values at which we evaluate the curve
    peak_stage_vals = 10^seq(-2, 1.3, by=0.025)

    for(nm_i in names_plot_order){

        # Construct 95% "analytical" confidence intervals for the Monte-Carlo exceedance-rates
        lower_CI_stratified = rep(NA, length(peak_stage_vals))
        upper_CI_stratified = rep(NA, length(peak_stage_vals))
        lower_CI_stratified_importance = rep(NA, length(peak_stage_vals))
        upper_CI_stratified_importance = rep(NA, length(peak_stage_vals))
        variance_stratified = rep(NA, length(peak_stage_vals))
        variance_stratified_importance = rep(NA, length(peak_stage_vals))

        exrate_mean = rep(NA, length(peak_stage_vals))
        event_peak_stage_local = event_peak_stage_list[[nm_i]]

        for(i in 1:length(peak_stage_vals)){

            # stratified/importance sampling, non-uniform samples in each bin
            tmp = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
                event_Mw, event_rates, event_peak_stage_local,
                stage_threshold=peak_stage_vals[i],
                samples_per_Mw = nonuniform_sampling_fun,
                # Importance sampling based on target-point
                event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)
            lower_CI_stratified_importance[i] = tmp[1] + qnorm(0.025)*sqrt(tmp[2])
            upper_CI_stratified_importance[i] = tmp[1] + qnorm(0.975)*sqrt(tmp[2])
            exrate_mean[i] = tmp[1]
            variance_stratified_importance[i] = tmp[2]


            # stratified sampling, uniform samples in each bin
            tmp = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
                event_Mw, event_rates, event_peak_stage_local,
                stage_threshold=peak_stage_vals[i],
                samples_per_Mw = uniform_sampling_fun,
                event_importance_weighted_sampling_probs=(event_rates))
            lower_CI_stratified[i] = tmp[1] + qnorm(0.025)*sqrt(tmp[2])
            upper_CI_stratified[i] = tmp[1] + qnorm(0.975)*sqrt(tmp[2])
            variance_stratified[i] = tmp[2]
        }

        YLIM = c(1.0e-04, 1.0e-02)
        XLIM = c( approx(lower_CI_stratified, peak_stage_vals, xout=YLIM[2])$y,
                  approx(upper_CI_stratified, peak_stage_vals, xout=YLIM[1])$y)
        # Fixes for corner cases when there are pretty much no waves
        if(!is.finite(XLIM[1])) XLIM[1] = 1.0e-02
        if(XLIM[2] <= XLIM[1]) XLIM[2] = XLIM[1] + 1.0

        # Add curves to the plot
        plot(peak_stage_vals, exrate_mean, t='l', log='xy', xlim=XLIM, ylim=YLIM,
             xlab="", ylab="",
             cex.lab=1.4, cex.axis=1.3, las=1)
        mtext(side=1, 'Tsunami maxima (m)', line=2.2, cex=1.1)
        mtext(side=2, 'Exceedance-rate', line=5, cex=1.2)
        points(peak_stage_vals, lower_CI_stratified, t='l', col='blue', lty='dashed')
        points(peak_stage_vals, upper_CI_stratified, t='l', col='blue', lty='dashed')
        points(peak_stage_vals, lower_CI_stratified_importance, t='l', col='red', lty='twodash')
        points(peak_stage_vals, upper_CI_stratified_importance, t='l', col='red', lty='twodash')
        #abline(v=c(2,4), lwd=2, col='purple')
        title(all_titles[[nm_i]], cex.main=1.8)
        add_log_axis_ticks(side=1)
        add_log_axis_ticks(side=2)
        abline(h=c(1.0e-04, 1.0e-03, 1.0e-02), lty='dotted', col='orange')
        abline(v=c(0.05, 0.1, 0.5, 1, 5, 10), lty='dotted', col='orange')

        # Estimate variance reductions (interpolate with log-transform)
        variance_reduction_at_2m = (exp(approx(peak_stage_vals, log(variance_stratified), xout=2)$y)/
                                    exp(approx(peak_stage_vals, log(variance_stratified_importance), xout=2)$y))
        variance_reduction_at_4m = (exp(approx(peak_stage_vals, log(variance_stratified), xout=4)$y)/
                                    exp(approx(peak_stage_vals, log(variance_stratified_importance), xout=4)$y))
        text(XLIM[1]*2.5, 3e-04, paste0('VR2 = ', signif(variance_reduction_at_2m, 3)), cex=1.5)
        text(XLIM[1]*2.5, 1.33e-04, paste0('VR4 = ', signif(variance_reduction_at_4m, 3)), cex=1.5)
    }

    # Add legend in the final panel
    par(mar=c(0,0,0,0))
    plot(c(0,1), c(0,1), ann=FALSE, col='white', axes=FALSE)
    legend(0.00, 1.0, c('95% interval (analytical) \nStratified with uniform\nMw-bin sampling'),
           lty='dashed', col='blue', lwd=2, cex=1.4, bty='n')
    legend(0.00, 0.7, c('95% interval (analytical) \nStratified/importance with\nnon-uniform Mw-bin sampling'),
           lty='twodash', col='red', lwd=2, cex=1.4, bty='n')
    text(0.5, 0.25, bquote(paste('VR2 = Variance-reduction @ ', Q^T, '=2m')), cex=1.4)
    text(0.5, 0.15, bquote(paste('VR4 = Variance-reduction @ ', Q^T, '=4m')), cex=1.4) 

}
# Compare the scheme used herein with 'stratified sampling and uniform N(M_w)', using the logic-tree mean
# kermadectonga2 Mw-frequency curve
png('Curve_comparison_multiple_sites.png', width=9, height=4.5, units='in', res=300)
plot_curve_comparison_multiple_sites(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    event_rates,
    event_importance_weighted_sampling_probs=(event_peak_stage_target_point*event_rates),
    nonuniform_sampling_fun = approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'),
    uniform_sampling_fun = function(Mw){ const_samples[1]*(Mw < 9.65) }
    )
dev.off()

#
# As above, but do the calculation for the unsegmented/segmented logic-tree mean, and assume that stratified
# sampling can spend 50% of the scenarios on the unsegmented branch. Later we also look at the segments, assuming
# we can spend 30% on Tonga, 20% on Kermadec, and 10% on Hikurangi. In a real application the latter numbers
# should add to 50%, but these plots are just exploratory - we never use the combined results - so there is no
# problem.
#

ptha18_source_rate_env = new.env()
#source('../../../../../AustPTHA/CODE/ptha/ptha_access/get_detailed_PTHA18_source_zone_info.R',
source('../../../ptha_access/get_detailed_PTHA18_source_zone_info.R',
       local=ptha18_source_rate_env, chdir=TRUE)

# Unsegmented case -- stratified sampling gets 50% of the total scenarios (other 50% on segments)
unsegmented_KT2 = ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
    source_zone='kermadectonga2', segment='')
png('Curve_comparison_multiple_sites_UNSEGMENTED.png', width=9, height=4.5, units='in', res=300)
plot_curve_comparison_multiple_sites(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    # Assume the true exceedance-rates correspond to the UNSEGMENTED logic-tree mean
    event_rates=unsegmented_KT2$HS_event_rates,
    event_importance_weighted_sampling_probs=(event_peak_stage_target_point*event_rates),
    nonuniform_sampling_fun = approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'),
    # Assume for stratified sampling that only half of the samples are spent on the unsegmented model.
    # The other half would have to be used to sample the segments
    uniform_sampling_fun = function(Mw){ 0.5*const_samples[1]*(Mw < 9.65) }
    )
dev.off()
# Check it works OK via the other plot (which shows MC samples)
plot_hazard_curve('stratified_importance', event_rates=unsegmented_KT2$HS_event_rates,
    Number_MC_reps=10000,
    event_importance_weighted_sampling_probs=event_rates*event_peak_stage_target_point,
    fig_title='target_point_stratified_importance_unequal_UNSEGMENTED.png',
    mw_sampling_fun=approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'))

# As above, but do the calculation for the Tonga segment logic-tree mean, and assume that stratified
# sampling can spend 30% of the scenarios on the Tonga segment (optimistic)

tonga_segmented_KT2 = ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
    source_zone='kermadectonga2', segment='tonga')
png('Curve_comparison_multiple_sites_TONGA_SEGMENT.png', width=9, height=4.5, units='in', res=300)
plot_curve_comparison_multiple_sites(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    # Assume the true exceedance-rates correspond to the tonga-segment logic-tree mean
    event_rates=tonga_segmented_KT2$HS_event_rates,
    event_importance_weighted_sampling_probs=(event_peak_stage_target_point*event_rates),
    nonuniform_sampling_fun = approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'),
    # Assume for stratified sampling that only half of the samples are spent on the unsegmented model.
    # The other half would have to be used to sample the segments
    uniform_sampling_fun = function(Mw){ 0.3*const_samples[1]*(Mw < 9.65) }
    )
dev.off()
# Check it works OK via the other plot (which shows MC samples)
plot_hazard_curve('stratified_importance', event_rates=tonga_segmented_KT2$HS_event_rates,
    Number_MC_reps=10000,
    event_importance_weighted_sampling_probs=event_rates*event_peak_stage_target_point,
    fig_title='target_point_stratified_importance_unequal_TONGA_SEGMENT.png',
    mw_sampling_fun=approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'))

# As above, but do the calculation for the Kermadec segment logic-tree mean, and assume that stratified
# sampling can spend 20% of the scenarios on the Kermadec segment (optimistic)

kermadec_segmented_KT2 = ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
    source_zone='kermadectonga2', segment='kermadec')
png('Curve_comparison_multiple_sites_KERMADEC_SEGMENT.png', width=9, height=4.5, units='in', res=300)
plot_curve_comparison_multiple_sites(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    # Assume the true exceedance-rates correspond to the tonga-segment logic-tree mean
    event_rates=kermadec_segmented_KT2$HS_event_rates,
    event_importance_weighted_sampling_probs=(event_peak_stage_target_point*event_rates),
    nonuniform_sampling_fun = approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'),
    # Assume for stratified sampling that only half of the samples are spent on the unsegmented model.
    # The other half would have to be used to sample the segments
    uniform_sampling_fun = function(Mw){ 0.2*const_samples[1]*(Mw < 9.65) }
    )
dev.off()
# Check it works OK via the other plot (which shows MC samples)
plot_hazard_curve('stratified_importance', event_rates=kermadec_segmented_KT2$HS_event_rates,
    Number_MC_reps=10000,
    event_importance_weighted_sampling_probs=event_rates*event_peak_stage_target_point,
    fig_title='target_point_stratified_importance_unequal_KERMADEC_SEGMENT.png',
    mw_sampling_fun=approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'))


# As above, but do the calculation for the Hikurangi segment logic-tree mean, and assume that stratified
# sampling can spend 10% of the scenarios on the Hikurangi segment (very optimistic - impossible if we
# just spent 50% on unsegmented and 30% on Tonga and 20% on Kermadec !!)

hikurangi_segmented_KT2 = ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
    source_zone='kermadectonga2', segment='hikurangi')
png('Curve_comparison_multiple_sites_HIKURANGI_SEGMENT.png', width=9, height=4.5, units='in', res=300)
plot_curve_comparison_multiple_sites(
    event_peak_stage_list,
    unique_mw,
    mean_optimal_IS_including_constant,
    event_Mw,
    # Assume the true exceedance-rates correspond to the tonga-segment logic-tree mean
    event_rates=hikurangi_segmented_KT2$HS_event_rates,
    event_importance_weighted_sampling_probs=(event_peak_stage_target_point*event_rates),
    nonuniform_sampling_fun = approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'),
    # Assume for stratified sampling that only half of the samples are spent on the unsegmented model.
    # The other half would have to be used to sample the segments
    uniform_sampling_fun = function(Mw){ 0.1*const_samples[1]*(Mw < 9.65) }
    )
dev.off()
# Check it works OK via the other plot (which shows MC samples)
plot_hazard_curve('stratified_importance', event_rates=hikurangi_segmented_KT2$HS_event_rates,
    Number_MC_reps=10000,
    event_importance_weighted_sampling_probs=event_rates*event_peak_stage_target_point,
    fig_title='target_point_stratified_importance_unequal_HIKURANGI_SEGMENT.png',
    mw_sampling_fun=approxfun(unique_mw, mean_optimal_IS_including_constant, method='constant'))


#
# Extra hazard curve plots
#
print('#')
print('# EXTRA HAZARD CURVE PLOTS AT OTHER POINTS BELOW HERE')
print('#')

for(nm_i in names(event_peak_stage_list)){
    print(paste0('Plotting ', nm_i))
    plot_hazard_curve('stratified', Number_MC_reps=DEFAULT_MC_REPS,
        event_peak_stage_local=event_peak_stage_list[[nm_i]],
        fig_title=paste0(nm_i, '_stratified.png'))
    plot_hazard_curve('stratified_importance', Number_MC_reps=DEFAULT_MC_REPS,
        event_peak_stage_local=event_peak_stage_list[[nm_i]],
        fig_title=paste0(nm_i, '_stratified_importance.png'))
}

print('#')
print('# EXTRA HAZARD CURVE PLOTS WITH UNEQUAL SAMPLING BELOW HERE')
print('#')

mw_sampling_fun_stratified = approxfun(unique_mw, mean_optimal_including_constant)
mw_sampling_fun_stratified_importance = approxfun(unique_mw, mean_optimal_IS_including_constant)

for(nm_i in names(event_peak_stage_list)){
    print(paste0('Plotting ', nm_i))
    plot_hazard_curve('stratified', Number_MC_reps=DEFAULT_MC_REPS,
        event_peak_stage_local=event_peak_stage_list[[nm_i]],
        mw_sampling_fun=mw_sampling_fun_stratified,
        fig_title=paste0(nm_i, '_stratified_unequal.png'))
    plot_hazard_curve('stratified_importance', Number_MC_reps=DEFAULT_MC_REPS,
        event_peak_stage_local=event_peak_stage_list[[nm_i]],
        mw_sampling_fun=mw_sampling_fun_stratified_importance,
        fig_title=paste0(nm_i, '_stratified_importance_unequal.png'))
}


#
# Add in a plot comparing MC variability at multiple sites
#
compute_MC_exceedances_multi_sites<-function(
    sampling_type='stratified',
    Number_MC_reps=DEFAULT_MC_REPS,
    event_importance_weighted_sampling_probs = NULL,
    mw_sampling_fun = NULL){

    if(sampling_type == 'stratified'){

        # No importance sampling
        if(!is.null(event_importance_weighted_sampling_probs)){
            stop('Do not specify event_importance_weighted_sampling_probs with stratified sampling')
        }
        # Weights equivalent to regular stratified sampling
        event_importance_weighted_sampling_probs = event_rates

    }else if(sampling_type == 'stratified_importance'){

        # By default use stage*rate at the target point
        if(is.null(event_importance_weighted_sampling_probs)){
            event_importance_weighted_sampling_probs = event_rates*event_peak_stage_target_point
        }

    }else{
        stop('unknown sampling_type')
    }

    if(is.null(mw_sampling_fun)){
        # How many samples in each Mw bin? By default use a constant for 'possible' scenarios.
        mw_sampling_fun = function(Mw){ const_samples[1]*(Mw < 9.65) }
    }

    target_stage = 2

    set.seed(1234) # Reproducible randomness

    empty_df = list(mean=rep(NA, Number_MC_reps), var=rep(NA, Number_MC_reps), analytical_mean_var=rep(NA,2))
    exrate_ts_store = list(
        target_point=empty_df,
        test_point_1=empty_df,
        test_point_2=empty_df,
        test_point_3=empty_df,
        test_point_NZ=empty_df)

    for(i in 1:Number_MC_reps){

        random_scenarios = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
            event_rates,
            event_Mw,
            samples_per_Mw = mw_sampling_fun,
            event_importance_weighted_sampling_probs = event_importance_weighted_sampling_probs)

        for(nm_i in names(event_peak_stage_list)){
            # For both kinds of sampling this will compute the correct mean/variance
            mean_and_variance = ptha18$estimate_exrate_uncertainty(
                random_scenarios, event_peak_stage_list[[nm_i]], target_stage)
            exrate_ts_store[[nm_i]]$mean[i] = mean_and_variance[1]
            exrate_ts_store[[nm_i]]$var[i] = mean_and_variance[2]

            if(i == 1){
                # Store the analytical mean/variance as well, for cross-checking. Only need to do this once
                exrate_ts_store[[nm_i]]$analytical_mean_var = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
                    event_Mw,
                    event_rates,
                    event_peak_stage_list[[nm_i]],
                    stage_threshold=target_stage,
                    samples_per_Mw=mw_sampling_fun,
                    event_importance_weighted_sampling_probs=event_importance_weighted_sampling_probs)
            }
        }

    }

    return(exrate_ts_store)
}

# Exceedance-rate variability, no importance sampling
exrate_ts_store_stratified = compute_MC_exceedances_multi_sites(sampling_type='stratified')
# ..........................., with importance based on site P
exrate_ts_store_stratified_importance = compute_MC_exceedances_multi_sites(sampling_type='stratified_importance')
# ..........................., with importance based on each sites local stage
exrate_ts_store_alternativeimportance = list()
for(nm_i in names(event_peak_stage_list)){
    exrate_ts_store_alternativeimportance[[nm_i]] = compute_MC_exceedances_multi_sites(
        sampling_type='stratified_importance',
        event_importance_weighted_sampling_probs=event_peak_stage_list[[nm_i]]*event_rates)
}

# Store the 'ideal' variance reduction -- i.e. the variance reduction if we used the 'local'
# tsunami maxima to represent the event importance, as compared with stratified sampling.
ideal_VR  = vector(mode='list', length=length(exrate_ts_store_stratified))
names(ideal_VR) = names(exrate_ts_store_stratified)
ideal_VR_analytical = ideal_VR
for(nm_i in names(event_peak_stage_list)){
    ideal_VR[[nm_i]] = var(exrate_ts_store_stratified[[nm_i]]$mean)/
        var(exrate_ts_store_alternativeimportance[[nm_i]][[nm_i]]$mean)
    ideal_VR_analytical[[nm_i]] = (
        exrate_ts_store_stratified[[nm_i]]$analytical_mean_var[2]/
        exrate_ts_store_alternativeimportance[[nm_i]][[nm_i]]$analytical_mean_var[2])
}


png('Density_comparison_multiple_sites.png', width=9, height=5, units='in', res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,2,1))
# Target_point
all_titles = vector(mode='list', length=length(exrate_ts_store_stratified))
names(all_titles) = names(exrate_ts_store_stratified)
all_titles$target_point = 'Site P'
all_titles$test_point_1 = 'Nearby site 1 '
all_titles$test_point_2 = 'Nearby site 2 '
all_titles$test_point_3 = 'Nearby site 3 '
all_titles$test_point_NZ = 'Distant site NZ'
names_plot_order = c('target_point', 'test_point_1', 'test_point_NZ', 'test_point_2', 'test_point_3')
store_VR = list()
for(nm_i in names_plot_order){

    ds_stratified = density(exrate_ts_store_stratified[[nm_i]]$mean)
    ds_stratified_importance = density(exrate_ts_store_stratified_importance[[nm_i]]$mean)

    var_stratified = var(exrate_ts_store_stratified[[nm_i]]$mean)
    var_stratified_importance = var(exrate_ts_store_stratified_importance[[nm_i]]$mean)

    # We could either compute the Variance ratio from the samples, or with the analytical approach.
    # They should be similar with enough MC samples, but perhaps not identical
    VR_from_sampling = var_stratified/var_stratified_importance
    VR_analytical = exrate_ts_store_stratified[[nm_i]]$analytical_mean_var[2]/
         exrate_ts_store_stratified_importance[[nm_i]]$analytical_mean_var[2]

    YLIM = max(c(max(ds_stratified_importance$y), max(ds_stratified$y)))
    plot(ds_stratified$x, ds_stratified$y, t='l', col='black', ylim=c(0, YLIM),
        xlab='', ylab='Density', cex.lab=1.4, cex.axis=1.4,
        main=all_titles[[nm_i]], cex.main=1.8, lwd=2)
    points(ds_stratified_importance, col='red', t='l', lty='dashed', lwd=2)
    grid(col='lightblue', lty='dotted')
    mtext(side=1, bquote(paste('Monte-Carlo exceedance-rate (', Q^T, '=2m)')), cex=0.95, line=3)
    text(max(ds_stratified$x)*0.97, YLIM*0.7,
         paste0('VR = ', signif(VR_analytical, 3), '\n',
                '     ', '(', signif(ideal_VR_analytical[[nm_i]], 3), ')'),
        adj=c(1,0), cex=2)

    # Useful cross-check
    store_VR[[nm_i]] = c(VR_analytical, VR_from_sampling)
}
par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1), ann=FALSE, col='white', axes=FALSE)
legend('left', c('Stratified', 'Stratified/Importance'),
       lty=c('solid', 'dashed'), col=c('black', 'red'), lwd=c(2, 2), cex=1.7, bty='n')
text(0.5, 0.3, 'VR = Variance-reduction factor', cex=1.5)
dev.off()

# Check that the analytical and sampling-based variance ratios agree well
vr_relative_err = unlist(lapply(store_VR, function(x) abs(1-x[1]/x[2])))
stopifnot(all(vr_relative_err < 0.031))


sink()

