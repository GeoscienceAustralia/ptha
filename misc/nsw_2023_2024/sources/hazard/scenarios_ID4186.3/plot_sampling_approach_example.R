# Configuration for these samples
sc = new.env()
source('sampling_config.R', local=sc)

# For getting PTHA18 results
ptha18 = new.env()
source(sc$get_ptha_results_script_file, local=ptha18, chdir=TRUE)

# For getting detailed PTHA18 results
#ptha18_detailed = new.env()
#source(sc$get_detailed_ptha18_source_zone_info_script_file, local=ptha18_detailed, chdir=TRUE)

# Importance sampling calculations
isu = new.env()
source('../importance_sampling_utilities.R', local=isu)

# Names of the source zones
source_zones = names(sc$source_zones_and_samples)

# Read existing data
source_zone_data = readRDS('source_zone_data_backup.RDS')

OUTDIR = './FIG_for_paper/'
dir.create(OUTDIR, showWarnings=FALSE)

plot_and_print_statistics<-function(sz, outdir=OUTDIR){
    #
    # Plots illustrating the sampling for a source zone
    #

    print(c('SZ: ', sz))

    MAX_STAGE_LIMIT = 5

    sz_bin_boundaries = sc$importance_function(
        peak_stage_target_pt = source_zone_data[[sz]]$event_peak_stage_at_targetpt[[sz]]$max_stage, 
        scenario_rate_logictreemean = source_zone_data[[sz]]$event_rate_logic_tree_mean,
        stage_vs_rate_all_source_zones_at_target_pt = sc$stage_vs_rate_all_source_zones_at_target_pt, 
        return_importance_function_bin_boundaries=TRUE)

    print(c(' length(sz_bin_boundaries) = ', length(sz_bin_boundaries)))

    S84_stage = sc$stage_vs_rate_all_source_zones_at_target_pt$stage
    S84_exrate = sc$stage_vs_rate_all_source_zones_at_target_pt$stochastic_slip_rate_84pc

    #
    # First panel -- illustrate the bins used to create I(e)
    #
    png(paste0(outdir, 'sampling_plot_', sz, '_bins.png'), width=9, height=6, units='in', res=300)
    options(scipen=5)
    par(mar=c(4.5, 7.6, 1.8, 1.8))
    plot(S84_stage, S84_exrate, log='xy', xlim=c(0.02, MAX_STAGE_LIMIT), ylim=c(1e-05, 1), t='l', col='orange', lwd=2,
        las=1, xlab="", ylab="", cex.axis=1.5, main=paste0('A) Bins used to define I(e), ', sz), cex.main=2.)
    mtext("PTHA18 maximum stage (m) @ reference site", side=1, cex=1.6, line=2.5)
    mtext('Exceedance rate (events/year)', side=2, cex=1.6, line=6)
    ptfit = approx(S84_exrate, S84_stage, xout=c(1/250, 1/2500, 1/10000))
    abline(v=sz_bin_boundaries, col='brown')
    abline(h=ptfit$x, lty='dotted')
    legend('topright', 
        c('PTHA18 84th percentile max stage\nexceedance rate (all source zones)', paste0('I(e) bin boundaries (', sz, ')')),
        lty=c('solid', 'solid'), col=c('orange', 'brown'), lwd=2, cex=1.4)

    xpt = 0.5*(sz_bin_boundaries[-1] + sz_bin_boundaries[-length(sz_bin_boundaries)])
    ypt = c(0.9, 0.4, 0.15)
    #text(xpt, ypt, paste0('B', seq(1, length(xpt))), adj=0.5, srt=90)
    text(xpt[1], 1/250, '1/250', pos=3, offset=0.25, cex=1.3)
    text(xpt[1], 1/2500, '1/2500', pos=3, offset=0.25, cex=1.3)
    text(xpt[1], 1/10000, '1/10000', pos=3, offset=0.25, cex=1.3)

    points(ptfit$y, ptfit$x)
    add_log_axis_ticks(side=1)
    add_log_axis_ticks(side=2)
    dev.off()

    #
    # Second panel: Source zone specific exrate curve (real, and Monte Carlo estimate)
    #
    for(plot_site in c('target_pt','tweed_offshore')){

        print(c(' Site = ', plot_site))

        if(plot_site == 'tweed_offshore'){
            MAX_STAGE_LIMIT = 2.0
        }else{
            MAX_STAGE_LIMIT = 5.0
        }

        stage_thresholds = 10**seq(log10(0.02), log10(MAX_STAGE_LIMIT), len=200)
        event_rate_logic_tree_mean = source_zone_data[[sz]]$event_rate_logic_tree_mean
        sampling_prob = source_zone_data[[sz]]$sampling_prob
        N_MC = source_zone_data[[sz]]$N_MC

        if(plot_site == 'target_pt'){
            peak_stage_plot_site = source_zone_data[[sz]]$event_peak_stage_at_targetpt[[sz]]$max_stage
        }else if(plot_site == 'tweed_offshore'){
            peak_stage_plot_site = source_zone_data[[sz]]$peak_stage_other_pts$Tweed_offshore[[sz]]$max_stage
        }else{
            stop('plot_site not supported')
        }

        analytical_mc_results = lapply(stage_thresholds, function(x){
            isu$analytical_exrate_and_MC_variance(
                event_rate=event_rate_logic_tree_mean, 
                event_stage=peak_stage_plot_site, 
                stage_threshold=x, 
                sampling_prob=sampling_prob, 
                N_MC=N_MC)})
        analytical_exrate = unlist(lapply(analytical_mc_results, function(x) x[1]))     
        mc_exrate_variance = unlist(lapply(analytical_mc_results, function(x) x[2]))

        png(paste0(outdir, 'sampling_plot_', sz, '_LTM_', plot_site, '.png'), width=9, height=6, units='in', res=300)
        options(scipen=5)
        par(mar=c(4.5, 7.6, 1.8, 1.8))
        plot(stage_thresholds, analytical_exrate, log='xy', xlim=c(0.02, MAX_STAGE_LIMIT), ylim=c(1e-05, 1), t='l', col='purple', lwd=2,
            las=1, xlab="", ylab="", cex.axis=1.5, main=paste0('B) ', sz, ' Monte Carlo exceedance rate estimate'), cex.main=1.7)
        if(plot_site == 'target_pt'){
            # The bins only really have meaning at the target pt
            abline(v=sz_bin_boundaries, col='brown')
            mtext("PTHA18 maximum stage (m) @ reference site", side=1, cex=1.6, line=2.5)
        }else if(plot_site == 'tweed_offshore'){
            mtext("PTHA18 maximum stage (m) @ northern NSW", side=1, cex=1.6, line=2.5)
        }
        abline(h=ptfit$x, lty='dotted')

        # Make test samples
        Ntest = 1000
        test_samples = vector(mode='list', length=Ntest)
        for(i in 1:Ntest){
            test_samples[[i]] = isu$sample_scenarios_monte_carlo(
                    event_rate=event_rate_logic_tree_mean, 
                    event_stage=peak_stage_plot_site, 
                    stage_thresholds=stage_thresholds, 
                    sampling_prob=sampling_prob, 
                    N_MC=N_MC)
            transparent_grey = rgb(0.5, 0.5, 0.5, alpha=0.2, maxColorValue=1)
            points(stage_thresholds, test_samples[[i]]$exrates, t='l', col=transparent_grey)
        }
        # Analytical confidence intervals
        points(stage_thresholds, analytical_exrate, col='black', lwd=2, t='l')
        points(stage_thresholds, analytical_exrate+qnorm(0.975)*sqrt(mc_exrate_variance), t='l', col='black', lwd=2, lty='dashed')
        points(stage_thresholds, analytical_exrate+qnorm(0.025)*sqrt(mc_exrate_variance), t='l', col='black', lwd=2, lty='dashed')
        legend('topright', c(paste0('PTHA18 exceedance rate \n(single source zone)'), 'Monte Carlo 95% CI', '1000 Monte Carlo samples'),
            lty=c('solid', 'dashed', 'solid'), lwd=c(2,2,1), col=c('black', 'black', 'grey'), cex=1.4, bg='white')
        mtext('Exceedance rate (events/year)', side=2, cex=1.6, line=6)
        dev.off()

        # Find a stage near this
        test_stage = 1.0
        nearby_ind = which.min(abs(stage_thresholds - test_stage))
        mc_samples_at_threshold = unlist(lapply(test_samples, function(x) x$exrates[nearby_ind]))
        analytical_exrate_var_at_threshold = mc_exrate_variance[nearby_ind]
        print(c(' At a stage threshold of ', stage_thresholds[nearby_ind]))
        print(c('  Mean(1000 MC samples)          : ', mean(mc_samples_at_threshold)))
        print(c('  sqrt(variance(1000 MC samples)): ', sqrt(var(mc_samples_at_threshold))))
        print(c('  Analytical mean                : ', analytical_exrate[nearby_ind]))
        print(c('  sqrt(analytical MC var       )): ', sqrt(analytical_exrate_var_at_threshold)))

        # Confidence interval coverage
        ci_lower = rep(NA, Ntest)
        ci_upper = rep(NA, Ntest)
        for(i in 1:Ntest){
            test_var = isu$estimate_exrate_and_variance_sampled_scenarios(
                rates_sampled_scenarios = event_rate_logic_tree_mean[test_samples[[i]]$si],
                sampling_prob_sampled_scenarios = sampling_prob[test_samples[[i]]$si],
                stage_sampled_scenarios = peak_stage_plot_site[test_samples[[i]]$si],
                desired_stage_thresholds = stage_thresholds[nearby_ind])
            ci_lower[i] = test_var$exrate[1] + qnorm(0.025)*sqrt(test_var$exrate_var[1])        
            ci_upper[i] = test_var$exrate[1] + qnorm(0.975)*sqrt(test_var$exrate_var[1])        
        }

        coverage_fraction = mean((ci_lower < analytical_exrate[nearby_ind]) & (ci_upper > analytical_exrate[nearby_ind]))
        print(c('  Sample 95% confidence interval coverage fraction: ', coverage_fraction))
    }

    return(invisible(0))
}

sz = 'kermadectonga2'
plot_and_print_statistics(sz)
sz = 'southamerica'
plot_and_print_statistics(sz)
sz = 'puysegur2'
plot_and_print_statistics(sz)
sz = 'newhebrides2'
plot_and_print_statistics(sz)
sz = 'solomon2'
plot_and_print_statistics(sz)
sz = 'outerrise_kermadectonga'
plot_and_print_statistics(sz)

