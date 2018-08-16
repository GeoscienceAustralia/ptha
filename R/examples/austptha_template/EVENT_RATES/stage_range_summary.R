#
# This script looks at the 'stage-range' at DART buoys for events in the
# 'corresponding family of model scenarios'. 
#
# It compares the distribution of model scenarios with the observations,
# using some simple statistical methods, which complement the graphical
# displays elsewhere.
# 
#
#
library(rptha)

variable_mu = TRUE

# Do not consider events with {peak_slip > peak_slip_limit_factor*mean-scaling-relation-slip}
# Note we use the 'reference Mw' (i.e. constant shear modulus) for doing this
source('peak_slip_limit_factor.R', local=TRUE)


#
# Get a vector with the Rdata files produced by running
# ../SOURCE_ZONES/sourcezone/TSUNAMI_EVENTS/plots/gauge_summary_statistics.R
# for each of the historical comparison events
#
if(variable_mu){
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*varyMu.Rdata')
}else{
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*[0-9].Rdata')
}

uniform_store = list()
stochastic_store = list()
variable_uniform_store = list()

# For each historical event, extract the stage range at all darts
for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    load(all_Rdata[i], envir=event_env)

    ngauges = length(event_env$uniform_slip_stats)

    uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges),
        reference_Mw = unlist(lapply(event_env$uniform_slip_stats[[1]], f<-function(x) x$events_with_Mw$Mw)),
        peak_slip = unlist(lapply(event_env$uniform_slip_stats[[1]], f<-function(x) x$peak_slip)),
        gauge_names = basename(names(event_env$uniform_slip_stats))
        )
    stochastic_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges),
        reference_Mw = unlist(lapply(event_env$stochastic_slip_stats[[1]], f<-function(x) x$events_with_Mw$Mw)),
        peak_slip = unlist(lapply(event_env$stochastic_slip_stats[[1]], f<-function(x) x$peak_slip)),
        gauge_names = basename(names(event_env$stochastic_slip_stats))
        )
    variable_uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges),
        reference_Mw = unlist(lapply(event_env$variable_uniform_slip_stats[[1]], f<-function(x) x$events_with_Mw$Mw)),
        peak_slip = unlist(lapply(event_env$variable_uniform_slip_stats[[1]], f<-function(x) x$peak_slip)),
        gauge_names = basename(names(event_env$variable_uniform_slip_stats))
        )

    # Loop over DART buoys and get the stage-range
    for(gauge_ind in 1:ngauges){

        # Uniform
        uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(uniform_store[[i]]$model[,gauge_ind]))

        # Stochastic
        stochastic_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        stochastic_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        stochastic_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(stochastic_store[[i]]$model[,gauge_ind]))

        # variable_uniform
        variable_uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        variable_uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        variable_uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(variable_uniform_store[[i]]$model[,gauge_ind]))
        
    }

    #
    # Remove events where peak slip is too high, if required
    # This is not adapted yet for events with mu != 30
    #
    print('Warning: Removing events with high slip assuming "reference mu" 30GPa and Strasser scaling relation')
    print('         This needs to be updated if we want to test normal faults')
    k = which(uniform_store[[i]]$peak_slip < (peak_slip_limit_factor * slip_from_Mw(uniform_store[[i]]$reference_Mw)))
    print(length(k))
    uniform_store[[i]]$data = uniform_store[[i]]$data[k,,drop=FALSE]
    uniform_store[[i]]$model = uniform_store[[i]]$model[k,,drop=FALSE]

    k = which(stochastic_store[[i]]$peak_slip < (peak_slip_limit_factor * slip_from_Mw(stochastic_store[[i]]$reference_Mw)))
    print(length(k))
    stochastic_store[[i]]$data = stochastic_store[[i]]$data[k,,drop=FALSE]
    stochastic_store[[i]]$model = stochastic_store[[i]]$model[k,,drop=FALSE]

    k = which(variable_uniform_store[[i]]$peak_slip < (peak_slip_limit_factor * slip_from_Mw(variable_uniform_store[[i]]$reference_Mw)))
    print(length(k))
    variable_uniform_store[[i]]$data = variable_uniform_store[[i]]$data[k,,drop=FALSE]
    variable_uniform_store[[i]]$model = variable_uniform_store[[i]]$model[k,,drop=FALSE]
}


if(variable_mu){
    save.image('model_data_envelope_summary_statistics_varyMu.Rdata')
}else{
    save.image('model_data_envelope_summary_statistics.Rdata')
}
#
# Info on number of unique models [since double-ups can occur, e.g.
# if 2 variable_uniform slip models have the same area and Mw. This
# even happens for single unit-source stochastic models
#
for(i in 1:length(all_Rdata)){
    print(paste0('################', i))
    print(basename(all_Rdata[i]))
    print(c(uniform_store[[i]]$nunique_models, '/', length(uniform_store[[i]]$model[,1])))
    print(c(stochastic_store[[i]]$nunique_models, '/', length(stochastic_store[[i]]$model[,1])))
    print(c(variable_uniform_store[[i]]$nunique_models, '/', length(variable_uniform_store[[i]]$model[,1])))
}


#
#
# Score each event, by comparing the 'gauge-averaged' non-exceedance rate of
# the observed stage range, to that predicted by the model
#
#
multi_gauge_summary_fun = median #mean 

#' Compute a statistic which 'scores' the data as high or low, based on the
#' distribution of events
#'
#' The idea is: For a single event, at each DART buoy, we can score the data
#' based on "the fraction of model stage ranges that are smaller". If we then
#' consider the 'median of this statistic over all DARTS', we have a 'score'
#' which indicates (to some extent) whether the data is 'high' or 'low' compared
#' with the models. Obviously much information is 'thrown away' using this statistic,
#' however, it has an advantage: We can compare this number with its sampling
#' distribution under the null hypothesis that a particular earthquake model
#' type is correct. If the NULL hypothesis is correct, then over many tsunami events,
#' the observation 'score' should look like random samples from the sampling
#' distribution. In other words, we can make a statistical test that considers whether
#' the model family is consistently high or low, over many events.
#' 
#' @param dta list with model and data results at all DARTS for a single events
#' @param fake_data_by_perturbing_random_model logical. If FALSE, get data from
#'   dta$data. If TRUE, suppose th data is a random model sample. The latter
#'   option is useful for testing if the statistic satisfies the null hypothesis
#'   on random data (i.e. which it should, but this permits "sanity testing" of the approach).
#' @param return_environment logical. If TRUE, return the function environment.
#'   Otherwise return a statistic giving a 'score' of the model compared with
#'   the data. This statistic is in (0,1)
#' @param weights optional vector of weights for each model. If present, it must
#'   have one entry for every column in dta$model. In that case, the model results
#'   are sampled randomly with replacement, with probability of being selected
#'   proportional to the weight. The sample size is some multiple of
#'   length(weights), controlled by sample_size_scale. This gives us a method to
#'   analyse the effect of data re-weighting on the results.
#' @param sample_size_scale integer. Only used if !is.na(weights). In that case,
#'   this controls the size of the random sample from the model. It should be
#'   large enough so that the statistic does not change in an important way
#'   between runs
#' @return a statistic between 0 and 1, with low values suggesting the
#'   observation is 'high' compared with the model family sampling distribution,
#'   and high values suggesting the observation is 'low' compared with the
#'   model family sampling distribution
score_gauge<-function(dta, fake_data_by_perturbing_random_model=FALSE, return_environment=FALSE, 
    weights=NA, sample_size_scale=10){

    # Use clear variable names!
    model_results = dta$model
    data_results = dta$data # Same size as model_results! Rows are repeated.

    # Optionally take a weighted sample from the data.
    # This will replace 'model_results' and 'data_results', so that the
    # remaining code does not have to change.
    if(!any(is.na(weights))){
        N = length(model_results[,1])
        stopifnot(length(weights) == N)
        model_sample_by_weight = sample.int(N, size=(N*sample_size_scale), replace=TRUE, prob=weights)

        new_model_results = model_results[model_sample_by_weight,]
        new_data_results = data_results[model_sample_by_weight,]

        # Weighted sample from the model and data results. Note that the
        # data_results are anyway effectively constant (occasional slight
        # variations due to time-window comparison shifting slightly for
        # models)
        model_results = new_model_results
        data_results = new_data_results
    }else{
        # Check for erronious input
        if(length(weights) > 1) stop('weights cannot contain NA values if it has length > 1')
    }

    # For each model run, get the rank of its stage-range at each gauge, and
    # then apply multi_gauge_summary_fun across all gauges. 
    nr = nrow(model_results)
    nc = ncol(model_results)

    ## APPROACH 1. 
    model_ranks = apply(model_results, 2, rank) # Rank each scenario at every dart
    model_summary = apply(model_ranks, 1, multi_gauge_summary_fun) # summary of rank over all gauges. (Each scenario) --> (a number)

    ## APPROACH 2 -- use a 'leave-one-out' approach to compute the ranks.
    ## Gives similar results. 
    #model_ranks = model_results * NA
    #for(i in 1:nr){
    #    # Compare (the model without row_i) to (row_i)
    #    for(j in 1:nc){
    #        model_ranks[i,j] = sum(model_results[-i,j] < model_results[i,j]) 
    #    }
    #}
    #model_summary = apply(model_ranks, 1, multi_gauge_summary_fun) # mean rank


    # For the data, get the rank of its stage-range at each gauge, and then apply
    # multi_gauge_summary_fun across all gauges
    if(!fake_data_by_perturbing_random_model){
        # Typical case
        data_mat = data_results
    }else{
        # Testing-only case.
        # Make up some fake data by perturbing the model data. This is only
        # useful for testing the performance of the method (e.g. that the results
        # suggest the null hypothesis is satisfied, when it really is!)
        #
        data_mat = data_results * NA
        fake_data_ind = sample(1:length(model_results[,1]), size=1)
        # Here is our random data -- a perturbation on a model result. Because
        # the statistic is rank based, the size of the perturbation is not an
        # issue, so long as it is 'small enough' to not push the random data 
        # outside the typically observed range. 
        fake_data = model_results[fake_data_ind,] + 1.0e-04*(0.5 - runif(nc))
        for(i in 1:nr) data_mat[i,] = fake_data
    }
    data_ranks = apply(model_results < data_mat, 2, f<-function(x) sum(x)) # equivalent of 'model_ranks' for data
    data_summary = multi_gauge_summary_fun(data_ranks)

    # This gives an empirical quantile of 'data statistic' relative to 'model statistic'
    output = (sum(model_summary < data_summary) + 0.5)/(length(model_summary)+1)

    if(return_environment){
        return(environment())
    }else{
        return(output)
    }
}

# Compute 'scores' for the data, for all observed events, with each earthquake
# generation type
stoc_fit = unlist(lapply(stochastic_store, score_gauge))
unif_fit = unlist(lapply(uniform_store, score_gauge))
vary_unif_fit = unlist(lapply(variable_uniform_store, score_gauge))

#
# Plot results for all cases
#
if(variable_mu){
    pdf('Gauge_statistic_sampling_distributions_varyMu.pdf', width=10, height=2.5)
}else{
    pdf('Gauge_statistic_sampling_distributions.pdf', width=10, height=2.5)
}
for(i in 1:length(stoc_fit)){
    par(mfrow=c(1,3))
   
    event_lab = gsub('gauge_summary_stats_session_', '', basename(all_Rdata[i]) )
    # Variable uniform 
    vary_unif_env = score_gauge(variable_uniform_store[[i]], return_environment=TRUE)
    n = length(vary_unif_env$model_summary)
    hist(vary_unif_env$model_summary/(n), col='green', main='Uniform slip stoc. size',
        xlab=paste0('S ', event_lab))
    abline(v=vary_unif_env$data_summary/(n), lwd=4)

    # Uniform
    unif_env = score_gauge(uniform_store[[i]], return_environment=TRUE)
    n = length(unif_env$model_summary)
    hist(unif_env$model_summary/(n), col='blue', main='Uniform slip fixed size',
        xlab=paste0('S ', event_lab))
    abline(v=unif_env$data_summary/(n), lwd=4)

    # Stochastic slip
    stoc_env = score_gauge(stochastic_store[[i]], return_environment=TRUE)
    n = length(stoc_env$model_summary)
    hist(stoc_env$model_summary/(n), col='red', main='Heterogeneous slip',
        xlab=paste0('S ', event_lab))
    abline(v=stoc_env$data_summary/(n), lwd=4)
}
dev.off()

# Record events to remove {GCMT Mw < 7.7}. Nice because the sample becomes 'the
# whole population of events > Mw 7.7 at DART 2007-2015' bar some aftershocks. 
# Otherwise, we include a few Mw GCMT 7.6's, which were opportunistically included
# because they were near Australia. Qualitative message is the same either way.
nk = c(3, 8)
stopifnot(grepl('kermadectonga_tonga_2009_03_19', all_Rdata[nk[1]]))
stopifnot(grepl('outerrise_kermadectonga_2011_07_06', all_Rdata[nk[2]]))
# Just looking at these small samples, we see for unif and vary_unif, much
# of the distribution is above 50% (corresponding to model underestimation)
# ks.test seems too weak to distinguish this -- recall it is well known to
# not have much tail influence. But the Anderson-Darling test
# (which involves tail weighting) does better. 
library(ADGofTest)
unif_ad = ad.test(unif_fit[-nk], punif)
unif_ad
vary_unif_ad = ad.test(vary_unif_fit[-nk], punif)
vary_unif_ad
stoc_ad = ad.test(stoc_fit[-nk], punif)
stoc_ad

# Plot the distribution of stoc_fit/unif_fit/vary_unif_fit should be very close
# to uniform if the NULL hypothesis is true, but in practice we see substantial
# deviations for models with uniform-slip
if(variable_mu){
    pdf('Null_hypothesis_test_varyMu.pdf', width=10, height=3)
}else{
    pdf('Null_hypothesis_test.pdf', width=10, height=3)
}
par(mfrow=c(1,3))
hist(unif_fit[-nk], col='blue', 
    main='F(S) over all events, uniform slip fixed size', 
    xlab='F(s)', xlim=c(0,1), ylim=c(0, 6))
text(0.3, 5, paste0('p = ', signif(unif_ad$p.value, 3)), cex=2.0) 
hist(vary_unif_fit[-nk], col='green', 
    main='F(S) over all events, uniform slip stoch. size', 
    xlab='F(s)', xlim=c(0,1), ylim=c(0,6))
text(0.3, 5, paste0('p = ', signif(vary_unif_ad$p.value, 3)), cex=2.0) 
hist(stoc_fit[-nk], col='red', 
    main='F(S) over all events, heterogeneous slip', 
    xlab='F(s)', xlim=c(0,1), ylim=c(0,6))
text(0.2, 5, paste0('p = ', signif(stoc_ad$p.value, 3)), cex=2.0)
dev.off()

#
# Now test that if data is 'random', then the ad.test returns non-significant
# results with the expected rate.
#
# This function computes a random value of the ad.test p.value, when the 'data'
# is actually a random draw from the model. If the logic is sound, then for
# each of stoch/unif/vary-unif, the p-values should have a near uniform
# distribution (with slight deviations tolerable due to the 'model
# distribution' actually being based on a finite set of samples.
#
run_random_test<-function(i){
    null_stoc_fit = unlist(lapply(stochastic_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_stoc = ad.test(null_stoc_fit[-nk], punif)$p.value

    null_unif_fit = unlist(lapply(uniform_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_unif = ad.test(null_unif_fit[-nk], punif)$p.value

    null_vary_unif_fit = unlist(lapply(variable_uniform_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_vary_unif = ad.test(null_vary_unif_fit[-nk], punif)$p.value

    return(list(p_store_stoc=p_store_stoc, p_store_unif=p_store_unif, 
        p_store_vary_unif=p_store_vary_unif))
}
# Reproducible randomness in parallel
set.seed(1234, "L'Ecuyer")
random_model_runs = mclapply(as.list(1:10000), run_random_test, mc.cores=12)
p_store_stoc = unlist(lapply(random_model_runs, f<-function(x) x$p_store_stoc))
p_store_unif = unlist(lapply(random_model_runs, f<-function(x) x$p_store_unif))
p_store_vary_unif = unlist(lapply(random_model_runs, f<-function(x) x$p_store_vary_unif))

# The 'p_store_xxx' variables should look similar to random samples from a
# uniform distribution! There might be deviations though due to the
# discreteness of the 'model distribution'
# However the uniform approximation seems very good.
mean(p_store_stoc < 0.05) # We want a 5% chance of being below 0.05 by chance #ad.test(p_store_stoc, punif)
mean(p_store_unif < 0.05) # We want a 5% chance of being below 0.05 by chance #ad.test(p_store_unif, punif)
mean(p_store_vary_unif < 0.05) # We want a 5% chance of being below 0.05 by chance #ad.test(p_store_vary_unif, punif)

#
# Make some plots. They use hard-coded Mw values and will break with
# file-changes, so wrap inside if(FALSE)
#
if(FALSE){

    #
    # How does the coverage statistic vary with Mw?
    #
    png('coverage_statistic_vs_mw.png', width=8, height=7, units='in', res=200)
    mws = c(8.1, 8, 7.7, 8.3, 9.1, 7.9, 7.8, 7.6, 7.8, 8.1, 7.8, 8.8, 8.0, 7.8, 8.2, 8.3, 7.8, 7.9, 8.5, 7.8)

    par(mfrow=c(2,2))
    plot(mws, stoc_fit, col='red', pch=19, ylim=c(0,1), main='Coverage statistic vs Mw \n for each model type')
    points(mws, vary_unif_fit, col='green', pch=19)
    points(mws, unif_fit, col='blue', pch=19)
    grid(); abline(h=0.5)

    plot(mws, stoc_fit - vary_unif_fit, main='Stochastic - variable uniform'); abline(h=0); grid()
    plot(mws, stoc_fit - unif_fit, main='Stochastic - uniform'); abline(h=0); grid()
    plot(mws, vary_unif_fit - unif_fit, main='Variable uniform - uniform'); abline(h=0); grid()

    dev.off()

    #
    # Quick look at 'high magnitude' events
    #
    k = which(mws >= 8.3)
    stoc_fit[k]
    vary_unif_fit[k]
    unif_fit[k]

    #
    # Question: Is the sd of log(stage-range) increasing with magnitude? 
    #
    # We might want this quantity to be 'fairly stable' with Mw, otherwise, if
    # it were increasing with Mw, it might suggest artefacts associated with
    # resolving more and more variability with Mw
    #
    sd_log_stagerange_stochastic = lapply(stochastic_store, f<-function(x) apply(log10(x$model), 2, sd)) 
    names(sd_log_stagerange_stochastic) = mws
    sd_log_stagerange_variable_uniform = lapply(variable_uniform_store, f<-function(x) apply(log10(x$model), 2, sd)) 
    names(sd_log_stagerange_variable_uniform) = mws
    sd_log_stagerange_uniform = lapply(uniform_store, f<-function(x) apply(log10(x$model), 2, sd)) 
    names(sd_log_stagerange_uniform) = mws
    # Doesn't seem to be a difference in the variance of variable_uniform or stochastic 
    summary(unlist(lapply(sd_log_stagerange_stochastic, median)))
    summary(unlist(lapply(sd_log_stagerange_variable_uniform, median)))
    summary(unlist(lapply(sd_log_stagerange_uniform, median)))

    # This figure illustrates the above points
    png('mean_sd_log_stage_range_vs_Mw.png', width=15, height=6, units='in', res=200)
    par(mfrow=c(1,2))
    plot(range(mws), c(0.05, 0.55), col=0, xlab='Mw', ylab='sd(log10(stage-range))',
        main='Mw vs sd(log10(stage-range)) at each DART \n Note similar variation of stochastic and variable-uniform')
    for(i in 1:length(sd_log_stagerange_stochastic)){
        x = sd_log_stagerange_stochastic[[i]]
        points(mws[i] + 0*x, x, col='red', cex=0.5)
        x = sd_log_stagerange_variable_uniform[[i]]
        points(mws[i] + 0*x, x, col='green', cex=0.5)
        x = sd_log_stagerange_uniform[[i]]
        points(mws[i] + 0*x, x, col='blue', cex=0.5)
    }

    #
    # Above, we find that the 'spread' of model results is similar for stochastic and variable uniform
    # However the differences in coverage imply the means clearly differ
    # Do they differ systematically with Mws? We might expect at low Mw, both approaches are similar 
    # due to resolution issues, while at higher Mw, stochastic starts to resolve peaks better (and they
    # might matter more as well, given the implicit smoothing of Okada/Kajiura at smaller spatial scales) 
    #
    if(variable_mu == TRUE) stop('The plot below will fail because variable mu cases might not have the same number of events')

    median_stochastic_vs_variable_uniform = mws*NA
    for(i in 1:length(mws)) median_stochastic_vs_variable_uniform[i] = median(stochastic_store[[i]]$model/variable_uniform_store[[i]]$model)
    plot(mws, median_stochastic_vs_variable_uniform, xlab='Mw', ylab='Median(stochastic/variable_uniform)', 
        main='Median(stochastic/variable_uniform) vs Mw. \n The increase with Mw seems consistent with resolution effects \n A factor 1.3 = 10^0.12, corresponds to a shift of about 0.5_x_the_logSD')
    abline(h=1.3, col='red'); grid()
    
    dev.off()

    #
    # Summary: 
    #          Stochastic has similar log10 variability to variable-uniform (log10 sd ~ 0.24), but has a higher median (about 1.3-1.4x, ignoring low Mw where resolution is low)
    #
    #          Fixed-size-uniform-slip has lower log10 variability than the other models, at least once we have enough resolution
    #
    #          The variability of log10-stage does not show increasing patterns with Mw in any models -- if anything it decreases, 
    #          presumably because geometric effects (like occurring under land) become less significant
    #
    #          Overall, these results DO NOT suggest strong artefacts associated with increases in Mw, such as strongly increasing variability or mean in the stochastic model.
    #          There is some increase in 'median(stochastic / variable-uniform)' with Mw, but it seems consistent with 'resolving stochastic slip', 
    #          and not so great as to be of strong concern (i.e. difference corresponds to a fraction of a log-standard-deviation)

    #
    # NEXT: 
    #        How do the statistical properties of 'good' solutions compare with the statistical properties of 'all' solutions?
}
