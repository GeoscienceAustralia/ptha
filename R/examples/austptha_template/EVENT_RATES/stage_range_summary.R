variable_mu = TRUE

if(variable_mu){
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*varyMu.Rdata')
}else{
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*[0-9].Rdata')
}

uniform_store = list()
stochastic_store = list()
variable_uniform_store = list()

for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    load(all_Rdata[i], envir=event_env)

    ngauges = length(event_env$uniform_slip_stats)

    uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    stochastic_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    variable_uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )

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

}

#
# Compute the fraction of model events that had stage-range exceeding the data.
# If multiple gauges exist, take the median value over all gauges.
#
meds_uniform          = unlist(lapply(uniform_store         , f<-function(x) median(colMeans(x$model > x$data))))
## PRELIMINARY
# [1] 0.4166667 0.8333333 0.0000000 0.0000000 0.9285714 0.1785714 0.3157895 0.8333333 0.3484848 0.0000000 0.5000000 0.5714286
## SECOND RUN
# [1] 0.3333333 0.8333333 0.0000000 0.0000000 0.6666667 0.2142857 0.2894737 0.8333333 0.3030303 0.0000000 0.5000000 0.5000000
## THIRD RUN, 15 CASES
# [1] 0.3333333 0.8333333 0.0000000 0.0000000 0.6666667 0.0000000 0.2142857
# [8] 0.0000000 0.2894737 0.0000000 0.8333333 0.3030303 0.0000000 0.5000000
# [15] 0.5000000
## 4TH RUN 16 CASES Mw-variability-0.15
# [1] 0.4722222 0.6388889 0.0000000 0.0000000 0.5263158 0.2500000 0.1481481
# [8] 0.0000000 0.2268041 0.0000000 0.7125000 0.7311828 0.4040404 0.0000000
# [15] 0.4901961 0.4375000
## 5TH RUN 16 CASES mu-variability & Mw-variability-0.15
# [1] 0.36666667 0.77419355 0.00000000 0.00000000 0.52631579 0.46153846
# [7] 0.10000000 0.00000000 0.15217391 0.00000000 0.56944444 0.49450549
# [13] 0.22666667 0.09090909 0.47706422 0.28409091


meds_stochastic       = unlist(lapply(stochastic_store      , f<-function(x) median(colMeans(x$model > x$data))))
## PRELIMINARY
# [1] 0.6444444 0.8944444 0.6176471 0.1111111 0.7952381 0.1761905 0.4438596 0.8244444 0.4939394 0.2825397 0.7595238 0.5833333
## SECOND RUN
# [1] 0.6666667 0.9444444 0.5843137 0.1333333 0.8888889 0.3095238 0.4561404 0.8633333 0.4595960 0.2603175 0.7809524 0.5142857
## THIRD RUN 15 CASES
# [1] 0.6222222 0.9277778 0.5549020 0.1000000 0.8555556 0.6333333 0.2714286
# [8] 0.1000000 0.5175439 0.1777778 0.8566667 0.5464646 0.2380952 0.8309524
#[15] 0.5190476
## 4TH RUN 16 CASES Mw-variability-0.15
# [1] 0.6314815 0.9185185 0.5725490 0.1500000 0.8175439 0.6000000 0.2827160
# [8] 0.1222222 0.4927835 0.2353535 0.7650000 0.8275986 0.5569024 0.2676768
# [15] 0.7823529 0.5533333
## 5TH RUN 16 CASES mu-variability & Mw-variability-0.15
# [1] 0.58880309 0.91990847 0.42026144 0.17880795 0.83082707 0.84615385
# [7] 0.29446064 0.07246377 0.40831629 0.19553073 0.55456172 0.63519471
# [13] 0.47129909 0.42986425 0.78581363 0.47242921

meds_variable_uniform = unlist(lapply(variable_uniform_store, f<-function(x) median(colMeans(x$model > x$data))))
## PRELIMINARY
# [1] 0.35555556 0.65555556 0.18823529 0.02222222 0.49047619 0.11904762 0.25964912 0.61666667 0.27474747 0.14920635 0.50714286 0.44047619
## SECOND RUN
# [1] 0.34444444 0.73888889 0.22745098 0.01666667 0.55555556 0.19047619 0.24561404 0.62000000 0.26868687 0.12698413 0.49761905 0.40952381
## THIRD RUN 15 CASES
# [1] 0.38888889 0.76111111 0.20588235 0.06666667 0.53333333 0.25000000
# [7] 0.15476190 0.03750000 0.27017544 0.04444444 0.66000000 0.27777778
#[13] 0.10476190 0.49047619 0.41190476
## 4TH RUN 16 CASES Mw-variability-0.15
# [1] 0.39444444 0.74074074 0.21895425 0.03888889 0.53684211 0.37777778
# [7] 0.16419753 0.04722222 0.26735395 0.07878788 0.80000000 0.62473118
# [13] 0.33905724 0.14141414 0.49477124 0.44166667
## 5TH RUN 16 CASES mu-variability Mw-variability-0.15
# [1] 0.37044146 0.72325581 0.15490196 0.07594937 0.52029520 0.61025641
# [7] 0.16806723 0.04196643 0.18490823 0.04454976 0.59275618 0.44222222
# [13] 0.26091954 0.24238876 0.52620690 0.34050445

## SUMMARY OF 4TH RUN -- fixed mu & mw variability -- compare with similar variable mu run below
## nk = c(2, 6) # Remove events with Mw < 7.8, so we can apply the analysis to 'all Mw>=7.8 events in 2007-2015'.
## > summary(meds_variable_uniform[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.03889 0.14710 0.30320 0.32770 0.48150 0.80000 
## > summary(meds_stochastic[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1222  0.2714  0.5551  0.5041  0.7316  0.8276 
## > summary(meds_uniform[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.0000  0.3154  0.2964  0.4857  0.7312 
#
#### NOTE: The ks.tests below are invalid theoretically, because we should not expect the
####       median (over gauges) of the 'proportion of model runs >' should have a 
####       uniform distribution (unless the gauges were co-monotonic -- which they
####       are not, although they are strongly positively correlated).
## > ks.test(meds_uniform[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_uniform[-nk]
## D = 0.35714, p-value = 0.05623
## alternative hypothesis: two-sided
## 
## Warning message:
## In ks.test(meds_uniform[-nk], "punif") :
##   ties should not be present for the Kolmogorov-Smirnov test
## > ks.test(meds_stochastic[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_stochastic[-nk]
## D = 0.1724, p-value = 0.7385
## alternative hypothesis: two-sided
## 
## > ks.test(meds_variable_uniform[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_variable_uniform[-nk]
## D = 0.3203, p-value = 0.08898
## alternative hypothesis: two-sided



#
#
## SUMMARY OF 5th RUN -- variable mu & mw variability -- compare with similar fixed mu run above
## > nk = c(2, 6) # Remove events with Mw < 7.8.
## > summary(meds_uniform[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00000 0.02273 0.18940 0.23480 0.44950 0.56940 
## > summary(meds_stochastic[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.07246 0.32290 0.45060 0.45280 0.58020 0.83080 
## > summary(meds_variable_uniform[-nk])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.04197 0.15820 0.25170 0.28330 0.42430 0.59280 
##
## > ks.test(meds_uniform[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_uniform[-nk]
## D = 0.43056, p-value = 0.01114
## alternative hypothesis: two-sided
## 
## Warning message:
## In ks.test(meds_uniform[-nk], "punif") :
##   ties should not be present for the Kolmogorov-Smirnov test
## > ks.test(meds_stochastic[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_stochastic[-nk]
## D = 0.22195, p-value = 0.4332
## alternative hypothesis: two-sided
## 
## > ks.test(meds_variable_uniform[-nk], 'punif')
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  meds_variable_uniform[-nk]
## D = 0.40724, p-value = 0.01302
## alternative hypothesis: two-sided

meds2_uniform          = (lapply(uniform_store         , f<-function(x) (colMeans(x$model > x$data))))
meds2_stochastic       = (lapply(stochastic_store      , f<-function(x) (colMeans(x$model > x$data))))
meds2_variable_uniform = (lapply(variable_uniform_store, f<-function(x) (colMeans(x$model > x$data))))

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
# Score each event, by comparing the 'gauge-averaged' non-exceedance rate of
# the observed stage range, to that predicted by the model
#
multi_gauge_summary_fun = median # mean 
score_gauge<-function(dta, fake_data_by_perturbing_random_model=FALSE, return_environment=FALSE){
    # For each model run, get the rank of its stage-range at each gauge, and then
    # apply multi_gauge_summary_fun across all gauges. 
    nr = nrow(dta$model)
    nc = ncol(dta$model)

    ## APPROACH 1. 
    model_ranks = apply(dta$model, 2, rank)
    model_summary = apply(model_ranks, 1, multi_gauge_summary_fun) # summary of rank over all gauges

    ## APPROACH 2 -- use a 'leave-one-out' approach to compute the ranks.
    ## Gives similar results. 
    #model_ranks = dta$model * NA
    #for(i in 1:nr){
    #    # Compare (the model without row_i) to (row_i)
    #    for(j in 1:nc){
    #        model_ranks[i,j] = sum(dta$model[-i,j] < dta$model[i,j]) 
    #    }
    #}
    #model_summary = apply(model_ranks, 1, multi_gauge_summary_fun) # mean rank


    # For the data, get the rank of its stage-range at each gauge, and then apply
    # multi_gauge_summary_fun across all gauges
    if(!fake_data_by_perturbing_random_model){
        # Typical case
        data_mat = dta$data
    }else{
        # Testing-only case.
        # Make up some fake data by perturbing the model data. This is only
        # useful for testing the performance of the method (e.g. that the results
        # suggest the null hypothesis is satisfied, when it really is!)
        #
        data_mat = dta$data * NA
        fake_data_ind = sample(1:length(dta$model[,1]), size=1)
        # Here is our random data -- a perturbation on a model result. Because
        # the statistic is rank based, the size of the perturbation is not an
        # issue, so long as it is 'small enough' to not push the random data 
        # outside the typically observed range. 
        fake_data = dta$model[fake_data_ind,] + 1.0e-04*(0.5 - runif(nc))
        for(i in 1:nr) data_mat[i,] = fake_data
    }
    data_ranks = apply(dta$model < data_mat, 2, f<-function(x) sum(x)) # equivalent of 'model_ranks' for data
    data_summary = multi_gauge_summary_fun(data_ranks)

    # This gives an empirical quantile of 'data statistic' relative to 'model statistic'
    output = (sum(model_summary < data_summary) + 0.5)/(length(model_summary)+1)

    if(return_environment){
        return(environment())
    }else{
        return(output)
    }
}

# Compute 'scores' for the data, for all observed events, with each earthquake generation type
stoc_fit = unlist(lapply(stochastic_store, score_gauge))
unif_fit = unlist(lapply(uniform_store, score_gauge))
vary_unif_fit = unlist(lapply(variable_uniform_store, score_gauge))

# Plot results for all cases
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
# Now test that if data is 'random', then the ad.test returns non-significant results
# with the expected rate.
#
NN = 300 # Large enough sample size to detect p-values having significant deviations from uniform distribution
p_store_stoc = rep(NA, NN)
p_store_unif = rep(NA, NN)
p_store_vary_unif = rep(NA, NN)
for(i in 1:NN){
    print(i)
    null_stoc_fit = unlist(lapply(stochastic_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_stoc[i] = ad.test(null_stoc_fit[-nk], punif)$p.value

    null_unif_fit = unlist(lapply(uniform_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_unif[i] = ad.test(null_unif_fit[-nk], punif)$p.value

    null_vary_unif_fit = unlist(lapply(variable_uniform_store, 
        f<-function(x) score_gauge(x, fake_data_by_perturbing_random_model=TRUE)))
    p_store_vary_unif[i] = ad.test(null_vary_unif_fit[-nk], punif)$p.value
}
# The 'p_store_xxx' variables should look like random samples from a uniform
# distribution! Indeed they do. This suggests the above results are OK.
ad.test(p_store_stoc, punif)
ad.test(p_store_unif, punif)
ad.test(p_store_vary_unif, punif)

#
# Make some plots. They use hard-coded Mw values and will break with file-changes, so wrap inside if(FALSE)
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
