# Get base datasets and functions for analysis -- these were split out to keep the
# file size managable
source('analysis_data_and_functions.R')

#
# Boxplots of various statistics, all in the one pdf.
#
pdf(paste0(FIG_OUTPUT_BASEDIR, '/model_vs_data_statistic_dense_downsampled.pdf'), width=9, height=9)
unique_event_names = unique(event_stats$event_name)
desired_stats = c('_max_24', '_max_36', '_max',
                  '_stgrng_24', '_stgrng_36', '_stgrng',
                  '_rms_24', '_rms_36', '_rms')
# Store the indices in case we want to do something to them
index_store = vector(mode='list', length=length(desired_stats))
names(index_store) = paste0('ratio', desired_stats)
counter = 0
for(stat_type in desired_stats){
    counter = counter+1
    local_stat = vector(mode='list', length=length(unique_event_names))
    names(local_stat) = unique_event_names
    for(i in 1:length(unique_event_names)){
        # Make the boxplot for every event, and store the indices that are
        # included
        local_stat[[i]] = boxplot_sites(
            event_stats_downsampled_model, 
            event_name = unique_event_names[i], 
            stat_type = stat_type, 
            #good_nearshore_and_close_and_highres=FALSE)
            good_nearshore_and_close_and_highres=TRUE)
    }
    # Note: index_store[[counter]] will not change with counter, we don't
    # actually need to store every time.
    index_store[[counter]] = local_stat
    # Add in a check to confirm the above
    stopifnot(isTRUE(all.equal(index_store[[counter]], index_store[[1]])))
}
dev.off()


#
# Scatterplot of model-vs-data, combining all events, for 'max' and
# 'stage-range'
#
png(paste0(FIG_OUTPUT_BASEDIR, '/Model-vs-data_maxima_stagerange.png'), width=12, height=12, units='in', res=300)
par(mfrow=c(3,3))
library(rptha)

sites_compared = lapply(unique(event_stats$event_name), function(event_name){
    get_table_indices_to_process(event_stats, event_name, 
        good_nearshore_and_close_and_highres=TRUE, #stagerange_threshold_15min=0.4, 
        expand_multi_counts=TRUE, run_type='random_like_historic',
        rigidity_type='constant', exclude_batch2_doubleups=FALSE)
    })
all_inds = unlist(sites_compared) # Collapse every event together
for(stat_type in c('_max', '_stgrng', '_rms_during_obs')){

    data_stat = paste0('data', stat_type)
    model_stat = paste0('model', stat_type)

    if(stat_type == '_max'){
        stat_type_label='Maxima'
    }else if(stat_type == '_stgrng'){
        stat_type_label='Stage-range'
    }else if(stat_type == '_rms_during_obs'){
        stat_type_label='Root-mean-square'
    }else{
        stop('unknown stat_type')
    }

    for(slip_type in c('FAUS', 'VAUS', 'HS')){
        
        base_col = list(FAUS='blue', VAUS='green', HS='red')

        is_slip_type = (1 + 
            (event_stats_downsampled_model$slip_type[all_inds] == slip_type))
        pt_col = c('white', base_col[[slip_type]])[is_slip_type]
        pt_size = c(0, 1)[is_slip_type]

        plot(event_stats_downsampled_model[[data_stat]][ all_inds], 
             event_stats_downsampled_model[[model_stat]][all_inds], 
             log='xy', cex=pt_size, col=pt_col, 
             xlab=paste0('Observed ', stat_type_label, ' (m)'), 
             ylab=paste0('Model ', stat_type_label, ' (m)'),
             cex.lab=1.4)
        title(paste0(stat_type_label, ' ', slip_type), cex.main=1.6)
        add_log_axis_ticks(side=1)
        add_log_axis_ticks(side=2)
        grid()
        x_seq = seq(0.001, 5, len=1000)
        points(x_seq, x_seq, t='l', col='purple')
        points(x_seq, 2*x_seq, t='l', col='purple', lty='dashed')
        points(x_seq, 0.5*x_seq, t='l', col='purple', lty='dashed')
    }
}
dev.off()

#stop('deliberate halt here')
if(FALSE){
    # Plot time-series of all models vs data

    pdf('site_time_series_good_nearshore_highres_domain.pdf', width=8, height=4)
    for(i in which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)){
        tmp = try(plot_model_and_data_at_event_row(event_stats, i))
    }
    dev.off()

    pdf('site_time_series_good_nearshore_highres_domain_DOWNSAMPLE.pdf', width=8, height=4)
    for(i in which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)){
        tmp = try(plot_model_and_data_at_event_row(
            event_stats_downsampled_model, i, downsample_model_like_data=TRUE))
    }
    dev.off()
}
#
# Correlation with energy? Yep. 
# NOTE: Here I am not accounting for 'double-up' scenarios. That's OK if we
# are simply trying to establish relationships, but not to make claims about
# model bias.
#
pdf(paste0(FIG_OUTPUT_BASEDIR, '/energy_correlations.pdf'), width=6, height=6)
for(stat_type in c('_max_24', '_max_36', '_max')){
    k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
    s1 = event_stats$site_and_event[k]
    m1 = event_stats[[paste0('model', stat_type)]][k]
    d1 = event_stats[[paste0('data' , stat_type)]][k]
    e1 = event_stats$energy_start[k]
    slip_type_int = match(event_stats$slip_type[k], c('FAUS', 'HS', 'VAUS'))
    slip_type_col = c('blue', 'red', 'green')[slip_type_int]
    run_type_int = match(event_stats$run_type[k], 
        c('nonrandom_like_historic', 'random_like_historic'))
    run_type_pch = c(19, 1)[run_type_int]
    for(site in unique(s1)){
        n = which(s1 == site)
        print(site)
        local_cor = cor(m1[n], e1[n], method='s')
        plot(e1[n], m1[n], log='xy', 
             main=paste0(site, '\n', stat_type, ' ', round(local_cor, 3)),
             col=slip_type_col[n], pch=run_type_pch[n])
        abline(h=d1[n], col='red')
        abline(h=d1[n]*1.5, col='red', lty='dotted')
        abline(h=d1[n]*1/1.5, col='red', lty='dotted')
        grid()
        sqrt_energy = sqrt(e1[n])
        max_stage = m1[n]
        linear_fit = lm(max_stage ~ sqrt_energy + 0)
        predicted_max_stage = coef(linear_fit)*sqrt_energy
        points(sqrt_energy**2, predicted_max_stage, pch=19, col='black', 
               cex=0.5)
    }
}
dev.off()

#  #
#  # How often are we too high, too low, whatever
#  # FIXME: Not accounting for counts, good nearshore, ...
#  fraction_smaller = aggregate(
#      event_stats_downsampled_model$model_max/event_stats_downsampled_model$data_max, 
#      by=list(event_stats_downsampled_model$site_and_event, 
#              event_stats_downsampled_model$slip_type, 
#              event_stats_downsampled_model$run_type),
#      function(x) mean(x < 1))
#  
#  # Compare noise in gauges before the tsunami arrival with the stage-range. Gives
#  # an indication of the tsunami prominence. 
#  k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
#  #data_noise_rms_indicator = aggregate(event_stats$data_rms_before_arrival[k]/event_stats$data_stgrng_noNA[k], 
#  #    by=list(site=event_stats$sites[k], event=event_stats$event_name[k]), function(x) c(mean(x), sd(x)) )
#  data_noise_rms_indicator = aggregate(event_stats$data_rms_before_arrival[k]/event_stats$data_rms_during_obs[k], 
#      by=list(site=event_stats$sites[k], event=event_stats$event_name[k]), function(x) c(mean(x), sd(x)) )
#  data_noise_rms_indicator2 = aggregate(sqrt(1-(event_stats$data_rms_before_arrival[k]/event_stats$data_rms_during_obs[k])**2), 
#      by=list(site=event_stats$sites[k], event=event_stats$event_name[k]), function(x) c(mean(x), sd(x)) )
#  data_noise_max_indicator = aggregate(event_stats$data_max_before_arrival[k]/event_stats$data_stgrng_noNA[k], 
#      by=list(site=event_stats$sites[k], event=event_stats$event_name[k]), function(x) c(mean(x), sd(x)) )
#  data_noise_stgrng_indicator = aggregate(event_stats$data_stgrng_before_arrival[k]/event_stats$data_stgrng_noNA[k], 
#      by=list(site=event_stats$sites[k], event=event_stats$event_name[k]), function(x) c(mean(x), sd(x)) )
#  
#
# More detailed analysis of tsunami period and spectra
#
# For each site with 1min data, get the energy-in_bands for ALL modelled time-series
# AND for the nearest offshore time-series [1000 m contour]. Then see if you can
# model the nearshore spectra using regression of the offshore spectra.
#
# Why does this matter? 
#    - The most obvious reason is comparison with data. But that requires
#      analysis of the data, and we expect some contamination due to
#      background atmospheric stuff.
#    - It may give us a way to distinguish models with different
#      slip types -- e.g. do the FAUS models consistently have smaller
#      high-frequency content both offshore and at gauges? Does this
#      lead to poor performance at sites that are sensitive to short
#      frequencies [e.g. Ulladullah], and good performance at sites that are
#      not? We already know there are differences in initial-potential-energy -- 
#      does this tell us any more?
#    - Also it will be interesting to see if the nearshore spectra are
#      simply related to the offshore spectra. And if the offshore spectra
#      are related to the slip-type.
#
#
## Here is a visualisation routine which can help to 'get a feel' for the spectra.
#get_model_and_observed_spectral_summary<-function(event_stats_row, PLOT_AMP=FALSE, PLOT_TS=FALSE){
#    # Call like:
#    #     get_model_and_observed_spectral_summary(event_stats_row=event_stats_downsampled_model[18866,], PLOT_AMP=TRUE, PLOT_TS=TRUE)
#    matching_RDS = match(event_stats_row$run_name, names(target_RDS_data))
#    site = event_stats_row$sites
#
#    # Observed data
#    obs = data.frame(
#        juliant = target_RDS_data[[matching_RDS]][[site]]$event_data$obs$juliant,
#        stage = target_RDS_data[[matching_RDS]][[site]]$event_data$obs$resid)
#    # Modelled data
#    model = data.frame(
#        juliant = (target_RDS_data[[matching_RDS]][[site]]$model_time/(3600*24) + 
#                   target_RDS_data[[matching_RDS]][[site]]$model_start_julian),
#        stage = target_RDS_data[[matching_RDS]][[site]]$model_stage)
#
#    # Observed spectra
#    kr = which(  obs$juliant > event_stats_row$model_arrival_time)
#    obs_amp = get_amplitudes_by_period(
#        time=as.numeric(obs$juliant[kr] - obs$juliant[kr[1]])*24*3600,
#        stage=obs$stage[kr], taper=0.1)
#    kr_obs = kr
#
#    # Modelled spectra
#    kr = which(model$juliant > event_stats_row$model_arrival_time)
#    model_amp = get_amplitudes_by_period(
#        time=as.numeric(model$juliant[kr] - model$juliant[kr[1]])*24*3600,
#        stage=model$stage[kr], taper=0.1)
#
#    # Make some plots below here
#
#    if(PLOT_AMP & PLOT_TS){
#        par(mfrow=c(2,1))
#    }else if(PLOT_AMP | PLOT_TS){
#        par(mfrow=c(1,1))
#    }
#
#    if(PLOT_AMP){
#        YLIM = c(0, max(c(max(obs_amp$amp), max(model_amp$amp))))
#        plot(obs_amp$period/3600, obs_amp$amp, t='l', log='x', ylim=YLIM)
#        points(model_amp$period/3600, model_amp$amp, t='l', col='red')
#        title(paste0(event_stats_row$slip_type, ' ', event_stats_row$event_name, 
#                     ' ', event_stats_row$sites, '\n', event_stats_row$run_name))
#    }
#
#    if(PLOT_TS){
#        plot(model$juliant[kr], model$stage[kr], t='l', col='red')
#        points(obs$juliant[kr_obs], obs$stage[kr_obs], t='l')
#        grid()
#        abline(h=0, col='orange')
#    }
#
#    output = list(model_max_period = model_amp$period[which.max(model_amp$amp)],
#                  obs_max_period = obs_amp$period[which.max(obs_amp$amp)])
#    return(output)
#}





#
# Boxplot and scatterplot of model/obs
#
source('boxplot_stats_and_scatterplots.R')

stop('deliberate stop here')
# ##
# # Given a single variable "data", convert it to a list of tables, with each table having sites as
# # columns, and scenarios as rows, for an invidual event and slip-type.
# # Must be a better way to do this?!! 
# # However once you have the tables, it is interesting to consider
# # - How closely do the gauges correlate with each other?
# # - How closely do they correlate with energy?
# convert_to_table_per_scenario_per_slip_type<-function(data, site, event, slip_type, sliptype_event_scenarioID){
# 
#     unique_event = unique(event)
#     unique_slip_type = unique(slip_type)
# 
#     list_of_tables = vector(mode='list', length=length(unique_event))
#     names(list_of_tables) = unique_event
#     for(event_ind in 1:length(unique_event)){
#         u_event = unique_event[event_ind]
# 
#         list_of_tables[[event_ind]] = vector(mode='list', length=length(unique_slip_type))
#         names(list_of_tables[[event_ind]]) = unique_slip_type
# 
#         for(slip_type_ind in 1:length(unique_slip_type)){
# 
#             u_slip_type = unique_slip_type[slip_type_ind]
# 
#             k = which(event == u_event & slip_type == u_slip_type)
#             unique_site_k = unique(site[k])
#             
#             numcol = length(unique_site_k)
#             numrow = length(k)/numcol
#             browser()
#             stopifnot(numrow == 60) # 60 scenarios per historical event
#             model_mat = matrix(NA, nrow=numrow, ncol=numcol, byrow=TRUE)
#             counter = 0
#             for(u_site in unique_site_k){
#                 counter = counter+1
#                 # Data for site i
#                 p = which(site[k] == u_site)
#                 stopifnot(length(p) == numrow)
# 
#                 # Check alignment of rows
#                 if(counter == 1){
#                     sc0 = sliptype_event_scenarioID[k[p]]
#                 }else{
#                     stopifnot(all(sliptype_event_scenarioID[k[p]] == sc0))
#                 }
# 
#                 model_mat[,counter] = data[k[p]]
#             }
#             colnames(model_mat) = unique_site_k
#             rownames(model_mat) = sc0
# 
#             # Store the matrix
#             list_of_tables[[event_ind]][[slip_type_ind]] = model_mat
#         }
#     }
#     return(list_of_tables)
# }
#  
# if(FALSE){
#     # Trying out some statistical modelling to understand/justify the structure of the model errors 
#     #
#     # Key issues:
#     #     - The box-plots suggest decomposing the variation into
#     #         + A constant slip-model effect
#     #         + An event effect
#     #         + A gauge effect
#     #         + The remaining variability is highly correlated among different gauges for the same scenario.
#     # 
#     # With only 6 observed events, we cannot really hope to reduce each event
#     # to binary variable per gauge AND show deviations from randomness.
#     #
#     # However, we have 2 events with clear under-prediction for FAUS [Chile
#     # 1960, Solomon 2007].  The other models are also 'on the low side' for
#     # these events, but their tendency for larger waves overall, and greater
#     # variability, means they can better produce tsunamis that look like
#     # scenarios.
#     #
#     # For the other events the FAUS model still has a tendency to be 'always
#     # too low' at individual gauges, while it is rarely too high. An exception
#     # is the Freemantle-gauge (Sumatra 2004) but this was also over-predicted
#     # by all of the source-inversions in Davies et al. (2020)
#     #
#     # For the Puysegur 2009 scenario, the FAUS model looks reasonable. For this
#     # event Davies and Griffin (2018) found the FAUS model predicted one DART
#     # high-ish, and one DART low-ish, so it's hard to make strong statements
#     # based on that. 
#     #
#     # For Sumatra2004, Chile 2010 and Tohoku 2011, the FAUS model is OK albeit
#     # on the low-side. Recall that Davies (2019) reports strong underestimation
#     # of Tohoku leading waves by the FAUS model -- however inspection of their 
#     # DART records indicate that was not the case for DARTs 55012 and 55023 in the
#     # Coral Sea, and other darts 'nearer our site' (52406) also don't show FAUS
#     # under-estimation. So it could be a regional effect.
#     #     - For the Tohoku event, visually it looks like the FAUS spectra may be
#     #       'often too long period' ?
#     #
# 
# 
#     # Collapse every event together
#     sites_compared = lapply(unique(event_stats$event_name), function(event_name){
#         get_table_indices_to_process(event_stats, event_name, 
#             good_nearshore_and_close_and_highres=FALSE, stagerange_threshold_15min=0.4, 
#             expand_multi_counts=TRUE, run_type='random_like_historic',
#             rigidity_type='constant', exclude_batch2_doubleups=FALSE)
#         })
#     all_inds = unlist(sites_compared) 
# 
#     #obs   = event_stats_downsampled_model$data_max_noNA[ all_inds]
#     #model = event_stats_downsampled_model$model_max_noNA[all_inds]
#     obs   = event_stats_downsampled_model$data_stgrng_noNA[ all_inds]
#     model = event_stats_downsampled_model$model_stgrng_noNA[all_inds]
# 
#     slip_type = as.factor(event_stats_downsampled_model$slip_type[all_inds])
#     site = as.factor(event_stats_downsampled_model$sites[all_inds])
#     event = as.factor(event_stats_downsampled_model$event_name[all_inds])
#     scenario_ID = as.factor(event_stats_downsampled_model$scenario_ID[all_inds])
#     energy_start_sqrt = sqrt(event_stats_downsampled_model$energy_start[all_inds])
#     log_energy_start_sqrt = log(energy_start_sqrt)
#     site_and_event = as.factor(event_stats_downsampled_model$site_and_event[all_inds])
# 
#     sliptype_site = as.factor(paste0(slip_type, '_', site))
#     sliptype_event = as.factor(paste0(slip_type, '_', event))
#     sliptype_site_event = as.factor(paste0(slip_type, '_', site, '_', event))
#     sliptype_event_scenarioID = as.factor(paste0(slip_type, '_', event, '_', scenario_ID))
# 
#     # Make a normalised log energy
#     event_median_energy = aggregate(log_energy_start_sqrt, by=list(event), median)
#     event_median_energy = event_median_energy$x[ match(event, event_median_energy[,1]) ]
#     norm_log_energy = log_energy_start_sqrt - event_median_energy
# 
#     scenario_event_sliptype = as.factor(paste0(scenario_ID, '_', event, '_', slip_type))
# 
#     err = log(model/obs)
# 
#     #convert_to_table_one_event_scenario_vs_site<-function(event_name, site, 
#     #    event, sliptype_event_scenarioID){
# 
#     #    k = which(event == event_name)
#     #    site = site[k]
# 
#     #}
# 
#     # This is getting toward what I want -- the modelled variation depends on the
#     # slip_type and the event[potentially varying by slip-type] and a site-specific
#     # effect that is not dependent on the slip type.
#     model_1 = lm(err ~ 0 + slip_type + event:slip_type + event:site)
#     # We can see that the HS and VAUS have lower overall error than the FAUS [note this
#     # should be combined with the 'event:slip_type' terms to interpret it]
#     # Investigation of the event:slip_type coefficients suggests that for a fixed
#     # event, they are not too different across slip-types. 
#     # Another interesting thing is that, for a given site, the "event:site" terms 
#     # bounce around all over the place -- so it's not like we see consistent site
#     # biases.
#     # A weakness of this model is that the variances are constant -- whereas I might
#     # want a model-specific variance.
#     #
#     # Given the insensitvity of event:slip_type to slip_type, here we fit another
#     # model with a single 'event' term. It has a higher AIC, but still maybe it's more
#     # parsimonious in some sense? 
#     model_2 = lm(err ~ 0 + slip_type + event + event:site)
#     # The above model still doesn't treat the unequal-variances issue. Furthermore, 
#     # I "feel like" the "event" and "event:site" terms should be random effects, as I
#     # conceptualise them like that (realisations of something random).
#     # 
#     # Here I try lme, which I understand can do variances changing by variable.
#     library(nlme)
#     (model_2.5 = lme(err ~ 0 + slip_type, random= ~1|event/site, 
#                     weights=varIdent(form=~1|slip_type)) )
#     model_2.5$coefficients # Here are the coefficients -- looks sensible
#     intervals(model_2.5) # Uncertainties on coefficients -- beware these might be based on simple approximations.
# 
#     # This one allows unequal variances for each sliptype/site combination
#     (model_2.5B = lme(err ~ 0 + slip_type, random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_site), method='ML',
#                     control=list(maxIter=200, msMaxIter=200)) )
#     # This one allows unequal variances for each sliptype/site/event combination,
#     # so there are many variance parameters, but the AIC is better than the previous
#     # model
#     (model_2.5C = lme(err ~ 0 + slip_type, random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_site_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     # This one has a simpler variance model than above -- varies with slip-type and
#     # event -- and it has a better AIC than the previous model, with fewer parameters.
#     (model_2.5D = lme(err ~ 0 + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
# 
#     # Something that's not depicted in this model is possible dependence of the residuals
#     # between gauges, within groups having the same "slip-type + event"
#     # For illustration
#     k = which(slip_type == 'HS' & event == 'sumatra2004')
#     # By inspection, scenario_ID[k] is in contiguous blocks of 11, reflecting that
#     # for every scenario we have 11 gauge observations. This is the part that I think
#     # can be well 'collapsed' based on the scenario potential energy.
#     m1 = matrix((model)[k], ncol=11, byrow=TRUE); colnames(m1) = as.character(site[k[1:11]])
#     pairs(m1) # Clearly strong correlations between observations at different sites.
#     #
#     # How to get that?
#     #(model_2.5E = lme(err ~ 0 + slip_type, 
#     #                  random= ~1|event/site, 
#     #                  #correlation=corSymm(form=~1|sliptype_event),
#     #                  weights=varIdent(form=~1|sliptype_event), 
#     #                  method='ML',
#     #                  control=list(maxIter=500, msMaxIter=500)) )
# 
#     (model_2.5E = lme(err ~ event + log_energy_start_sqrt:event, 
#     #(model_2.5E = lme(err ~ 0 + energy_start_sqrt*event, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     # This has a far better AIC than previous. 
# 
#     (model_2.5Eb = lme(err ~ event + log_energy_start_sqrt:event, 
#                     random= ~1|event/site, 
#                     #random= ~site|event, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     (model_2.5Ec = lme(err ~ log_energy_start_sqrt:event, 
#                     random= ~1|event/site, 
#                     #random= ~site|event, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     # What about if we make the energy coefficient depend on the site_and_event, rather
#     # than just the event? This leads to many more coefficients, and worse AIC
#     (model_2.5F = lme(err ~ event + log_energy_start_sqrt:site_and_event, 
#                     random= ~1|event/site, 
#                     method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     # This one has slightly worse AIC and BIC, than the one that just had
#     # energy.  Coefficients are hard to interpret due to the scaling of
#     # log_energy_start_sqrt. I like the earlier model that had clear
#     # FAUS/VAUS/HS coefficients to understand the bias. FIXME: Understand how
#     # to interpret this stuff!
#     (model_2.5G = lme(err ~ event + log_energy_start_sqrt:event + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
# 
#     #
#     # CONCEPT: 
#     # -- A) Show that FAUS, VAUS and HS have different biases -- this model does that,
#     #       and it emphasises that HS/VAUS have less bias then FAUS, but more variance.
#     library(nlme)
#     (model_2.5D = lme(err ~ 0 + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     # -- B) Show that this collapses (largely) if we use energy as a predictor.
#     # Perhaps we can do this if we 'normalise' the energy for each event. 
#     (model_2.5D2 = lme(err ~ norm_log_energy, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     #        So here, adding 'slip_type' as a predictor no longer helps much.
#     #        AIC and BIC are worse -- and the slip-type coefficients are barely
#     #        significant.
#     (model_2.5D3 = lme(err ~ norm_log_energy + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     #
#     #
#     # I REALLY NEED TO PLOT THE MODEL vs DATA TO SEE IF THIS IS WORKING.
#     # - Also good to do some normality checks
# 
#     # What happens to the 'categorical' model if we use a different residual
#     # standard-deviation for each gauge? -- AIC/BIC is a bit worse. Little change
#     # to other model parameters. 
#     (model_2.5D4 = lme(err ~ 0 + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_site_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     # As above for the 'energy' model -- leads to AIC-better, BIC-worse as compared
#     # to the model that treated the variance as related to sliptype and site.
#     # Again the coefficients are fairly similar. 
#     (model_2.5D5 = lme(err ~ 1 + norm_log_energy, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_site_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     # Here we add the predictor of slip-type. It turns out that for these
#     # coefficients, the uncertainties are much greater than the differences with
#     # slip-types. Also the coefficients themselves are of 'borderline statistical
#     # significance'.
#     (model_2.5D5_alternate = lme(err ~ 1 + norm_log_energy + slip_type, 
#                     random= ~1|event/site, 
#                     weights=varIdent(form=~1|sliptype_site_event), method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     (TESTER = lme(err ~ 0 + slip_type, 
#                     random= ~1|sliptype_event_scenarioID, 
#                     #weights=varIdent(form=~1|sliptype_site_event), 
#                     method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
#     (TESTERb = lme(err ~ 1, 
#                     random= ~norm_log_energy|event/site, 
#                     #weights=varIdent(form=~1|sliptype_site_event), 
#                     method='ML',
#                     control=list(maxIter=500, msMaxIter=500)) )
# 
#     ##
#     ## If we want to account for correlation in the residuals, then I guess
#     ## "generalised least squares" is the way to go? Actually 'lme' supports this too.
#     ##
#     ### This one is similar to the lme models.
#     #test = gls(err ~ 0 + slip_type + site_and_event, 
#     #                weights=varIdent(form=~1|sliptype_site_event), 
#     #                method='ML',
#     #                control=list(maxIter=500, msMaxIter=500))
#     ## This one adds correlation in the residuals for each event
#     ## Taking ages to fit!
#     #test2 = gls(err ~ 0 + slip_type + site_and_event, 
#     #                correlation=corSymm(form=~1|site_and_event), #
#     #                #weights=varIdent(form=~1|sliptype_site_event), 
#     #                method='ML',
#     #                control=list(maxIter=500, msMaxIter=500))
# 
# 
#     # Above we have been assuming normality -- good to check if that is really the case.
#     # A simple graphical approach is to employ qqnorm()
#     # Visually log-normalitity is usually a good approximation. We do see a few
#     # possibly heavy-tailed sites. Recall for the DARTs paper analysis, I found this
#     # would happen at sites where some events could have most deformation 'on-land'
#     # -- i.e. some obstructions to wave generation.
#     pdf('qqnorm_plots.pdf', width=10, height=10)
#     unique_site_event = unique(site_and_event)
#     unique_slip_type = unique(slip_type)
#     for(j in 1:length(unique_slip_type)){
#         for(i in 1:length(unique_site_event)){
#             k = which(site_and_event == unique_site_event[i] & 
#                       slip_type == unique_slip_type[j])
#             err_k = err[k]
#             qqnorm(err_k, main=""); qqline(err_k)
#             test_stat = shapiro.test(err_k)
#             title(paste0(unique_slip_type[j], ' --- ', unique_site_event[i], 
#                          '\n p=', round(test_stat$p.value, 3)))
#         }
#     }
#     dev.off()
# 
# 
#     # A simple approach -- plotting the median value of 'model/obs' at each gauge
#     # for each slip-type.
#     meds = aggregate(err, by=list(site=site, slip_type=slip_type, event=event), median)
#     plot(site_medians$x, pch=as.integer(meds$slip_type), col=as.integer(meds$slip_type))
#     plot(meds$x, pch=as.integer(meds$slip_type), col=as.integer(meds$slip_type))
# 
# 
#     ### EXPERIMENTS WITH LME4 -- Actually the notation is similar to 'lme' except
#     ### the "random=" of lme instead appears in the equation, in parentheses.
#     #library(lme4)
#     ## Fixed effect is slip-type, intercept varying among "event" and "site" within
#     ## "event"
#     #model_4 = lmer(err ~ 0 + slip_type + (1|event/site), REML=FALSE) 
#     ## This is almost identical to model_2.5 above, but doesn't have the different
#     ## variances.
#     #print(model_4)
#     ## Random effects -- seems sensible
#     #ranef(model_4)
#     ## Notice how the 'events' constants sum to 1
#     #sum(ranef(model_4)$events)
#     ## The fixed effects are intuitive.
#     #fixef(model_4)
# 
# 
# 
# }
# 
# if(FALSE){
#     # Scatterplot of model-vs-data statistic
# 
#     for(i in 1:length(sites_compared)){
#         inds = sites_compared[[i]]
#         plot(event_stats_downsampled_model$model_max_noNA[inds], 
#              event_stats_downsampled_model$data_max_noNA[inds], 
#              log='xy')
#         add_log_axis_ticks(side=1)
#         add_log_axis_ticks(side=2)
#         grid()
#         abline(0, 1, col='red')
#     }
# 
# }

