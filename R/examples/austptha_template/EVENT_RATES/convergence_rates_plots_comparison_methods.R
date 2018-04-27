#
# Plot integrated convergence for different cases
# 1) The 'real' analysis
# 2) The 'real' analysis without edge adjustment
# 3) Giving the same rate to all events with the same Mw
# We should also show the Bird convergence rates.
# This will illustrate how our adjustments are a 'good thing'
# 
# How to do this?
# 1) Load compute_rates_all_sources environment
#    - This gives us the 'uniform slip' event rates = 
#
#    event_rates = event_conditional_probabilities * 
#        (mw_rate_function(event_table$Mw - dMw/2) - 
#         mw_rate_function(event_table$Mw + dMw/2) )
#
# 2) The above event_conditional_probabilities involves edge correction,
#     but it is easy to compute the one without edge correction
#
#    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
#        event_table, 
#        conditional_probability_model = conditional_probability_model)    
#
# 3) The computation of event rates where 'equal Mw --> equal rate' is straightforward.
#
#
# For each of the above cases, I should be able to 'fairly easily' re-compute the integrated
# fault slip, and make a plot. 
#
# Hmm, what about segmentation? This complicates the above reasoning.
#
# Alternative: We know edge correction is 'just' applied to corresponding uniform
# slip events touching edges of the full sourcezone. We can easily detect these in the files,
# and do a crude back-adjustment (e.g. give them the same rate as their neighbouring events).
# This may be an easier way to make the plot?

load('compute_rates_all_sources_session.RData')

## Choose the source zone
src = 'sunda'
plot_xlim_default = c(100, 125) #NULL
plot_ylim_default = c(-15, -5)  #NULL
plot_arrow_scale_default = 60

# Unit source statistics
uss = bird2003_env$unit_source_tables[[src]]
#
# Note the following properties only apply to uniform slip cases. The variable slip cases are derived
# by partitioning the uniform slip cases over the corresponding variable slip cases.
#
# Event rates: version used in the PTHA (with edge correction + spatially variable convergence)
event_rates = source_envs[[src]]$event_rates
# Magnitude
event_Mw = source_envs[[src]]$event_table$Mw
# Unit source indices involved
event_index_string = source_envs[[src]]$event_table$event_index_string
# Slip (uniform)
slip = source_envs[[src]]$event_table$slip
# The edge inflation factor we used for the final PTHA
edge_multiplier_final = source_envs[[src]]$sourcepar$best_edge_multiplier$minimum

#
#
#

# Flag events which have a unit-source at the edge of the rupture
uss_alongstrike_range = range(uss$alongstrike_number)
eis = sapply(event_index_string, f<-function(x) as.numeric(strsplit(x, split="-")[[1]]), simplify=FALSE)
is_on_edge = unlist(lapply(eis, f<-function(x) any(uss$alongstrike_number[x] %in% uss_alongstrike_range)))

# Make the 'alternative' event rates
event_rates_no_edge_inflation = event_rates
event_rates_uniformly_spread = event_rates
unique_Mw = unique(event_Mw)
for(mm in unique_Mw){

    k = which(event_Mw == mm)
    total_rate_mm = sum(event_rates[k])
    if(total_rate_mm > 0){
        conditional_probability_event = event_rates[k]/total_rate_mm
    }else{
        conditional_probability_event = 0*event_rates[k]
    }

    # This would be the rate if we assigned all events with the same Mw the
    # same rate
    event_rates_uniformly_spread[k] = sum(total_rate_mm)/length(k)
 
    # This would be the rate if we didn't use the edge multiplier
    tmp = conditional_probability_event * (1/(1+is_on_edge[k]*edge_multiplier_final))
    if(sum(tmp) > 0) tmp = tmp/sum(tmp)
    event_rates_no_edge_inflation[k] = tmp * total_rate_mm

}

# Here is the 'data based' convergence rate, in m/year
bird_convergence_div = uss$bird_vel_div * 1/1000
bird_convergence_rl = uss$bird_vel_rl   * 1/1000
data_convergence = sqrt(bird_convergence_div**2 + bird_convergence_rl**2)

# Make some 'model based' convergence rate
conv_slip_final = uss[,1] * 0
conv_slip_no_edge_inflation = uss[,1] * 0
conv_slip_uniformly_spread = uss[,1] * 0
for(i in 1:length(event_rates)){
   k = as.numeric(strsplit(event_index_string[[i]], '-')[[1]])
   conv_slip_final[k] = conv_slip_final[k] + slip[i] * event_rates[i]
   conv_slip_no_edge_inflation[k] = conv_slip_no_edge_inflation[k] + slip[i]*event_rates_no_edge_inflation[i]
   conv_slip_uniformly_spread[k] = conv_slip_uniformly_spread[k] + slip[i]*event_rates_uniformly_spread[i]
}


#
# Make a plot comparing all the appproaches 
#
panel_plot<-function(uss, point_size_scale, arrow_size_scale, extra_arrow_size=plot_arrow_scale_default,
    plot_xlim = plot_xlim_default, plot_ylim = plot_ylim_default){

      plot(uss$lon_c, uss$lat_c, xlim=plot_xlim, ylim=plot_ylim,
          cex=point_size_scale, asp=1); grid()
      v = cbind(-cos(uss$strike/180*pi)*arrow_size_scale,
                 sin(uss$strike/180*pi)*arrow_size_scale)*extra_arrow_size + 1.0e-06
      arrows(uss$lon_c, uss$lat_c,
          uss$lon_c + v[,1], uss$lat_c + v[,2],
          length=0, col='red')
}

pdf(paste0('convergence_plot_comparison_methods_', src, '.pdf'), width=12, height=9)
par(mfrow=c(2,2))
panel_plot(uss, point_size_scale=1, arrow_size_scale=data_convergence)
title('Input convergence ')
panel_plot(uss, point_size_scale=1, arrow_size_scale=conv_slip_final)
title('Modelled convergence (coupled fraction) \n Spatially variable rate + edge rate inflation ')
panel_plot(uss, point_size_scale=1, arrow_size_scale=conv_slip_no_edge_inflation)
title('Modelled convergence (coupled fraction) \n Spatially variable rate')
panel_plot(uss, point_size_scale=1, arrow_size_scale=conv_slip_uniformly_spread)
title('Modelled convergence (coupled fraction) \n event_rate = f(Mw)')
dev.off()

