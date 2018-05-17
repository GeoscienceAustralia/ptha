#
# See if there are any 'general' patterns between DART-buoy GOF statistics and
# event properties. 
#
# For instance, we might imagine the 'good' events consistently have small area
# compared to the corresponding family, or low peak slip compared to the
# corresponding family, etc. If we find relationships like these, it
# might help us design better models, or filter out unrealistic events.
#
# If we develop a bias-adjustment method, key questions are:
# -- Does it still lead to good 'coverage statistics' for the 'corresponding family of model scenarios'?
#    -- Can test this by applying the method within the statistical tests code
# -- How does the hazard from the 'corrected' stochastic and variable-uniform slip models compare?
#    -- Can test this by applying the method at a few sites, before applying it fully
#
#
#

#' Extract statistics from the "corresponding family of model scenarios"
#'
#' This computes some summary statistics from objects made in
#' gauge_summary_statistics.R, that are lists (per dart buoy) containing lists
#' (per model scenario) with summary statistics
#'
#' Returns a data.frame with the goodness of fit statistic and other columns summarising
#' the rupture properties
#' 
#' @param gauge_stats object like 'stochastic_slip_stats' or 'uniform_slip_stats', etc, as created
#' by the script gauge_summary_statistics.R {in e.g. SOURCE_ZONES/TEMPLATE/TSUNAMI_EVENTS/plots/ }
#' @param unit_source_statistics the unit source statistics
#'
family_stats<-function(gauge_stats, unit_source_statistics){

    # Get time goodness-of-fit statistic for each model scenario
    gf_mat = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) x$model_data_similarity_time))
    # Convert from list of lists to matrix
    for(i in 1:length(gf_mat)) gf_mat[[i]] = unlist(gf_mat[[i]])
    gf_mat = matrix(unlist(gf_mat), ncol=length(gf_mat))
    # Use the 'median GF over all dart buoys' as our GOF value
    gf_median = apply(gf_mat, 1, median)

    # Get the peak wave height (median over all darts). This is a crude indicator 
    # of the tsunami size
    stage_range_mat = lapply(gauge_stats, 
        f<-function(x) lapply(x, f<-function(x) diff(x$model_range)))
    # Convert from list of lists to matrix
    for(i in 1:length(stage_range_mat)) stage_range_mat[[i]] = unlist(stage_range_mat[[i]])
    stage_range_mat = matrix(unlist(stage_range_mat), ncol=length(stage_range_mat))
    stage_range_median = apply(stage_range_mat, 1, median)

    # Get the peak slip for each model scenario
    if('event_slip_string' %in% names(gauge_stats[[1]][[1]]$events_with_Mw)){
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) max(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
        # Mean slip 
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) mean(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
    }else{
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
    }

    # Nonzero slip area 
    area = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            area = sum(unit_source_statistics$length[inds] * unit_source_statistics$width[inds])
            return(area)
        }))

    # Length -- we find the alongstrike range of the event
    length = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            alongstrike_range = range(unit_source_statistics$alongstrike_number[inds])
            #inds2 = which(unit_source_statistics$alongstrike_number >= alongstrike_range[1] & 
            #              unit_source_statistics$alongstrike_number <= alongstrike_range[2] & 
            #              unit_source_statistics$downdip_number == 1)
            #length = sum(unit_source_statistics$length[inds2])
            length =  diff(alongstrike_range) + 1
            return(length)
            }))
            
    # width -- we find the downdip range of the event
    width = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            downdip_range = range(unit_source_statistics$downdip_number[inds])
            #inds2 = which(unit_source_statistics$downdip_number >= downdip_range[1] & 
            #              unit_source_statistics$downdip_number <= downdip_range[2] & 
            #              unit_source_statistics$alongstrike_number == unit_source_statistics$alongstrike_number[inds[1]])
            #width = sum(unit_source_statistics$width[inds2])
            width = diff(downdip_range) + 1 
            return(width)
            }))


    area_including_zeros  = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            downdip_range = range(unit_source_statistics$downdip_number[inds])
            alongstrike_range = range(unit_source_statistics$alongstrike_number[inds])
            k = which(
                    unit_source_statistics$downdip_number >= downdip_range[1] &
                    unit_source_statistics$downdip_number <= downdip_range[2] &
                    unit_source_statistics$alongstrike_number >= alongstrike_range[1] &
                    unit_source_statistics$alongstrike_number <= alongstrike_range[2] )
            area = sum(unit_source_statistics$length[k] * unit_source_statistics$width[k])
            return(area)
        }))

    # Magnitude.
    if( 'Mw_variable_mu' %in% names(gauge_stats[[1]][[1]]$events_with_Mw) ){
        # Variable shear modulus case
        Mw = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$Mw_variable_mu))
    }else{
        # Constant shear modulus case
        Mw = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$Mw))
    }

    # The following variables need different treatments for the 'heterogeneous slip' case
    # vs both uniform slip cases
    if('physical_corner_wavenumber_x' %in% names(gauge_stats[[1]][[1]]$events_with_Mw)){
        # Variable slip

        # Corner wavenumbers
        corner_wavenumber_x = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) x$events_with_Mw$physical_corner_wavenumber_x))
        corner_wavenumber_y = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) x$events_with_Mw$physical_corner_wavenumber_y))

        # Peak slip location
        peak_slip_downdip = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) x$events_with_Mw$peak_slip_downdip_ind))
        peak_slip_alongstrike = unlist(lapply(gauge_stats[[1]], 
            f<-function(x) x$events_with_Mw$peak_slip_alongstrike_ind))

    }else{
        #
        # Both uniform slip models
        #
        corner_wavenumber_x = width*NA
        corner_wavenumber_y = width*NA

        # Peak slip location
        peak_slip_downdip = unlist(lapply(gauge_stats[[1]], f<-function(x){
                inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
                downdip_mean = mean(unit_source_statistics$downdip_number[inds])                 
                return(downdip_mean) 
        }))
        # Peak slip location
        peak_slip_alongstrike = unlist(lapply(gauge_stats[[1]], f<-function(x){
                inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
                alongstrike_mean = mean(unit_source_statistics$alongstrike_number[inds])                 
                return(alongstrike_mean) 
        }))

    }

    output = data.frame(gf = gf_median, 
        stage_range_median = stage_range_median,
        peak_slip = peak_slip_sum, 
        mean_slip = mean_slip_sum, 
        area = area,
        length=length, 
        width=width, 
        area_including_zeros = area_including_zeros,
        Mw = Mw,
        corner_wavenumber_x = corner_wavenumber_x, 
        corner_wavenumber_y = corner_wavenumber_y,
        peak_slip_downdip = peak_slip_downdip, 
        peak_slip_alongstrike = peak_slip_alongstrike)

    return(output)
}

#
#
# Main script here
#
#

variable_mu = FALSE
if(variable_mu){
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*varyMu.Rdata')
}else{
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*[0-9].Rdata')
}
print(paste0(' Remove the normal fault Mw 7.7 from Kermadec, and the 7.65, so that \n',
    ' we consistently examine thrust-like events that are above Mw 7.7 by GCMT. '))
all_Rdata = all_Rdata[-c(3,8)]

# Store the statistics in a list (one entry per historical event)
uniform_stat = vector(mode='list', length=length(all_Rdata))
stochastic_stat = vector(mode='list', length=length(all_Rdata))
variable_uniform_stat = vector(mode='list', length=length(all_Rdata))

for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    # Load the R session associated with the gauge_summary_statistics.R script
    # for the ith event
    load(all_Rdata[i], envir=event_env)

    # Main computation here
    stochastic_stat[[i]] = family_stats(event_env$stochastic_slip_stats, 
        event_env$unit_source_statistics)
    uniform_stat[[i]] = family_stats(event_env$uniform_slip_stats, 
        event_env$unit_source_statistics)
    variable_uniform_stat[[i]] = family_stats(event_env$variable_uniform_slip_stats, 
        event_env$unit_source_statistics)
}
event_env = new.env() # Clear the memory
# Put informative names on the lists
names(uniform_stat) = basename(all_Rdata)
names(stochastic_stat) = basename(all_Rdata)
names(variable_uniform_stat) = basename(all_Rdata)

# Store the results
if(variable_mu){
    save.image('event_properties_and_GOF_session_varyMu.Rdata')
}else{
    save.image('event_properties_and_GOF_session.Rdata')
}
#stop('Deliberate stop here')

#
# Useful plotting code below
#

# For a single event, plot the variables for all models (by value, or by rank),
# and highlight the 'top few (ng=...)' with 'best fit' in red. This helps us 'eyeball'
# relations, e.g. to notice if the good models always have a 'low' or 'high'
# value of the statistic
event_statistics_plot<-function(stat, rank_transform=FALSE){

    # Number of 'best fit' scenarios to highlight in red
    ng = 5

    stat = stat[rev(order(stat$gf)),] # Order to put the 'good' ones on top in the plot
    if(rank_transform){
        pairs(apply(stat, 2, rank), col=(rank(stat$gf) <= ng) + 1)
    }else{
        pairs(stat, col=(rank(stat$gf) <= ng) + 1)
    }
}

# Do a 'event_statistics_plot' type plot, for all events
pdf_events_statistics_plot<-function(stats, type=''){
    file_name = paste0('Rank_statistic_plots_', type, '.pdf')
    pdf(file_name, width=12, height=10)
    for(i in 1:length(stats)){
        par(oma=c(0,0,3,0))
        # Plot the raw data. Often harder to see the 'spread' because
        # points cluster around the distribution
        event_statistics_plot(stats[[i]])
        title(names(stats)[i], outer=TRUE)
        # Plot the data ranks. Easier to see if the good models are 'broadly uniformly'
        # drawn from the data, or not
        event_statistics_plot(stats[[i]], rank_transform=TRUE)
        title(names(stats)[i], outer=TRUE)
    }
    dev.off()
}

# Alternative plot showing statistics of 'good' events in comparison with
# 'all simulated' events
events_scaling_plot<-function(stats, title_extra=""){
   
    # Number of events we will call 'good' 
    ng = 5

    all_mw   = unlist(lapply(stats, f<-function(x) x$Mw))
    good_mw   = unlist(lapply(stats, f<-function(x) x$Mw[order(x$gf)[1:ng]]))

    all_vs_good<-function(var){
        all_var = unlist(lapply(stats, f<-function(x) x[[var]]))
        good_var = unlist(lapply(stats, f<-function(x) x[[var]][order(x$gf)[1:ng]]))

        plot(all_mw, all_var, log='y', main=paste0(var, title_extra), 
            col=rgb(0,0,0,alpha=0.05), pch=19, xlab="Mw", ylab=var)
        points(good_mw, good_var, col='red', pch=19)

        # Add quantiles of 'ambient' result, for visual reference
        quants = c(0.05, 0.2, 0.5, 0.8, 0.95)
        cols = c('purple', 'blue', 'green', 'blue', 'purple')
        labels = c('Middle 90%', 'Middle 60%', 'median')
        for(i in 1:length(quants)){
            mw_ambient = aggregate(all_var, by=list(round(all_mw,1)), 
                f<-function(x) quantile(x, p=quants[i]))
            LWD = 1#-2*abs(0.5-quants[i])**1.2 # line width
            points(mw_ambient[,1], mw_ambient[,2], t='l', col=cols[i], lwd=LWD)
        }
        legend('bottomright', c('Models', paste0('Top ', ng, ' Models'), labels), 
            col=c('black', 'red', cols[1:3]), lty=c(NA, NA, rep('solid', 3)), 
            pch=c(19, 19, rep(NA, 3)))
        
        return(invisible(environment()))
    }

    par(mfrow=c(2,3))
    area = all_vs_good('area')
    peak_slip = all_vs_good('peak_slip')
    mean_slip = all_vs_good('mean_slip')
    length = all_vs_good('length')
    width = all_vs_good('width')
    stage_range_median = all_vs_good('stage_range_median')

    return(invisible(environment()))
}

#
# Make the plot
#
if(variable_mu){
    pdf('event_properties_with_good_fitting_events_varyMu.pdf', width=12, height=8)
}else{
    pdf('event_properties_with_good_fitting_events.pdf', width=12, height=8)
}
events_scaling_plot(stochastic_stat, title_extra=' heterogeneous slip')
events_scaling_plot(variable_uniform_stat, title_extra=' variable_uniform slip')
events_scaling_plot(uniform_stat, title_extra=' fixed_uniform slip')
dev.off()

#
# Get the rank of the statistic for the n-'best' events in terms of the other
# events in the corresponding family of model scenarios
#
best_event_quantiles<-function(stats, nbest=5){

    # Store a list (one entry per variable), each containing
    # an array with one row for each event, and one column for
    # each value of nbest
    myvar = names(stats[[1]])
    nvar = length(myvar) # number of variables we look at
    output = vector(mode='list', length=nvar)
    names(output) = myvar
    # Prepare data structure
    for(i in 1:length(myvar)){
        output[[i]] = matrix(NA, ncol=nbest, nrow=length(stats))
        rownames(output[[i]]) = names(stats)
    }

    # Populate data structure
    for(j in 1:length(stats)){ # Every event
        for(i in 1:nvar){ # Every variable
            for(k in 1:nbest){ # 1st best, 2nd best, ... nbest best.
                # Get the i'th variable for the j'th event
                var_of_interest = stats[[j]][[myvar[i]]]
                # Find the one with goodness-of-fit rank = k (rank=1 is best-fit)
                gof_value = stats[[j]][['gf']]
                eoi = which(rank(gof_value, ties='first') == k)# Beware ties treatment
                # Find the fraction of var_of_interest that are < the value
                # associated with 'eoi'
                empirical_fraction_less_than = 
                    sum(var_of_interest <= var_of_interest[eoi])/(length(var_of_interest)+1)
                #if(is.na(empirical_fraction_less_than)) browser()
                output[[i]][j,k] = empirical_fraction_less_than
            }
        }
    }

    return(output)
}

#
# Find how the statistics of 'good' events are distributed, compared
# with their 'corresponding family of model scenarios'
#
stochastic_best_event_quantiles = best_event_quantiles(stochastic_stat, nbest=5)
variable_uniform_best_event_quantiles = best_event_quantiles(variable_uniform_stat, nbest=5)
uniform_best_event_quantiles = best_event_quantiles(uniform_stat, nbest=5)

#
# Estimate a possible bias correction for the models, based on
# the summary_statistic 'var'
#
quantile_adjuster<-function(stochastic_best_event_quantiles, var, colind = 'mean', title_start=""){
    library(ADGofTest)

    # Sort the quantiles for best, 2nd best, 3rd best, ...
    sorted_var = apply(stochastic_best_event_quantiles[[var]], 2, sort)

    # Get mean of sorted rows. This gives a sense of the quantile-match between
    # the 'good' events and 'all events'
    if(colind == 'mean'){
        mean_of_sorted = apply(sorted_var, 1, mean)
    }else{
        mean_of_sorted = sorted_var[,colind]
    }
    #stopifnot(all(diff(mean_of_sorted) > 0))

    # For an unbiased model, mean_of_sorted should be reasonably close to the
    # following
    ideal = seq(1,length(mean_of_sorted))/(length(mean_of_sorted)+1) 

    # We can fit a cubic that follows this curve. Note that the cubic
    # f(x) = a*x + b*x^2 + (1-a-b)*x^3
    # will always pass through 0 and 1. We can compute a/b by minimisation.
    #
    x = c(0, ideal, 1)
    y = c(0, mean_of_sorted, 1)
    cubic_01<-function(a,x){
        if(a[1] < 0) return(rep(9999, length(x)))
        a[1]*x + a[2]*x^2 + (1-a[1]-a[2])*x^3
    }
    cubic_01_deriv<-function(a, x) a[1] + 2*a[2]*x + 3*(1-a[1]-a[2])*x^2

    # We can also fit a '2-piece linear' function, as an alternative to the cubic
    bi_linear<-function(a, x){ 
        if(any(a>1 | a<0)) return(rep(99999, length(x)))
        # Two-piece linear fit joining c(0,0), c(a[1], a[2]), c(1,1)
        g1 = (a[2]-0)/(a[1]-0)
        g2 = (1 - a[2])/(1 - a[1])
        output = g1 * (x - 0)*(x <= a[1]) + (g2 * (x-1) + 1) * (x > a[1])
        return(output)
    }
    bi_linear_deriv<-function(a,x){
        if(any(a>1 | a<0)) return(rep(99999, length(x)))
        # Two-piece linear fit joining c(0,0), c(a[1], a[2]), c(1,1)
        g1 = (a[2]-0)/(a[1]-0)
        g2 = (1 - a[2])/(1 - a[1])
        
        g1 * (x<=a[1]) + g2*(x>a[1])
    }
   
    # Find the best parameters 
    best_parameters = optim(c(1, 0), f<-function(a) sum((cubic_01(a,x) - y)^2))
    best_parameters_bilinear = optim(c(0.5, 0.5), f<-function(a) sum((bi_linear(a,x) - y)^2))

    plot(x, y, xlim=c(0,1), ylim=c(0,1), 
        ylab='Q good relative to all ',
        xlab='Q good relative to good')
    extra_title=""
    # If we are only using a single column, add a p-value for uniformity
    if(is.numeric(colind)) extra_title = signif(ad.test(mean_of_sorted)$p.value, 3)
    title(paste0(title_start, var, ' ', extra_title))
    xl = seq(-0.001, 1.001, len=1001)
    cubic_vals = cubic_01(best_parameters$par, xl) 
    points(xl, cubic_vals, t='l', col='red')
    grid(); abline(0,1,col='green')
    bilin_vals = bi_linear(best_parameters_bilinear$par, xl) 
    points(xl, bilin_vals, t='l', col='purple')

    #
    # If the old event rates are uniform (in a group with events having quantile E_q),
    # then the new rates are bias_adjustment_factor(E_q)/sum(bias_adjustment_factor(E_q))
    # 
    bias_adjustment_factor<-function(x){
        (approx(cubic_vals, xl, xout=x)$y - approx(cubic_vals, xl, xout=x-1.0e-06)$y)/1.0e-06
    }
    bias_adjustment_factor_bilinear<-function(x){
        (approx(bilin_vals, xl, xout=x)$y - approx(bilin_vals, xl, xout=x-1.0e-06)$y)/1.0e-06
    }

    return(invisible(environment()))
}

#
# Plot of the quantile adjustments
#
if(variable_mu){
    pdf('quantile_adjustment_variable_mu.pdf', width=10, height=9)
}else{
    pdf('quantile_adjustment.pdf', width=10, height=9)
}

stochastic_quantile_adjust = list()
variable_uniform_quantile_adjust = list()
uniform_quantile_adjust = list()
for(var in c('mean_slip', 'area', 'area_including_zeros', 'peak_slip')){
    par(mfrow=c(3,4))
    stochastic_quantile_adjust[[var]] = quantile_adjuster(stochastic_best_event_quantiles, var, title_start='stoc ')
    for(i in 1:3) quantile_adjuster(stochastic_best_event_quantiles, var, colind=i, title_start='stoc ')
    variable_uniform_quantile_adjust[[var]] = quantile_adjuster(variable_uniform_best_event_quantiles, var, title_start='VU ')
    for(i in 1:3) quantile_adjuster(variable_uniform_best_event_quantiles, var, colind=i, title_start='VU ')
    uniform_quantile_adjust[[var]] = quantile_adjuster(uniform_best_event_quantiles, var, title_start='U ')
    for(i in 1:3) quantile_adjuster(uniform_best_event_quantiles, var, colind=i, title_start='U ')
}
dev.off()

if(variable_mu){
    save.image('event_properties_and_GOF_session_varyMu_end.Rdata')
}else{
    save.image('event_properties_and_GOF_session_end.Rdata')
}
