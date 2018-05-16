library(rptha)

# Read key user variables
config = new.env()
source('config.R', local=config)

# Read key sourcezone parameters
sourcezone_parameter_file = config$sourcezone_parameter_file 
sourcezone_parameters = read.csv(sourcezone_parameter_file, stringsAsFactors=FALSE)

#
# Utility function to search a vector of file paths ('files') for
# a single file with a particular base-name ('file')
#
match_file<-function(files, file){
    ind = grep(file, basename(files), fixed=TRUE)
    if(length(ind) != 1) stop(paste0('Could not match unique file like ', file))
    return(files[ind])
}

##
## Variable shear modulus function
##
#mu_fun<-function(dd){
#    # Point values based on plot of Lay and Bilek (2007)
#    # See
#    #/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DATA/EARTHQUAKE/Shear_modulus
#    depths = c(0, 7.5, 15, 35, 9999)
#    mu = c(10, 10, 30, 67, 67)*1e+09
#    output = 10**(approx(depths, log10(mu), xout=dd)$y)
#    return(output)
#}

#
# Plot function
#
plot_sourcezone_rate_curve_with_fixed_and_variable_mu<-function(sourcezone){

    sourcezone_events = list()

    # 
    # Get the stochastic slip events (without tsunami) 
    #
    stochastic_file = match_file(
        gsub('_tsunami', '', config$all_source_stochastic_slip_tsunami, fixed=TRUE),
        paste0('all_stochastic_slip_earthquake_events_', sourcezone, '.nc'))
    stochastic_eq_file = gsub('_tsunami', '', stochastic_file, fixed=TRUE)
    sourcezone_events$events = read_table_from_netcdf(stochastic_eq_file)

    # Get the unit source statistics
    uss_file = match_file(config$unit_source_statistics_netcdf_files,
        paste0('unit_source_statistics_', sourcezone, '.nc'))
    sourcezone_events$unit_source_statistics = read_table_from_netcdf(uss_file)

    # Event indices
    event_inds = lapply(sourcezone_events$events$event_index_string, 
        f<-function(x) as.numeric(strsplit(x, '-')[[1]]) )
    event_slip = lapply(sourcezone_events$events$event_slip_string, 
        f<-function(x) as.numeric(strsplit(x, '_')[[1]]) )


    #
    # Extract event properties from the file
    #
    depth = sourcezone_events$unit_source_statistics$depth
    mu_constant = depth * 0 + 3e+10
    #mu_variable = mu_fun(depth)
    mu_variable = shear_modulus_depth(depth)
    area = sourcezone_events$unit_source_statistics$length * 
        sourcezone_events$unit_source_statistics$width
    # Pre-allocate memory and populate in loop
    moment1 = rep(0, length(event_inds)) 
    moment2 = rep(0, length(event_slip))
    nonzero_rupture_area = rep(0, length(event_slip))
    full_rupture_area = rep(0, length(event_slip))
    area_80pc = rep(0, length(event_slip))
    area_50pc= rep(0, length(event_slip))
    peak_slip = rep(0, length(event_slip))
    for(i in 1:length(moment1)){
        inds = event_inds[[i]]
        slip = event_slip[[i]]
        moment1[i] = sum(slip * area[inds] * mu_constant[inds] * 1e+06)
        moment2[i] = sum(slip * area[inds] * mu_variable[inds] * 1e+06)
        nonzero_rupture_area[i] = sum(area[inds])

        ddi = sourcezone_events$unit_source_statistics$downdip_number[inds]
        asi = sourcezone_events$unit_source_statistics$alongstrike_number[inds]
        full_inds = which(
            (sourcezone_events$unit_source_statistics$downdip_number >= min(ddi) ) &
            (sourcezone_events$unit_source_statistics$downdip_number <= max(ddi) ) &
            (sourcezone_events$unit_source_statistics$alongstrike_number >= min(asi) ) &
            (sourcezone_events$unit_source_statistics$alongstrike_number <= max(asi) ) )

        full_rupture_area[i] = sum(area[full_inds])

        # Area with 80% of the moment
        slip_order = order(slip, decreasing=TRUE)
        slip_x_area_cumprop = cumsum( (slip * area[inds])[slip_order])/sum(slip * area[inds])
        n = sum(slip_x_area_cumprop < 0.8) + 1
        area_80pc[i] = sum((area[inds][slip_order])[1:n])
        # Area with 50% 
        n = sum(slip_x_area_cumprop < 0.5) + 1
        area_50pc[i] = sum((area[inds][slip_order])[1:n])

        peak_slip[i] = max(slip)
    }


    # If mu is variable, we can consider it a relabeling of magnitude.
    Mw_mu_fixed = round(M0_2_Mw(moment1), 1) # Rounding is 'almost' not needed, just accounting for numerical imprecision in file storage 
    Mw_mu_vary = M0_2_Mw(moment2)

    # Here we compute the exceedance rate
    Mw_bin_size = 0.1
    # Note our fixed-mu Mws range 7.2, 7.3, ...
    # These rates can be interpreted as representing a 'bin' of events
    mws = seq(7.15, 9.95, by=Mw_bin_size) 
    rate_mu_fixed = sapply(mws, f<-function(x) sum(sourcezone_events$events$rate_annual * (Mw_mu_fixed >= x)))
    rate_mu_fixed_lower = sapply(mws, f<-function(x) sum(sourcezone_events$events$rate_annual_lower_ci * (Mw_mu_fixed >= x)))
    rate_mu_fixed_upper = sapply(mws, f<-function(x) sum(sourcezone_events$events$rate_annual_upper_ci * (Mw_mu_fixed >= x)))
    rate_mu_vary = sapply(mws, f<-function(x) sum(sourcezone_events$events$rate_annual * (Mw_mu_vary >= x)))

    plot(mws, pmax(rate_mu_fixed, 1.0e-100), t='o', log='y', ylim=c(1.0e-06, 1), xlab="Mw", ylab='Exceedance rate')
    points(mws, pmax(rate_mu_vary, 1.0e-100), t='l', col='red')
    points(mws, rate_mu_fixed_lower, t='l', col='orange', lty='dashed')
    points(mws, rate_mu_fixed_upper, t='l', col='orange', lty='dashed')
    grid()
    title(paste0('Rates with fixed and variable mu, ', sourcezone))

    abline(h=c(1/500, 1/1000, 1/2500), col='brown', lty='dotted')

    # Return the key data as output
    output = data.frame(mws=mws, rate_mu_fixed=rate_mu_fixed, rate_mu_vary=rate_mu_vary)
    return(output)
}


#
# Main run code
#
all_source_zones = sort(unique(sourcezone_parameters$sourcename))
rate_curves_store = vector(mode='list', length=length(all_source_zones))
names(rate_curves_store) = all_source_zones

pdf('Mw_rate_curves_with_and_without_mu_variation.pdf', width=10, height=8)
for(i in 1:length(all_source_zones)){

    # Only work on subduction and thrust type zones, because we use a different
    # shear modulus treatment elsewhere.
    k = min(which(sourcezone_parameters$sourcename == all_source_zones[i]))
    if(sourcezone_parameters$scaling_relation[k] != 'Strasser') next

    print(all_source_zones[i])

    # Make the panel plot, and return the mw-rate results
    tmp = plot_sourcezone_rate_curve_with_fixed_and_variable_mu(all_source_zones[i])

    # Store a globally integrated rate
    if(i == 1){
        global_rates = tmp
    }else{
        global_rates[,2:3] = global_rates[,2:3] + tmp[,2:3]
    }

    # Store the results
    rate_curves_store[[i]] = tmp
}
#
# Plot the globally integrated rate
#
plot(global_rates[,1], pmax(global_rates[,2], 1e-100), t='o', log='y', 
    ylim=c(1e-04,10), xlab="Mw", ylab='Exceedance rate')
points(global_rates[,1], pmax(global_rates[,3], 1e-100), t='l', col='red')
grid()
title('Globally integrated rates')

dev.off()

save.image('variable_mu_checks_session.Rdata')

# Look at how much mw changes at different return periods as a result of varying mu
# We generally see a small shift that doesn't change much with return period, but
# varies by source-zone.
mw_change = lapply(rate_curves_store, 
    f<-function(x){ try(
        # Note that x[,2] and x[,3] are rates, which are decreasing but might
        # finish with a collection of zeros. In that case, we want to use the 
        # minimum magnitude corresponding to rate=0 for interpolation. The ties='min'
        # argument takes care of that.
        approx(x[,2], x[,1], xout=c(1/100, 1/500, 1/1000, 1/2500), ties='min')$y - 
        approx(x[,3], x[,1], xout=c(1/100, 1/500, 1/1000, 1/2500), ties='min')$y
    )} )

