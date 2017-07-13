library(rptha)

# 
rate_files = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/tsunami_*.nc')

# Get all the rates
rates = lapply(as.list(rate_files), f<-function(x) nc_open(x))
names(rates) = basename(dirname(dirname(rate_files))) 


# Other key variables
lon = ncvar_get(rates[[1]], 'lon')
lat = ncvar_get(rates[[1]], 'lat')
elev = ncvar_get(rates[[1]], 'elev')
stage_seq = ncvar_get(rates[[1]], 'stage')

#
#
#
plot_rates_at_a_station<-function(lon_p, lat_p, greens_law_adjust=FALSE){

    # Find index of nearest gauge. FIXME: Do spherical coordinates here
    site = which.min(abs(lon-lon_p) + abs(lat - lat_p))

    #
    # Get uniform and stochastic rates in a list
    #
    # Access will be like, e.g.
    #     site_rates$uniform$puysegur, site_rates$stochastic_lower_ci$southamerica

    # Names of variables in netcdf file
    rate_var_nc_names = c(
        'uniform_slip_rate', 
        'stochastic_slip_rate', 
        'uniform_slip_rate_lower_ci', 
        'stochastic_slip_rate_lower_ci', 
        'uniform_slip_rate_upper_ci', 
        'stochastic_slip_rate_upper_ci') 
    # Names corresponding to rate_var_nc_names, which we will use in the site_rates list
    rate_var_local_names = c(
        'uniform', 
        'stochastic', 
        'uniform_lower_ci', 
        'stochastic_lower_ci', 
        'uniform_upper_ci', 
        'stochastic_upper_ci') 
    # Titles corresponding to plots of the above variables 
    rate_var_titles = c(
        'Uniform slip tsunami peak stage exceedance rates', 
        'Stochastic slip tsunami peak stage exceedance rates', 
        'Uniform slip tsunami peak stage exceedance rates, lower CI', 
        'Stochastic slip tsunami peak stage exceedance rates, lower CI', 
        'Uniform slip tsunami peak stage exceedance rates, upper CI', 
        'Stochastic slip tsunami peak stage exceedance rates, upper CI')

    # Read into the data structure
    site_rates = vector(mode='list', length=length(rate_var_nc_names))

    for(i in 1:length(rate_var_nc_names)){

        site_rates[[i]] = vector(mode='list', length=length(rates))
        names(site_rates)[i] = rate_var_local_names[i]

        for(j in 1:length(rates)){
            x = rates[[j]]
            site_rates[[i]][[j]] = ncvar_get(x, rate_var_nc_names[i], start=c(1,site), count=c(-1,1))
            names(site_rates[[i]])[j] = names(rates)[j]
        }
    }

    # Plot parameters
    par(mfrow=c(3,2))
    rate_min = 1e-07
    rate_max = 3*max(unlist(lapply(site_rates$uniform, max)))
    rate_max = 3*max(rate_max, max(unlist(lapply(site_rates$stochastic, max))))

    # Scaling factor to allow for greens-law adjustment to wave height
    greens_adjust = 1 * (1-greens_law_adjust) + ( max(-elev[site],0)**0.25 )*(greens_law_adjust)
    greens_adjust_title = c('', '\n translated to 1m depth with greens law')[greens_law_adjust+1]

    # Convenience function to plot a panel of rates
    panel_rate_plot<-function(stage_seq, site_rates_uniform, rate_min, rate_max, greens_adjust, site, titlep, greens_adjust_title){

        plot(range(stage_seq)*greens_adjust, c(rate_min, rate_max), col=0, log='xy',
            main=paste0(titlep, ' @(',
                round(lon[site], 2), ',', round(lat[site], 2), ',', round(elev[site], 2), ')', 
                greens_adjust_title),
            xlab='Peak stage (m)', ylab = 'Exceedance Rate (events/year)' )
        site_rates_uniform_sum = site_rates_uniform[[1]]*0
        for(i in 1:length(site_rates_uniform)){
            points(stage_seq*greens_adjust, site_rates_uniform[[i]], pch=19, t='o', col=i, cex=0.3)
            site_rates_uniform_sum = site_rates_uniform_sum + site_rates_uniform[[i]]
        }
        points(stage_seq*greens_adjust, site_rates_uniform_sum, t='l', lwd=2)
        grid()
        abline(h=c(1e-02, 1e-04, 1e-06), lty='dotted', col='grey')
        legend('topright', names(site_rates_uniform), col=1:8, pch=19)

    }

    for(i in 1:length(rate_var_local_names)){

        panel_rate_plot(stage_seq, site_rates[[i]], rate_min, rate_max, greens_adjust, site, 
            titlep=rate_var_titles[i], greens_adjust_title)
    }
        
    return(invisible())

}

