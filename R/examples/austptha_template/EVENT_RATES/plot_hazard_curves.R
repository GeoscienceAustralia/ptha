library(rptha)
config_env = new.env()
source('config.R', local=config_env)

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
plot_rates_at_a_station<-function(lon_p, lat_p, greens_law_adjust=FALSE, verbose=FALSE){

    # Find index of nearest gauge. FIXME: Do spherical coordinates here
    lon_p_x = rep(lon_p, length.out=length(lon))
    lat_p_x = rep(lat_p, length.out=length(lon))

    #site = which.min(abs(lon-lon_p) + abs(lat - lat_p))
    site = which.min(distHaversine(cbind(lon_p_x, lat_p_x), cbind(lon, lat)))

    if(verbose){
        print(c('Coordinates: ', round(lon[site], 4), round(lat[site], 4)))
    }

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

    #
    # Convenience function to plot a panel of rates
    #
    panel_rate_plot<-function(stage_seq, site_rates_uniform, rate_min, rate_max, 
        greens_adjust, site, titlep, greens_adjust_title){
    
        # Set up plot
        plot(range(stage_seq)*greens_adjust, c(rate_min, rate_max), col=0, log='xy',
            main=paste0(titlep, ' @(',
                round(lon[site], 2), ',', round(lat[site], 2), ',', round(elev[site], 2), ')', 
                greens_adjust_title),
            xlab='Peak stage (m)', ylab = 'Exceedance Rate (events/year)' )

        site_rates_uniform_sum = site_rates_uniform[[1]]*0

        # Add rate curves, and keep running sum of total rate
        for(i in 1:length(site_rates_uniform)){
            points(stage_seq*greens_adjust, site_rates_uniform[[i]], pch=19, t='o', col=i, cex=0.3)
            site_rates_uniform_sum = site_rates_uniform_sum + site_rates_uniform[[i]]
        }

        # Add total rate
        points(stage_seq*greens_adjust, site_rates_uniform_sum, t='l', lwd=2)

        # Extras
        grid()
        abline(h=c(1e-02, 1e-04, 1e-06), lty='dotted', col='grey')
        legend('topright', names(site_rates_uniform), col=1:length(site_rates_uniform), pch=19)

    }

    for(i in 1:length(rate_var_local_names)){

        panel_rate_plot(stage_seq, site_rates[[i]], rate_min, rate_max, greens_adjust, site, 
            titlep=rate_var_titles[i], greens_adjust_title)
    }
        
    return(invisible())

}

#' Boxplot of wave-heights for all events from a given source-zone, for all Mw,
#' at a station.
#'
#' @param lon_p, lat_p station coordinates
#' @param source_zone name of source_zone
#' @param slip_type either 'uniform' or 'stochastic'
#' @param plot_y_range y range of plot
#' @param boxwex Bar thickness 
#' @param ... further arguments to plot
#'
plot_wave_heights_at_a_station<-function(lon_p, lat_p, source_zone, slip_type = 'uniform',
    plot_y_range=c(1e-04, 1e+02), boxwex=1, ...){

    # Find index of nearest gauge. FIXME: Do spherical coordinates here
    lon_p_x = rep(lon_p, length.out=length(lon))
    lat_p_x = rep(lat_p, length.out=length(lon))
    site = which.min(distHaversine(cbind(lon_p_x, lat_p_x), cbind(lon, lat)))

    # Get the filename with max_stage
    if(slip_type == 'uniform'){

        nc_file_ind = grep(source_zone, config_env$all_source_uniform_slip_tsunami)

        if(length(nc_file_ind) != 1){
            stop(paste0('Could not find unique uniform slip file matching source_zone = ', source_zone))
        }

        nc_file = config_env$all_source_uniform_slip_tsunami[nc_file_ind]

    }else if(slip_type == 'stochastic'){

        nc_file_ind = grep(source_zone, config_env$all_source_stochastic_slip_tsunami)

        if(length(nc_file_ind) != 1){
            stop(paste0('Could not find unique stochastic slip file matching source_zone = ', source_zone))
        }

        nc_file = config_env$all_source_stochastic_slip_tsunami[nc_file_ind]

    }else{
        stop('unrecognized slip type')
    }

    fid = nc_open(nc_file)
    gauge_max_stage = ncvar_get(fid, 'max_stage', start=c(1,site), count=c(-1,1))
    event_Mw = ncvar_get(fid, 'event_Mw')
    # Deal with finite-precision netcdf limitations
    event_Mw = round(event_Mw, 3)

    stopifnot(length(gauge_max_stage) == length(event_Mw))

    # Make the plot
    plot(range(event_Mw), plot_y_range, col=0, log='y', ...)
    boxplot(gauge_max_stage ~ event_Mw, at=unique(event_Mw), boxwex=boxwex, add=TRUE)

    nc_close(fid)

}

