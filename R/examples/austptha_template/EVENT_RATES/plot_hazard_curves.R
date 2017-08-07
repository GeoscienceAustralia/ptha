library(rptha)
config_env = new.env()
source('config.R', local=config_env)

# Files with stage-vs-exceedance-rate for every hazard point.
# (We can't put these locations in config.R because they are created
# by scripts in this folder) 
rate_files = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/tsunami_*.nc')
# Get all the rates
rates = lapply(as.list(rate_files), f<-function(x) nc_open(x))
names(rates) = basename(dirname(dirname(rate_files))) 

# Check rate files are ordered by source-zone in the same way as key files in
# 'config'
stopifnot(all(names(rates) == config_env$source_names_1))

# Other key variables
# Gauge location information
lon = ncvar_get(rates[[1]], 'lon')
lat = ncvar_get(rates[[1]], 'lat')
elev = ncvar_get(rates[[1]], 'elev')
# Sequence of stages at which exceedance rates are calculated,
# for each station
stage_seq = ncvar_get(rates[[1]], 'stage') 

#' Find index of gauge nearest a given point
#'
get_station_index<-function(lon_p, lat_p){

    #lon_p_x = rep(lon_p, length.out=length(lon))
    #lat_p_x = rep(lat_p, length.out=length(lon))
    #site = which.min(distHaversine(cbind(lon_p_x, lat_p_x), cbind(lon, lat)))

    site = lonlat_nearest_neigbours(cbind(lon_p, lat_p), cbind(lon, lat))

    return(site)
}

#
# Make multi-panel stage-vs-exceedance-rate plots for a station (for both
# uniform and stochastic slip, and upper and lower credible intervals)
#
plot_rates_at_a_station<-function(lon_p, lat_p, greens_law_adjust=FALSE, verbose=FALSE){

    # Find index of nearest gauge.
    site = get_station_index(lon_p, lat_p)

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
            site_rates[[i]][[j]] = ncvar_get(x, rate_var_nc_names[i], 
                start=c(1,site), count=c(-1,1))
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
                round(lon[site], 2), ',', round(lat[site], 2), ',', 
                round(elev[site], 2), ')', greens_adjust_title),
            xlab='Peak stage (m)', ylab = 'Exceedance Rate (events/year)' )

        site_rates_uniform_sum = site_rates_uniform[[1]]*0

        # Add rate curves, and keep running sum of total rate
        for(i in 1:length(site_rates_uniform)){
            points(stage_seq*greens_adjust, site_rates_uniform[[i]], pch=19, 
                t='o', col=i, cex=0.3)
            site_rates_uniform_sum = site_rates_uniform_sum + site_rates_uniform[[i]]
        }

        # Add total rate
        points(stage_seq*greens_adjust, site_rates_uniform_sum, t='l', lwd=2)

        # Extras
        grid()
        abline(h=c(1e-02, 1e-04, 1e-06), lty='dotted', col='grey')
        legend('topright', names(site_rates_uniform), 
            col=1:length(site_rates_uniform), pch=19)

    }

    for(i in 1:length(rate_var_local_names)){

        panel_rate_plot(stage_seq, site_rates[[i]], rate_min, rate_max, 
            greens_adjust, site, titlep=rate_var_titles[i], greens_adjust_title)
    }
        
    return(invisible(site))

}

#' Plot of wave-heights for all events from a given source-zone, for all Mw,
#' at a station.
#'
#' Note that the plot made by this function changes greatly if
#' split_into_subsets is not NULL! If the latter is NULL, then we box-plot stage
#' by Mw, and add a plot of the Mw-vs-rate curve. However, if split_into_subsets
#' is an integer, we instead plot stage-vs-exceedance rate the given number of
#' times, each time with a distinct subset of the data. The idea is that if we
#' have made 'enough' stochastic slip scenarios, then these rate curves should
#' not significantly differ [say with 2 subsets]
#'
#' @param lon_p, lat_p station coordinates
#' @param source_zone name of source_zone
#' @param slip_type either 'uniform' or 'stochastic'
#' @param plot_y_range y range of plot
#' @param boxwex Bar thickness 
#' @param site_index provide the index of lon_p/lat_p in lon/lat, thereby
#' avoiding the nearest-neighbour search
#' @param split_into_subsets If not NULL, split the events into the given number
#' of subsets before plotting the stage-vs-exceedance-rate curve -- to help
#' assess whether we have enough events for convergence of the rate curve.
#' @param ... further arguments to plot
#'
plot_wave_heights_at_a_station<-function(lon_p, lat_p, source_zone, 
    slip_type = 'uniform', plot_y_range=c(1e-04, 1e+02), boxwex=0.1, 
    site_index = NULL, split_into_subsets=NULL, ...){

    if(is.null(site_index)){
        # Find index of nearest gauge.
        site = get_station_index(lon_p, lat_p)
    }else{
        site = site_index
    }

    site_lat = lat[site]
    site_lon = lon[site]
    site_elev = elev[site]
  
    print(paste0('Station info: ', site_lon, ', ', site_lat, ', ', site_elev, 
        ', ', site))

    # Get the filename with max_stage
    if(slip_type == 'uniform'){

        nc_file_ind = grep(source_zone, 
            config_env$all_source_uniform_slip_tsunami)

        if(length(nc_file_ind) != 1){
            stop(paste0(
                'Could not find unique uniform slip file matching source_zone = ', 
                source_zone))
        }

        nc_file = config_env$all_source_uniform_slip_tsunami[nc_file_ind]

    }else if(slip_type == 'stochastic'){

        nc_file_ind = grep(source_zone, 
            config_env$all_source_stochastic_slip_tsunami)

        if(length(nc_file_ind) != 1){
            stop(paste0(
                'Could not find unique stochastic slip file matching source_zone = ', 
                source_zone))
        }

        nc_file = config_env$all_source_stochastic_slip_tsunami[nc_file_ind]

    }else{
        stop('unrecognized slip type')
    }

    fid = nc_open(nc_file)
    gauge_max_stage = ncvar_get(fid, 'max_stage', start=c(1,site), count=c(-1,1))
    event_Mw = ncvar_get(fid, 'event_Mw')
    event_nominal_rate = ncvar_get(fid, 'event_rate_annual')
    event_nominal_rate_upper = ncvar_get(fid, 'event_rate_annual_upper_ci')
    event_nominal_rate_lower = ncvar_get(fid, 'event_rate_annual_lower_ci')
    nc_close(fid)
    # Deal with finite-precision netcdf limitations
    event_Mw = round(event_Mw, 3)

    # When plotting rates, multiply by this constant beforehand. 
    rate_rescale = 100
    if(is.null(split_into_subsets)){
        #
        # Boxplot stage by Mw
        #

        unique_Mws = sort(unique(event_Mw))

        # Exceedance rates. For plotting, set rate to a very small number, so
        # the line dips to zero
        Mw_exceedance_rate = pmax(1e-12, sapply(unique_Mws, 
            f<-function(x) sum(event_nominal_rate*(event_Mw >= x)) ))
        Mw_exceedance_rate_upper = pmax(1e-12, sapply(unique_Mws, 
            f<-function(x) sum(event_nominal_rate_upper*(event_Mw >= x)) ))
        Mw_exceedance_rate_lower = pmax(1e-12, sapply(unique_Mws, 
            f<-function(x) sum(event_nominal_rate_lower*(event_Mw >= x)) ))

        stopifnot(length(gauge_max_stage) == length(event_Mw))

        # Make the plot
        plot(range(event_Mw), plot_y_range, col=0, log='y', axes=FALSE, 
            frame.plot=TRUE, xlab="Mw", ylab="", ...)
        boxplot(gauge_max_stage ~ event_Mw, at=unique(event_Mw), boxwex=boxwex, 
            add=TRUE)

        points(unique_Mws, Mw_exceedance_rate * rate_rescale, t='l', col='red')
        points(unique_Mws, Mw_exceedance_rate_upper * rate_rescale, t='l', 
            col='blue')
        points(unique_Mws, Mw_exceedance_rate_lower * rate_rescale, t='l', 
            col='blue')

        grid()
        abline(h=5*10**(seq(-4, 2)), lty='dotted', col='grey')
        abline(h=10**(seq(-4, 2)), lty='dotted', col='orange')
        mtext('Stage (m) and "Mw exceedance rate per century"', side=2, line=1.8)
        title(paste0(source_zone,', ', slip_type, 
            ' slip: Mw vs stage (m) AND Mw vs exceedance rates \n @ ', 
            round(site_lon, 3), ', ', round(site_lat, 3), ', ', 
            round(site_elev, 3)))
    }else{
        #
        # Plot stage vs exceedance rate, in subsets
        #
        for(i in 1:split_into_subsets){

            subset = seq(i, length(event_Mw), by=split_into_subsets)

            # Subset stages and nominal rates
            gauge_max_stage_subset = gauge_max_stage[subset]
            event_nominal_rate_subset = event_nominal_rate[subset] * 
                split_into_subsets * rate_rescale
            event_nominal_rate_upper_subset = event_nominal_rate_upper[subset] * 
                split_into_subsets * rate_rescale
            event_nominal_rate_lower_subset = event_nominal_rate_lower[subset] * 
                split_into_subsets * rate_rescale

            # Get stage-vs-exceedance rate, for the subset
            stage_seq = 10**(seq(-2,1, len=100)) #seq(0.01, max(gauge_max_stage_subset), len=100)
            stage_nominal_exceed = sapply(stage_seq, 
                f<-function(x) sum(event_nominal_rate_subset * (gauge_max_stage_subset > x)))
            stage_nominal_upper_exceed = sapply(stage_seq, 
                f<-function(x) sum(event_nominal_rate_upper_subset * (gauge_max_stage_subset > x)))
            stage_nominal_lower_exceed = sapply(stage_seq, 
                f<-function(x) sum(event_nominal_rate_lower_subset * (gauge_max_stage_subset > x)))

            # Plotting 
            if(i == 1){           
                plot(stage_seq, stage_nominal_exceed, t='l', ylim=plot_y_range, 
                    log='xy', xlab='Stage (m)', ylab='Exceedance per century', 
                    col=i)
                points(stage_seq, stage_nominal_upper_exceed, t='l', col = i, 
                    lty='dotted')
                points(stage_seq, stage_nominal_lower_exceed, t='l', col = i, 
                    lty='dotted')
            }else{
                points(stage_seq, stage_nominal_exceed, t='l', col = i)
                points(stage_seq, stage_nominal_upper_exceed, t='l', col = i, 
                    lty='dotted')
                points(stage_seq, stage_nominal_lower_exceed, t='l', col = i, 
                    lty='dotted')
            }
            
        }

    }


    return(invisible(site))

}
# plot_wave_heights_at_a_station(lon_p, lat_p, source_zone='southamerica', slip_type = 'stochastic', boxwex=0.1)
# plot_wave_heights_at_a_station(lon_p, lat_p, source_zone='southamerica', slip_type = 'uniform', boxwex=0.1)


#
# Hazard deaggregation plot
#
get_station_deaggregated_hazard<-function(lon_p, lat_p, station_name = "", 
    slip_type = 'uniform', exceedance_rate = NULL, stage = NULL){

    # Check input args
    if(is.null(exceedance_rate) & is.null(stage)){
        stop('Must provide either stage or exceedance rate')
    }

    site_index = get_station_index(lon_p, lat_p)
    source_names = names(rates)

    # Get tsunami files
    if(slip_type == 'uniform'){
        all_source_tsunami = config_env$all_source_uniform_slip_tsunami
    }else if(slip_type == 'stochastic'){
        all_source_tsunami = config_env$all_source_stochastic_slip_tsunami
    }else{
        stop('slip_type has invalid value -- must be either "uniform" or "stochastic"')
    }


    # Only one of stage/exceedance rate can be provided 
    # If exceedance rate is provided, then use it to compute the stage
    if(!is.null(exceedance_rate)){
        if(!is.null(stage)) stop('Cannot provide both exceedance rate and stage')

        # Sum the exceedance rates over all source-zones for the chosen site
        if(slip_type == 'uniform'){
            varname = 'uniform_slip_rate'
        }else{
            varname = 'stochastic_slip_rate'
        }
        rate_sum = ncvar_get(rates[[1]], varname, 
            start=c(1,site_index), count=c(-1,1))
        for(i in 2:length(rates)){
            rate_sum = rate_sum + ncvar_get(rates[[i]], varname, 
                start=c(1,site_index), count=c(-1,1))
        }

        stage_exceed = approx(rate_sum, stage_seq, xout=exceedance_rate)$y
    }
    # If exceedance rate is not provided, then stage must be
    if(!is.null(stage)){
        stage_exceed = stage
    }

    # Get unit-source locations & summary statistics
    sources_list = list()
    for(i in 1:length(source_names)){
        sources_list[[i]] = list()
        sources_list[[i]]$unit_source_statistics = read_table_from_netcdf(
            config_env$unit_source_statistics_netcdf_files[i])
        # Check that subfault numbers are sorted from 1 to max, since our
        # algorithm relies on this
        nr = nrow(sources_list[[i]]$unit_source_statistics)
        stopifnot(all(sources_list[[i]]$unit_source_statistics$us_subfault_number == 1:nr))
    }
    names(sources_list) = source_names

    # Get the contribution for each source-zone
    for(i in 1:length(source_names)){
        # One value for every unit source
        nr = nrow(sources_list[[i]]$unit_source_statistics)
        sources_list[[i]]$contribution = rep(0, length=nr)
        sources_list[[i]]$stage_exceed = stage_exceed
        sources_list[[i]]$station_location = c(lon[site_index], lat[site_index], elev[site_index])

        # Extract required info from the netcdf files
        fid = nc_open(all_source_tsunami[i], readunlim=FALSE)
        event_index_string = ncvar_get(fid, 'event_index_string')
        event_rates = ncvar_get(fid, 'event_rate_annual')
        event_stage = ncvar_get(fid, 'max_stage', start=c(1,site_index), 
            count=c(-1, 1))
        if(slip_type == 'stochastic'){
            event_slip = ncvar_get(fid, 'event_slip_string')
        }
        nc_close(fid)

        # Find events with stage > value
        events_exceeding = which(event_stage >= stage_exceed)

        # Quick exit
        if(length(events_exceeding) == 0) next

        for(j in events_exceeding){
            # Find unit-sources in event
            us_id = get_unit_source_indices_in_event(data.frame(
                event_index_string = event_index_string[j]))
   
            # Get slip on unit source [for weighting] 
            if(slip_type == 'stochastic'){
                unit_source_weights = as.numeric(strsplit(event_slip[j], '_')[[1]])
            }else{
                unit_source_weights = rep(1, length(us_id))
            }

            # Spread the rates over the active unit sources
            local_rates = event_rates[j] * unit_source_weights/sum(unit_source_weights)

            # Add to the store for the unit source
            sources_list[[i]]$contribution[us_id] = sources_list[[i]]$contribution[us_id] + local_rates 
    
        }

        # Logical check
        stopifnot(all.equal(
            sum(sources_list[[i]]$contribution), 
            sum(event_rates * (event_stage >=stage_exceed) )) )
    }

    return(sources_list)

}

#
#
#
plot_station_deaggregated_hazard<-function(deaggregated_hazard, scale = 1){

    plot(c(-40, 320), c(-80, 80), col=0, asp=1, xlab='lon', ylab='lat')

    contrib_range = range(unlist(lapply(deaggregated_hazard, f<-function(x) range(x$contribution))))

    scale = scale * 1/(contrib_range[2]) # Make make bar have size = 'scale' in degrees longitude

    ncol = 200
    mycol = rev(rainbow(255)[1:ncol])
    

    for(i in 1:length(deaggregated_hazard)){

        xx = deaggregated_hazard[[i]]

        colz = floor( (xx$contribution - contrib_range[1])/(contrib_range[2] - contrib_range[1]) * ncol)
        colz = pmax(colz, 1)

        arrows(xx$unit_source_statistics$lon_c, xx$unit_source_statistics$lat_c, 
            xx$unit_source_statistics$lon_c, xx$unit_source_statistics$lat_c + xx$contribution * scale,
            col = mycol[colz], length=0)
    }
    points(xx$station[1], xx$station[2], col='red', pch=19)

}

# lon_p = 151.42
# lat_p = -34.05

plot_station_exceedance_rate_pdf<-function(lon_p, lat_p, station_name = ""){

    if(station_name == ""){
        station_name = 'unnamed_station'
    }


    pdf(paste0(station_name, '_stage_exceedance_rate.pdf'), width=12, height=10)

    site_index = plot_rates_at_a_station(lon_p, lat_p)

    for(i in 1:length(rates)){
        par(mfrow=c(2,1))
        plot_wave_heights_at_a_station(lon_p, lat_p, source_zone = names(rates)[i], 
            slip_type = 'uniform', site_index=site_index)
        gc()
        plot_wave_heights_at_a_station(lon_p, lat_p, source_zone = names(rates)[i], 
            slip_type = 'stochastic', site_index=site_index)
        gc()
    }


    for(exceedance_rate in c(1/50, 1/100, 1/500, 1/1000, 1/2500)){
        for(slip_type in c('uniform', 'stochastic')){
            # Deaggregated hazard 
            par(mfrow=c(1,1))
            site_deagg = get_station_deaggregated_hazard(lon_p, lat_p, 
                slip_type = slip_type, exceedance_rate=exceedance_rate)
            plot_station_deaggregated_hazard(site_deagg, scale=8)
            points(lon, lat, col=rgb(0.5, 0.5, 0.5, alpha=0.5), pch='.')
            points(lon_p, lat_p, col='red', pch=19, cex=2)
            title(paste0('Deaggregated hazard: ', slip_type, ' slip, AEP = ', 
                signif(exceedance_rate,4), ',\n Peak stage = ', 
                signif(site_deagg[[1]]$stage_exceed,4)))
            rm(site_deagg); gc()
        }
    }

    dev.off()
}


# For a given source-zone, plot the peak stage at every gauge associated with a
# given exceedance rate.
#
plot_source_zone_stage_vs_exceedance_rate<-function(
    source_zone, 
    desired_rates=c(1/1000),
    scale = 5,
    max_stage_color = 10, 
    greens_law_to_1m = FALSE,
    shapefile_output_dir=NULL,
    shapefile_rate_names = NULL,
    ...){

    # Find the site in rates
    site = match(source_zone, names(rates))
    if(length(site) != 1){
        stop(paste0('Could not uniquely match ', source_zone, 
            ' among all source-zones'))
    }

    site_stages_uniform = matrix(NA, nrow=length(desired_rates), ncol=length(lon))
    site_stages_stochastic = matrix(NA, nrow=length(desired_rates), ncol=length(lon))

    # Stage sequence that rates have been computed for
    stages = ncvar_get(rates[[site]], 'stage')

    # Get rates at every gauge, for every value in the stage sequence
    site_rates_uniform = ncvar_get(rates[[site]], 'uniform_slip_rate')
    site_rates_stochastic = ncvar_get(rates[[site]], 'stochastic_slip_rate')


    # Compute stages for every desired rate, at every gauge
    # Uniform slip
    for(i in 1:length(lon)){

        sr = site_rates_uniform[,i]
        if(sr[1] > 0 & all(!is.na(sr))){
            # Check sorting
            dsr = diff(sr)
            if(any(dsr) < 0) stop('Unsorted rates')
            # Only use points with rate>0 to interpolate -- since rate=0 goes
            # outside our model range -- and further, we might have multiple
            # points with rate=0, which will confuse the approx function (needs
            # the rates to be monotonic). We also append a 'stage=0' value to the start,
            # and ensure it has a rate slightly above the peak rate, to ensure we have >= 2
            # interpolation points in the approx function
            kk = which(sr > 0)
            site_stages_uniform[,i] = approx(c(sr[1]+0.001, sr[kk]), c(0,stages[kk]), 
                xout=desired_rates, rule=2)$y
        }
    }
    # Stochastic slip
    for(i in 1:length(lon)){

        sr = site_rates_stochastic[,i]
        if(sr[1] > 0 & all(!is.na(sr))){
            # Check sorting
            dsr = diff(sr)
            if(any(dsr) < 0) stop('Unsorted rates')

            # Only use points with rate>0 to interpolate -- since rate=0 goes
            # outside our model range -- and further, we might have multiple
            # points with rate=0, which will confuse the approx function (needs
            # the rates to be monotonic). We also append a 'stage=0' value to the start,
            # and ensure it has a rate slightly above the peak rate, to ensure we have >= 2
            # interpolation points in the approx function

            kk = which(sr > 0)
            site_stages_stochastic[,i] = approx(c(sr[1]+0.001,sr[kk]), c(0,stages[kk]), 
                xout=desired_rates, rule=2)$y
        }
    }

    # Write out as shapefiles
    if(!is.null(shapefile_output_dir)){

        dir.create(shapefile_output_dir, recursive=TRUE, showWarnings=FALSE)

        # Make uniform slip shapefile, with rates and elevation
        s1 = SpatialPoints(cbind(lon, lat), proj4string=CRS('+init=epsg:4326'))
        local_data_frame = as.data.frame(t(site_stages_uniform))
        if(!is.null(shapefile_rate_names)){
            names(local_data_frame) = shapefile_rate_names
        }
        local_data_frame = cbind(local_data_frame, data.frame('elev'=as.numeric(elev)))
        s2 = SpatialPointsDataFrame(s1, data=local_data_frame)
        writeOGR(s2, dsn=paste0(shapefile_output_dir, '/uniform_', source_zone),
            layer=source_zone, driver='ESRI Shapefile', overwrite=TRUE)

        # Make stochastic slip shapefile, with rates and elevation
        s1 = SpatialPoints(cbind(lon, lat), proj4string=CRS('+init=epsg:4326'))
        local_data_frame = as.data.frame(t(site_stages_stochastic))
        if(!is.null(shapefile_rate_names)){
            names(local_data_frame) = shapefile_rate_names
        }
        local_data_frame = cbind(local_data_frame, data.frame('elev'=as.numeric(elev)))
        s2 = SpatialPointsDataFrame(s1, data=local_data_frame)
        writeOGR(s2, dsn=paste0(shapefile_output_dir, '/stochastic_', source_zone),
            layer=source_zone, driver='ESRI Shapefile', overwrite=TRUE)

        # Make shapefile with 'ratio of uniform to stochastic'
        # uniform_stage/(stochastic_stage+eps)
        s1 = SpatialPoints(cbind(lon, lat), proj4string=CRS('+init=epsg:4326'))
        eps = 1.0e-06
        local_data_frame = as.data.frame(t(site_stages_uniform)/(t(site_stages_stochastic)+eps))
        if(!is.null(shapefile_rate_names)){
            names(local_data_frame) = shapefile_rate_names
        }
        local_data_frame = cbind(local_data_frame, data.frame('elev'=as.numeric(elev)))
        s2 = SpatialPointsDataFrame(s1, data=local_data_frame)
        writeOGR(s2, dsn=paste0(shapefile_output_dir, '/uniform_on_stochastic_', source_zone),
            layer=source_zone, driver='ESRI Shapefile', overwrite=TRUE)

    }

    # Rescale based on greens law, if desired
    if(greens_law_to_1m){
        rescale = pmax(0, -elev)**0.25
        greens_title = ", Green's law rescaled to 1m depth"
    }else{
        rescale = lon*0 + 1
        greens_title = ''
    }

    # Make a plot
    nc = 200
    colscheme = rev(rainbow(255)[1:nc])

    # Function to convert stage values to colors
    col_value<-function(x){
        colval = sqrt(pmin(x, max_stage_color)/max_stage_color) * nc
        return(colscheme[pmin(pmax(1, round(colval)), nc)])
    }
    # Legend for wave heights, inside Australia!
    arrow_legend<-function(){
        lon_c = seq(130, 140, by=1)
        lat_c = rep(-28, length=length(lon_c))
        arrow_ht = seq(0.001, max_stage_color, len=length(lon_c))
        arrows(lon_c, lat_c, lon_c, lat_c+arrow_ht*scale, col=col_value(arrow_ht),
            length=0)
        #rect(129, -28 - 1, 141, -28+11)
        text(129+5, -28+7, paste0('Wave height \n of 0-', max_stage_color, 'm'))
    }

        
    # Make the plot, for all desired_rates
    for(i in 1:length(desired_rates)){
        par(mfrow=c(2,1))
        # Uniform
        plot(lon, lat, pch='.', col='grey', asp=1, 
            main=paste0(source_zone, ': Uniform slip, rate = ', 
                signif(desired_rates[i]), greens_title),
            ...)
        # arrow height
        ht = site_stages_uniform[i,]*rescale
        # Order arrows so long ones are plotted last -- makes it visually clearer
        o1 = order(ht) 
        arrows(lon[o1], lat[o1], lon[o1], lat[o1] + ht[o1]*scale, 
            col=col_value(ht[o1]), length=0)
        grid()
        arrow_legend()

        # Stochastic
        plot(lon, lat, pch='.', col='grey', asp=1, 
            main=paste0(source_zone, ': Stochastic slip, rate = ', 
                signif(desired_rates[i]), greens_title),
            ...)
        # arrow height
        ht = site_stages_stochastic[i,]*rescale 
        # Order arrows so long ones are plotted last, for visual clarity
        o1 = order(ht) 
        arrows(lon[o1], lat[o1], lon[o1], lat[o1] + ht[o1]*scale,
            col=col_value(ht[o1]), length=0)
        grid()
        arrow_legend()
    }

    return(invisible())
}


#
# Plot global peak-stage vs return period, for each source-zone
#
make_global_stage_return_period_plots<-function(){
    # Make a pdf for each source-zone, containing a number of different return periods
    for(source_name in names(rates)){
        pdf(paste0(source_name, '_stage_vs_exceedance_rate.pdf'), width=20, height=15)
        plot_source_zone_stage_vs_exceedance_rate(source_name, 
            desired_rates=c(1/100, 1/500, 1/1000, 1/5000, 1/10000, 1/50000, 1/100000),
            greens_law_to_1m=FALSE, 
            shapefile_output_dir=paste0('peakstage_shapefiles/', source_name),
            shapefile_rate_names = c('R_100', 'R_500', 'R_1000', 'R_5000', 'R_10000', 
                'R_50000', 'R_100000'),
            xlim=c(40,320), ylim=c(-60,60), scale=1)
        dev.off()
    }
}

#
# Make plots of stage-vs-exceedance rate, and Mw-vs-stage for each source, at a
# bunch of sites
#
make_standard_site_exceedance_rate_plots<-function(){

    sites = list(
        'Cairns_100m' = c(146.34, -16.74),
        'Bundaberg_100m' = c(152.833, -24.137), 
        'Brisbane_100m' = c(153.72, -27.56),
        'Sydney_100m' = c(151.42, -34.05), 
        'Eden_100m' = c(150.188, -37.141),
        'Hobart_100m' = c(148.207, -42.897),
        'Adelaide_100m' = c(135.614, -35.328),
        'Perth_100m' = c(115.250, -31.832),
        'SteepPoint_100m' = c(112.958, -26.093),
        'Broome_100m' = c(120.74, -18.25),
        'Darwin_100m' = c(129.03, -11.88) 
        )
    
    
    for(i in 1:length(sites)){
        plot_station_exceedance_rate_pdf(sites[[i]][1], sites[[i]][2], names(sites)[i])
    } 
 
}

#
# A useful plot for examining if we have 'enough' stochastic slip events
#
make_stage_exceedance_rate_convergence_plot<-function(){

    sites = list(
        'Cairns_100m' = c(146.34, -16.74),
        'Bundaberg_100m' = c(152.833, -24.137), 
        'Brisbane_100m' = c(153.72, -27.56),
        'Sydney_100m' = c(151.42, -34.05), 
        'Eden_100m' = c(150.188, -37.141),
        'Hobart_100m' = c(148.207, -42.897),
        'Adelaide_100m' = c(135.614, -35.328),
        'Perth_100m' = c(115.250, -31.832),
        'SteepPoint_100m' = c(112.958, -26.093),
        'Broome_100m' = c(120.74, -18.25),
        'Darwin_100m' = c(129.03, -11.88) 
        )

    pdf('stage_convergence_plot.pdf', width=8, height=6)
    for(i in 1:length(sites)){
        for(source_name in names(rates)){
            plot_wave_heights_at_a_station(
                sites[[i]][1], sites[[i]][2], 
                source_name, slip_type='stochastic', split_into_subsets=2)
            title(paste0(names(sites)[i], ' : ', source_name))
            grid()
        }
    }
    dev.off()
}

