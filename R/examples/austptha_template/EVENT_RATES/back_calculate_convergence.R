library(rptha)
if(!exists('config')){
    config = new.env()
    source('config.R', local=config)
}
if(class(config) != 'environment'){
    stop('Apparent name clash with "config" variable')
}

#' Routine to back-calculate the convergence rate on each unit source, based on
#' the earthquake events and their rates
#' 
#' It can also be used to examine the effect of inflating the rate of edge
#' events (i.e. earthquakes with an edge at one of the source-zone edges) FOR
#' UNIFORM SLIP ONLY
#' 
#' Use of an appropriate 'edge-multiplier' can help to better satisfy seismic moment
#' conservation. This is because with (edge_multiplier = 0), seismic slip tends to
#' be concentrated towards the centre of the fault, even if the desired tectonic slip
#' is uniform, because unit-sources near the edge are covered by less events. Inflating
#' the rate at the edge can greatly improve this.
#'
#' @param sourcename name of source, e.g. 'sunda'
#' @param slip_type either 'uniform' or 'variable_uniform' or 'stochastic'. 
#' @param edge_multiplier If slip_type == 'uniform', inflate the conditional
#'   probabilities of edge events (the ones touching extreme points of the unit
#'   source geometry) by the factor edge_multiplier. This can help reduce edge
#'   effects in the spatial patterns of moment conservation
#' @param uss If (slip_type == 'uniform'), you may provide the
#'   unit_source_summary_statistics as a data.frame. Otherwise they are read from a
#'   file. If provided, then all subsequent optional arguments must be provided too.
#' @param event_rates Rate for each event (see help for uss, requires provision of ALL event data)
#' @param event_index_string String containing unit source indices for each event (see help for uss,
#'   requires provision of ALL event data)
#' @param slip vector with slip values for each event (see help for uss,
#'   requires provision of ALL event data)
#' @return the function environment
back_calculate_convergence<-function(sourcename, slip_type='uniform', edge_multiplier=0,
    uss = NULL, event_rates=NULL, event_Mw = NULL, event_index_string = NULL, 
    slip = NULL){

    # This routine can be used to experiment with inflated rates of edge events
    # (to enhance seismic moment conservation), but it only makes sense for
    # uniform slip
    if(edge_multiplier != 0){
        stopifnot(slip_type == 'uniform')
    }

    # Check input arguments
    k = is.null(uss) + is.null(event_rates) + is.null(event_Mw) + is.null(event_index_string) + is.null(slip)
    if(! (k %in% c(0, 5))){
        msg = 'You must either provide NO event properties, or ALL of them'
        stop(msg)
    }

    if(k == 0 & (slip_type != 'uniform')){
        stop('Function only works with non-null event properties as arguments if slip_type=="uniform"')
    }

    # Useful flag, to know if we get values from a file or not
    read_data_from_file = (k==5)

    # Get the unit source summary statistics
    if(read_data_from_file){

        # Find the unit source summary statistics file
        basename_uss_files = basename(config$unit_source_statistics_netcdf_files)
        which_uss_file = which(basename_uss_files ==
            paste0('unit_source_statistics_', sourcename, '.nc'))

        if(length(which_uss_file) != 1){
            stop(paste0('Could not find a unique unit_source_statistics nc file matching ',
                sourcename))
        }

        uss_file = config$unit_source_statistics_netcdf_files[which_uss_file]
        uss = read_table_from_netcdf(uss_file)

    }else{
        # Bring function argument into local environment
        uss = uss
    }

    # Get the event data
    if(read_data_from_file){

        # Find the relevant netcdf file
        if(slip_type == 'uniform'){
            tsunami_events_file = config$all_source_uniform_slip_tsunami[which_uss_file]
        }else if(slip_type == 'stochastic'){
            tsunami_events_file = config$all_source_stochastic_slip_tsunami[which_uss_file]
        }else if(slip_type == 'variable_uniform'){
            tsunami_events_file = config$all_source_variable_uniform_slip_tsunami[which_uss_file]
        }
        # For speed of reading, use the file that only contains earthquake events
        tsunami_events_file = gsub('_tsunami', '', tsunami_events_file)


        # Get the relevant data
        fid = nc_open(tsunami_events_file, readunlim=FALSE)

        # Event rates -- may be modified if edge_multiplier != 0
        event_rates = ncvar_get(fid, 'rate_annual')

        # Event magnitudes
        event_Mw = ncvar_get(fid, 'Mw')
        unique_Mw = sort(unique(event_Mw))

        # Unit source indices in each event
        event_index_string = ncvar_get(fid, 'event_index_string')

    }else{

        # Bring function arguments into this environment
        event_rates = event_rates
        event_Mw = event_Mw
        unique_Mw = sort(unique(event_Mw))
        event_index_string = event_index_string

    }

    # Convert event_index_string to an easier-to-use form
    eis = sapply(event_index_string,
        f<-function(x) as.numeric(strsplit(x, split="-")[[1]]), simplify=FALSE)
    rm(event_index_string); gc()

    # Flag events which have a unit-source at the edge of the rupture
    uss_alongstrike_range = range(uss$alongstrike_number)
    is_on_edge = unlist(lapply(eis, 
        f<-function(x) any(uss$alongstrike_number[x] %in% uss_alongstrike_range)))

    if(slip_type == 'uniform'){

        if(read_data_from_file){
            slip = ncvar_get(fid, 'slip')
        }else{
            # Bring argument into function environment
            slip = slip
        }
        # Broadcast the slip to a vector for each event, so the format is
        # identical as for the stochastic slip case 
        ess = eis
        for(i in 1:length(eis)){
            ess[[i]] = eis[[i]]*0 + slip[i]
        }

    }else{
        # Non-uniform slip -- the data MUST be read from a file, as enforced above
        event_slip_string = ncvar_get(fid, 'event_slip_string')
        # List containing a vector of slip values for each event
        ess = sapply(event_slip_string,
            f<-function(x) as.numeric(strsplit(x, split="_")[[1]]),
            simplify=FALSE)
        rm(event_slip_string); gc()
    }
    if(read_data_from_file) nc_close(fid)


    if(edge_multiplier != 0){

        stopifnot(edge_multiplier > 0)
        new_event_rates = event_rates*0
        new_conditional_probability = event_rates*0

        for(i in 1:length(unique_Mw)){
            k = which(event_Mw == unique_Mw[i])
            if(all(event_rates[k] == 0)){
                # Dummy values for impossible event
                new_conditional_probability[k] = 1/length(k)
                next
            }
            # Compute new rate, with appropriate normalisation so that
            #   sum(event_rates[k]) = sum(new_event_rates[k])
            new_event_rates[k] = (event_rates[k]*(1 + edge_multiplier*(is_on_edge[k])))
            new_event_rates[k] = new_event_rates[k] * sum(event_rates[k])/sum(new_event_rates[k])
            new_conditional_probability[k] = new_event_rates[k]/sum(new_event_rates[k])
        }

    }else{
        # Do not change anything, but make 'new_' variables containing the OLD
        # data -- so that the code below works
        new_event_rates = event_rates
        new_conditional_probability = event_rates*0
        for(i in 1:length(unique_Mw)){
            k = which(event_Mw == unique_Mw[i])
            if(all(event_rates[k] == 0)){
                new_conditional_probability[k] = 1/length(k)
            }else{
                new_conditional_probability[k] = new_event_rates[k]/sum(new_event_rates[k])
            }
        }
    }
    
   
    # Integrate slip rate. We would like this to look similar to the desired
    # convergence rate 
    us_slip_rate = rep( 0, length(uss[,1]) )
    for(i in 1:length(ess)){
        ind = eis[[i]]
        slp = ess[[i]]
        us_slip_rate[ind] = us_slip_rate[ind] + slp*new_event_rates[i]
    }

    output = cbind(uss, data.frame(integrated_slip = us_slip_rate))

    #
    # Now get the integrated slip for each unique mw value
    # 

    for(j in 1:length(unique_Mw)){

        mwj = unique_Mw[j]
        us_slip_rate_mwj = us_slip_rate*0
        #us_slip_rate_mwj_edge_only = us_slip_rate*0
        events_with_Mw = which(event_Mw == mwj)


        # Integrate slip over all events
        for(i in events_with_Mw){

            ind = eis[[i]]
            slp = ess[[i]]

            # Integrate slip for all events
            us_slip_rate_mwj[ind] = us_slip_rate_mwj[ind] +
                slp * new_event_rates[i]

        }
        output = cbind(output, data.frame(us_slip_rate_mwj))
        names(output)[ncol(output)] = paste0('int_slip_', round(mwj*10))

    }

    return(environment())
}

##
# yy2 = back_calculate_convergence('sunda', slip_type='stochastic')
# plot(yy2$lon_c, yy2$integrated_slip, col=yy2$downdip_number, pch=19)
# points(yy$lon_c, yy$integrated_slip, col=yy$downdip_number)

batch_plot_convergence<-function(){
    
    for(i in 1:length(config$source_names_1)){

        sname = config$source_names_1[i]
        print(sname)

        # Try to reduce memory
        source_uniform = NULL
        source_stochastic = NULL
        source_variable_uniform = NULL
        gc()

        # Calculate convergence
        source_uniform = back_calculate_convergence(sname, 'uniform')
        source_stochastic = back_calculate_convergence(sname, 'stochastic')
        source_variable_uniform = back_calculate_convergence(sname, 'variable_uniform')

        normalizer = max(source_uniform$output$integrated_slip)

        panel_plot<-function(source_X, normalizer, type, size_scale=NULL){

            if(is.null(size_scale)){
                size_scale = source_X$output$integrated_slip
            }
            plot(source_X$output$lon_c, source_X$output$lat_c, 
                cex=size_scale/normalizer, asp=1, 
                main=paste0(sname, ' ', type)); grid()
            v = cbind(-cos(source_X$uss$strike/180*pi)*size_scale, 
                      sin(source_X$uss$strike/180*pi)*size_scale)*100 + 1.0e-06
            arrows(source_X$output$lon_c, source_X$output$lat_c,
                source_X$output$lon_c + v[,1], source_X$output$lat_c + v[,2], 
                length=0, col='red')
        }

        par(mfrow=c(1,3))
        panel_plot(source_uniform, normalizer, 'Uniform')
        panel_plot(source_stochastic, normalizer, 'Stochastic')
        panel_plot(source_variable_uniform, normalizer, 'Variable-uniform')

        # Look at rate of > 9.0
        u1 = rowSums(as.matrix(source_uniform$output[,35:43]))
        s1 = rowSums(as.matrix(source_stochastic$output[,35:43]))
        vu1 = rowSums(as.matrix(source_variable_uniform$output[,35:43]))
        normalizer = max(u1)

        par(mfrow=c(1,3))
        panel_plot(source_uniform, normalizer, 'Uniform >= 9', u1)
        panel_plot(source_stochastic, normalizer, 'Stochastic >= 9', s1)
        panel_plot(source_variable_uniform, normalizer, 'Variable-uniform >= 9', vu1)

        # Look at rate of < 8.5
        u1 = rowSums(as.matrix(source_uniform$output[,17:30]))
        s1 = rowSums(as.matrix(source_stochastic$output[,17:30]))
        vu1 = rowSums(as.matrix(source_variable_uniform$output[,17:30]))
        normalizer = max(u1)

        par(mfrow=c(1,3))
        panel_plot(source_uniform, normalizer, 'Uniform <= 8.5', u1)
        panel_plot(source_stochastic, normalizer, 'Stochastic <= 8.5', s1)
        panel_plot(source_variable_uniform, normalizer, 'Variable-uniform <= 8.5', vu1)
        #plot(c(0,1), col='white')
    }

}

make_plot<-function(){
    pdf('convergence_integrated.pdf', width=18, height=6)
    batch_plot_convergence()
    dev.off()
}


