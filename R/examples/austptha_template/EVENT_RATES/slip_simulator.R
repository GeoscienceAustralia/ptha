#
# For seismic moment conservation, it is clear that we cannot assign equal
# weight to all uniform slip events ( or use the modified approach like in
# Davies et al. which adjusts for tectonic slip) (2017) since it is quite prone
# to edge effects, and pushes the moment release into the centre of
# source-zones. 
#
# Here, we simulate 'greedy' earthquake behaviour, confirming that a complex
# pattern of earthquake event rates results. The idea is:
#   - From the fitted event rate curve, we simulate timings of the next 'nevents=1e+06' events, and
#     assign them random Mw values.
#   - We use the Bird/Jono convergence rates for spatially variable convergence
#   - We only treat uniform slip events with a fixed number of unit-sources (based on Mw)
#   - Given the earthquake magnitude, we choose the scenario which occurs in the place with maximum 'slip
#     deficit', i.e. difference between the long-term-tectonic-displacement
#     (slip_rate * time_of_event) and the historic earthquake slip (which tracks the
#     cumulative effect of the events, and is initially zero).
#
# Under these assumptions, then with a fixed slip rate, we see a 'spike' in the
# rate of events on the edge. This makes sense, since there have to be 'higher-rates'
# of the edge event, to keep the edge cells converging at the same rate as other points.
#

config_env = new.env()                                                          
source('config.R', local=config_env)      

#
# Get slip rates
#
source('make_spatially_variable_source_zone_convergence_rates.R')               
bird2003_env = event_conditional_probability_bird2003_factory(                  
    return_environment=TRUE)                      

#
# Simple uniform slip simulator
#
# Try this to deal with 'edge effects' in spatially variable subduction zone rates
#
# @param sourcename name of source-zone, e.g. 'sunda' or 'sunda_java'
# @param nevents Simulate this many earthquakes above MW_MIN, from which
#        individual event rates can be derived
# @param mu Shear modulus, used for computing tectonic moment (and thus coupling)
# @param fix_tectonic_slip If not NULL, then use this as the long-term slip rate (m/year) everywhere
# @return the function environment
#
simulate_event_rates<-function(sourcename, nevents=1e+06, mu=3e+010, fix_tectonic_slip=NULL){

    # Read the unit source summary statistics                                   
    basename_uss_files = basename(config_env$unit_source_statistics_netcdf_files)
    which_uss_file = which(basename_uss_files ==                                
        paste0('unit_source_statistics_', sourcename, '.nc'))                   
                                                                                
    if(length(which_uss_file) != 1){                                            
        stop(paste0('Could not find a unique unit_source_statistics nc file matching ',
            sourcename))                                                        
    }                    
    tsunami_events_file = config_env$all_source_uniform_slip_tsunami[which_uss_file]
    uss_file = config_env$unit_source_statistics_netcdf_files[which_uss_file]   
    uss = read_table_from_netcdf(uss_file)

    # Get the relevant data for events, including:
    # - The unit-sources they involve
    # - Their slip
    # - The 'overall' rate of events in each magnitude bin
    #
    fid = nc_open(tsunami_events_file, readunlim=FALSE)                         
    # Only use event_rates 'integrated over all events with the same Mw', since
    # the whole point of this routine is to distribute event rates in a new way
    event_rates = ncvar_get(fid, 'event_rate_annual')
    event_Mw = ncvar_get(fid, 'event_Mw')
    event_index_string = ncvar_get(fid, 'event_index_string')
    # Indices of unit-sources in event
    eis = sapply(event_index_string,
        f<-function(x) as.numeric(strsplit(x, split="-")[[1]]), simplify=FALSE) 
    # Slip on those indices. Set it up like stochastic slip, even though slip
    # is constant in the uniform case treated here
    slip = ncvar_get(fid, 'event_slip')
    # Broadcast the slip to a vector for each event, so the format is
    # identical as for the stochastic slip case
    ess = eis                                                               
    for(i in 1:length(eis)){                                                
        ess[[i]] = eis[[i]]*0 + slip[i]                                     
    }                                                                       

    # rate in m/year -- always positive
    if(is.null(fix_tectonic_slip)){
        source_slip_rate = abs(bird2003_env$unit_source_tables[[sourcename]]$bird_vel_div*1/1000)
    }else{
        source_slip_rate = rep(fix_tectonic_slip, length(uss[,1]))
    }
    # Smoothing -- probably change rupture info
    


    historical_slip = source_slip_rate*0 # we will accumulate this

    # Times and magnitudes for next 'nevents' events
    nevents = nevents
    mu = mu
    sim_rate = sum(event_rates) # Rate of events above min(event_Mw)
    # event times
    sim_event_t = cumsum(rexp(nevents, rate=sim_rate))
    # event mw -- sample with weights proportional to rates
    sim_rate_mw = aggregate(event_rates, by=list(event_Mw), f<-function(x) sum(x))
    sim_event_t_mw = sample(sim_rate_mw[,1], size=length(sim_event_t), 
        replace=TRUE, prob=sim_rate_mw[,2]/sum(sim_rate_mw[,2]))
    
    # Back-compute source-zone coupling coefficient
    coupling_coef = sum(M0_2_Mw(sim_event_t_mw, inverse=TRUE))/
        sum(max(sim_event_t)*mu*uss$length*uss$width*1e+06*source_slip_rate)

    # Store when an event occurred
    event_occurred = event_Mw * 0

    # Fast lookup of events with Mw = particular value
    unique_mw = unique(event_Mw)
    unique_mw_start = sapply(unique_mw, f<-function(x) min(which(event_Mw == x)))
    unique_mw_end = sapply(unique_mw, f<-function(x) max(which(event_Mw == x)))
    eis_matrix_list = list()
    for(i in 1:length(unique_mw)){
        n1 = unique_mw_start[i]
        n2 = unique_mw_end[i]
        eis_matrix_list[[i]] = matrix(unlist(eis[n1:n2]), ncol=(n2-n1+1))
    }

    sim_event_t_mw_ind = match(sim_event_t_mw, unique_mw)
    sim_event_t_event_ind = rep(NA, nevents)
    

    for(i in 1:length(sim_event_t)){
        # Find events with mw = sim_event_t_mw[i]        
        mw_ind = sim_event_t_mw_ind[i]
        n1 = unique_mw_start[mw_ind]
        n2 = unique_mw_end[mw_ind]

        t = sim_event_t[i]

        slip_deficit = historical_slip - t*source_slip_rate*coupling_coef

        # For each event, find the difference between the tectonic motion
        # [source_slip_rate*sim_event_t_mw[i]] and the historical motion on unit-sources that
        # are involved in the event. Choose the one where the difference is as
        # small (i.e. large negative) as possible.
        ## event_ind = which.min( (lapply(eis[n1:n2], f<-function(x) sum(slip_deficit[x]))) ) + (n1-1)
        # This approach is faster
        slip_deficit_mat = slip_deficit[eis_matrix_list[[mw_ind]]]
        dim(slip_deficit_mat) = dim(eis_matrix_list[[mw_ind]])
        # Break ties in ranks at random [probably doesn't matter]
        #event_ind = which.min(rank(colSums(slip_deficit_mat), ties.method='random')) + (n1-1)
        event_ind = which.min(colSums(slip_deficit_mat)) + (n1-1)

        sim_event_t_event_ind[i] = event_ind

        # Record that the event occurred
        event_occurred[event_ind] = event_occurred[event_ind] + 1
        active_sources = eis[[event_ind]]
        slip_sources = ess[[event_ind]]
        # Update the motion on the unit source
        historical_slip[active_sources] = historical_slip[active_sources] + slip_sources
        if(i%%1e+04 == 0){
            print(c(i, range(historical_slip/t - source_slip_rate*coupling_coef)))
            gc()
        }
    }

    return(environment())
}
