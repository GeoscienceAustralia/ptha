#
# Extract statistics from the models-vs-observations.
#
# Note:
#     - For some gauges, the observation is 'just noise' -- for instance, the Cocos Island gauge
#       won't well detect small tsunamis in the Pacific.
#     - The observations often include long-times prior to the tsunami. Also, they retain the Rayleigh-wave
#       perturbations caused by the earthquake. These can be removed with various strategies, for instance,
#       only consider the data after the MODEL arrival (using some estimator thereof).
#     - These issues are NOT accounted for here -- later analyses need to apply their own -hoc filtering
#

# Work with the RDS files produced with the ../swals/plots/plot_all.R script
target_RDS = c(Sys.glob('../gauge_RDS_files/OUTPUTS/*-full-*-highres_NSW/*/gauge*.RDS'),
               Sys.glob('../gauge_RDS_files/OUTPUTS/*-full-*-highres_australia/*/gauge*.RDS'))
# Avoid ptha scenario runs that include some other events -- only get the inversions
to_keep = unique(unlist(lapply(c('Chile1960', 'Sumatra2004', 'Tohoku2011', 'Chile2010', 'Chile2015'), 
                        f<-function(x) grep(x, target_RDS))))
target_RDS = target_RDS[to_keep]
# Here we only want to work with 'standard' models for each inversion.
# Do not want to work with the highres model, or the leapfrog-nonlinear model -- their purpose
# is to test the other models, but including them in this analysis would be double-counting / inconsistent.
# Also don't want to work with ptha scenarios for now.
# Also skip the australia-wide Chile simulations which double-up in a few cases -- the ones we need are 
# 'NSW' Chile 1960 simulations.
to_remove = unique( c(grep('HIGH', target_RDS), grep('leapfrog_nonlinear', target_RDS), 
                      which(grepl("Chile1960", target_RDS) & grepl('australia', target_RDS)) ) )
target_RDS = target_RDS[-to_remove]

# Check I have 12 events for nkinds_of_model
nkinds_of_model = 5
if(length(target_RDS) != nkinds_of_model*12) stop('Did not get the right number of RDS files, beware false matches')

ONE_HOUR=3600.0

WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD = 5e-04
    # When 
    #      abs(modelled-stage) > ( max(abs(modelled-stage))*WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD )
    # we say the tsunami has 'arrived'. This begins the time-period over which
    # we might compare model and data (although that is also affected by other
    # model parameters below).

DATA_SIZE_THRESHOLD = 1e-02
    # Skip gauges where the data max is below the threshold, OR, the "negative data min" is below the threshold
MODEL_SIZE_THRESHOLD = 5e-03
    # Skip gauges where the model max is below the threshold, OR, the "negative model min" is below the threshold
    # The model is not contaminated by noise


# Interrogate the target_RDS for various stats
get_statistics_at_site<-function(site_data){

    # The observed gauge time-series
    obs = site_data$event_data$obs

    # The model time-series at this gauge
    model = data.frame(
        juliant = as.difftime(site_data$model_time/(ONE_HOUR*24), units='days') + site_data$model_start,
        resid = site_data$model_stage)

    n = length(model$juliant)
    model_range = diff(range(model$resid))

    if(n == 0 | model_range == 0){
        # Fail gracefully in the case where we didn't extract a model time-series.
        # This corresponds to one DART in the middle of the Pacific (fix this
        # by adding a hazard point nearby).
        model_range = c(NA, NA) 
        data_range = c(NA, NA) 
        model_range8hrs = c(NA, NA) 
        data_range8hrs = c(NA, NA) 
        model_range12hrs = c(NA, NA) 
        data_range12hrs = c(NA, NA) 
        model_range24 = c(NA, NA) 
        data_range24 = c(NA, NA) 
        model_range36 = c(NA, NA) 
        data_range36 = c(NA, NA) 
        model_range_lastday = c(NA, NA) 
        data_range_lastday = c(NA, NA) 
        model_range_24to36 = c(NA, NA) 
        data_range_24to36 = c(NA, NA) 

        model_arrival_time = NA
        model_time_to_arrive = NA
        obs_max_time_to_arrive = NA

        model_time_of_max = NA
    }else{
        
        wave_threshold = WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD *max(abs(model$resid))
        k = min(which(abs(model$resid) > wave_threshold))
        model_arrival_time = model$juliant[k]
        model_time_to_arrive = site_data$model_time[k]

        #
        # Look at observations within the model time range. This can be problematic because of
        # Rayleigh waves in the DART records -- these need to be removed somehow, e.g. based on
        # travel-time approaches.
        #
        obs_keep = which( (obs$juliant >= model_arrival_time) & (obs$juliant <= model$juliant[n]))
        data_range = range(obs$resid[obs_keep], na.rm=TRUE)
        model_range = range(model$resid)
      
        # Find the time of the observed maxima 
        tmp = which( obs$resid[obs_keep] == max(data_range) )[1]
        obs_max_time_to_arrive = (obs$juliant[obs_keep[tmp]] - model$juliant[1])*ONE_HOUR*24

        #
        # Now look at results when we cut off the final 1 day (i.e. run 36 hours instead of 60)
        #
        obs_keep   = which( (obs$juliant   >= model_arrival_time) & (obs$juliant   <= (model$juliant[n]-1)) )
        model_keep = which( (model$juliant >= model_arrival_time) & (model$juliant <= (model$juliant[n]-1)) )
        model_range36 = range(model$resid[model_keep])
        data_range36 = range(obs$resid[obs_keep], na.rm=TRUE)

        #
        # Now look at results when we cut off the final 1.5 day (i.e. 24 hours instead of 60)
        #
        obs_keep   = which( (obs$juliant   >= model_arrival_time) & (obs$juliant   <= (model$juliant[n]-1.5)) )
        model_keep = which( (model$juliant >= model_arrival_time) & (model$juliant <= (model$juliant[n]-1.5)) )
        model_range24 = range(model$resid[model_keep])
        data_range24 = range(obs$resid[obs_keep], na.rm=TRUE)

        #
        # Now look at results "within 8 hours of tsunami arrival"
        #
        obs_keep   = which( (obs$juliant   >= model_arrival_time) & (obs$juliant   <= (model_arrival_time+8/24)) )
        model_keep = which( (model$juliant >= model_arrival_time) & (model$juliant <= (model_arrival_time+8/24)) )
        model_range8hrs = range(model$resid[model_keep])
        data_range8hrs = range(obs$resid[obs_keep], na.rm=TRUE)


        #
        # Now look at results "within 12 hours of tsunami arrival"
        #
        obs_keep   = which( (obs$juliant   >= model_arrival_time) & (obs$juliant   <= (model_arrival_time+12/24)) )
        model_keep = which( (model$juliant >= model_arrival_time) & (model$juliant <= (model_arrival_time+12/24)) )
        model_range12hrs = range(model$resid[model_keep])
        data_range12hrs = range(obs$resid[obs_keep], na.rm=TRUE)


        #
        # Now look at results RESTRICTED TO the final 1 day (i.e. from 36-60 hours)
        # This can help to illuminate late-time biases
        #
        obs_keep   = which( (obs$juliant   >= model$juliant[n]-1) & (obs$juliant   <= (model$juliant[n])) )
        model_keep = which( (model$juliant >= model$juliant[n]-1) & (model$juliant <= (model$juliant[n])) )
        model_range_lastday = range(model$resid[model_keep])
        data_range_lastday = range(obs$resid[obs_keep], na.rm=TRUE)

        #
        # Now look at results RESTRICTED TO the 24-36 hours, i.e. 1-1.5 days.
        # This may help illuminate the transition from 'few differences' to 'some differences'
        #
        obs_keep   = which( (obs$juliant   >= model$juliant[n]-1.5) & (obs$juliant   <= (model$juliant[n]-1)) )
        model_keep = which( (model$juliant >= model$juliant[n]-1.5) & (model$juliant <= (model$juliant[n]-1)) )
        model_range_24to36 = range(model$resid[model_keep])
        data_range_24to36 = range(obs$resid[obs_keep], na.rm=TRUE)

        # Time of modelled maxima
        k = which.max(model$resid)
        model_time_of_max = site_data$model_time[k]
    }

    # For the output, make sensible names
    output = list(#model_range = model_range, 
                  model_min = model_range[1],
                  model_max = model_range[2],
                  model_stgrng = diff(model_range[1:2]),
                  #data_range=data_range, 
                  data_min = data_range[1],
                  data_max = data_range[2],
                  data_stgrng = diff(data_range[1:2]),
                  #model_range8hrs=model_range8hrs, 
                  model_min_8hrs = model_range8hrs[1],
                  model_max_8hrs = model_range8hrs[2],
                  model_stgrng_8hrs = diff(model_range8hrs[1:2]),
                  #data_range8hrs=data_range8hrs,
                  data_min_8hrs = data_range8hrs[1],
                  data_max_8hrs = data_range8hrs[2],
                  data_stgrng_8hrs = diff(data_range8hrs[1:2]),
                  #model_range12hrs=model_range12hrs, 
                  model_min_12hrs = model_range12hrs[1],
                  model_max_12hrs = model_range12hrs[2],
                  model_stgrng_12hrs = diff(model_range12hrs[1:2]),
                  #data_range12hrs=data_range12hrs,
                  data_min_12hrs = data_range12hrs[1],
                  data_max_12hrs = data_range12hrs[2],
                  data_stgrng_12hrs = diff(data_range12hrs[1:2]),
                  #model_range24=model_range24, 
                  model_min_24 = model_range24[1],
                  model_max_24 = model_range24[2],
                  model_stgrng_24 = diff(model_range24[1:2]),
                  #data_range24=data_range24,
                  data_min_24 = data_range24[1],
                  data_max_24 = data_range24[2],
                  data_stgrng_24 = diff(data_range24[1:2]),
                  #model_range36=model_range36, 
                  model_min_36 = model_range36[1],
                  model_max_36 = model_range36[2],
                  model_stgrng_36 = diff(model_range36[1:2]),
                  #data_range36=data_range36,
                  data_min_36 = data_range36[1],
                  data_max_36 = data_range36[2],
                  data_stgrng_36 = diff(data_range36[1:2]),
                  #model_range_lastday=model_range_lastday, 
                  model_min_lastday = model_range_lastday[1],
                  model_max_lastday = model_range_lastday[2],
                  model_stgrng_lastday = diff(model_range_lastday[1:2]),
                  #data_range_lastday=data_range_lastday
                  data_min_lastday = data_range_lastday[1],
                  data_max_lastday = data_range_lastday[2],
                  data_stgrng_lastday = diff(data_range_lastday[1:2]),
                  #
                  #model_range_lastday=model_range_lastday, 
                  model_min_24to36 = model_range_24to36[1],
                  model_max_24to36 = model_range_24to36[2],
                  model_stgrng_24to36 = diff(model_range_24to36[1:2]),
                  #data_range_lastday=data_range_lastday, 
                  data_min_24to36 = data_range_24to36[1],
                  data_max_24to36 = data_range_24to36[2],
                  data_stgrng_24to36 = diff(data_range_24to36[1:2]),
                  #
                  model_arrival_time = model_arrival_time,
                  model_time_to_arrive = model_time_to_arrive,
                  obs_max_time_to_arrive = obs_max_time_to_arrive,
                  #
                  model_time_of_max = model_time_of_max

                  )

    # We could have "Inf" or "-Inf" values above if the data didn't cover the
    # time period (e.g. Cronulla gauge for Chile 1960, which only lasts a day
    # or so in the record we have). Convert those values to NA
    for(i in 1:length(output)){
        if(any(is.infinite(output[[i]]))) output[[i]] = NA
    }

    return(output)
}

#
# Run the statistic for each event/site
#
stat_store = vector(mode='list', length=length(target_RDS))
names(stat_store) = target_RDS
target_RDS_data = vector(mode='list', length=length(target_RDS))
names(target_RDS_data) = target_RDS
for(i in 1:length(target_RDS)){
    print(i)
    x = readRDS(target_RDS[i])
    all_site_data = lapply(x, get_statistics_at_site)
    names(all_site_data) = names(x)
    stat_store[[i]] = all_site_data

    target_RDS_data[[i]] = x

}

#
# Convert to a table. We want
#   - event
#   - site
#   - model_range1
#   - model_range2
#   - data_range1
#   - data_range2
#   - inversion
#   - hydrodynamic_model_type
#

# Helper function 1
unpack_site<-function(site_data) unlist(lapply(site_data, unlist))
# Helper function 2
unpack_event<-function(event){
    all_sites_stats = lapply(event, unpack_site)
    out_df = do.call(rbind, all_sites_stats)
    out_df = cbind(data.frame(sites=names(event), stringsAsFactors=FALSE), out_df)
    rownames(out_df) = NULL
    return(out_df)
}
# Helper function 3
unpack_all_events<-function(stat_store){ 
    event_dfs = lapply(stat_store, unpack_event)
    for(i in 1:length(event_dfs)){
        n = nrow(event_dfs[[i]])
        event_dfs[[i]] = cbind(
            data.frame(run_name=rep(names(stat_store)[i], n)), 
            event_dfs[[i]])
    }
    return(do.call(rbind, event_dfs))
}

# Here is the table
event_stats = unpack_all_events(stat_store)
rownames(event_stats) = NULL


# Flag to note the INVERSION simulations.
is_inversion = rep(TRUE, length(event_stats$run_name)) #grepl('INVERSIONS', event_stats$run_name)

# Denote the model-type, and append to event_stats
model_type_from_run_name<-function(run_name){
    run_name2 = tolower(run_name)

    model_type = NA
    if(grepl('-linear_with_manning-0.035', run_name2, fixed=TRUE)){
        model_type = 'Manning0.035'
    }else if(grepl('-linear_with_no_friction-', run_name2, fixed=TRUE)){
        model_type = 'Frictionless'
    }else if(grepl('-linear_with_linear_friction-', run_name2, fixed=TRUE)){
        model_type = 'LinearFriction'
    }else if(grepl('-linear_with_reduced_linear_friction-', run_name2, fixed=TRUE)){
        model_type = 'LinearReducedFriction'
    }else if(grepl('-linear_with_delayed_linear_friction-', run_name2, fixed=TRUE)){
        model_type = 'LinearDelayedFriction'
    }
    return(model_type)
}
model_type = sapply(event_stats$run_name, model_type_from_run_name)
model_type_int = match(model_type, c('Frictionless', 'Manning0.035', 'LinearFriction', 
    'LinearReducedFriction', 'LinearDelayedFriction'))

# Denote the event name
event_name_and_source_from_run_name<-function(run_name){
    # Dates are in names, separated with either '-' or '_'. Simplify the search by
    # replacing '-' with '_'
    run_name2 = gsub('_', '-', run_name)

    event_name = NA
    source_name = NA
    #if(grepl('ptha18', run_name2)) source_name = 'ptha18'

    # Do the match for the inversions
    if(grepl('Chile1960', run_name2)){

        event_name = 'Chile-1960'
        if(grepl('FujiSatake2013', run_name2)){
            source_name = 'Fujii13'
        }else if(grepl('HoEtAl2019', run_name2)){
            source_name = 'Ho2019'
        }

    }else if(grepl('Sumatra2004', run_name2)){

        event_name = 'Sumatra-2004'
        if(grepl('PiatanesiLorito2007', run_name2)){
            source_name = 'Piatanesi08'
        }else if(grepl('LoritoEtAl2010', run_name2)){
            source_name = 'Lorito10'
        }else if(grepl('FujiSatake2007', run_name2)){
            source_name = 'Fujii06'
        }

    }else if(grepl('Chile2010', run_name2)){

        event_name = 'Chile-2010'
        if(grepl('Fuji', run_name2)){
            source_name = 'Fujii13'
        }else if(grepl('Lorito', run_name2)){
            source_name = 'Lorito11'
        }

    }else if(grepl('Tohoku2011', run_name2)){

        event_name = 'Tohoku-2011'

        if(grepl('SatakeEtAl2013', run_name2)){
            source_name = 'Satake13'
        }else if(grepl('YamakaziEtAl2018', run_name2)){
            source_name = 'Yamakazi18'
        }else if(grepl('RomanoEtAl2015', run_name2)){
            # In the raster/model naming I often call this Romano2015, but the paper is from 2014
            source_name = 'Romano14'
        }

    }else if(grepl('Chile2015', run_name2)){

        event_name = 'Chile-2015'

        if(grepl('WilliamsonEtAl2017', run_name2)){
            source_name = 'Williamson17'
        }else if(grepl('RomanoEtAl2016', run_name2)){
            source_name = 'Romano16'
        }

    }

    return(c(event_name, source_name))
}
event_name_and_source = sapply(event_stats$run_name, 
    event_name_and_source_from_run_name, USE.NAMES=FALSE)
event_name = event_name_and_source[1,]
source_name = event_name_and_source[2,]
# Flag with identical values corresponding to a single observed record
site_and_event = paste0(event_name, '_', event_stats$sites)

#
# We won't want to use all of the data. We should skip;
# - "double-ups (2 gauges @ 1 site)"
# - "clearly erronious data"
# - "very small data" (where the signal is more likely noise than tsunami)
# - data with infrequent sampling (e.g. 15 minute in NSW)
#
# For the nearshore data, use 1-minute data in the nearshore -- or 5 minute data (for WA). Except
# for Chile 1960 where we take whatever we have (digitized scans of old tidal
# gauges). We always avoid Cocos Island always because the model doesn't have
# good resolution (because I don't have good bathymetry there).
good_nearshore_data = 
    ( grepl('_1min', event_stats$sites) | grepl('_5min', event_stats$sites) | (event_name == 'Chile-1960')) & 
    (pmin(event_stats$data_max, -event_stats$data_min) > DATA_SIZE_THRESHOLD) &
    (pmin(event_stats$model_max, -event_stats$model_min) > MODEL_SIZE_THRESHOLD)
#
# For Chile 1960, we need to give special treatment to the Cronulla gauge that doesn't cover the full
# event time. From visual inspection all stats are OK except for those post 29 hours: 36hr, last-day, 24to36, and "full"
k = which(site_and_event == "Chile-1960_Cronulla_CSIRO_Fisheries_1960")
event_stats$data_max[k] = NA
event_stats$data_min[k] = NA
event_stats$data_min_36[k] = NA
event_stats$data_max_36[k] = NA
event_stats$data_min_lastday[k] = NA
event_stats$data_max_lastday[k] = NA
event_stats$data_min_24to36[k] = NA
event_stats$data_max_24to36[k] = NA

event_name_int = as.numeric(as.factor(event_name))
event_and_source_int = as.numeric(as.factor(paste0(event_name, source_name)))

# Append to Table
event_stats = cbind( 
                    data.frame(is_inversion = is_inversion, 
                               model_type   = model_type, 
                               event_name = event_name,
                               event_name_int = event_name_int,
                               source_name = source_name,
                               event_and_source_int = event_and_source_int,
                               event_name_and_source = paste0(event_name, ' ', source_name),
                               good_nearshore = good_nearshore_data,
                               site_and_event = site_and_event,
                               stringsAsFactors=FALSE),
                    event_stats)

#
# Make some plots of the 'good' timeseries for QC purposes. 
#
for(use_model in c(TRUE, FALSE)){
    #for(data_type in c('nearshore', 'dart')){
    for(data_type in 'nearshore'){

        # Make a plot of each data time-series, identifying the maxima in various time-ranges.
        if(use_model){
            pdf(paste0('QC_plot_', data_type, '_data_with_model.pdf'), width=10, height=5)
        }else{
            pdf(paste0('QC_plot_', data_type, '_data.pdf'), width=10, height=5)
        }
        for(i in 1:nrow(event_stats)){

            if(data_type == 'nearshore'){
                # Only look at 'good-nearshore' series which we use for comparison
                if(!event_stats$good_nearshore[i]) next
            }else if(data_type == 'dart'){
                # Only look at 'good-dart' series which we use for comparison
                if(!event_stats$good_dart[i]) next
            }else{
                stop('unknown data_type, must be a bug')
            }
            # No need to do this for every model.
            if(!(event_stats$model_type[i] == 'Manning0.035')) next
            # We will still have repetition by inverion

            # Find the RDS data containing this time-series
            k = match(event_stats$run_name[i], names(target_RDS_data))
            j = match(event_stats$sites[i], names(target_RDS_data[[k]]))

            plot_ylim = c(event_stats$data_min[i], event_stats$data_max[i])*1.1
            if(event_stats$site_and_event[i] == "Chile-1960_Cronulla_CSIRO_Fisheries_1960"){
                # Work-around for the Chile 1960 Cronulla series, which has NA for
                # max/min because the observations only last ~ 30hrs.
                plot_ylim = c(event_stats$data_min_24[i], event_stats$data_max_24[i])* 1.1
            }

            plot_xlim = c(target_RDS_data[[k]][[j]]$model_start_julian,
                          target_RDS_data[[k]][[j]]$model_start_julian + 2.5)
            dat_x = target_RDS_data[[k]][[j]]$event_data$obs$juliant
            dat_y = target_RDS_data[[k]][[j]]$event_data$obs$resid
            
            plot(dat_x, dat_y, t='l', xlim=plot_xlim, ylim=plot_ylim, xlab='Day', ylab='Stage (m)')
            title(event_stats$site_and_event[i])
            title(sub=event_stats$source_name[i])
            all_col = blues9[4:9] 
            # Add computed maxima for various times corresponding to max/min statistics in event_stats
            line_times = plot_xlim[1] + 
                c(0, 
                  event_stats$model_time_to_arrive[i]/(3600*24) + c(0, 8, 12)/24, 
                  1, 1.5, 2.5)
            points(line_times[c(2,3)], rep(event_stats$data_max_8hrs[i], 2), t='l', col=all_col[1], lwd=0.5)
            points(line_times[c(2,4)], rep(event_stats$data_max_12hrs[i], 2), t='l', col=all_col[2], lwd=0.5)
            points(line_times[c(2,5)], rep(event_stats$data_max_24[i], 2), t='l', col=all_col[3], lwd=0.5)
            points(line_times[c(2,6)], rep(event_stats$data_max_36[i], 2), t='l', col=all_col[4], lwd=0.5)
            points(line_times[c(2,7)], rep(event_stats$data_max[i], 2), t='l', col=all_col[5], lwd=0.5)
            points(line_times[c(6,7)], rep(event_stats$data_max_lastday[i], 2), t='l', col=all_col[6], lwd=0.5)

            # Add times that are significant to interpret the above
            abline(v=line_times[2:length(line_times)],
                   col='orange', lty=c('dotted', 'dotted', 'dotted', 'solid', 'solid', 'solid'))
            abline(h=0, col='orange')

            if(use_model){
                points(target_RDS_data[[k]][[j]]$model_start_julian + target_RDS_data[[k]][[j]]$model_time/(3600*24),
                       target_RDS_data[[k]][[j]]$model_stage, t='l', col='green')
            }

        }
        dev.off()
    }
}
