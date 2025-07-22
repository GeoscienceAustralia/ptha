#
# Extract statistics comparing models and observations
#

parse_gauge_outputs<-function(DOWNSAMPLE_MODEL_TO_DATA){
    # The purpose of this function is to extract various statistics about gauges vs models: 
    # - if(DOWNSAMPLE_MODEL_TO_DATA == FALSE) 
    #       Compare the 'raw' model time series (stored at 30s time interval)
    #       with data (recorded at a range of time intervals).
    # - if(DOWNSAMPLE_MODEL_TO_DATA == TRUE) 
    #      Re-sample the model time series before computing the statistics.
    #      This involves a 1-min smooth of the model, followed by downsampling to the tide-gauge
    #      resolution. The idea is to check if our results might be affected by
    #      the tide-gauge smoothing and reduced sampling frequency.
    #      For example, when using 5 minute tide-gauge data, we down-sample the
    #      model to 5 minutes, to represent the under-sampling of the observations (Rabinovich 2011).
    #      All nearshore tide-gauge data is assumed to represent a 1-min smooth (at
    #      minimum) because the nearshore gauges generally have smoothing [read the
    #      metadata!]. In reality the situation is more complex (e.g. Davies et al 2024 Hunga-Tonga paper)
    #      but this will offer some idea as to whether it's important.

    # Work with the RDS files containing gauge data and model results 
    target_RDS = Sys.glob('../REDUCED_OUTPUTS/*/RUN*/gau*.RDS')

    ONE_HOUR=3600.0

    # When 
    #      abs(modelled-stage) > ( max(abs(modelled-stage))*WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD )
    # we say the tsunami has 'arrived'. This begins the time-period over which
    # we might compare model and data (although that is also affected by other
    # model parameters below).
    WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD = 5e-04

    source('time_domain_hybrid_norm.R')
        # Function giving goodness-of-fit of 2 time-series, like used 
        # in Davies (2019) at DARTs. 

    source('downsample_and_smooth_model_at_nearshore_tide_gauge.R')
        # Function that can downsample a model to be more consistent with
        # the data at nearshore tide-gauges (particularly important for the 15
        # minute records).
        # This function will only be called when (DOWNSAMPLE_MODEL_TO_DATA=TRUE)

    source('manually_despike_data.R')
        # Remove spikes from a couple of the time-series.

    is_gauge_in_domain_env = new.env()
    source('is_gauge_in_highres_domain.R', local=is_gauge_in_domain_env)
        # Determine whether a gauge is in the high-res domain

    source('good_nearshore_data_definition.R')
        # Classify data as "good" and "nearshore", which is usually the data we
        # want to analyse.

    source('create_times_to_start_tide_gauge_comparison_with_models.R')
        # A routine that defines a consistent start time for comparing observed
        # and modelled gauge time-series.
        # It relies on the script already having been run once -- if it does
        # not exist then the code will do separate calculations for each
        # modelled time series, which can be used to create consistent start
        # times (see the above script for details)
    get_comparison_start_times_from_function = exists('get_julian_time_at_which_to_start_comparison_with_tide_gauge')
    
    #' Extract model and data statistics at a single site
    #'
    #' This will be applied to all sites.
    #'
    #' @param site_data a list with the model and data at a site.
    #' @param gauge_name the name of the file that we read site_data from
    #' @param filter_model_like_data If TRUE, then pre-process the model result
    #' to add some smoothing/down-sampling like we expect from the data. This
    #' is important for 15minute sampled data, and for some other cases where
    #' the tsunami wave frequency is very high. In this case one must also pass
    #' gauge_name, which will contain information about the data source.
    #' @param target_RDS_file associated RDS filename. Used to extract information
    #' on the domain high-resolution areas and whether the gauge coordinate is inside.
    #' @return a bunch of statistics
    #'
    get_statistics_at_site<-function(site_data, gauge_name, filter_model_like_data, target_RDS_file){

        # The observed gauge time-series
        obs = site_data$event_data$obs
        obs_coordinate = site_data$event_data$coord
        model_coordinate = site_data$model_gaugeloc[1:2]

        # Store the distance between the tide gauge and the nearest model gauge.
        # These can be further than we want (e.g. due to corrections of tide-gauge
        # data coordinates, or gauges we obtained after setting up the model)
        library(geosphere)
        distance_to_gauge = distHaversine(
            cbind(obs_coordinate[1], obs_coordinate[2]), 
            cbind(model_coordinate[1], model_coordinate[2]))

        # Object with the "raw" model time-series
        model = data.frame(
            juliant = as.difftime(site_data$model_time/(ONE_HOUR*24), units='days') + site_data$model_start_julian,
            resid = site_data$model_stage)

        if(filter_model_like_data){
            # Try to undersample the model like the tide-gauge data 
            # At DARTs, we do not adjust the model
            if(gauge_name == ''){
                stop('Must provide a non-empty gauge_name if filter_model_like_data=TRUE')
            }

            if(!grepl('DART', gauge_name)){
                # Apply processing only to nearshore tidal gauges, not DART
                model = downsample_and_smooth_model_at_nearshore_tide_gauge(model, gauge_name)
            }
        }

        # Extract info on high-res modelled regions via a flag from the filename
        tmp = strsplit(basename(target_RDS_file), split='-')[[1]]
        model_resolution_tag = gsub('.RDS', '', tmp[length(tmp)], fixed=TRUE)
        stopifnot(model_resolution_tag %in% c(
            'highres_australia', 'highres_NSW', 'highres_perth', 'highres_australiaSWWA', 
            'highres_SWWA', 'highres_NWWA', 'highres_australiaWA', 'highres_WA'))

        # Determine if the gauge is in a highres domain on this model.
        is_gauge_in_highres_domain = is_gauge_in_domain_env$is_gauge_in_highres_domain(
            obs_coordinate[1], obs_coordinate[2], model_resolution_tag)

        n = length(model$juliant)
        model_range = diff(range(model$resid))

        # Last non-missing observation in data. 
        final_obs_juliant = max(obs$juliant[!is.na(obs$resid)])

        if(n == 0 | model_range == 0){
            # Fail gracefully in the case where we didn't extract a model time-series.
            # This has happened previously when a gauge is "on land" in the
            # model due to the prescence of coarse cells.

            model_range = c(NA, NA) 
            data_range = c(NA, NA) 
            model_range_during_obs = c(NA, NA)
            data_range_during_obs = c(NA, NA)

            model_rms = NA
            data_rms = NA
            model_rms_during_obs = NA
            data_rms_during_obs = NA

            model_range8hrs = c(NA, NA) 
            data_range8hrs = c(NA, NA) 
            model_range8hrs_during_obs = c(NA, NA)
            data_range8hrs_during_obs = c(NA, NA)

            model_rms8hrs = NA
            data_rms8hrs = NA
            model_rms8hrs_during_obs = NA
            data_rms8hrs_during_obs = NA

            model_range12hrs = c(NA, NA) 
            data_range12hrs = c(NA, NA) 
            model_range12hrs_during_obs = c(NA, NA) 
            data_range12hrs_during_obs = c(NA, NA) 

            model_rms12hrs = NA
            data_rms12hrs = NA
            model_rms12hrs_during_obs = NA
            data_rms12hrs_during_obs = NA

            model_range24hrs = c(NA, NA) 
            data_range24hrs = c(NA, NA) 
            model_range24hrs_during_obs = c(NA, NA) 
            data_range24hrs_during_obs = c(NA, NA) 

            model_rms24hrs = NA
            data_rms24hrs = NA
            model_rms24hrs_during_obs = NA
            data_rms24hrs_during_obs = NA

            model_range24 = c(NA, NA) 
            data_range24 = c(NA, NA) 
            model_range24_during_obs = c(NA, NA) 
            data_range24_during_obs = c(NA, NA) 

            model_rms24 = NA
            data_rms24 = NA
            model_rms24_during_obs = NA
            data_rms24_during_obs = NA

            model_range36 = c(NA, NA) 
            data_range36 = c(NA, NA) 
            model_range36_during_obs = c(NA, NA) 
            data_range36_during_obs = c(NA, NA) 

            model_rms36 = NA
            data_rms36 = NA
            model_rms36_during_obs = NA
            data_rms36_during_obs = NA

            model_range_lastday = c(NA, NA) 
            data_range_lastday = c(NA, NA) 
            model_range_lastday_during_obs = c(NA, NA) 
            data_range_lastday_during_obs = c(NA, NA) 

            model_rms_lastday = NA
            data_rms_lastday = NA
            model_rms_lastday_during_obs = NA
            data_rms_lastday_during_obs = NA

            model_range_24to36 = c(NA, NA) 
            data_range_24to36 = c(NA, NA) 
            model_range_24to36_during_obs = c(NA, NA) 
            data_range_24to36_during_obs = c(NA, NA) 

            model_rms_24to36 = NA
            data_rms_24to36 = NA
            model_rms_24to36_during_obs = NA
            data_rms_24to36_during_obs = NA

            model_arrival_time = NA
            model_time_to_arrive = NA
            obs_max_time_to_arrive = NA

            model_time_of_max = NA
            model_time_of_max_during_obs = NA
            data_time_of_max = NA

            time_domain_hybrid_norm_stat = NA
            time_domain_hybrid_norm_time_offset = NA
            time_domain_hybrid_norm_stat2 = NA
            time_domain_hybrid_norm_time_offset2 = NA
            time_domain_hybrid_norm_stat3 = NA
            time_domain_hybrid_norm_time_offset3 = NA
            time_domain_hybrid_norm_stat4 = NA
            time_domain_hybrid_norm_time_offset4 = NA

        }else{
            # Regular case -- compute all the statistics we need
            #
            # Look at observations within the model time range.  
            # For some time-series care is needed to focus on the tsunami time, not other waves
            # before arrival, and we account for this here (and use plotting to see that it is
            # making sense)
            wave_threshold = WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD *max(abs(model$resid))
            k = min(which(abs(model$resid) > wave_threshold))
            model_arrival_time = model$juliant[k]
            model_time_to_arrive = site_data$model_time[k]

            if(get_comparison_start_times_from_function){
                # Get a consistent time to start comparing the model with data at each tide gauge.
                # This is important because it means the "data" statistics are
                # the same for all model runs corresponding to the the
                # event/tide-gauge.
                # This uses the minimum of all modelled waveforms for each event/tide-gauge pair.
                model_vs_gauge_start_time = julian(strptime('1970-01-01', format='%Y-%m-%d', tz='Etc/GMT')) + # Zero julian days (trick to get the right tz)
                    as.difftime(
                        get_julian_time_at_which_to_start_comparison_with_tide_gauge(
                            gauge_name, as.numeric(site_data$model_start_julian)), 
                        units='days')
    
                # Fallback when new gauges have been introduced
                #if(is.na(model_vs_gauge_start_time)) model_vs_gauge_start_time = model$juliant[k]
            }else{
                # Get a time to start comparing the model with data which
                # varies for each simulation. Consistent times can be created
                # from these initial estimates using
                # "create_times_to_start_tide_gauge_comparison_with_models.R"
                model_vs_gauge_start_time = model$juliant[k]
            }

            # Indices of data which can be compared with model
            obs_keep = which( (obs$juliant >= model_vs_gauge_start_time) & 
                              (obs$juliant <= model$juliant[n]))

            # Range of model and de-tided data
            data_range = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range = range(model$resid)
            data_range_during_obs = data_range
            model_range_during_obs = range(model$resid[model$juliant <= final_obs_juliant])

            # Root mean square of the model and the data. 
            data_rms  = sqrt( mean(obs$resid[obs_keep]**2, na.rm=TRUE) )
            model_rms = sqrt( mean(model$resid[k:n]**2   , na.rm=TRUE) )
            data_rms_during_obs  = data_rms
            n_before_model_finished = max(which(model$juliant <= final_obs_juliant))
            model_rms_during_obs = sqrt( mean(model$resid[k:n_before_model_finished]**2   , na.rm=TRUE) )

            # Time of the observed maxima 
            tmp = which( obs$resid[obs_keep] == max(data_range) )[1]
            data_time_of_max = obs$juliant[obs_keep[tmp]]
            obs_max_time_to_arrive = (data_time_of_max - model$juliant[1])*ONE_HOUR*24

            #
            # Results when we cut off the final 1 day (i.e. run 36 hours instead of 60)
            #
            obs_keep   = which( (obs$juliant   >= model_vs_gauge_start_time) & (obs$juliant   <= (model$juliant[n]-1)) )
            model_keep = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model$juliant[n]-1)) )
            model_keep_during_obs = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model$juliant[n]-1)) & (model$juliant <= final_obs_juliant ) )
            model_range36 = range(model$resid[model_keep])
            data_range36 = range(obs$resid[obs_keep], na.rm=TRUE)
            data_range36_during_obs = data_range
            model_range36_during_obs = range(model$resid[model_keep_during_obs])
    
            data_rms36  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms36 = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms36_during_obs = data_rms36
            model_rms36_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  

            #
            # Results when we cut off the final 1.5 day (i.e. 24 hours instead of 60)
            #
            obs_keep   = which( (obs$juliant   >= model_vs_gauge_start_time) & (obs$juliant   <= (model$juliant[n]-1.5)) )
            model_keep = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model$juliant[n]-1.5)) )
            model_keep_during_obs = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model$juliant[n]-1.5)) & (model$juliant <= final_obs_juliant) )
            model_range24 = range(model$resid[model_keep])
            data_range24 = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range24_during_obs = range(model$resid[model_keep_during_obs])
            data_range24_during_obs = data_range24

            data_rms24  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms24 = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms24_during_obs = data_rms24
            model_rms24_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  

            #
            # Results within 8 hours of the time at which we start comparing the model and data
            #
            obs_keep   = which( (obs$juliant   >= model_vs_gauge_start_time) & (obs$juliant   <= (model_vs_gauge_start_time+8/24)) )
            model_keep = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+8/24)) )
            model_keep_during_obs = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+8/24)) & (model$juliant <= final_obs_juliant) )
            model_range8hrs = range(model$resid[model_keep])
            data_range8hrs = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range8hrs_during_obs = range(model$resid[model_keep_during_obs])
            data_range8hrs_during_obs = data_range8hrs

            data_rms8hrs  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms8hrs = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms8hrs_during_obs = data_rms8hrs
            model_rms8hrs_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )

            #
            # Results within 12 hours of the time at which we start comparing the model and data
            #
            obs_keep   = which( (obs$juliant   >= model_vs_gauge_start_time) & (obs$juliant   <= (model_vs_gauge_start_time+12/24)) )
            model_keep = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+12/24)) )
            model_keep_during_obs = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+12/24)) & (model$juliant <= final_obs_juliant))
            model_range12hrs = range(model$resid[model_keep])
            data_range12hrs = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range12hrs_during_obs = range(model$resid[model_keep_during_obs])
            data_range12hrs_during_obs = data_range12hrs

            data_rms12hrs  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms12hrs = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms12hrs_during_obs = data_rms12hrs
            model_rms12hrs_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  

            #
            # Results within 24 hours of the time at which we start comparing the model and data
            #
            obs_keep   = which( (obs$juliant   >= model_vs_gauge_start_time) & (obs$juliant   <= (model_vs_gauge_start_time+24/24)) )
            model_keep = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+24/24)) )
            model_keep_during_obs = which( (model$juliant >= model_vs_gauge_start_time) & (model$juliant <= (model_vs_gauge_start_time+24/24)) & (model$juliant <= final_obs_juliant) )
            model_range24hrs = range(model$resid[model_keep])
            data_range24hrs = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range24hrs_during_obs = range(model$resid[model_keep_during_obs])
            data_range24hrs_during_obs = data_range24hrs

            data_rms24hrs  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms24hrs = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms24hrs_during_obs = data_rms24hrs
            model_rms24hrs_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  


            #
            # Results RESTRICTED TO the final 1 day (i.e. from 36-60 hours)
            #
            obs_keep   = which( (obs$juliant   >= model$juliant[n]-1) & (obs$juliant   <= (model$juliant[n])) )
            model_keep = which( (model$juliant >= model$juliant[n]-1) & (model$juliant <= (model$juliant[n])) )
            model_keep_during_obs = which( (model$juliant >= model$juliant[n]-1) & (model$juliant <= (model$juliant[n])) & (model$juliant <= final_obs_juliant) )
            model_range_lastday = range(model$resid[model_keep])
            data_range_lastday = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range_lastday_during_obs = range(model$resid[model_keep_during_obs])
            data_range_lastday_during_obs = data_range_lastday
            
            data_rms_lastday  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms_lastday = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms_lastday_during_obs = data_rms_lastday
            model_rms_lastday_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  

            #
            # Results RESTRICTED TO 24-36 hours, i.e. 1-1.5 days.
            #
            obs_keep   = which( (obs$juliant   >= model$juliant[n]-1.5) & (obs$juliant   <= (model$juliant[n]-1)) )
            model_keep = which( (model$juliant >= model$juliant[n]-1.5) & (model$juliant <= (model$juliant[n]-1)) )
            model_keep_during_obs = which( (model$juliant >= model$juliant[n]-1.5) & (model$juliant <= (model$juliant[n]-1)) & (model$juliant <= final_obs_juliant) )
            model_range_24to36 = range(model$resid[model_keep])
            data_range_24to36 = range(obs$resid[obs_keep], na.rm=TRUE)
            model_range_24to36_during_obs = range(model$resid[model_keep_during_obs])
            data_range_24to36_during_obs = data_range_24to36

            data_rms_24to36  = sqrt( mean(obs$resid[obs_keep]**2    , na.rm=TRUE) )
            model_rms_24to36 = sqrt( mean(model$resid[model_keep]**2, na.rm=TRUE) )  
            data_rms_24to36_during_obs = data_rms_24to36
            model_rms_24to36_during_obs = sqrt( mean(model$resid[model_keep_during_obs]**2, na.rm=TRUE) )  

            # Time of modelled maxima
            k = which.max(model$resid)
            model_time_of_max = site_data$model_time[k]
            model_time_of_max_during_obs = which.max(model$resid * (model$juliant <= final_obs_juliant))

            #
            # Goodness of fit statistic.
            #
            # Comparison up to 10 hrs after model arrival. Note 'model_vs_gauge_start_time'
            # is an early start time (when derived from all model scenarios), so no need to shorten it further.
            trng = c(model_vs_gauge_start_time, model_vs_gauge_start_time + 10/24)
            # Allow for truncated data
            if(trng[2] > final_obs_juliant) trng[2] = final_obs_juliant
            # Force the model/observed times to overlap
            obs_keep = which( (obs$juliant >= trng[1]) & 
                              (obs$juliant <= trng[2]))
            model_keep = which( (model$juliant >= trng[1]) &
                                (model$juliant <= trng[2]))
            # Get model/data timeseries, with times in seconds for the GOF code.
            model_time_seconds = site_data$model_time[model_keep]
            model_stage = site_data$model_stage[model_keep]
            data_time_seconds = as.numeric(obs$juliant[obs_keep] - site_data$model_start_julian)*3600*24
            data_stage = obs$resid[obs_keep]
            # Optimal time-shift of +- 20min applied to GOF stat
            tdhn = try(gauge_similarity_time_domain(
                data_time_seconds,
                data_stage,
                model_time_seconds,
                model_stage,
                interp_dt = 15, # seconds
                allowed_lag_minutes = c(-20, 20),
                detailed=TRUE))
            tdhn_2 = try(gauge_similarity_time_domain(
                data_time_seconds,
                data_stage,
                model_time_seconds,
                model_stage,
                interp_dt = 15, # seconds
                allowed_lag_minutes = c(-20, 20),
                statistic_variant = 'no_model_in_denominator',
                detailed=TRUE))
    
            #
            # Repeat the GOF calculation, but limit the time period to at most
            # 3 hours after the observed maxima.
            # This will matter for tsunamis with leading large waves --
            # focussing the comparison time on the more important part.
            #
            if(trng[2] > data_time_of_max + 3/24) trng[2] = data_time_of_max + 3/24
            obs_keep = which( (obs$juliant >= trng[1]) & 
                              (obs$juliant <= trng[2]))
            model_keep = which( (model$juliant >= trng[1]) &
                                (model$juliant <= trng[2]))
            # Get model/data timeseries, with times in seconds for the GOF code.
            model_time_seconds = site_data$model_time[model_keep]
            model_stage = site_data$model_stage[model_keep]
            data_time_seconds = as.numeric(obs$juliant[obs_keep] - site_data$model_start_julian)*3600*24
            data_stage = obs$resid[obs_keep]
            # Optimal time-shift of +- 20min applied to GOF stat
            tdhn_3 = try(gauge_similarity_time_domain(
                data_time_seconds,
                data_stage,
                model_time_seconds,
                model_stage,
                interp_dt = 15, # seconds
                allowed_lag_minutes = c(-20, 20),
                detailed=TRUE))
            tdhn_4 = try(gauge_similarity_time_domain(
                data_time_seconds,
                data_stage,
                model_time_seconds,
                model_stage,
                interp_dt = 15, # seconds
                allowed_lag_minutes = c(-20, 20),
                statistic_variant = 'no_model_in_denominator',
                detailed=TRUE))

            if(is(tdhn, 'try-error') | is(tdhn_2, 'try-error')){
                # Something failed
                print(paste0('Failure in time_domain_hybrid_norm at ', gauge_name, ', with file ', target_RDS_file))
                #browser()
                time_domain_hybrid_norm_stat = NA
                time_domain_hybrid_norm_time_offset = NA
                time_domain_hybrid_norm_stat2 = NA
                time_domain_hybrid_norm_time_offset2 = NA
                time_domain_hybrid_norm_stat3 = NA
                time_domain_hybrid_norm_time_offset3 = NA
                time_domain_hybrid_norm_stat4 = NA
                time_domain_hybrid_norm_time_offset4 = NA
            }else{
                time_domain_hybrid_norm_stat = tdhn$objective
                time_domain_hybrid_norm_time_offset = tdhn$minimum
                time_domain_hybrid_norm_stat2 = tdhn_2$objective
                time_domain_hybrid_norm_time_offset2 = tdhn_2$minimum
                time_domain_hybrid_norm_stat3 = tdhn_3$objective
                time_domain_hybrid_norm_time_offset3 = tdhn_3$minimum
                time_domain_hybrid_norm_stat4 = tdhn_4$objective
                time_domain_hybrid_norm_time_offset4 = tdhn_4$minimum
            }
        }

        # For the output, make sensible names
        output = list(
            model_min = model_range[1],
            model_max = model_range[2],
            model_min_during_obs = model_range_during_obs[1],
            model_max_during_obs = model_range_during_obs[2],

            model_stgrng = diff(model_range[1:2]),
            model_stgrng_during_obs = diff(model_range_during_obs[1:2]),

            model_rms = model_rms,
            model_rms_during_obs = model_rms_during_obs,

            data_min = data_range[1],
            data_max = data_range[2],
            data_min_during_obs = data_range_during_obs[1],
            data_max_during_obs = data_range_during_obs[2],

            data_stgrng = diff(data_range[1:2]),
            data_stgrng_during_obs = diff(data_range_during_obs[1:2]),

            data_rms = data_rms,
            data_rms_during_obs = data_rms_during_obs,

            model_min_8hrs = model_range8hrs[1],
            model_max_8hrs = model_range8hrs[2],
            model_min_8hrs_during_obs = model_range8hrs_during_obs[1],
            model_max_8hrs_during_obs = model_range8hrs_during_obs[2],

            model_stgrng_8hrs = diff(model_range8hrs[1:2]),
            model_stgrng_8hrs_during_obs = diff(model_range8hrs_during_obs[1:2]),

            model_rms_8hrs = model_rms8hrs,
            model_rms_8hrs_during_obs = model_rms8hrs_during_obs,
                                           
            data_min_8hrs = data_range8hrs[1],
            data_max_8hrs = data_range8hrs[2],
            data_min_8hrs_during_obs = data_range8hrs_during_obs[1],
            data_max_8hrs_during_obs = data_range8hrs_during_obs[2],

            data_stgrng_8hrs = diff(data_range8hrs[1:2]),
            data_stgrng_8hrs_during_obs = diff(data_range8hrs_during_obs[1:2]),

            data_rms_8hrs = data_rms8hrs,
            data_rms_8hrs_during_obs = data_rms8hrs_during_obs,
                                                
            model_min_12hrs = model_range12hrs[1],
            model_max_12hrs = model_range12hrs[2],
            model_min_12hrs_during_obs = model_range12hrs_during_obs[1],
            model_max_12hrs_during_obs = model_range12hrs_during_obs[2],

            model_stgrng_12hrs = diff(model_range12hrs[1:2]),
            model_stgrng_12hrs_during_obs = diff(model_range12hrs_during_obs[1:2]),

            model_rms_12hrs = model_rms12hrs, 
            model_rms_12hrs_during_obs = model_rms12hrs_during_obs, 
                                             
            data_min_12hrs = data_range12hrs[1],
            data_max_12hrs = data_range12hrs[2],
            data_min_12hrs_during_obs = data_range12hrs_during_obs[1],
            data_max_12hrs_during_obs = data_range12hrs_during_obs[2],

            data_stgrng_12hrs = diff(data_range12hrs[1:2]),
            data_stgrng_12hrs_during_obs = diff(data_range12hrs_during_obs[1:2]),
            data_rms_12hrs = data_rms12hrs, 
            data_rms_12hrs_during_obs = data_rms12hrs_during_obs, 

            model_min_24hrs = model_range24hrs[1],
            model_max_24hrs = model_range24hrs[2],
            model_min_24hrs_during_obs = model_range24hrs_during_obs[1],
            model_max_24hrs_during_obs = model_range24hrs_during_obs[2],

            model_stgrng_24hrs = diff(model_range24hrs[1:2]),
            model_stgrng_24hrs_during_obs = diff(model_range24hrs_during_obs[1:2]),

            model_rms_24hrs = model_rms24hrs, 
            model_rms_24hrs_during_obs = model_rms24hrs_during_obs, 
                                             
            data_min_24hrs = data_range24hrs[1],
            data_max_24hrs = data_range24hrs[2],
            data_min_24hrs_during_obs = data_range24hrs_during_obs[1],
            data_max_24hrs_during_obs = data_range24hrs_during_obs[2],

            data_stgrng_24hrs = diff(data_range24hrs[1:2]),
            data_stgrng_24hrs_during_obs = diff(data_range24hrs_during_obs[1:2]),

            data_rms_24hrs = data_rms24hrs, 
            data_rms_24hrs_during_obs = data_rms24hrs_during_obs, 

            model_min_24 = model_range24[1],
            model_max_24 = model_range24[2],
            model_min_24_during_obs = model_range24_during_obs[1],
            model_max_24_during_obs = model_range24_during_obs[2],

            model_stgrng_24 = diff(model_range24[1:2]),
            model_stgrng_24_during_obs = diff(model_range24_during_obs[1:2]),

            model_rms_24 = model_rms24, 
            model_rms_24_during_obs = model_rms24_during_obs, 
                                       
            data_min_24 = data_range24[1],
            data_max_24 = data_range24[2],
            data_min_24_during_obs = data_range24_during_obs[1],
            data_max_24_during_obs = data_range24_during_obs[2],

            data_stgrng_24 = diff(data_range24[1:2]),
            data_stgrng_24_during_obs = diff(data_range24_during_obs[1:2]),

            data_rms_24 = data_rms24,
            data_rms_24_during_obs = data_rms24_during_obs,
                                          
            model_min_36 = model_range36[1],
            model_max_36 = model_range36[2],
            model_min_36_during_obs = model_range36_during_obs[1],
            model_max_36_during_obs = model_range36_during_obs[2],

            model_stgrng_36 = diff(model_range36[1:2]),
            model_stgrng_36_during_obs = diff(model_range36_during_obs[1:2]),

            model_rms_36 = model_rms36,
            model_rms_36_during_obs = model_rms36_during_obs,
                                       
            data_min_36 = data_range36[1],
            data_max_36 = data_range36[2],
            data_min_36_during_obs = data_range36_during_obs[1],
            data_max_36_during_obs = data_range36_during_obs[2],

            data_stgrng_36 = diff(data_range36[1:2]),
            data_stgrng_36_during_obs = diff(data_range36_during_obs[1:2]),

            data_rms_36 = data_rms36,
            data_rms_36_during_obs = data_rms36_during_obs,
                                                      
            model_min_lastday = model_range_lastday[1],
            model_max_lastday = model_range_lastday[2],
            model_min_lastday_during_obs = model_range_lastday_during_obs[1],
            model_max_lastday_during_obs = model_range_lastday_during_obs[2],

            model_stgrng_lastday = diff(model_range_lastday[1:2]),
            model_stgrng_lastday_during_obs = diff(model_range_lastday_during_obs[1:2]),

            model_rms_lastday = model_rms_lastday,
            model_rms_lastday_during_obs = model_rms_lastday_during_obs,
                                                  
            data_min_lastday = data_range_lastday[1],
            data_max_lastday = data_range_lastday[2],
            data_min_lastday_during_obs = data_range_lastday_during_obs[1],
            data_max_lastday_during_obs = data_range_lastday_during_obs[2],

            data_stgrng_lastday = diff(data_range_lastday[1:2]),
            data_stgrng_lastday_during_obs = diff(data_range_lastday_during_obs[1:2]),

            data_rms_lastday = data_rms_lastday,
            data_rms_lastday_during_obs = data_rms_lastday_during_obs,
             
            model_min_24to36 = model_range_24to36[1],
            model_max_24to36 = model_range_24to36[2],
            model_min_24to36_during_obs = model_range_24to36_during_obs[1],
            model_max_24to36_during_obs = model_range_24to36_during_obs[2],

            model_stgrng_24to36 = diff(model_range_24to36[1:2]),
            model_stgrng_24to36_during_obs = diff(model_range_24to36_during_obs[1:2]),

            model_rms_24to36 = model_rms_24to36,
            model_rms_24to36_during_obs = model_rms_24to36_during_obs,

            data_min_24to36 = data_range_24to36[1],
            data_max_24to36 = data_range_24to36[2],
            data_min_24to36_during_obs = data_range_24to36_during_obs[1],
            data_max_24to36_during_obs = data_range_24to36_during_obs[2],

            data_stgrng_24to36 = diff(data_range_24to36[1:2]),
            data_stgrng_24to36_during_obs = diff(data_range_24to36_during_obs[1:2]),

            data_rms_24to36 = data_rms_24to36,
            data_rms_24to36_during_obs = data_rms_24to36_during_obs,

            model_arrival_time = model_arrival_time,
            model_vs_gauge_start_time = model_vs_gauge_start_time,
            model_time_to_arrive = model_time_to_arrive,
            obs_max_time_to_arrive = obs_max_time_to_arrive,
             
            model_time_of_max = model_time_of_max,
            model_time_of_max_during_obs = model_time_of_max_during_obs,
            data_time_of_max = data_time_of_max,

            # Gauge location
            obs_lon = obs_coordinate[1],
            obs_lat = obs_coordinate[2],
            model_lon = model_coordinate[1],
            model_lat = model_coordinate[2],
            distance_to_gauge = distance_to_gauge,

            # This string has a 1:1 mapping with the high-res domain locations.
            model_resolution_tag = model_resolution_tag,
            is_gauge_in_highres_domain = is_gauge_in_highres_domain,

            # GOF statistic and time offset for data that gives best match to model
            # Similar approach to Davies (2019) at DARTs
            time_domain_hybrid_norm_stat = time_domain_hybrid_norm_stat,
            time_domain_hybrid_norm_time_offset = time_domain_hybrid_norm_time_offset,
            time_domain_hybrid_norm_stat2 = time_domain_hybrid_norm_stat2,
            time_domain_hybrid_norm_time_offset2 = time_domain_hybrid_norm_time_offset2,
            time_domain_hybrid_norm_stat3 = time_domain_hybrid_norm_stat3,
            time_domain_hybrid_norm_time_offset3 = time_domain_hybrid_norm_time_offset3,
            time_domain_hybrid_norm_stat4 = time_domain_hybrid_norm_stat4,
            time_domain_hybrid_norm_time_offset4 = time_domain_hybrid_norm_time_offset4,
            
            model_start = site_data$model_start_julian,
            model_end = max(model$juliant),
            data_end = final_obs_juliant
            )

        # Make sure there are not NULL values in the output
        for(i in 1:length(output)){
            if(any(is.null(output[[i]]))){
                print(gauge_name)
                print(target_RDS_file)
                print(names(output)[i])
                print(output[[i]])
                stop('NULL value in output - halting as this would cause problems when unlisting later')
            }
        }
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
    parse_parallel<-function(i){
        # Read the gauges
        x = readRDS(target_RDS[i])
        # In a couple of cases we clean up spikes in the record
        x = manually_despike_data(x, target_RDS[i])
        # Compute the statistics
        all_site_data = mapply(get_statistics_at_site, 
           x, names(x), 
           MoreArgs=list(filter_model_like_data=DOWNSAMPLE_MODEL_TO_DATA, target_RDS_file=target_RDS[i]), 
           SIMPLIFY=FALSE)
        names(all_site_data) = names(x)
        # Store what we need
        stat_store = all_site_data
        target_RDS_data = x
        output = list(target_RDS_data=target_RDS_data, 
                      stat_store=stat_store) 
        # Be conservative with memory
        rm(x, all_site_data, target_RDS_data, stat_store)
        gc()

        return(output)
    }

    #library(parallel)
    parse_all = lapply(1:length(target_RDS), parse_parallel)
    target_RDS_data = lapply(parse_all, function(x) x$target_RDS_data)
    names(target_RDS_data) = target_RDS
    stat_store = lapply(parse_all, function(x) x$stat_store)
    names(stat_store) = target_RDS
    rm(parse_all)
    gc()
    # The above 'stat_store' object is a nested list
    #
    # Convert to a table.  
    #
    # Helper function
    unpack_event<-function(event){
        all_sites_stats = lapply(event, as.data.frame)
        out_df = do.call(rbind, all_sites_stats)
        # Include the site name
        out_df = cbind(data.frame(sites=names(event), stringsAsFactors=FALSE), 
                       out_df)
        rownames(out_df) = NULL
        return(out_df)
    }
    # Helper function 2
    unpack_all_events<-function(stat_store){ 
        event_dfs = lapply(stat_store, unpack_event)
        # Include the "run_name", i.e. file where we read the model results.
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

    # Extract a bunch of key variables from the model folder name. If I had been
    # more disciplined about the naming convention, this could have been easier!
    model_type_from_run_name<-function(run_name){

        if(length(run_name) != 1){
            stop('run_name should have length(1), the code is not vectorized')
        }

        # Convert the file path to lower case -- the folder/file names contain
        # information we need, and we will use this to do the extraction.
        # Because the name structure does not include a consistent separator,
        # it is easiest to edit the run_name2 as we go [removing information
        # we've already looked at].
        run_name2 = tolower(run_name)

        # Extract the slip type (HS / VAUS / FAUS) and the run type
        # (random_like_historic / nonrandom_like_historic)
        slip_model_type = NA
        run_type = NA
        batch_number = 1 # Maybe updated below

        if(grepl('_random_like_historic_', run_name2)){
            #
            # Treat random PTHA18 scenarios [every case for some analyses]
            #

            run_type = 'random_like_historic'

            # Deal with other batches of random scenarios, which mirror the
            # first and have '_batch[2-6]' in the filename. 
            BATCHES_THAT_ARE_NUMBERED=seq(2,6) # the first batch didn't use this naming convention
            for(bn in BATCHES_THAT_ARE_NUMBERED){
                batch_match = paste0('-batch', bn) #  = "-batch2", "-batch3", ...
                if(grepl(batch_match, run_name2, fixed=TRUE)){
                    # We are in batch bn
                    batch_number = bn
                    # Remove -batchN from run_name2 (helps with processing)
                    run_name2 = gsub(batch_match, '', run_name2, fixed=TRUE)
                }
            }

            # Strip first part of run_name2 (helps with processing)
            run_name2 = gsub('ptha18_random_like_historic_', '', run_name2, 
                             fixed=TRUE)

            # Identify the slip type
            if(grepl('_heterogeneous_slip_', run_name2, fixed=TRUE)){

                slip_model_type = 'HS'
                run_name2 = gsub('_heterogeneous_slip', '', run_name2, 
                                 fixed=TRUE) # Remove slip type from run_name2

            }else if(grepl('_variable_area_uniform_slip_', run_name2, fixed=TRUE)){

                slip_model_type = 'VAUS'
                run_name2 = gsub('_variable_area_uniform_slip', '', run_name2, 
                                 fixed=TRUE)# Remove slip type from run_name2

            }else if(grepl('_uniform_slip_', run_name2, fixed=TRUE)){

                # Important that this comes AFTER _variable_area_uniform_slip_
                if(grepl('_variable_area_uniform_slip', run_name2)){
                    stop('Must search for VAUS before searching for FAUS, to avoid spurious matches')
                }

                slip_model_type = 'FAUS'
                run_name2 = gsub('_uniform_slip', '', run_name2, fixed=TRUE)# Remove slip type from run_name2

            }else{
                stop('Could not match the slip type')

            }

            # Get the rigidity type, and remove the term from run_name2
            mu_type = ifelse(grepl('_varymu_', run_name2), 'varyMu', 'constant')
            run_name2 = gsub('_varymu', '', run_name2, fixed=TRUE)

            split_name = strsplit(run_name2, '_')[[1]]

            # Identify the PTHA18 scenario ID
            scenario_ID = split_name[2]
            historic_event = split_name[1]

            # Identify the source-name as used in PTHA18
            if(historic_event %in% c('chile1960', 'chile2010', 'southamerica2014', 'southamerica2015')){
                source_zone = 'southamerica'
            }else if(historic_event == 'puysegur2009'){
                source_zone = 'puysegur2'
            }else if(historic_event %in% c('sumatra2004', 'sumatra2005', 'sumatra2007', 'java2006')){
                source_zone = 'sunda2'
            }else if(historic_event == 'solomon2007'){
                source_zone = 'solomon2'
            }else if(historic_event == 'tohoku2011'){
                source_zone = 'kurilsjapan'
            }else if(historic_event == 'newhebrides2021'){
                source_zone = 'newhebrides2'
            }else if(historic_event == 'kermadec2021'){
                source_zone = 'kermadectonga2'
            }else if(historic_event == 'sandwich2021'){
                source_zone = 'sandwich'
            }else{
                stop(paste0('no source-zone match', run_name2))
            }

            # Value for the scenario count [which only matters for random scenarios]
            ind = grep('count', split_name)[1] + 1
            scenario_count = as.numeric(gsub('-risetime', '', split_name[ind], 
                                             fixed=TRUE))
        }else if(grepl('_nonrandom_like_historic_', run_name2)){
            #
            # Treat non-random PTHA18 scenarios.
            # Not used in the analysis of random ptha18 scenarios
            #
            # These have a somewhat different naming convention, vs the random
            # scenarios above.

            run_type = 'nonrandom_like_historic'
            run_name2 = gsub('ptha18_nonrandom_like_historic_', '', run_name2, fixed=TRUE) # Edit this as we extract the info

            # Identify the slip type
            if(grepl('_stochastic_', run_name2, fixed=TRUE)){

                slip_model_type = 'HS'
                run_name2 = gsub('_stochastic', '', run_name2, fixed=TRUE)

            }else if(grepl('_variable_uniform_', run_name2, fixed=TRUE)){

                slip_model_type = 'VAUS'
                run_name2 = gsub('_variable_uniform', '', run_name2, fixed=TRUE)

            }else{
                # Currently we did not run uniform-slip models for "non-random"
                # PTHA scenarios.
                stop('Could not match the slip type')
            }

            # The rigidity type
            mu_type = ifelse(grepl('_varymu_', run_name2), 'varyMu', 'constant')
            run_name2 = gsub('_varymu', '', run_name2, fixed=TRUE)

            split_name = strsplit(run_name2, '_')[[1]]

            # Identify the PTHA18 scenario ID
            scenario_ID = split_name[2]
            source_zone = split_name[1]
            historic_event = paste0(split_name[8], split_name[9])

            # Nominal value for the count [which only matters for random scenarios]
            scenario_count = 1

        }else{
            #
            # Potentially a source inversion. We need to do hacky things to
            # extract the fields of interest. In future the naming conventions
            # should be established ahead of time!
            #
            batch_number = NA
            run_type = 'source_inversion'
            slip_model_type = 'source_inversion'
            mu_type = 'source_inversion'
            scenario_ID = NA

            split_name = strsplit(run_name2, '_')[[1]]

            # The source inversions always start with the source name, but conventions
            # are different than from PTHA18.
            # Force results similar to what we get from PTHA random scenarios
            event_translation = list(
                'Chile2014_AnEtAl2014-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'southamerica2014',
                'Chile2015_WilliamsonEtAl2017-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'southamerica2015',
                #'Chile2015_RomanoEtAl2016-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'southamerica2015',
                #'Tohoku2011_RomanoEtAl2015-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'tohoku2011', 
                #'Tohoku2011_SatakeEtAl2013-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'tohoku2011',
                'Tohoku2011_YamazakiEtAl2018Fixed-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'tohoku2011',
                'Chile2010_LoritoEtAl2011-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'chile2010',
                #'Chile2010_FujiSatake2013-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'chile2010',
                'Sandwich2021_Rogeretal2024-risetime_0-full-linear_with_manning-0.035-highres_perth' = 'sandwich2021',
                #'Sumatra2005_Fujiietal2020-risetime_0-full-linear_with_manning-0.035-highres_SWWA' = 'sumatra2005',
                'Solomon2007_Weietal2015-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'solomon2007', 
                'Kermadec2021_Romano_source-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'kermadec2021',
                'Puysegur2009_Bevanetal2010-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'puysegur2009',
                'NewHebrides2021_GusmanEtAl2022-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'newhebrides2021', 
                #'Sumatra2004_FujiiandSatake2007-risetime_0-full-linear_with_manning-0.035-highres_australiaSWWA' = 'sumatra2004',
                #"Sumatra2004_LoritoEtAl2010-risetime_0-full-linear_with_manning-0.035-highres_australiaSWWA" = 'sumatra2004',     
                #"Sumatra2004_PiatanesiLorito2007-risetime_0-full-linear_with_manning-0.035-highres_australiaSWWA" = 'sumatra2004', 
                "Sumatra2004_FujiiandSatake2007_time_varying-risetime_0.0-full-linear_with_manning-0.035-highres_australiaWA" = 'sumatra2004',
                "Sumatra2005_Fujiietal2020-risetime_0-full-linear_with_manning-0.035-highres_WA" = 'sumatra2005',
                'Chile1960_HoEtAl2019-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'chile1960',
                #'Chile1960_FujiSatake2013-risetime_0-full-linear_with_manning-0.035-highres_NSW' = 'chile1960',
                'Java2006_FujiiandSatake2006_time_varying-risetime_0.0-full-linear_with_manning-0.035-highres_NWWA' = 'java2006',
                'Sumatra2007_Fujiietal2008-risetime_0-full-linear_with_manning-0.035-highres_NWWA' = 'sumatra2007'
                )
            ki = match(run_name2, tolower(names(event_translation)))
            if(length(ki) != 1) stop(paste0('Could not match the run type for ', run_name2, 
                '-- if this is a source inversion, maybe you need to update event_translation above?'))

            historic_event = event_translation[[ki]]

            # Identify the source-name as used in PTHA18
            if(historic_event %in% c('chile1960', 'chile2010', 'southamerica2014', 'southamerica2015')){
                source_zone = 'southamerica'
            }else if(historic_event == 'puysegur2009'){
                source_zone = 'puysegur2'
            }else if(historic_event %in% c('sumatra2004', 'sumatra2005', 'java2006', 'sumatra2007')){
                source_zone = 'sunda2'
            }else if(historic_event == 'solomon2007'){
                source_zone = 'solomon2'
            }else if(historic_event == 'tohoku2011'){
                source_zone = 'kurilsjapan'
            }else if(historic_event == 'newhebrides2021'){
                source_zone = 'newhebrides2'
            }else if(historic_event == 'kermadec2021'){
                source_zone = 'kermadectonga2'
            }else if(historic_event == 'sandwich2021'){
                source_zone = 'sandwich'
            }else{
                stop(paste0('no source-zone match', run_name2))
            }

            scenario_count = 1    
        }

        # Output a handful of statistics we want
        return(c(run_type, slip_model_type, scenario_ID, source_zone, historic_event, mu_type, scenario_count, batch_number))
    }
    run_type_and_slip_type = sapply(
        basename(dirname(dirname(event_stats$run_name))), 
        model_type_from_run_name)

    # Pack the variables (extracted above) into the event statistics
    event_stats = cbind(event_stats, 
        data.frame(
            run_type=run_type_and_slip_type[1,], 
            slip_type = run_type_and_slip_type[2,],
            scenario_ID = as.numeric(run_type_and_slip_type[3,]),
            source_zone = run_type_and_slip_type[4,],
            event_name = run_type_and_slip_type[5,],
            site_and_event = paste0(run_type_and_slip_type[5,], '_', event_stats$sites),
            rigidity = run_type_and_slip_type[6,],
            scenario_count = as.numeric(run_type_and_slip_type[7,]),
            batch_number = as.numeric(run_type_and_slip_type[8,]))
        )


    ## Get statistics about the tide gauge residual, before the tsunami arrives. 
    ## This is useful to understand the "non-tsunami variations" in the tide gauge.
    #get_obs_stats_before_tsunami_arrives<-function(){

    #    # No model "arrives" before this time at the site
    #    minimum_model_arrival_times = aggregate(event_stats$model_time_to_arrive, 
    #        by=list(site_and_event=event_stats$site_and_event), min)
    #    min_arrival_times_rowind = match(event_stats$site_and_event, minimum_model_arrival_times$site_and_event)
    #    stopifnot(all(event_stats$site_and_event == minimum_model_arrival_times$site_and_event[min_arrival_times_rowind]))

    #    # Compute some statistics once per site and event pair
    #    obsrms_store = rep(NA, nrow(event_stats))
    #    obsmax_store = rep(NA, nrow(event_stats))
    #    obsrng_store = rep(NA, nrow(event_stats))
    #    for(i in 1:length(minimum_model_arrival_times[,1])){
    #        site_and_event = minimum_model_arrival_times[i,1]
    #        event_stats_rows = which(min_arrival_times_rowind == i) # Rows with this site_and_event
    #        site_and_event_table_rowind = min(event_stats_rows) # First row index
    #        # Index of the time-series data object that corresponds to this site and event
    #        target_RDS_ind = min(which(names(target_RDS_data) == event_stats$run_name[site_and_event_table_rowind]))

    #        tgs = target_RDS_data[[target_RDS_ind]] # Contains multiple sites
    #        tgs_name_matches = sapply(names(tgs), function(x) endsWith(site_and_event, x), USE.NAMES=FALSE) # Find the one we want
    #        gauge_ind = which(tgs_name_matches) # index of gauge we want
    #        if(length(gauge_ind) != 1){ # Sanity check
    #            print(names(tgs))
    #            print(site_and_event)
    #            print(gauge_ind)
    #            print(tgs_name_matches)
    #            stop('error: length(gauge_ind) is not equal to 1')
    #        }
    #        obs_data = tgs[[gauge_ind]]$event_data$obs # Data for the event/gauge of interest
    #        # Time window over which we compute the RMS -- 2 hrs beforehand, to 30min before first tsunami model arrives.
    #        win_start = tgs[[gauge_ind]]$model_start_julian - (2/24)
    #        win_end   = tgs[[gauge_ind]]$model_start_julian + minimum_model_arrival_times[i,2]/(3600*24) - (0.5/24) 
    #        obsinds = which(obs_data$juliant > win_start & obs_data$juliant < win_end)
    #        # Root mean square, removing NA values
    #        obsrms = sqrt( mean(obs_data$resid[obsinds]**2, na.rm=TRUE) )
    #        # Max
    #        obsmax = max(obs_data$resid[obsinds], na.rm=TRUE)
    #        # Stage range
    #        obsrng = diff(range(obs_data$resid[obsinds], na.rm=TRUE))
    #        # Put this value in all rows
    #        obsrms_store[event_stats_rows] = obsrms
    #        obsmax_store[event_stats_rows] = obsmax
    #        obsrng_store[event_stats_rows] = obsrng
    #    }
    #    return(data.frame(data_rms_before_arrival = obsrms_store, 
    #                      data_max_before_arrival = obsmax_store, 
    #                      data_stgrng_before_arrival=obsrng_store))
    #}
    # # Store these above statistics
    #event_stats = cbind(event_stats, get_obs_stats_before_tsunami_arrives())

    good_nearshore_data = is_good_nearshore_data(event_stats)


    # Append more useful statistics to the event_stats. 
    # For ease of coding later, it is good to have a 'model-statistic'
    # corresponding to each new 'data-statistic', even if the model-statistic
    # is just a repeat of some other column.
    event_stats = cbind( 
        data.frame(
            event_name_int = as.numeric(as.factor(event_stats$event_name)),
            event_and_source_int = as.numeric(as.factor(event_stats$site_and_event)), 
            good_nearshore = good_nearshore_data,
            stringsAsFactors=FALSE),
        event_stats)

    #
    # Get the scenario energy (stored in the zipped multidomain log file)
    # and append to the event stats. Also run some logical checks [e.g.
    # that energy is decreasing].
    #

    md_log_files_zipped = paste0(dirname(target_RDS), 
        '/multidomain_log_image_00000000000000000001.log.zip')
    # Requires a bit of work to get the log from the zipped file cleanly
    read_zipped_log<-function(x){
        # Do NOT DO THIS IN PARALLEL!
        temp_dir = tempdir()
        unzip(x, junkpaths=TRUE, exdir=temp_dir)
        output = readLines(paste0(temp_dir, 
            '/multidomain_log_image_00000000000000000001.log'))
        return(output)
        }
    all_md_logs = lapply(md_log_files_zipped, read_zipped_log)
    # Now read the energy time-series -- note this is "energy/rho"
    all_md_energy = lapply(all_md_logs, function(x){
        inds = grep('Global energy-total', x)
        energy_on_rho = as.numeric(x[inds+1])
        return(energy_on_rho)
        })

    # Check the energy stats -- do they decrease over time as expected?
    energy_start = unlist(lapply(all_md_energy, function(x) x[1]))
    energy_end = unlist(lapply(all_md_energy, function(x) x[length(x)]))
    energy_max = unlist(lapply(all_md_energy, function(x) max(x)))
    energy_min = unlist(lapply(all_md_energy, function(x) min(x)))
    #
    # Now match these energy statistics to the event_stats, and put them in the table.
    inds = match(dirname(event_stats$run_name), dirname(md_log_files_zipped))
    RHO = 1024 # SWALS gives energy/rho, so convert to raw-energy below
    event_stats$energy_start = energy_start[inds] * RHO
    event_stats$energy_end = energy_end[inds] * RHO
    event_stats$energy_max = energy_max[inds] * RHO
    event_stats$energy_min = energy_min[inds] * RHO

    # Check for non-negligible increases in energy. Tiny increases in the first
    # few timesteps are normal with the leapfrog solver (Davies et al 2020).
    # Only look at random scenarios (source inversions with a finite rise time will have
    # zero initial energy which sends these statistics infinite).
    i0 = which(event_stats$run_type == 'random_like_historic')
    print('SUMMARY (random_like_historic): (energy_max - energy_start) / energy_start')
    print(summary((event_stats$energy_max[i0] - event_stats$energy_start[i0])/event_stats$energy_start[i0])) 
    print('SUMMARY (random_like_historic): (energy_min - energy_end) / energy_end')
    print(summary((event_stats$energy_min[i0] - event_stats$energy_end[i0])/event_stats$energy_end[i0]))
    # Check that energy decreases from start to end. 
    energy_final_fraction = event_stats$energy_end/event_stats$energy_start
    print('SUMMARY (random_like_historic): energy_final_fraction')
    print(summary(energy_final_fraction[i0]))


    ##
    ## Save to RDS, and make QC plot of data [the latter only needed once]
    ##

    if(DOWNSAMPLE_MODEL_TO_DATA){
        # Here the only thing that changes, compared to not downsampling the model,
        # is the event_stats. So we avoid re-saving everything.
        saveRDS(event_stats, 'event_stats_DOWNSAMPLED.RDS')
    }else{
        saveRDS(event_stats, 'event_stats.RDS')
        saveRDS(target_RDS_data, 'target_RDS_data.RDS')

        ##
        ## Plot the nearshore data to check it is OK
        ##
        pdf('QC_good_nearshore_obs_in_highres_domains.pdf', width=8, height=5)
        k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
        unique_station_records = unique(event_stats$site_and_event[k])
        for(i in 1:length(unique_station_records)){

            # Find the data
            k = min(which(event_stats$site_and_event == unique_station_records[i]))
            site_name = event_stats$sites[k]
            tgi = match(event_stats$run_name[k], names(target_RDS_data)) 
            obs = target_RDS_data[[tgi]][[site_name]]$event_data$obs
            model_vs_gauge_start_time = event_stats$model_vs_gauge_start_time[k]

            # Plot it
            plot(obs$juliant, obs$resid, t='l', main=event_stats$site_and_event[k])
            abline(h=0, col='red')
            grid()

            # Also plot "first 15 hours after the comparison starts". 
            k = which( (obs$juliant - model_vs_gauge_start_time < 15/24) &
                       (obs$juliant >= (model_vs_gauge_start_time - 2/24) ) )
            plot(obs$juliant[k], obs$resid[k], t='o', cex=0.3, 
                main='From (-2) to (+15) hours around the model-vs-gauge-comparison start time')
            abline(h=0, col='red')
            grid()
        }
        dev.off()

        ##
        ## Repeat the above plot, but add in some summary stats
        ##
        pdf('QC_good_nearshore_obs_in_highres_domains_with_statistics_that_vary_by_model_time.pdf', width=8, height=5)
        k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
        unique_station_records = unique(event_stats$site_and_event[k])
        for(i in 1:length(unique_station_records)){

            # Find the data
            m_r= which(event_stats$site_and_event == unique_station_records[i]) # matching rows
            k = min(m_r)
            site_name = event_stats$sites[k]
            tgi = match(event_stats$run_name[k], names(target_RDS_data)) 
            obs = target_RDS_data[[tgi]][[site_name]]$event_data$obs
            model_vs_gauge_start_time = event_stats$model_vs_gauge_start_time[k]

            # Plot it
            plot(obs$juliant, obs$resid, t='l', main=event_stats$site_and_event[k])
            abline(h=0, col='red')
            grid()

            # Add various data ranges
            abline(h=unique(event_stats$data_max[m_r]), col='black')
            abline(h=unique(event_stats$data_min[m_r]), col='black')
            abline(h=unique(event_stats$data_max_24hrs[m_r]), col='red', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_min_24hrs[m_r]), col='red', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_max_12hrs[m_r]), col='green', lty='dashed')
            abline(h=unique(event_stats$data_min_12hrs[m_r]), col='green', lty='dashed')
            abline(h=unique(event_stats$data_max_8hrs[m_r]), col='purple', lty='dotted')
            abline(h=unique(event_stats$data_min_8hrs[m_r]), col='purple', lty='dotted')
            abline(v=event_stats$model_arrival_time[m_r], col='brown')

            points(obs$juliant, obs$resid, t='l') # Avoid obscuring this

            legend('bottom', 
                c('Max', '24hrs', '12hrs', '8hrs'), 
                lty=c('solid', 'dotted', 'dashed', 'dotted'),
                col=c('black', 'red', 'green', 'purple'), 
                lwd=c(1,2,1,1), horiz=TRUE, bty='n') 

            # Also plot the "fist 15 hours after arrival". Here the arrival time is
            # defined FROM THE MODEL in the usual way.
            k = which( (obs$juliant - model_vs_gauge_start_time < 15/24) &
                       (obs$juliant >= (model_vs_gauge_start_time - 2/24) ) )
            plot(obs$juliant[k], obs$resid[k], t='o', cex=0.3, 
                main='From (-2) to (+15) hours around a model-vs-gauge-comparison-start time')
            abline(h=0, col='red')

            # Add various data ranges
            abline(h=unique(event_stats$data_max[m_r]), col='black')
            abline(h=unique(event_stats$data_min[m_r]), col='black')
            abline(h=unique(event_stats$data_max_24hrs[m_r]), col='red', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_min_24hrs[m_r]), col='red', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_max_12hrs[m_r]), col='green', lty='dashed')
            abline(h=unique(event_stats$data_min_12hrs[m_r]), col='green', lty='dashed')
            abline(h=unique(event_stats$data_max_8hrs[m_r]), col='purple', lty='dotted')
            abline(h=unique(event_stats$data_min_8hrs[m_r]), col='purple', lty='dotted')
            abline(v=event_stats$model_arrival_time[m_r], col='brown')

            points(obs$juliant[k], obs$resid[k], t='o', cex=0.3) # Avoid obscuring this

            grid()
        }
        dev.off()

        pdf('QC_good_nearshore_obs_in_highres_domains_with_statistics_that_depend_on_earthquake_time.pdf', width=8, height=5)
        k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain)
        unique_station_records = unique(event_stats$site_and_event[k])
        for(i in 1:length(unique_station_records)){

            # Find the data
            m_r= which(event_stats$site_and_event == unique_station_records[i]) # matching rows
            k = min(m_r)
            site_name = event_stats$sites[k]
            tgi = match(event_stats$run_name[k], names(target_RDS_data)) 
            obs = target_RDS_data[[tgi]][[site_name]]$event_data$obs
            model_vs_gauge_start_time = event_stats$model_vs_gauge_start_time[k]

            # Plot it
            plot(obs$juliant, obs$resid, t='l', main=event_stats$site_and_event[k])
            abline(h=0, col='red')
            grid()

            # Add various data ranges
            abline(h=unique(event_stats$data_max[m_r]), col='black')
            abline(h=unique(event_stats$data_min[m_r]), col='black')
            abline(h=unique(event_stats$data_max_24[m_r]), col='blue', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_min_24[m_r]), col='blue', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_max_36[m_r]), col='purple', lty='dashed')
            abline(h=unique(event_stats$data_min_36[m_r]), col='purple', lty='dashed')
            abline(h=unique(event_stats$data_max_lastday[m_r]), col='orange', lty='dotted')
            abline(h=unique(event_stats$data_min_lastday[m_r]), col='orange', lty='dotted')
            abline(v=event_stats$model_arrival_time[m_r], col='brown')

            points(obs$juliant, obs$resid, t='l') # Avoid obscuring this

            legend('bottom', 
                c('Max', '24', '36', 'lastday'), 
                lty=c('solid', 'dotted', 'dashed', 'dotted'),
                col=c('black', 'blue', 'purple', 'orange'), 
                lwd=c(1,2,1,1), horiz=TRUE, bty='n') 

            # Also plot the "first 15 hours after start of comparison with gauges". 
            k = which( (obs$juliant - model_vs_gauge_start_time < 15/24) &
                       (obs$juliant >= (model_vs_gauge_start_time - 2/24) ) )
            plot(obs$juliant[k], obs$resid[k], t='o', cex=0.3, 
                main='From (-2) to (+15) hours around a model-vs-gauge-comparison start time')
            abline(h=0, col='red')

            # Add various data ranges
            abline(h=unique(event_stats$data_max[m_r]), col='black')
            abline(h=unique(event_stats$data_min[m_r]), col='black')
            abline(h=unique(event_stats$data_max_24[m_r]), col='blue', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_min_24[m_r]), col='blue', lty='dotted', lwd=2)
            abline(h=unique(event_stats$data_max_36[m_r]), col='purple', lty='dashed')
            abline(h=unique(event_stats$data_min_36[m_r]), col='purple', lty='dashed')
            abline(h=unique(event_stats$data_max_lastday[m_r]), col='orange', lty='dotted')
            abline(h=unique(event_stats$data_min_lastday[m_r]), col='orange', lty='dotted')
            abline(v=event_stats$model_arrival_time[m_r], col='brown')

            points(obs$juliant[k], obs$resid[k], t='o', cex=0.3) # Avoid obscuring this

            grid()
        }
        dev.off()
    }
}

# Parse the data without any smoothing of the model
parse_gauge_outputs(DOWNSAMPLE_MODEL_TO_DATA = FALSE)
                      
# Parse the data with smoothing of the model
parse_gauge_outputs(DOWNSAMPLE_MODEL_TO_DATA = TRUE)
