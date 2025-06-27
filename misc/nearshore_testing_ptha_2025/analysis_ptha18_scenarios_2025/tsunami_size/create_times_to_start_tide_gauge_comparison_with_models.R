#
# For consistency of the statistics, it's nice if for each tsunami event and tide-gauge we have a consistent 'start time' to begin comparison with models.
# Construct this here based on the minimum of the algorithmically defined times.
#
.times_to_start_comparison_RDS_file =  'times_to_start_comparison_with_tide_gauge.RDS'

.make_times_to_start_comparison<-function(){
    # Define the time at which to start comparing the model and data (a single
    # value for each gauge/event combination), by using the minimum of the
    # "modelled" arrival times from stochastic scenarios.
    #
    # The calculation here assumes the analysis has already been done! This
    # reflects our iterative development, but is just for convenience.
    #
    # In event_stats.RDS, we only use the 'model_start' julian time (which was
    # calculated from the earthquake start time) and the "model_arrival_time"
    # (which was computed from the modelled gauge time-series and the chosen
    # WAVE_ARRIVAL_TIME_RELATIVE_THRESHOLD), as well as the gauge and event names. 
    #
    # In principle we could compute this separately (but it would be slower
    # due to all the data reading)

    event_stats = readRDS('event_stats.RDS')
    stopifnot(length(unique(event_stats$model_start)) == 14) # 14 events, 14 unique values (even though they are reals!)

    # Restrict the calculation to the stochastic scenarios
    #event_stats = event_stats[event_stats$run_type == "random_like_historic",]

    # Find the minimum model arrival time for each site and event start time
    # (the latter is a proxy for the event).
    model_start_times = aggregate(event_stats$model_arrival_time, 
        by=list(site=event_stats$sites, event_start_juliant = event_stats$model_start), 
        FUN <- function(x) min(x, na.rm=TRUE))

    model_start_times$event_start_juliant = as.numeric(model_start_times$event_start_juliant) 

    k = which(names(model_start_times) == "x")
    names(model_start_times)[k] = "gauge_comparison_start_juliant"
   
    model_start_times$gauge_comparison_start_juliant = as.numeric(model_start_times$gauge_comparison_start_juliant)

    if(!file.exists(.times_to_start_comparison_RDS_file)){
        saveRDS(model_start_times, .times_to_start_comparison_RDS_file)
    }else{
        print(paste0('Will not overwrite existing RDS file ', .times_to_start_comparison_RDS_file, 
            ' \n You must delete it manually if you want to replace it'))
    }

    return(NULL)
}

#
# Create a function to return the julian time at which we start comparing models and gauges for each event.
# We make this a function of the gauge name and the model start time, because they are conveniently available
# at the point in the code at which we need this.
#

if(file.exists(.times_to_start_comparison_RDS_file)){

    .times_to_start_comparison = readRDS(.times_to_start_comparison_RDS_file)

    get_julian_time_at_which_to_start_comparison_with_tide_gauge<-function(gauge_name, model_start_time){

        k1 = which(.times_to_start_comparison$event_start_juliant == model_start_time)
        k2 = which(.times_to_start_comparison$site == gauge_name)

        k3 = intersect(k1, k2)
        if(length(k3) != 1){
            stop(paste0('Could not find unique row matching ', gauge_name, ' ', model_start_time))
        }
        return(.times_to_start_comparison$gauge_comparison_start_juliant[k3])
    }

}else{

    print(paste0('Warning: To get consistent times for comparing models and gauges, you need to create the file \n', 
        .times_to_start_comparison_RDS_file, ', \n for instance using .make_times_to_start_comparison()'))
}
