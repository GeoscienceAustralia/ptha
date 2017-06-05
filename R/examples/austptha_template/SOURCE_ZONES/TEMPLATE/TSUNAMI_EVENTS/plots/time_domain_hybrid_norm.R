

gauge_similarity_time_domain<-function(data1_t, data1_s, data2_t, data2_s, interp_dt = NULL,
    allowed_lag_minutes=15, time_range=range(data1_t), detailed=FALSE){
    
    if(is.null(interp_dt)){
        interp_dt = min(c(diff(data1_t), diff(data2_t)))
    }

    # Interpolate the data
    interp1_t = approx(data1_t, data1_s, xout = seq(time_range[1], 
        time_range[2], by = interp_dt), rule=2)

    # Find min/max of data, and associated times
    #peak_stage_info = gauge_range_filtered(data1_t, data1_s, interp_dt = interp_dt, 
    #    filter_freq = peak_stage_high_pass_cutoff, detailed=TRUE)
    #mint = min(c(peak_stage_info$min_stage_time, peak_stage_info$max_stage_time))
    #maxt = max(c(peak_stage_info$min_stage_time, peak_stage_info$max_stage_time))
    newt = time_range
    mint = min(newt)
    maxt = max(newt)
    
   
    # Extract the model and data between the min/max times 
    data1_interp = list()
    data1_interp$x = interp1_t$x[which(interp1_t$x >= mint & interp1_t$x <= maxt)]
    data1_interp$y = interp1_t$y[which(interp1_t$x >= mint & interp1_t$x <= maxt)]
   
    # Compute statistic from Lorito et al., 2008, eqn 2 
    f<-function(lag){ 
        data2_interp = approx(data2_t - lag, data2_s, xout = data1_interp$x, rule=2)
        Em = 1 - 2*sum(data1_interp$y * data2_interp$y)/( sum(data1_interp$y^2) + sum(data2_interp$y^2))
        return(Em)
    }

    best_lag = optimize(f, interval=c(-60, 60)*allowed_lag_minutes)

    if(!detailed){
        return(best_lag$objective)
    }else{
        return(best_lag)
    }

}
