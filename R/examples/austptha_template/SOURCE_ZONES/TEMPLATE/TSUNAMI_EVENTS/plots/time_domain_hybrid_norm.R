
#
# Function to measure the 'similarity' between 2 time-series [stage-gauges]
#
# gauge1: data1_t, data1_s ; vectors of time/stage values
# gauge2: data2_t, data2_s ; vectors of time/stage values [note times do not
#   need to correspond exactly to data1_t, since interpolation is used to create
#   a common time scale]
#
gauge_similarity_time_domain<-function(
    data1_t, 
    data1_s, 
    data2_t, 
    data2_s, 
    interp_dt = NULL,
    allowed_lag_minutes=c(-15, 0), 
    time_range=range(data1_t), 
    detailed=FALSE){
    
    if(is.null(interp_dt)){
        interp_dt = min(c(diff(data1_t), diff(data2_t)))
    }

    # Interpolate the data
    interp1_t = approx(data1_t, data1_s, xout = seq(time_range[1], 
        time_range[2], by = interp_dt), rule=2)

    # Find min/max of data, and associated times
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

    # Test lags between min/max, with spacing of about 10s
    # Implement our own brute force minimization, since optimize appears like it might not hit the minimum
    sec_in_min = 60
    lag_vals = seq(sec_in_min * allowed_lag_minutes[1], sec_in_min*allowed_lag_minutes[2], 
        length=max(1, ceiling(diff(allowed_lag_minutes)*60/10)) )
    lag_f = lag_vals*0
    for(i in 1:length(lag_vals)){
        lag_f[i] = f(lag_vals[i])
    }

    #best_lag = optimize(f, interval=sec_in_min*allowed_lag_minutes)

    best_lag = list(objective=min(lag_f), minimum = lag_vals[which.min(lag_f)])

    if(!detailed){
        return(best_lag$objective)
    }else{
        return(best_lag)
    }

}
