
#'
#' Function to measure the 'similarity' between 2 time-series [stage-gauges]
#'
#' @param data1_t Times (seconds) of data 
#' @param data1_s Stage values of data at times in data1_t. 
#' @param model2_t Time (seconds) of model. This should overlap the data time
#' period, but they don't need to have exactly the same times, as interpolation
#' will be used to align the model to the data.
#' @param model2_s Stage values of model at times in model2_t
#' @param interp_dt The data is first interpolated to this time interval (seconds) and then
#' the model is interpolated to the same times.
#' @param allowed_lag_minutes To compute the GOF, we test time shifts in the range of allowed_lag_minites. 
#' @param statistic_variant either 'model_in_denominator' {e.g. Lorito et al 2008, Davies 2019, Romano et al. 2020) or experimental alternative 'no_model_in_denominator'
#' The model times are shifted by (-1)*lag, where lag (in seconds) ranges over allowed_lag_minutes in 10s intervals. 
#' The GOF statistic uses the best of the test lag values. 
#' @param time_range The time period over which to compute the statistic. Both the model and data should
#' be operating
#' @param detailed If TRUE, return both the GOF and the best lag value. If FALSE, just return the GOF.
gauge_similarity_time_domain<-function(
    data1_t, 
    data1_s, 
    model2_t, 
    model2_s, 
    interp_dt = NULL,
    allowed_lag_minutes=c(-15, 0), 
    time_range=range(data1_t), 
    statistic_variant = 'model_in_denominator',
    detailed=FALSE){
    
    if(is.null(interp_dt)){
        interp_dt = min(c(diff(data1_t), diff(model2_t)))
    }

    # Interpolate the data
    interp1_t = approx(data1_t, data1_s, xout = seq(time_range[1], 
        time_range[2], by = interp_dt), rule=2)
   
    # Extract data between the min/max times 
    newt = time_range
    mint = min(newt)
    maxt = max(newt)
    data1_interp = list()
    data1_interp$x = interp1_t$x[which(interp1_t$x >= mint & interp1_t$x <= maxt)]
    data1_interp$y = interp1_t$y[which(interp1_t$x >= mint & interp1_t$x <= maxt)]
   
    ## Statistic similar to that from Lorito et al., 2008, eqn 2 
    ##     = 1 - 2*(a.b)/(a.a + b.b)
    ## Note this is mathematically identical to (a-b).(a-b) / (a.a + b.b), i.e. normalised least squares
    ## Here we optionally add weights to emphasise the points with significant waves.
    f<-function(lag, weighted=TRUE){ 
        # Interpolate the model at the data times
        model2_interp = approx(model2_t - lag, model2_s, xout = data1_interp$x, rule=2)

        if(weighted){
            # Weight large-abs-value points more. Minimum point weight is not
            # less than 1/3 maximum point weight.
            w = pmax(abs(data1_interp$y),  max(abs(data1_interp$y))/3)
        }else{
            # Weight all points equally
            w = 1
        }

        if(statistic_variant == 'model_in_denominator'){
            Em = 1 - 2*sum(w*w*data1_interp$y * model2_interp$y)/
                ( sum(w*w*data1_interp$y^2) + sum(w*w*model2_interp$y^2) )
        }else if(statistic_variant == 'no_model_in_denominator'){
            Em = sum((w*data1_interp$y - w*model2_interp$y)^2)/(sum(w*w*data1_interp$y^2))
        }else{
            stop('unknown statistic_variant')
        }
        return(Em)
    }

    # Test lags between min/max, with spacing of about 10s
    # Implement our own brute force minimization, since the R function
    # "optimize" appears like it might not hit the minimum reliably for this
    # problem.
    sec_in_min = 60
    lag_vals = seq( sec_in_min * allowed_lag_minutes[1], sec_in_min*allowed_lag_minutes[2], 
        length=max(1, ceiling(diff(allowed_lag_minutes)*sec_in_min/10)) )
    lag_f = lag_vals*0
    for(i in 1:length(lag_vals)){
        lag_f[i] = f(lag_vals[i])
    }

    # Make the output look like output from a call to optimize
    best_lag = list(objective=min(lag_f), minimum = lag_vals[which.min(lag_f)])

    if(!detailed){
        return(best_lag$objective)
    }else{
        return(best_lag)
    }

}
