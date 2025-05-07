#' Workhorse function for downsampling model results at nearshore tide-gauges. 
#' This will be called if (DOWNSAMPLE_MODEL_TO_DATA==TRUE), at nearshore gauges only
#'
#' Why do this? Nearshore tide-gauges generally record a time-smooth of the data 
#' (often 1min), and some of them only report the value at 15 min intervals. Thus
#' the tsunami is 'under-sampled', and the model will tend to have greater maxima etc. 
#' To compare with a model with such data, one approach is to under-sample the
#' model in a similar way.
#'
#' @param model A list with entries model$juliant and model$resid giving the julian
#' time and tsunami FROM THE MODEL
#' @param gauge_name the name of the gauge. From this we can infer information on the gauge
#' sampling.
#' @return a list model with model$juliant and model$resid, where the model results have been
#' smoothed and downsampled for greater consistency with the gauge.
#'
downsample_and_smooth_model_at_nearshore_tide_gauge<-function(model, gauge_name){

    # - At most good nearshore gauges we know that the tide-gauge data is a running mean over 1min. 
    # - At the 15min gauges, that data is sub-sampled to 1/15 of the samples.
    # - At the BOM 6min gauges, there is some additional smoothing going on. For now I won't worry
    #   about this because anyway all that 6min data is redundant [I have 1min data]

    if(grepl('DART', gauge_name)) stop('downsample_and_smooth_...() is not designed for DART buoy data')

    # The nearshore tide-gauges have different sampling frequencies, but in all
    # cases, we start with 1 min smooth. The metadata for all nearshore gauges
    # I've seen suggests this would be the minimal smoothing.
    ONE_HOUR = 3600.0
    HOURS_IN_DAY = 24
    model_secs = as.numeric((model$juliant - model$juliant[1])*ONE_HOUR*HOURS_IN_DAY)
    model_approxfun = approxfun(model_secs, model$resid, rule=2)
    model_1min = seq(0, max(model_secs), by=60)
    model_smooth = model_1min * 0
    count = 0
    # Smooth over 1 minute [from -30s to +30s]
    for(i in seq(-30, 30, by=10)){ 
        model_smooth = model_smooth + model_approxfun(model_1min + i)
        count = count + 1
    }
    model_smooth = model_smooth/count

    #
    # The gauge data filenames include info on the sampling
    # Use that to determine the output times.
    #
    if(grepl('_1min_', gauge_name) | endsWith(gauge_name, '_1min') | gauge_name == "ROSSLYN_BAY_BOM" | gauge_name == "GoldCoastSandBypass"){
        # 1min data 
        output_times = seq(0, max(model_secs), by=60)
        
    }else if(grepl('_6min_', gauge_name) | endsWith(gauge_name, '_6min')){
        # The 6minute gauges are from BOM. There is DEFINITELY
        # extra smoothing applied, see their metadata. For this 
        # study the 6min data is redundant [as I always have
        # corresponding 1min data] so I don't implement smoothing, but
        # in reality there is some.
        output_times = seq(0, max(model_secs), by=60*6)

    }else if(grepl('_5min_', gauge_name) | endsWith(gauge_name, '_5min')){
        # The 5minute gauges are from WA. I don't know what extra
        # smoothing is applied, if any.
        output_times = seq(0, max(model_secs), by=60*5)

    }else if(grepl('_10min_', gauge_name) | endsWith(gauge_name, '_10min') | endsWith(gauge_name, '_MSQ') | endsWith(gauge_name, '_DESI')){
        # This data isn't used herein
        output_times = seq(0, max(model_secs), by=60*15)

    }else if(grepl('_15min_', gauge_name) | endsWith(gauge_name, '_15min')){
        # The 15 min data from Manly Hydraulics Lab.
        # They explicitly do not smooth -- it is just decimiated 1
        # minute data -- which I have confirmed in several cases.
        output_times = seq(0, max(model_secs), by=60*15)

    }else if(endsWith(gauge_name, '1960')){
        # Special case -- chile 1960 manually digitized data
        # In this case let's assume sampling at 2min
        output_times = seq(0, max(model_secs), by=60*2)

    }else{
        stop(paste0('Could not match gauge_name= ', gauge_name))
    }

    # Here are the values to pack into the 'model' list
    model_values_out = approx(model_1min, model_smooth, xout=output_times)$y
    model_times_out = as.difftime(output_times/(ONE_HOUR*HOURS_IN_DAY), units='days') + model$juliant[1]

    # HERE WE OVERWRITE 'model' with the smoothed and down-sampled data.
    model = data.frame(
        juliant = model_times_out,
        resid = model_values_out)
    return(model)
}


