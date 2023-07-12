#
# Functions to remove tides and make plots
#
source('global_variables.R')
source('spectral_highpass_filter.R')

#
# Get tidal predictions from tpxo, using the 'stormwavecluster' interface to tpxo.
#
predict_tide_location = './TPXO9_interface/predict_tide.R'
tpxo = new.env()
source(predict_tide_location, chdir=TRUE, local=tpxo)


#' Extract the tsunami signal from the provided stage timeseries
#'
#' @param tidal_stage Tidal stage data.frame with columns 'time' (strptime), 'juliant' (julian time), and 'stage'
#' @param gauge_location numeric vector with c(lon,lat) coordinate of the gauge
#' @return a list with the times, the original data, and the low-frequency and high-frequency components
#'
extract_tsunami_from_stage<-function(tidal_stage, gauge_location, low_frequency_limit=1/(3*3600)){


    if(!exists('tpxo')) stop('You must create the tpxo object before running this function')

    # Time increment used to slightly extended the requested interval (for interpolation)
    sixmin = as.difftime(6, units='mins') 

    # Tidal prediction @ gauge location. This is a global model.
    tidal_stage_tidal_pred = tpxo$get_tidal_prediction(
        'GAUGELOCATION', gauge_location, 
        start_time=(tidal_stage$time[1]-sixmin), 
        end_time=(tidal_stage$time[length(tidal_stage$time)]+sixmin),
        time_interval='6 min') # If this interval is too small it can be slow
    tidal_stage_tidal_pred$juliant = julian(tidal_stage_tidal_pred$time)
    tidal_stagefun = approxfun(tidal_stage_tidal_pred$juliant, tidal_stage_tidal_pred$tide)

    # Tidal residual
    tidal_stages = tidal_stagefun(tidal_stage$juliant)
    tidal_resid = tidal_stage$stage - tidal_stages
    mean_resid = mean(tidal_resid, na.rm=TRUE)
    # We want to keep the 'tsunami' and remove the long-period component.
    tidal_resid_nomean = tidal_resid - mean_resid

    # Now, remove the remaining long-period component (defined by low_frequency_limit)
    tidal_resid_nomean_filtered = spectral_highpass_filter(
        data_t = as.numeric(tidal_stage$juliant - tidal_stage$juliant[1])*60*60*24, # Must be in seconds, not days
        data_s = tidal_resid_nomean,
        interp_dt = 15,
        cutoff_frequency = low_frequency_limit)

    # Alternative de-tiding: Just apply the high-pass filter to the original data, without considering tidal predictions
    residual_filtered = spectral_highpass_filter(
        data_t = as.numeric(tidal_stage$juliant - tidal_stage$juliant[1])*60*60*24, # Must be in seconds, not days
        data_s = tidal_stage$stage,
        interp_dt = 15,
        cutoff_frequency = low_frequency_limit)

    # If the stage is NA, then ensure the residual is also NA (because spectral_highpass_filter fills these parts with interpolation).
    k = which(is.na(tidal_stage$stage))
    if(length(k) > 0){
        tidal_resid_nomean_filtered$highfreq[k] = NA        
        residual_filtered$highfreq[k] = NA        
    }

    # Return the original times, the original record, the predicted tide, the high_frequency
    # and low_frequency components
    output = list(time = tidal_stage$time, 
                  juliant = tidal_stage$juliant, 
                  stage = tidal_stage$stage,
                  tidal_prediction = tidal_stages,
                  tide = tidal_stages + tidal_resid_nomean_filtered$lowfreq + mean_resid, 
                  tsunami = tidal_resid_nomean_filtered$highfreq,
                  tsunami_highpass = residual_filtered$highfreq,
                  mean_resid = mean_resid)

    return(output)
}

# Quick plot of the stage results after de-tiding.
plot_detided_results<-function(tide_gauge_filtered){

    for(i in 1:length(tide_gauge_filtered)){
        par(mfrow=c(4,1))
        par(mar=c(3,2,1,1))
        tgf = tide_gauge_filtered[[i]]
        if(is(tgf, 'try-error')) next
        tgf_name = names(tide_gauge_filtered)[i]

        #
        plot(tgf$time, tgf$stage, t='l', main=paste0(tgf_name, 
            ': Raw stage [black] & tide [red, prediction + low-frequency residual]'))
        points(tgf$time, tgf$tide, t='l', col='red')
        grid(col='orange')
        #
        plot(tgf$time, tgf$tsunami, t='l', main=paste0(tgf_name, 
            ': Tsunami [= Stage - (tide + low-frequency residual)]'))
        grid(col='orange')
        # 
        plot(tgf$time, tgf$stage, t='l', main=paste0(tgf_name, 
            ': Stage [black] & tidal_prediction without low-frequency residual [red]'))
        points(tgf$time, tgf$tidal_prediction + mean(tgf$tide, na.rm=TRUE), t='l', col='red')

        # Zoom around 15th-16th January (after volcanic explosion) 
        t0 = explosion_start_time
        t1 = explosion_start_time + as.difftime(48, units='hours')
        k = which(tgf$time > t0 & tgf$time < t1)
        plot(tgf$time[k], tgf$tsunami[k], t='l', main=paste0(tgf_name, 
            ': Two alternative de-tiding techniques [red = spectral; black = tide-pred & spectral]'))
        points(tgf$time[k], tgf$tsunami_highpass[k], t='l', col='red')

    }
}

