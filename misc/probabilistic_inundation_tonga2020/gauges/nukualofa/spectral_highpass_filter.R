
#' Remove low-frequencies (below cutoff_frequency) from a time-series
#'
#' This offers a reasonable method to remove non-tsunami components from a
#' signal -- just remove everything with frequency less than a cutoff of (say) 3
#' hours or similar. Beware that unless the signal is periodic, this function
#' lead to some artefacts around the start/end of the time-series (so it is good
#' to provide sufficent buffer around the start/end of the data you are actually
#' interested in).
#'
#' @param data_t times (in seconds) of data samples. Can be unevenly spaced
#' @param data_s data values at each data_t
#' @param interp_dt Interpolate the data to this spacing to estimate the low-frequency component (for FFT).
#' @param cutoff_frequency Frequency (1/s) below which we consider the data 'low-frequency'
#' @return a list containing the "highfreq" filtered data, as well as other useful stuff
#'
spectral_highpass_filter<-function(
    data_t,
    data_s,
    interp_dt = 15,
    cutoff_frequency = 1/(3*3600)){

    # Straight line between first and last data point. If we subtract this, then we ensure
    # the transformed data begins and ends with zero, so is periodic, nice for fft.
    linear_t<-function(t){
        (data_s[length(data_s)] - data_s[1])*(t - data_t[1])/(data_t[length(data_t)] - data_t[1]) + data_s[1]
    }

    # Interpolate the data at the given dt
    interp_t = approx(data_t, data_s - linear_t(data_t), xout=seq(min(data_t), max(data_t), by=interp_dt), rule=2)

    # These two 'energies' are equal
    #signal_energy = sum(interp_t$y^2)
    #spectral_energy = sum(Mod(fft(interp_t$y))**2)/length(interp_t$y)

    # Frequencies of specral components
    interp_f = pmin(0:(length(interp_t$x)-1),
        length(interp_t$x) - 0:(length(interp_t$x)-1))/length(interp_t$x) * 1/interp_dt

    # Get a 'low-frequency' component
    data_fft = fft(interp_t$y)
    low_frequency = Re(fft(data_fft*(interp_f < cutoff_frequency), inverse=TRUE))/length(data_fft)
    # Interpolate back to the original times, and also add-back the linear trend
    low_freq = approx(interp_t$x, low_frequency, xout=data_t, rule=2)$y + linear_t(data_t)
    # Subtract the low_freq from the original data
    high_frequency = data_s - low_freq

    # Several outputs might be of interest
    return(list(data_t = data_t, data_s = data_s, highfreq=high_frequency,
                lowfreq=low_freq, cutoff=cutoff_frequency))
}

