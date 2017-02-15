#' Gauge time-series filtered range
#'
#' Compute the range of a gauge stage time-series after spectral filtering of high frequencies.
#' This is practically useful for avoiding seismic wave signals in some DART buoys.
#' The results should always be checked visually though -- since it cannot automatically
#' treat all 'bad-behaviour' in the observations.
#'
#' @param data_t vector of numeric times (usually in seconds)
#' @param data_s vector of numeric stages corresponding to data_t (usually in metres)
#' @param filter_freq We remove frequencies above filter_freq from the series
#' before computing the range.
#' @param interp_dt numeric time-step (same units as data_t). We interpolate
#' the data to a fixed time-step given by interp_dt, before performing the
#' discrete fourier transform. If NULL, the minimum spacing in data_t is used.
#' @return vector giving the min and max of the filtered stage
#'
#' @examples
#' # Make a wave train with 2 spectral components, having frequencies 1 and 1/50
#' t = seq(1,1000,by=0.1)
#' stage = 0.2*sin(2*pi*t) + 1.0*sin(2*pi*t/50)
#'
#' range(stage) # Should be about +-1.2
#' # Get the stage range after filtering out frequencies > 1/10. By construction it
#' # should be about +- 1.0
#' longperiod_range = gauge_range_filtered(t, stage, filter_freq = 1/10)
#' stopifnot(all(abs(longperiod_range - c(-1,1)) < 1.0e-02))
#' # Check that if we filter frequencies > 1, then it still works
#' full_range = gauge_range_filtered(t, stage, filter_freq = 1/0.5)
#' stopifnot(all(abs(full_range - range(stage)) < 1.0e-02))
#'
#' @export
gauge_range_filtered<-function(data_t, data_s, filter_freq = 1/(2*60), interp_dt = NULL){

    if(is.null(interp_dt)) interp_dt = min(diff(data_t))

    # Interpolate the data at the given dt
    interp_t = approx(data_t, data_s, xout=seq(min(data_t), max(data_t), 
        by=interp_dt))

    # Frequencies of specral components
    interp_f = pmin(0:(length(interp_t$x)-1), 
        length(interp_t$x) - 0:(length(interp_t$x)-1))/length(interp_t$x) * 1/interp_dt

    # Fourier transform
    series_fft = fft(interp_t$y)
    # ... zero the high frequencies
    series_fft[interp_f > filter_freq] = 0
    # ... back-transform
    filtered_series = fft(series_fft, inverse=TRUE)/length(interp_t$y)

    # Make sure we didn't accidently create complex numbers! [i.e. code bug].
    if(max(Im(filtered_series)) > 1.0e-06) stop('Non-negligable complex values in filter')
    filtered_series = Re(filtered_series)

    return(range(filtered_series))
}

#' Zero-crossing-period of a gauge time-series
#'
#' This function computes both the up-crossing and down-crossing periods and
#' returns their average
#'
#' @param data_t vector of numeric times (usually in seconds)
#' @param data_s vector of numeric stages corresponding to data_t (usually in metres)
#' @param interp_dt numeric time-step (same units as data_t). We interpolate
#' the data to a fixed time-step given by interp_dt, before computing the
#' zero crossing period. If NULL, the minimum spacing in data_t is used.
#' @return The zero crossing period of the gauge data
#' @export
#' @examples
#' x = seq(1,2000, by=2)
#' period = 30
#' y = sin(x*2*pi/period)
#' # Check we can compute this
#' y_period = gauge_zero_crossing_period(x, y)
#' stopifnot(abs(y_period/period - 1) < 1.0e-03)
#'
gauge_zero_crossing_period <-function(data_t, data_s, interp_dt=NULL){

    if(is.null(interp_dt)) interp_dt = min(diff(data_t))

    # Interpolate the data at the given dt
    interp_t = approx(data_t, data_s, xout=seq(min(data_t), max(data_t), 
        by=interp_dt))
 
    # Below here, code based on earlier zero_crossing_period function 
    dt = interp_dt
    x = interp_t$y  

    sg_x = sign(x)

    # Get 'positive' zero crossings
    up_cross = which(diff(sg_x) < 0)
    # Get 'negative' zero crossings
    down_cross = which(diff(sg_x) > 0)

    lu = length(up_cross)
    ld = length(down_cross) 

    if( lu < 2 | ld < 2){
        #warning('To few zero crossings for period computation')
        return(NA)
    }

    # Compute periods
    up_period = (up_cross[lu] - up_cross[1])*dt/(lu - 1)
    down_period = (down_cross[ld] - down_cross[1])*dt/(ld - 1)

    return(0.5*(up_period + down_period))
}


#' Spectral energy in frequency bands of gauge time-series
#'
#' Get the spectral energy associated with different frequency bands in a wave
#' time-series. Suppose data_s is a vector of stages recorded at an even time
#' interval interp_dt. From the properties of the fourier transform, the
#' following two quantities are equal: \cr
#' \code{signal_energy = sum(data_s^2)*interp_dt} \cr
#' \code{spectral_energy = sum(Mod(fft(data_s))**2)/length(data_s) * interp_dt} \cr
#' The current code splits the 'spectral_energy' quantity above into a number of bins
#' containing different frequencies. The divisions between the bins are defined
#' by \code{bin_divisors}. By default the mean of data_s is excluded.
#'
#' @param data_t vector of numeric times at which stages are recorded (usually in seconds)
#' @param data_s vector of numeric stages corresponding to data_t (usually in metres)
#' @param bin_divisors compute the spectral energy in bands which are separated
#' by these divisors. They should be ordered highest-to-lowest. 
#' @param interp_dt numeric time-step (same units as data_t). We interpolate
#' the data to a fixed time-step given by interp_dt, before performing the
#' discrete fourier transform. If NULL, the minimum spacing in data_t is used.
#' @param exclude_mean logical -- do not include the mean of data_s in any energy band
#' @return vector of length = length(bin_divisors)+1 giving the spectral energy
#' in a range of frequency bands. The first entry corresponds to the highest
#' frequency bin, containing frequencies ranging from (Infty, bin_divisors[1]).
#' The next value includes frequencies from (bin_divisors[1], bin_divisors[2]).
#' The last includes frequencies from (bin_divisors[n], eps), where n is the
#' length of bin_divisors, and eps is a number slightly larger than zero (which
#' ensures the mean is not included in the lowest bin).
#'
#' @export
#' @examples
#' dt = 0.5
#' t = seq(0, 2000*2*pi, by=dt)
#' # Make some frequencies to add together
#' s1 = 0.3*sin(2*pi*t/5) # Frequency of 1/5
#' s2 = 1.2*sin(2*pi*t/30) # Frequency of 1/30
#' s3 = 2.0*sin(2*pi*t/100) # Frequency of 1/100
#' stage = s1 + s2 + s3 
#' 
#' # gauge_energy_banding can be used to pick out the energy in s1, s2, s3, with 
#' # a suitable choice of the bin divisors.
#' energy_bands = gauge_energy_banding(t, stage, bin_divisors=c(1/10, 1/50))
#' 
#' # Excluding the mean, all the energy should be accounted for
#' stopifnot( abs(sum(energy_bands) - sum((stage-mean(stage))^2)*dt) < 1.0e-04)
#' 
#' # Each bin should have almost the same energy as the corresponding component,
#' # with some error caused by spectral leakage in the dft
#' stopifnot( abs( energy_bands[1] - sum(s1^2)*dt) < 0.01*energy_bands[1])
#' stopifnot( abs( energy_bands[2] - sum(s2^2)*dt) < 0.01*energy_bands[2])
#' stopifnot( abs( energy_bands[3] - sum(s3^2)*dt) < 0.01*energy_bands[3])
#'
gauge_energy_banding<-function(data_t, data_s, bin_divisors, interp_dt = NULL, exclude_mean=TRUE){

    stopifnot(max(diff(bin_divisors)) < 0)

    if(is.null(interp_dt)) interp_dt = min(diff(data_t))

    # Interpolate the data at the given dt
    interp_t = approx(data_t, data_s, xout=seq(min(data_t), max(data_t), 
        by=interp_dt))

    ## NOTE: These two 'energies' are equal
    ## signal_energy = sum(interp_t$y^2)
    ## spectral_energy = sum(Mod(fft(interp_t$y))**2)/length(interp_t$y)

    # Frequencies of specral components
    interp_f = pmin(0:(length(interp_t$x)-1), 
        length(interp_t$x) - 0:(length(interp_t$x)-1))/length(interp_t$x) * 1/interp_dt

    # Power of each component. By multiplying by interp_dt, we make the results
    # roughly independent of the choice of interp_dt. Otherwise mod_sq would
    # scale inversely with interp_dt. The latter still has some influence on
    # the result though, since it choice will change the interpolation. 
    mod_sq = Mod(fft(interp_t$y))**2/length(interp_t$y) * interp_dt

    if(exclude_mean){
        # Make bins -- do not include zero (mean)
        bin_divisors = c(Inf, bin_divisors, min(bin_divisors) * 1.0e-20)
    }else{
        bin_divisors = c(Inf, bin_divisors, 0)
    }
    energy_band = rep(NA, len=length(bin_divisors)-1)
    for(i in 1:(length(energy_band))){
        energy_band[i] = sum(mod_sq[which(interp_f <= bin_divisors[i] & interp_f > bin_divisors[i+1])])
    }

    return(energy_band)
}

