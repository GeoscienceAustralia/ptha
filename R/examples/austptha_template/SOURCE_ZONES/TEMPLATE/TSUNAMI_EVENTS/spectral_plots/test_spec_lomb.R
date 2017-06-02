#
# Get spectral energy associated with different frequency bands
#
energy_banding<-function(data_t, data_s, interp_dt = 15, bin_divisors = c(1/1200, 1/3600, 1/7200), 
    exclude_mean=TRUE){

    # Interpolate the data at the given dt
    interp_t = approx(data_t, data_s, xout=seq(min(data_t), max(data_t), 
        by=interp_dt))

    # These two 'energies' are equal
    #signal_energy = sum(interp_t$y^2)
    #spectral_energy = sum(Mod(fft(interp_t$y))**2)/length(interp_t$y)

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

test_spec_lomb<-function(){
    # Test that my understanding of spec.lomb is correct
    library(spectral)

    # Make data with spectral peaks at freq = 1/3600 and 1/500
    data_t = c(seq(0, 2040, by=60), seq(2100, 7200, by=120), 
        7200 + 60 * 15, 7200 + 2 * 60 * 15, 
        7200 + 2 * 60 * 16 + seq(0, 7200, by=60))
    data_s = 0.0 * sin(2*pi*data_t/500) + 1.2 * sin(2*pi*data_t/3600)

    data_spec = spec.lomb(data_s, data_t, ofac=1)

    # Plot shows that we get the correct frequencies
    #
    # But, notice how the 'p-value' at the 1/500 frequency is not significant. This
    # will depend on the test being used.
    #
    plot(data_spec)
    abline(v=1/3600, col='red')
    abline(h=1, col='red')
    abline(v=1/500, col='green')
    abline(h=0.3, col='green')
}
