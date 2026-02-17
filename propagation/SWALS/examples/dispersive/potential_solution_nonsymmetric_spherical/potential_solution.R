#
# Solve the potential flow equation for the free surface
# Equation 22 from 
#     Saito, T. Dynamic tsunami generation due to sea-bottom deformation: 
#     Analytical representation based on linear potential theory Earth Planets Space, 2013, 65, 1411-1423
#

# Domain = Box with side-length L
L = 200000 
g = 9.8
h0 = 4000 # Depth -- uniform
N = 801

xs = seq(-L/2, L/2, len=N)
ys = seq(-L/2, L/2, len=N)

dxs = diff(range(xs))/(N-1)
dys = diff(range(ys))/(N-1)

# Lon/lat corresponding to origin of xs/ys
central_lon = 0
central_lat = 45
radius_earth = 6371000
ys_lat = ys/(2*pi*radius_earth)*360 + central_lat
xs_lon = xs/(2*pi*radius_earth*cos(central_lat/180*pi))*360 + central_lon
dxs_sph = diff(range(xs_lon))/(N-1)
dys_sph = diff(range(ys_lat))/(N-1)

# Instantaneous displacement - NON-SYMMETRIC (elliptical/asymmetric form)
# Use different scales in x and y directions to break radial symmetry
d = matrix(0, ncol=length(xs), nrow=length(ys))
for(j in 1:ncol(d)){
    # Elliptical Gaussian with different standard deviations in x and y
    # Plus an asymmetric oscillatory component
    x_dist = xs - 0
    y_dist = ys[j] - 0

    # Elliptical component: different scales in x and y directions
    elliptical_term = exp(-0.5 * ((x_dist/(3*h0))**2 + (y_dist/(5*h0))**2)) * 2.0

    # Asymmetric oscillatory term (depends on x and y independently)
    r = sqrt(x_dist**2 + y_dist**2)
    oscillatory_term = 0.3 * sin(2*pi*x_dist/(8*h0)) * cos(2*pi*y_dist/(10*h0)) * (r < 16*h0)

    d[,j] = elliptical_term + oscillatory_term
}

# Wavenumbers
kx = 2*pi*pmin(0:(N-1), N - (0:(N-1)) )/L
ky = 2*pi*pmin(0:(N-1), N - (0:(N-1)) )/L

# k = sqrt(kx^2 + ky^2)
k = d*0
for(j in 1:ncol(k)){
    k[,j] = sqrt(kx**2 + ky[j]**2)
}

w0 = sqrt(g * k * tanh(k * h0))

free_surface_solution<-function(t){
    # Eqn 25 of Saito (2013) when H(t) = 1 and d(delta(t))/dt = 0
    # (i.e. away from the generation zone, after the rise time)

    if(t > 0){
        d_fft = fft(d)
        inv_fft = 1/prod(dim(d)) * Re( fft(d_fft/cosh(k*h0) * cos(w0*t), inverse=TRUE) )
    }else{
        stop('Not implemented for t <= 0')
    }
    return(inv_fft)
}

pressure_residual_bed_solution<-function(t, rho0 = 1024){
    # Eqn 25 of Saito (2013) when H(t) = 1 and d(delta(t))/dt = 0
    # (i.e. away from the generation zone, after the rise time)
    #
    # Note that to get the wave-height from the pressure oscillation, we should
    # multiply by 'cosh(k * h0) / (rho * g)'. This approaches the usual adjustment
    # (1/(rho * g)) for sufficiently small (k*h0), but is nearly 2x greater for a
    # wave-period of 2min in a depth of 4000m. 
    # library(rptha)
    # H0 = 4000
    # cosh(2*pi * H0/airy_wavelength(period=seq(120, 600, by=60), h=H0, g=9.81))
    #[1] 1.967542 1.310377 1.157754 1.096538 1.065464 1.047419 1.035976 1.028249
    #[9] 1.022780
    ## The effect is weaker in shallower water, as expected
    # H0 = 2000
    # cosh(2*pi * H0/airy_wavelength(period=seq(120, 600, by=60), h=H0, g=9.81))
    # [1] 1.359862 1.138280 1.074144 1.046435 1.031873 1.023254 1.017724 1.013961
    # [9] 1.011284
    #
    # At sites like NZJ which show clear late dispersive waves from Hunga-Tonga, and
    # around 14 waves per hour in 2 km depth, the effect is likely small:
    # H0 = 2000
    # cosh(2*pi*H0/airy_wavelength(period=3600/14, h=H0, g=9.81))
    # # [1] 1.064086
    #
    # But at a site like NZG, in about 6000m depth, and also about 14 waves/hr, the effect
    # is becoming significant.
    # H0 = 6000
    # cosh(2*pi*H0/airy_wavelength(period=3600/14, h=H0, g=9.81))
    #[1] 1.214376




    if(t > 0){
        d_fft = fft(d)
        inv_fft = rho0 * g/prod(dim(d)) * Re( fft(d_fft/(cosh(k*h0)**2) * cos(w0*t), inverse=TRUE) )
    }else{
        stop('Not implemented for t <= 0')
    }
    return(inv_fft)
}


# Solution for very small time
sol0 = free_surface_solution(0.0001)

# # Later times -- do not let the wave exit the boundaries
# maxT = 1500
# dT = 30
# Nt = round(maxT/dT)
# sols = array(NA, dim=c(dim(sol0), Nt))
# for(i in 0:(Nt-1)){
#     if(i == 0){
#         sols[,,1] = sol0
#     }else{
#         sols[,,i+1] = solution(i*dT)
#     }
# }

#
suppressMessages(library(raster))
r1 = raster(sol0, xmn=min(xs_lon)-dxs_sph/2, xmx=max(xs_lon)+dxs_sph/2, ymn=min(ys_lat)-dys_sph/2, ymx=max(ys_lat)+dys_sph/2)
writeRaster(r1, file='initial_condition_file_spherical.tif', overwrite=TRUE, options=c('COMPRESS=DEFLATE'))

make_final_solution<-function(comptime){
    # Solution to compare against
    # Ensure the solution at the end time is small enough that the flow has not reached
    # the boundary -- because these analytical calcuations are using periodic boundary
    # conditions, whereas SWALS will not be.
    #comptime = 100 
    solend = free_surface_solution(comptime)

    #r1 = raster(solend, xmn=min(xs), xmx=max(xs), ymn=min(ys), ymx=max(ys))
    #writeRaster(r1, file=paste0('solution_time_', comptime, '.tif'), overwrite=TRUE, options=c('COMPRESS=DEFLATE'))
    return(solend)
}

# The 'raster' function interprets the orientation differently to R's image
# function, and the following can account for that.
flip_like_raster<-function(mat){ a = t(mat); b = a[,ncol(a):1]; b}
