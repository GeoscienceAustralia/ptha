# Solution of Zhang and Zhu (1994) New solutions for the
# propagation of long waves over variable depth. Journal of fluid mechanics
# 278: 391-406
#
compute_zhang_solution<-function(){
    g = 9.8 # Should not matter, but conceptually useful
    h0 = 300 # Offshore water depth -- should not matter, but conceptually useful

    # Figure 2 panel B, Zhang
    # Only 2 parameters actually matter.
    b_on_a = 4.
    lambda0_on_b = 2
    a = 40000   # Radius of island
    b = b_on_a * a   # Radius of zone where depth < h0
    lambda0 = lambda0_on_b*b # Deep water wavelength
    period = lambda0/(sqrt(g*h0)) # Wave length / wave celerity
    omega = 2 * pi/period # Angular frequency


    # Numerical parameters
    max_m = 300 # Zhang and Zhu suggest this was enough for the most difficult problem they tried
    max_n = 90

    #
    island_slope = h0/(b-a)
    Ls = a # Length scale
    mu_sq = (omega^2 * Ls / g)
    nu_sq = mu_sq / island_slope 
    k0 = sqrt(nu_sq/(b/a - 1)) # Dimensionless wavenumber in 'constant depth region'


    # Eqn 10
    # R_{n} = sum_{m=0}^{Inf} alpha_{m,n} t^{m} 
    # (since c = 0)
    #
    # Here, a_{m,n} is stored in amn_coef[m+1, n+1]
    amn_coef = matrix(NA, nrow=max_m + 1, ncol=max_n + 1)
    amn_coef[1,] = 1 # a_{0, n} By definition, text after eqn 10 
    amn_coef[2,] = -nu_sq # a_{1, n} eqn 12
    amn_coef[3,] = 0.25*(nu_sq**2 - 3*nu_sq + (0:max_n)**2) # a_{2,n}, eqn 13
    for(m in 0:(max_m-3)){
        # Numerator. We offset indices of amn_coef by +1 because R vectors are base 1
        amn_coef[m+3 +1,] = ( 3*(m+2)**2 - nu_sq)*amn_coef[m+2+1,] - 
            (3*(m+1)**2 - (0:max_n)**2)*amn_coef[m+1+1,] + 
            (m^2-(0:max_n)**2)*amn_coef[m+1,]
        # Denominator
        amn_coef[m+3+1,] = amn_coef[m+3+1,]/(m+3)**2
    }

    # Define various special functions and derivatives we need
    # Note there seem to be neat analytic versions of the derivatives
    ## This definition gives a complex number
    hankle_second_kind<-function(x, n) (besselJ(x, n) - (0+1i)*besselY(x,n))
    deriv_hankle_second_kind<-function(x,n, eps=1.0e-04){
        (hankle_second_kind(x+eps, n) - hankle_second_kind(x-eps, n))/(2*eps)
    }
    # Analytical version
    deriv_hankle_second_kindB<-function(x,n){
         -0.5*(hankle_second_kind(x, n+1) - hankle_second_kind(x, n-1))
    }

    bessel_first_kind<-function(x, n) besselJ(x,n)
    deriv_bessel_first_kind<-function(x,n,eps=1.0e-05){
        (bessel_first_kind(x+eps,n) - bessel_first_kind(x-eps,n))/(2*eps)
    }

    # Eqn 15
    Rbar_n<-function(r, n){
        stopifnot(length(r) == 1)

        output = rep(NA, length(n))
        for(i in 1:length(n)){
            # a_{m,n} is in amn_coef[m+1, n+1]
            output[i] = sum(amn_coef[,n[i]+1] * (1-1/r)^(0:(max_m)))
        }
        return(output)
    }
    # Numerical derivative function
    deriv_Rbar_n<-function(r, n, eps=1.0e-05) (Rbar_n(r+eps, n) - Rbar_n(r-eps, n))/(2*eps)

    eps_n = 1 + ((0:max_n) > 0) # 1 when n=0, 2 otherwise

    R0 = b/a # Defined just before eqn 18

    triangle = (
        -k0 * Rbar_n(R0, 0:max_n) * deriv_hankle_second_kind(k0*R0, 0:max_n) + 
        deriv_Rbar_n(R0, 0:max_n) * hankle_second_kind(k0*R0, 0:max_n))
    triangle1 = ((0-1i)**(0:max_n)) * eps_n * (
        Rbar_n(R0, 0:max_n) * k0 * deriv_bessel_first_kind(k0*R0, 0:max_n) - 
        bessel_first_kind(k0*R0, 0:max_n) * deriv_Rbar_n(R0, 0:max_n))

    An = -2*eps_n * ((0-1i)**((0:max_n) + 1)) / (pi * R0 * triangle)
    Cn = triangle1/triangle

    R_n<-function(r, n) An[n+1] * Rbar_n(r, n)

    # Island height vs theta
    # @param theta is an angular coordiniate around the island
    # @param wt is the "time" in [0, 2pi]
    # @param r is the radial coordinate
    island_height<-function(theta, wt, r=1) {
        output = rep(NA, length(theta))
        for(i in 1:length(output)){
            output[i] = sum(An * Rbar_n(r,0:max_n) * cos(theta[i]*(0:max_n))) * exp(-(0+1i)*wt)
        }
        return(output)
    }

    #
    # Figure 2B from Zhang and Zhu paper
    # 'max wave' around the island.
    #
    thetas = seq(0, pi, len=101)
    wts = seq(0, 2*pi, len=100)
    hmat = matrix(NA, nrow=length(thetas), ncol=length(wts))
    for(i in 1:length(thetas)){
        for(j in 1:length(wts)){
            hmat[i,j] = Re(island_height(thetas[i], wts[j]))
        }
    }
    max_wave = apply(hmat, 1, max)
    #plot(180 - thetas/pi*180, max_wave, 
    #    main='Analytical solution: reproduce peak wave in Figure 2B of Zhang and Zhu 1994')
    #grid()



    get_peak_height_at_point<-function(x,y){
        # Note: My theta is shifted by 'pi' from Zhang and Zhu's for some reason,
        # suggesting a sign error in An above. Not sure where this is coming
        # from (?either error in the paper, or above?). Anyway, once shifted the
        # result is obviously correct.
        thetas = pi - atan2(y, x)
        rs = sqrt(x^2+y^2)/a
        wts = seq(0, 2*pi, len=100)

        hmat = matrix(NA, nrow=length(thetas), ncol=length(wts))
        for(i in 1:length(thetas)){
            for(j in 1:length(wts)){
                hmat[i,j] = Re(island_height(thetas[i], wts[j], r=rs[i]))
            }
        }

        max_wave = apply(hmat, 1, max)
        return(max_wave)
    }

    return(get_peak_height_at_point)

}
