#' Incomplete gamma function 
#'
#' The incomplete gamma function is "integral [ t^{a-1} exp(-t) dt] from t=x to t=Infinity".
#' Actually it is defined differently in different literature. This definition is
#' 'unnormalised' and uses the upper tail. The implementation assumes positive arguments. 
#'
#' @param a the shape parameter
#' @param x the lower limit on the integral in the definition
#' @importFrom stats pgamma
#' @export
.incomplete_gamma<-function(a, x) exp( pgamma(x, a, lower=FALSE, log.p=TRUE) + lgamma(a) )

#' Tapered Gutenberg Richter exceedance-rate
#'
#' Compute the rate of events with Moment > M according to the Tapered Gutenberg
#' Richter model
#'
#' @param M the moment of interest
#' @param N_Mt the number of events with moment exceeding the threshold moment Mt
#' @param Mt the threshold moment
#' @param Mc the corner moment
#' @param beta the beta parameter
#' @return The rate of events with moment > M
#' @export
#' @examples
#' # See example in ?taperedGR_exceedance_rate_derivative
#'
taperedGR_exceedance_rate<-function(M, N_Mt, Mt, Mc, beta){
    out = ifelse(M < Mt,
        N_Mt,
        N_Mt * (Mt/M)^beta * exp((Mt - M)/Mc)
    )
    return(out)
}

#' Derivative of Tapered Gutenberg Richter exceedance-rate
#'
#' This gives the derivative of the exceedance-rate with respect to the moment,
#' at the specified moment M
#'
#' @param M the moment of interest
#' @param N_Mt the number of events with moment exceeding the threshold moment Mt
#' @param Mt the threshold moment
#' @param Mc the corner moment
#' @param beta the beta parameter
#' @return The derivative
#' @export
#' @examples
#'  # Threshold moment
#'  Mt = M0_2_Mw(5.5, inverse=TRUE)
#'  # Corner moment
#'  Mc = M0_2_Mw(9.2, inverse=TRUE)
#'
#'  # How many events with M > Mt each year on average?
#'  N_Mt = 10.47
#'  # The beta parameter
#'  beta = 0.54
#'
#'  #
#'  # Double check that dN_dM is right by comparison with numerical derivative
#'  #
#'  x = M0_2_Mw(7.5, inverse=TRUE)
#'  # Numerical derivative -- needs a large 'delta h' because x is very large
#'  dh = x/1e+06 # A suitable numerical derivative increment
#'  dN_dM_estimate = (taperedGR_exceedance_rate(x+dh, N_Mt, Mt, Mc, beta)-
#'                    taperedGR_exceedance_rate(x-dh, N_Mt, Mt, Mc, beta)
#'                   )/(2*dh)
#'  dN_dM_exact = taperedGR_exceedance_rate_derivative(x, N_Mt, Mt, Mc, beta)
#'  test = (abs(dN_dM_estimate - dN_dM_exact) < 1e-06*abs(dN_dM_exact))
#'  stopifnot(test)
#'
#'  # Quick check of the vectorized version
#'  dN_dM_check = taperedGR_exceedance_rate_derivative(c(0, x, x, 1), N_Mt, Mt, Mc, beta)
#'  stopifnot(all(dN_dM_check[c(1, 4)] == 0) & 
#'            all(dN_dM_check[2:3] == dN_dM_exact))
#'
taperedGR_exceedance_rate_derivative<-function(M, N_Mt, Mt, Mc, beta) {
    out = ifelse(M < Mt,
        0,
        N_Mt * ( -( (Mt/M)^beta * exp( (Mt-M)/Mc)*(beta * Mc + M)/(Mc*M)   ))
    )
    return(out)
}

#' Long-term average moment release of tapered Gutenberg Richter model
#'
#' This is equal to the integral of [ M *
#' (-taperedGR_exceedance_rate_derivative) ] dM from Mt to Infinity
#'
#' @param N_Mt the number of events with moment > threshold moment Mt
#' @param Mt the threshold moment
#' @param Mc the corner moment
#' @param beta the beta parameter
#' @return The long-term moment release rate (moment/year)
#' @export
#' @examples
#'    # Threshold moment
#'    Mt = M0_2_Mw(5.5, inverse=TRUE)
#'    # Corner moment
#'    Mc = M0_2_Mw(9.2, inverse=TRUE)
#'
#'    # How many events with M > Mt each year on average?
#'    N_Mt = 10.47
#'    # The beta parameter
#'    beta = 0.54
#'
#'    # Numerically integrate the moment release rate, and compare with the analytical formula
#'    f<-function(M) - M * taperedGR_exceedance_rate_derivative(M, N_Mt, Mt, Mc, beta)
#'
#'    # Use "Mc *10000" as an approximation of infinity
#'    numerical_integral = integrate(f, Mt, Mc*10000, rel.tol=1e-08)
#'
#'    exact_integral = taperedGR_moment_release_rate(N_Mt, Mt, Mc, beta)
#'
#'    error = abs(exact_integral - numerical_integral$value)
#'    stopifnot(error < numerical_integral$value*1e-07)
taperedGR_moment_release_rate<-function(N_Mt, Mt, Mc, beta){
    N_Mt * ( 
        Mt^beta * Mc^(1-beta) * exp(Mt/Mc) * 
        ( beta*.incomplete_gamma(1-beta, Mt/Mc) + .incomplete_gamma(2-beta, Mt/Mc)) 
    )
}

#' Corner moment for tapered Gutenberg Richter model
#'
#' Given a specified long-term moment release rate and other taperedGR
#' parameters, compute the corner moment Mc.
#'
#' @param L - specified long-term moment release rate [e.g. from coupled fraction of tectonic convergence]
#' @param N_Mt - rate of events with moment exceeding Mt
#' @param Mt - threshold moment Mt
#' @param beta - taperedGR beta parameter
#' @param simple if TRUE then implement the approximation in Kagan (2002),
#' which is simple and inexact, but quite accurate when Mt << Mc. If FALSE
#' then solve for Mc using root finding
#' @return corner moment
#' @importFrom stats uniroot
#' @export
#' @examples
#'    # Threshold moment
#'    Mt = M0_2_Mw(5.5, inverse=TRUE)
#'    # Corner moment
#'    Mc = M0_2_Mw(9.2, inverse=TRUE)
#'
#'    # How many events with M > Mt each year on average?
#'    N_Mt = 10.47
#'    # The beta parameter
#'    beta = 0.54
#'
#'    L = taperedGR_moment_release_rate(N_Mt, Mt, Mc, beta)
#'
#'    back_computed_Mc = taperedGR_Mc_from_moment_release_rate(L, N_Mt, Mt, beta, simple=FALSE)
#'
#'    stopifnot(abs(back_computed_Mc - Mc) < Mc*1.0e-07)
#'
#'    back_computed_Mc_approx = taperedGR_Mc_from_moment_release_rate(L, N_Mt, Mt, beta, simple=TRUE)
#'    stopifnot(abs(back_computed_Mc - back_computed_Mc_approx) < 0.01*back_computed_Mc)
#'
taperedGR_Mc_from_moment_release_rate<-function(L, N_Mt, Mt, beta, simple=FALSE){
    
    if(simple){
        # This implements the approximations in Kagan (2002),
        # basically we set (Mt/Mc $\simeq$ 0) in a number of terms.
        # It's a pretty good approximation if Mt << Mc
        Mc =( (L * (1-beta)) / 
             (N_Mt * Mt^beta * .incomplete_gamma(2-beta, 0)))**(1/(1-beta))
        return(Mc)
    }else{

        if(length(L) > 1 | length(N_Mt) > 1 | length(Mt) > 1 | length(beta) > 1){
            stop('taperedGR_Mc_from_moment_release_rate cannot yet take vector arguments unless simple=TRUE -- you need to loop')
        }

        # This is the exact solution
        # Make a guess based on the simple approximation
        lower = taperedGR_Mc_from_moment_release_rate(L*0.1, N_Mt, Mt, beta, simple=TRUE)
        upper = taperedGR_Mc_from_moment_release_rate(L*10,  N_Mt, Mt, beta, simple=TRUE)
        # Function for root-finding: 
        getL<-function(Mc_local){
            L - taperedGR_moment_release_rate(N_Mt, Mt, Mc_local, beta)
        }

        # Find the value of Mc that gives long term moment release rate L
        fit_mc = uniroot(getL, c(lower, upper), tol=1e-08)
        return(fit_mc$root)
    }

}

