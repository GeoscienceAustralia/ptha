#' Compute Moment Magnitude Mw from Seismic Moment M0, and vice-versa
#'
#' By default we use the equation Mw = 2/3 * (log10(M0) - 9.05), see Hanks and Kanamori
#' (1977) and Bird and Kagan (2004).  Note that slightly different values of the
#' final constant (9.05) are sometimes used, so we provide the option to override this.
#'
#' @param M0 Seismic Moment (units of Nm) if inverse = FALSE, else Moment Magnitude
#' @param inverse logical. If FALSE, return Mw given M0. If TRUE, return M0 given Mw
#' @param constant Optionally override the value "9.05" used by default in the equation.
#' @return Mw if inverse = FALSE, or M0 (units of Nm) if inverse=TRUE
#' @export
#' @examples
#' Mw = M0_2_Mw(4e+17)
#' M0 = M0_2_Mw(Mw, inverse=TRUE)
#' stopifnot(isTRUE(all.equal(M0, 4e+17)))
#'
M0_2_Mw<-function(M0, inverse=FALSE, constant=9.05){

    if(inverse){
        # Input was actually Mw
        Mw = M0

        M0 = 10^(Mw*3/2 + constant)

        return(M0)
    }else{

        Mw = 2/3* ( log10(M0) - constant)

        return(Mw)
    }

}

#' Compute earthquake rupture area, width and length from Mw based on an
#' empirical scaling relation
#'
#' Units are km and km^2
#'
#' @param Mw Moment Magnitude 
#' @param relation Name for the scaling relation ('Strasser' uses the interface
#' event relation for Strasser et al 2010.)
#' @param detailed logical. If False return a vector with area/width/length,
#' otherwise provide a list with the latter as well as information on
#' log10-standard-deviations
#' @param CI_sd Logical. If detailed = TRUE, the output includes a positive
#' and negative confidence interval threshold, both of which are CI_sd
#' standard deviations away from the mean (in log space where the regression is
#' computed)
#' @return A numeric vector
#' @export
#' @examples
#' rupture_statistics1 = Mw_2_rupture_size(9.0)
#' rupture_statistics2 = Mw_2_rupture_size(9.0, detailed=TRUE)
Mw_2_rupture_size<-function(Mw, relation='Strasser', detailed=FALSE,
    CI_sd=1){

    if(relation == 'Strasser'){

        # Area
        area_absigma = c(-3.476, 0.952, 0.304)
        area = 10**(area_absigma[1] + Mw*area_absigma[2])

        # Width
        width_absigma = c(-0.882, 0.351, 0.173)
        width = 10**(width_absigma[1] + Mw*width_absigma[2])

        # Length
        length_absigma = c(-2.477, 0.585, 0.180)
        length = 10**(length_absigma[1] + Mw*length_absigma[2])

    }else{
        stop(paste0('Relation value ', relation, ' not recognized'))
    }

    output = c(area, width, length)
    names(output) = c('area', 'width', 'length')

    if(detailed){
        output = list(values = output)
        output$log10_sigmas = c(area_absigma[3], width_absigma[3], length_absigma[3])
        names(output$log10_sigmas) = c('area', 'width', 'length')
        output$plus_CI = 10**(log10(output$values) + CI_sd*output$log10_sigmas)
        output$minus_CI = 10**(log10(output$values) - CI_sd*output$log10_sigmas)
    }

    return(output)
}

#' Compute mean slip on a rupture of a given area, moment magnitude Mw, and
#'
#' @param Mw Moment magnitude
#' @param area area of rupture (km^2)
#' @param mu Shear Modulus (Pascals)
#' @param constant value of constant passed to \code{M0_2_Mw}
#' @return slip in m
#' @export
#' @examples
#' s0 = slip_from_Mw_area_mu(9.0, 100e+03) # Should be close to 10m
#'
slip_from_Mw_area_mu<-function(Mw, area, mu=3e+10, constant=9.05){

    M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

    area = area*1e+06 # m^2

    slip = M0/(area*mu)

    return(slip)
}

#' Compute slip given Mw (and optionally mu, and a function to compute the area from Mw)
#'
#' This gives an alternative interface to 'slip_from_Mw_area_mu', where a
#' function to compute area from Mw is provided instead of the area. The default
#' interface serves as a basic example of a function to compute mean slip given
#' Mw alone. Functions like this are required in the event_probability calculations
#'
#' @param Mw Moment magnitude
#' @param mu Shear Modulus (Pascals)
#' @param area_function function which returns the area (in km^2) given Mw
#' @return slip in m
#' @examples
#' x = slip_from_Mw(9.0) # Should be roughly 10m
#' @export
slip_from_Mw<-function(Mw, mu=3e+10, 
    area_function=function(Mw){Mw_2_rupture_size(Mw)[1]}){

        # Area in m^2
        area = area_function(Mw)*1e+06

        M0 = M0_2_Mw(Mw, inverse=TRUE)

        slip = M0/(mu*area)

        # Get rid of 'area' name
        names(slip) = NULL

        return(slip)
}

