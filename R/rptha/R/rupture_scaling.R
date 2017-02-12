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
#' @param relation Name for the scaling relation ('Strasser' uses the subduction
#' interface event relation for Strasser et al 2010.)
#' @param detailed logical. If False return a vector with area/width/length,
#' otherwise provide a list with the latter as well as information on
#' log10-standard-deviations
#' @param CI_sd Logical. If detailed = TRUE, the output includes a positive
#' and negative confidence interval threshold, both of which are CI_sd
#' standard deviations away from the mean (in log space where the regression is
#' computed)
#' @return A numeric vector with the area/width/length (if detailed = FALSE), otherwise
#' a list with the rupture size statistics as well as upper and lower bounds
#' for a confidence interval, and information on the log10 standard deviation
#' of each variable (which can be used to compute any other confidence interval).
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
        # Store coefficients too
        output$area_absigma = area_absigma
        output$width_absigma = width_absigma
        output$length_absigma = length_absigma
    }

    return(output)
}

#' Compute the inverse of Mw_2_rupture_size, given area as an input
#'
#' Given an area, this function computes the Mw value such that
#'\code{Mw_2_rupture_size(Mw, relation=relation, detailed=TRUE, CI_sd = CI_sd) = area}. 
#' It currently does not give information on length or width.
#' 
#' @param area numeric area
#' @param relation Type of scaling relation used ('Strasser')
#' @param CI_sd numeric (can be positive or negative). Positive values correspond to
#' lower Mw, negative values to higher Mw.
#' @return values of Mw
#' @export
#' @examples
#'    Mw = 8.0
#'    # Get detailed information on the expected rupture size range
#'    area0 = Mw_2_rupture_size(Mw, detailed=TRUE, CI_sd = 2)
#'    # Find Mw such that area0$values[1] is a lower 2-sigma area
#'    Mw_squeezed = Mw_2_rupture_size_inverse(area0$values[1], CI_sd = -2)
#'    # Confirm that it worked
#'    area1 = Mw_2_rupture_size(Mw_squeezed, detailed=TRUE, CI_sd = 2)
#'    # The minus_CI component of area1 should equal area0
#'    stopifnot(abs(area1$minus_CI[1] - area0$values[1]) < 1.0e-04)
Mw_2_rupture_size_inverse<-function(area, relation='Strasser', CI_sd = 0){

    if(length(CI_sd) > 1) stop('length(CI_sd) must = 1')
  
    # log10(area) = Mw*area_coef[2] + area_coef[1] + CI_sd*area_coef[3] 
    area_coef = Mw_2_rupture_size(6.0, detailed=TRUE)$area_absigma 

    Mw = (log10(area) - area_coef[1] - CI_sd*area_coef[3])/area_coef[2]

    names(Mw) = NULL

    return(Mw)
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
#' @param constant The value of constant passed to \code{M0_2_Mw}
#' @return slip in m
#' @examples
#' x = slip_from_Mw(9.0) # Should be roughly 10m
#' @export
slip_from_Mw<-function(Mw, mu=3e+10, 
    area_function=function(Mw){Mw_2_rupture_size(Mw)[1]}, constant=9.05){

        # Area in m^2
        area = sapply(Mw, area_function)*1e+06

        M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

        slip = M0/(mu*area)

        # Get rid of 'area' name
        names(slip) = NULL

        return(slip)
}

