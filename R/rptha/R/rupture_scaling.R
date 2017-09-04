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
#' Output units are km and km^2
#'
#' @param Mw Moment Magnitude (must have length(Mw) == 1)
#' @param relation Name for the scaling relation. 'Strasser' (default) uses the subduction
#' interface event relation for Strasser et al 2010; 'Strasser-intraslab' uses the
#' subduction intraslab relations of Strasser et al 2010 [for this case, 
#' Strasser et al 2010 suggest the width sigma might be too small, but we make
#' no effort to correct that]; 'AllenHayes' uses the interface relations of
#' Allen and Hayes (2017, Table 2), with sigma for prediction based on the sigma value
#' for the log10(L / or W / or A) of the orthogonal regression. Note this case 
#' has Area and Width being multi-segment linear; 'AllenHayes-inslab' gives the 
#' inslab relations of Allen and Hayes (2017, Table 5)
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
#' # Using the Strasser et al subduction scaling relation by default
#' rupture_statistics1 = Mw_2_rupture_size(9.0)
#' rupture_statistics2 = Mw_2_rupture_size(9.0, detailed=TRUE)
#' # Try Allen and Hayes relation
#' rupture_statistics3 = Mw_2_rupture_size(9.0, relation='AllenHayes', 
#'     detailed=TRUE, CI_sd=2)
Mw_2_rupture_size<-function(Mw, relation='Strasser', detailed=FALSE,
    CI_sd=1){

    # Assumes Mw is not a vector
    stopifnot(length(Mw) == 1)

    if(relation == 'Strasser'){

        # Strasser's subduction relations
        area_absigma = c(-3.476, 0.952, 0.304)
        width_absigma = c(-0.882, 0.351, 0.173)
        length_absigma = c(-2.477, 0.585, 0.180)

    }else if(relation == 'Strasser-intraslab'){

        # Strasser's intraslab relations 
        area_absigma = c(-3.225, 0.890, 0.184)
        width_absigma = c(-1.058, 0.356, 0.067)
        length_absigma = c(-2.350, 0.562, 0.146)

    }else if(relation == 'AllenHayes'){
        # Allen and Hayes, 2017 Interface rupture scaling coef (table 2)
        # 
        #
        # These are based on orthogonal regression,
        # but for prediction we only include the 'x' sigma value from the paper 
        # (which I have confirmed is the error on 'the LHS of the equations as
        # written in the paper')
        #
        # Some relations change above a magnitude threshold
        Mw_above_thresh = (Mw > 8.63)
        if(!Mw_above_thresh){
            area_absigma = c(-5.62, 1.22, 0.256)
        }else{
            area_absigma = c(2.23, 0.31, 0.256)
        }

        length_absigma = c(-2.90, 0.63, 0.182)

        Mw_above_thresh = (Mw > 8.67)
        if(!Mw_above_thresh){
            width_absigma = c(-1.91, 0.48, 0.137)
        }else{
            width_absigma = c(2.29, 0.0, 0.137)
        }

    }else if(relation == 'AllenHayes-inslab'){
        # Allen and Hayes, 2017 Inslab rupture scaling coef (table 5)

        length_absigma = c(-3.03, 0.63, 0.14)
        width_absigma = c(-1.01, 0.35, 0.15)
        area_absigma = c(-3.89,  0.96, 0.19)

    }else if(relation == 'AllenHayes-outer-rise'){
        # Allen and Hayes, 2017 Outer Rise rupture scaling coef (table 5)
        # 

        length_absigma = c(-2.87, 0.63, 0.08)
        width_absigma = c(-1.18, 0.35, 0.08)
        area_absigma = c(-3.89, 0.96, 0.11)

    }else{

        stop(paste0('Relation value ', relation, ' not recognized'))
    }

    area = 10**(area_absigma[1] + Mw*area_absigma[2])
    width = 10**(width_absigma[1] + Mw*width_absigma[2])
    length = 10**(length_absigma[1] + Mw*length_absigma[2])

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
        output$relation = relation
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
#' @param relation Type of scaling relation used (e.g. 'Strasser', see ?Mw_2_rupture_size)
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
    area_coef = Mw_2_rupture_size(6.0, relation=relation, detailed=TRUE)$area_absigma 

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
#' @param relation Scaling relation type. Value of relation passed to \code{Mw_2_rupture_size}
#' @param area_function function which returns the area (in km^2) given Mw
#' @param constant The value of constant passed to \code{M0_2_Mw}
#' @return slip in m
#' @examples
#' x = slip_from_Mw(9.0) # Should be roughly 10m
#' @export
slip_from_Mw<-function(Mw, 
    mu=3e+10,
    relation='Strasser',
    area_function=function(Mw){Mw_2_rupture_size(Mw, relation=relation)[1]}, 
    constant=9.05){

        # Area in m^2
        area = sapply(Mw, area_function)*1e+06

        M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

        slip = M0/(mu*area)

        # Get rid of 'area' name
        names(slip) = NULL

        return(slip)
}

