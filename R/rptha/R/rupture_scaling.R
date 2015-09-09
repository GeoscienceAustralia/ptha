#' Compute Moment Magnitude Mw from Seismic Moment M0, and vice-versa
#'
#' We use the equation Mw = 2/3 * (log10(M0) - 9.05), see Hanks and Kanamori
#' (1977) and Bird and Kagan (2004).  Note that slightly different values of the
#' final constant are sometimes used.
#'
#' @param M0 Seismic Moment (units of Nm) if inverse = FALSE, else Moment Magnitude
#' @param inverse. If FALSE, return Mw given M0. If TRUE, return M0 given Mw
#' @return Mw
#' @export
M0_2_Mw<-function(M0, inverse=FALSE){

    if(inverse){
        # Input was actually Mw
        Mw = M0

        M0 = 10^(Mw*3/2 + 9.05)

        return(M0)
    }else{

        Mw = 2/3* ( log10(M0) - 9.05)

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
#' @return A numeric vector
#' @export
Mw_2_rupture_size<-function(Mw, relation='Strasser', detailed=FALSE){

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
    }

    return(output)
}

#' Compute mean slip on a rupture of a given area, moment magnitude Mw, and
#'
#' @param Mw Moment magnitude
#' @param area area of rupture (km^2)
#' @param mu Shear Modulus (Pascals)
#' @return slip in m
#' @export
slip_from_Mw_area_mu<-function(Mw, area, mu=3e+10){

    M0 = M0_2_Mw(Mw, inverse=TRUE)

    area = area*1e+06 # m^2

    slip = M0/(area*mu)

    return(slip)
}

