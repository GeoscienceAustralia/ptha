#' Call the fortran kajiura_convolution routine via C
#'
#' Internal routine to speed up kajiura filter. Tried to efficiently
#' implement the inner loop for the kajiura filter. See the code in
#' kajuria_filter to understand in detail the arguments.
#'
#' @param depth_inv matrix with inverse of depth (no NaN's)
#' @param initial_deformation_padded matrix with initial deformations, padded
#' with rows and columns so the filter can be applied to all internal points (dimension lny+2lfy, lnx + 2lfx)
#' @param filterXYr matrix with dimensions of filter, containing distances from the central point
#' @param kajiuraGmax maximum input value to the kajiuraG function
#' @param lfx (see lfy) 
#' @param lfy filterXYr has dimension = c(lfy, lfx)
#' @param lnx (see lny)
#' @param lny depth_inv has dimensions = c(lny, lnx) 
#' @param output matrix with dimension = c(lny, lnx). Output goes here.
#' @return Nothing but modify outputs
#'
kajiura_convolution<-function(depth_inv, initial_deformation_padded, filterXYr, 
    kajiuraGmax, lfx, lfy, lnx, lny, output){
    # SUBROUTINE kajiura_convolution(depth_inv, initial_deformation_padded, filterXYr, &
    #     kajiuraGmax, lfx, lfy, lnx, lny, output) BIND(C, name='kajiura_convolution')

    .Call('kajiura_convolution_c', depth_inv, initial_deformation_padded, filterXYr,
        kajiuraGmax, as.integer(lfx), as.integer(lfy), as.integer(lnx), as.integer(lny), output)

}

