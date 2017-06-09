#' Interface to BLAS axpy
#'
#' Replaces y with a*x + y. Both x and y must be double and have the same
#' length
#'
#' @param y matrix or array
#' @param a constant of length 1
#' @param x matrix or array with same length as y
#' @return Nothing, but updates the value of y in-place
#' 
#' @export
axpy_local<-function(y, a, x){
    n = as.integer(length(y))
    a = as.double(a)
    if(!is.double(y)) stop('y must be double')
    if(!is.double(x)) stop('x must be double')
    if(length(x) != n) stop('x and y must have same length')

    .Call('axpy_c', y, a, x, n, PACKAGE='rptha')
}

