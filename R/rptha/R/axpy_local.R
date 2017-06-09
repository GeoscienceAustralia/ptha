#' Interface to BLAS axpy
#'
#' Replaces y with a*x + y. Both x and y must be double and have the same
#' length. On large arrays, this is faster than using R to do "y = a*x+y"
#'
#' @param y matrix or array
#' @param a constant of length 1
#' @param x matrix or array with same length as y
#' @return Nothing, but updates the value of y in-place
#' 
#' @export
#' @examples
#' y = 1.0*(1:10) # Make double
#' x = 1.0*(21:30) # Make double
#' axpy_local(y, 2.0, x)
#' stopifnot(all(abs(y - (2*(21:30) + 1:10)) < 1.0e-10))
#'
axpy_local<-function(y, a, x){
    n = as.integer(length(y))
    a = as.double(a)
    if(!is.double(y)) stop('y must be double')
    if(!is.double(x)) stop('x must be double')
    if(length(x) != n) stop('x and y must have same length')

    .Call('axpy_c', y, a, x, n, PACKAGE='rptha')
}

