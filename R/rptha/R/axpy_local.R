#' Interface to BLAS axpy
#'
#' Replaces y with a*x + y. Both x and y must be double and have the same
#' length. On large arrays, this is faster than using R to do "y = a*x+y".\cr
#' BEWARE THAT THIS PERFORMS IN-PLACE UPDATE WHICH PROBABLY VIOLATES
#' YOUR EXPECTATIONS OF R CODE. For instance, the code \cr
#' \code{x = runif(10); z = runif(10); y = z; axpy_local(y, 2.0, x)} will cause
#' BOTH y and z to be updated by the axpy call! \cr
#' This occurs because R's reference counting scheme makes
#' z and y point to the same memory. To suppress this behaviour, replace the \cr
#' \code{y=z} line above with an operation, e.g. \cr
#' \code{y=z*1} \cr
#'  This will make R create new memory for y, so z will not be updated. \cr
#' While use of this routine requires some care, it can be 3x faster than doing
#' the same calculation in native R, so it is worthwhile for performance bottlenecks.
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

