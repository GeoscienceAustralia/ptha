#' Add ticks for a log-axis to a plot. 
#'
#' The tick locations will be similar to
#' ...,0.1,0.2,0.3,...,0.9,1,2,3,4,...,9,10,20,30,40,.... 
#'
#' @param side integer in 1,2,3,4, giving the side to put the ticks on (see ?axis for details)
#' @param lower lower tick limit
#' @param upper upper tick limit
#' @param ... further arguments to \code{axis}
#' @return nothing, but adds ticks to the axis
#' @export
#' @examples
#' plot(10**runif(1000), log='y')
#' add_log_axis_ticks(side=2)
#'
add_log_axis_ticks<-function(side, lower=1e-30, upper=1e+30, ...){

    row_vec = 10**(seq(log10(lower), log10(upper)))
    col_vec = 1:9
    tick_places = c(t(outer(row_vec, col_vec)))

    axis(side=side, at=tick_places, labels=FALSE, tick=TRUE, ...)
}

