# Given wave period T and depth h, solve for the wavelength lambda
# The solution satisfies the implicit equation:
#  lambda = g/(2 pi) * T^2 * tanh ( 2 * pi * h / lambda)
#
# Special cases include the deep water limit (h/lambda)--> Inf
# where lambda = T^2 * g/(2*pi)
# (phase speed = T * (g/(2*pi)) ~ 1.56 T)
#
# and the shallow water limit (h/lambda) --> 0
# where lambda = sqrt(g * h) * T
# (phase speed = sqrt(g*h))

#' Compute wavelength given period and depth
#'
#' Assumes linear (Airy) wave theory
#'
#' @param period wave period 
#' @param h depth
#' @param g gravity
#' @param full_output if true return full optimization info
#' @return the wavelength, along with optimization information if full_output=TRUE
#' @export
#' @examples
#'    ## Shallow water limit
#'    wl = 20000
#'    h = 20
#'    period0 = airy_period(wl, h)
#'    wavelength0 = airy_wavelength(period0, h)
#'
#'    stopifnot(abs(wl - wavelength0) < 1.0e-05*wavelength0)
#'    stopifnot(abs(period0 - wl/sqrt(9.81 * h)) < 1.0e-05 * period0)
#'
#'    ## Deep water limit
#'    wl = 2
#'    h = 200
#'    
#'    period0 = airy_period(wl, h)
#'    wavelength0 = airy_wavelength(period0, h)
#'
#'    stopifnot(abs(wl - wavelength0) < 1.0e-05)
#'    stopifnot(abs(wl - period0**2 * 9.81/(2*pi)) < 1.0e-05 * wl)
#'
#'    ## Intermediate water limit
#'    
#'    wl = 200
#'    h = 20
#'    
#'    period0 = airy_period(wl, h)
#'    wavelength0 = airy_wavelength(period0, h)
#'
#'    stopifnot(abs(wl - wavelength0) < 1.0e-05)
#'
#'    # Check that the dispersion relation is satisfied to within an expected tolerance
#'    resid = wavelength0 - period0**2 * 9.81 /(2*pi) * tanh( h * 2 * pi / wavelength0)
#'    stopifnot(abs(resid) < 1.0e-06)
#'
airy_wavelength<-function(period, h, g=9.81, full_output=FALSE){

    lenout = max(length(period), length(h))
    period = rep(period, length.out=lenout)
    h = rep(h, length.out=lenout) 

    to_optim<-function(lambda) lambda - g/(2*pi) * period**2  * tanh(2*pi*h/lambda)

    guess = sqrt(g * h) * period

    output = bisection(to_optim, lower = rep(0, length.out=lenout), upper=1.1*guess)

    if(!full_output) output = output$root

    return(output)
}

#' Compute wave period given wavelength and depth
#' 
#' Assumes linear (Airy) wave theory
#'
#' @param lambda wavelength
#' @param h depth
#' @param g gravity
#' @param full_output if true return full optimization info
#' @return the wave period, along with optimization details if full_output=TRUE
#' @export
#' @examples
#' ## See further examples in the help for airy_wavelength
#' period = airy_period(10000, 100)
#' 
airy_period<-function(lambda, h, g=9.81, full_output=FALSE){
    
    lenout = max(length(lambda), length(h))
    lambda = rep(lambda, length.out=lenout)
    h = rep(h, length.out=lenout) 
    
    to_optim<-function(period) lambda - g/(2*pi) * period**2  * tanh(2*pi*h/lambda)

    guess = lambda / sqrt(g * h)

    output = bisection(to_optim, lower = 0.9*guess, 
        upper=rep(1e+20, length.out=lenout))

    if(!full_output) output = output$root

    return(output)
}


#' Function for vectorized 1d root finding with bisection. 
#'
#' Solves f(x) = 0 where f(x) is a vector of the same length as x,
#' and f(1) is only affected by x(1), f(2) is only affected by
#' x(2), etc. In R, to do this efficiently we need vectorization, hence the
#' current function
#'
#' @param f the function
#' @param lower a vector of lower bounds for x
#' @param upper a vector of upper bounds for x
#' @param ... further arguments to f
#' @param numiter maximum allowed number of iterations
#' @param tolerance = Allowed error in the root
#' @return a vector x with f(x) = 0 to within the tolerance
#' @export
#' @examples
#' #
#' # Here is a function satisfying the assumptions described above
#' #
#' f<-function(x) x^2 - c(2,3,4)
#' 
#' solution = bisection(f, lower=c(0,0,0), upper=c(10,10,10))
#' stopifnot(all(abs(solution$root - sqrt(c(2,3,4))) < 1.0e-06))
#'
bisection <- function(f, lower, upper, ..., numiter=20000, tolerance =
.Machine$double.eps^0.5){

  stopifnot(length(lower) == length(upper))

  stopifnot(all(lower<= upper))

  flower <- f(lower, ...)
  fupper <- f(upper, ...)

  for (n in 1:numiter) {
    mid <- (lower+upper)/2
    fmid <- f(mid, ...)
    #if (all(abs(fmid) < tolerance)) break
    if(all(upper - lower < tolerance)) break
    samesign <- ((fmid<0)&(flower<0))|((fmid>=0)&(flower>=0))
    lower <- mid*samesign + lower*(1-samesign) #ifelse(samesign, mid, lower )
    flower <- fmid*samesign + flower*(1-samesign) #ifelse(samesign, fmid, flower )
    upper <- mid*(1-samesign) + upper*samesign #ifelse(!samesign, mid, upper )
    fupper <- fmid*(1-samesign) + fupper*samesign #ifelse(!samesign, fmid, fupper )

    if(n==numiter) print('Bisection hit maximum number of iterations')
  }

  return(list( root=mid, fmid=fmid, lower=lower, upper=upper,
        flower=flower, fupper=fupper, n=n ))
} 

