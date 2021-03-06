% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wave_dispersion_relation.R
\name{bisection}
\alias{bisection}
\title{Function for vectorized 1d root finding with bisection.}
\usage{
bisection(
  f,
  lower,
  upper,
  ...,
  numiter = 20000,
  tolerance = .Machine$double.eps^0.5
)
}
\arguments{
\item{f}{the function}

\item{lower}{a vector of lower bounds for x}

\item{upper}{a vector of upper bounds for x}

\item{...}{further arguments to f}

\item{numiter}{maximum allowed number of iterations}

\item{tolerance}{= Allowed error in the root}
}
\value{
a vector x with f(x) = 0 to within the tolerance
}
\description{
Solves f(x) = 0 where f(x) is a vector of the same length as x,
and f(1) is only affected by x(1), f(2) is only affected by
x(2), etc. In R, to do this efficiently we need vectorization, hence the
current function
}
\examples{
#
# Here is a function satisfying the assumptions described above
#
f<-function(x) x^2 - c(2,3,4)

solution = bisection(f, lower=c(0,0,0), upper=c(10,10,10))
stopifnot(all(abs(solution$root - sqrt(c(2,3,4))) < 1.0e-06))

}
