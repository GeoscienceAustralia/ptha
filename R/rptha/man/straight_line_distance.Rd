% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometric_util.R
\name{straight_line_distance}
\alias{straight_line_distance}
\title{Find the 'straight line' distance of point p0 to p1, i.e. the Euclidean
distance between the points in 3D space. Both points are assumed to be on a sphere}
\usage{
straight_line_distance(p0, p1, r = 6378137)
}
\arguments{
\item{p0}{vector of lon/lat, or matrix with 1st column lon and 2nd column lat}

\item{p1}{vector of lon/lat, or matrix with 1st column lon and 2nd column lat}

\item{r}{radius of the earth}
}
\description{
This will be shorter than the great circle distance (but similar for small
differences between the points)
}
