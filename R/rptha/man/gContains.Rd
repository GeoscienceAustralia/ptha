% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alternatives_rgeos.R
\name{gContains}
\alias{gContains}
\title{limited replacement for rgeos::gContains using sf functionality}
\usage{
gContains(spgeom1, spgeom2, byid = FALSE, prepared = TRUE, returnDense = TRUE)
}
\description{
Like rgeos this assumes all coordinates are Cartesian. Suggest to use
sf::st_contains for non-legacy applications.
}
