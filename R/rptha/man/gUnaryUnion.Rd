% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alternatives_rgeos.R
\name{gUnaryUnion}
\alias{gUnaryUnion}
\title{limited replacement for rgeos::gUnaryUnion using sf functionality}
\usage{
gUnaryUnion(spgeom, id = NULL)
}
\description{
Suggest to use sf::aggregate(..., dissolve=TRUE) and/or sf::st_union for non-legacy applications.
}
