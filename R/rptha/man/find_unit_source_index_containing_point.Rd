% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unit_sources.R
\name{find_unit_source_index_containing_point}
\alias{find_unit_source_index_containing_point}
\title{Find the unit source containing a given xy coordinate}
\usage{
find_unit_source_index_containing_point(
  point_xy,
  unit_source_geometry,
  unit_source_statistics
)
}
\arguments{
\item{point_xy}{numeric vector of length 2 giving c(lon, lat)}

\item{unit_source_geometry}{SpatialPolygonsDataFrame with unit source
geometry, e.g. output of \code{unit_source_grid_to_SpatialPolygonsDataFrame}}

\item{unit_source_statistics}{data.frame with unit source statistics, e.g.
output of \code{discretized_source_summary_statistics}}
}
\value{
integer giving the row index in unit_source_statistics of the unit
source containing point_xy
}
\description{
Find the unit source containing a given xy coordinate
}
\examples{

# Get source contours
puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
# Get downdip lines
puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
   layer='puysegur_downdip')
# Make discretized_source with 50km x 50km unit-sources (approximately)
puysegur_discretized_source = discretized_source_from_source_contours(
    source_shapefile=puysegur,
   desired_subfault_length=50,
   desired_subfault_width=50,
   downdip_lines=puysegur_downdip)

# Get geometry
puysegur_geometry = unit_source_grid_to_SpatialPolygonsDataFrame(
    puysegur_discretized_source$unit_source_grid)

# Get summary statistics
puysegur_summary_statistics_approx = discretized_source_approximate_summary_statistics(
   puysegur_discretized_source)

# Find the unit source containing this lon/lat point
pp = c(166.5, -45.9)

# Plot the situation
plot(puysegur_geometry, axes=TRUE, asp=1)
points(pp[1], pp[2], col='red', pch=19)
# Find the index of the unit source containing pp [i.e. index in the summary statistics]
unit_source_index = find_unit_source_index_containing_point(
    pp, puysegur_geometry, puysegur_summary_statistics_approx)
# Here are the summary statistics
target_unit_source = puysegur_summary_statistics_approx[unit_source_index,]

# Plotting
points(target_unit_source$lon_c, target_unit_source$lat_c, col='green', pch=19, cex=0.3)
geo_ind = which(target_unit_source$downdip_number == puysegur_geometry$downdip_number &
    target_unit_source$alongstrike_number == puysegur_geometry$alongstrike_number)
plot(puysegur_geometry[geo_ind,], border='green', add=TRUE)

}
