% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kajiura_filter.R
\name{kajiura_filter}
\alias{kajiura_filter}
\title{Implementation of a Kajiura filter}
\usage{
kajiura_filter(
  xyDef,
  depth,
  grid_dx = max(depth)/2,
  grid_dy = max(depth)/2,
  edge_buffer_value = 0,
  edge_effect_correction_scale = 1.5,
  kajiuraGmax = 9,
  interpolator = "linear",
  interpolator_categories = function(xy) {
     xy[, 1] * 0
 },
  volume_change_error_threshold = 0.02,
  volume_change_MSL = 0,
  verbose = FALSE,
  fortran_inner_loop = TRUE
)
}
\arguments{
\item{xyDef}{3 column matrix with x,y, deformation. x and y can be unstructured.
The x,y coordinate system must be cartesian and in m}

\item{depth}{vector with the depth at each x,y point in xyDef (in m)}

\item{grid_dx}{Numeric (m). See grid_dy}

\item{grid_dy}{Numeric (m). To apply the filter, we regrid xyDef on a grid with
point spacing grid_dx, grid_dy, then smooth, then transform back from the grid
to our original xy points. This should ideally be a fraction of the depth in
areas of significant deformation.}

\item{edge_buffer_value}{Numeric. Outside the domain edges we assume this is
the value of xyDef[,3], when the filter is applied. Without further correction, this
would make edges tend towards edge_buffer_value}

\item{edge_effect_correction_scale}{To reduce edge effects, after
filtering, we compute the distance of every point to the edge of the domain
'd', and then return the solution \cr
WT = max( 1-d/(reference_depth*edge_effect_correction_scale), 0)**0.5 \cr
OUTPUT = WT*(ORIGINAL OUTPUT) + (1-WT)*FILTERED OUTPUT \cr}

\item{kajiuraGmax}{When empirically approximating kajiuraG, we fit it from
x=[0, kajiuraGmax]. Values above this are evaluated to zero. Be cautious about changing this.}

\item{interpolator}{'linear' or 'nearest'. Linear is better, but may be slow for
large point clouds. Don't use nearest unless grid_x, grid_y are 'small enough'}

\item{interpolator_categories}{Function of xy (matrix of coordinates)
returning a value for each point, with distinct values corresponding to
distinct 'groups' of points. Interpolation will be performed separately on
each group. It is passed as \code{category_function} to \code{interpolation_discontinuous}}

\item{volume_change_error_threshold}{If the difference in the positive or
negative or total volume before and after filtering, relative to the original
'absolute' volume, is more than this, then throw an error. This might not indicate a mistake,
but it does indicate a large change in the deformation, which is worth investigating.
Possible causes include strong depth variation, or regridding with coarse grid_dx, grid_dy.}

\item{volume_change_MSL}{The mean-sea-level used to define positive/negative deformation
for the volume check mentioned above.}

\item{verbose}{Print lots of information about the fit}

\item{fortran_inner_loop}{Use compiled fortran for the inner loop instead of R}
}
\value{
replacement version of xyDef, with smoothing applied to xyDef[,3]
}
\description{
Implement the filter similar to that of Glimsdal et al (2013, Eqns 10-11), based on Kajiura (1963). \cr
\deqn{ newDeformation(x,y) = depth(x,y)^{-2} \int \int \big{[}oldDeformation(x',y')G( \sqrt{(x'-x)^2+(y'-y)^2} / depth(x, y) )\big{]} dx' dy'}
This is a 2D generalisation of the cosh filter, justified for a 'temporally short'
earthquake with ocean governed by linear wave equations in constant depth
water. The definition of G is explained further in the help for \code{kajiura_g}\cr
Essentially: \cr
xyDef[,3] <-- convolution of the "old deformation" and "a smoothing kernel which integrates to 1"\cr
Numerically, we compute the new deformation as: \cr
\deqn{ numerator_{nm} = sum_{i} sum_{j} oldDeformation_{ij} * G( \sqrt{ (x_{nm} - x_{ij})^2 + (y_{nm} - y_{ij})^2} / depth_{nm} ) }
\deqn{ denominator_{nm} = sum_{i} sum_{j} G( \sqrt{ (x_{nm} - x_{ij})^2 + (y_{nm} - y_{ij})^2} / depth_{nm} ) }
\deqn{ newDeformation_{nm} = numerator_{nm} / denominator_{nm} }
where nm and ij denote pixel coordinates on a regular (cartesian) grid. So at any point
nm, the new deformation is a weighted average of nearby values, with weights
coming from Kajiuras function.\cr
This is slightly different to the approach used in Glimsdal et al (2013), but is
equivalent for constant depth. With non-constant depth the underlying theory is
not exactly valid, but should provide a reasonable approximation for slowly varying depths.\cr
We attempt to reduce edge effects by linearly weighting original and filtered values at edges,
since we cannot efficiently deal with edge effects in a better way. Therefore
it is best to have unimportant features around the edge of the input points. \cr
We allow xyDef to be unstructured, and start by gridding the results
on a grid with spacing approximately grid_dx,grid_dy.
The grid spacing is not exactly grid_dx,grid_dy, because it is forced to 
exactly divide the maximum input depth. \cr
For deformations with discontinuities, there can be artefacts due to
regridding, and it may be numerically beneficial to rotate the x,y input
coordinates so that the discontinuity is aligned with one of the coordinate
axes. (This is done in make_tsunami_unit_source).
}
