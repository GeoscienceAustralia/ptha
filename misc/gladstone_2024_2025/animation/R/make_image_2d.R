library(fields)

source("../../../../propagation/SWALS/plot.R")


#' Create a 2D image plot of water stage data with optional backdrop.
#'
#' This function generates a PNG image showing water level (stage) data for a
#' given time index, overlaid on an optional OpenStreetMap backdrop. It supports
#' adding city labels, a scale bar, and custom plot configurations such as
#' colour scales and dimensions.
#'
#' @param t Integer. Time index for the stage data to plot.
#' @param dir_save Character. Directory where the output PNG file will be saved.
#' @param backdrop Spatial object. Backdrop imagery (e.g., OpenStreetMap).
#' @param multidomain_dir Character. Directory containing multidomain model output files.
#' @param ll Numeric vector of length 2. Lower-left (lon, lat) of the plot extent.
#' @param ur Numeric vector of length 2. Upper-right (lon, lat) of the plot extent.
#' @param plot_config List. Plot configuration options including:
#'   \itemize{
#'     \item width, height: Dimensions of the output image in pixels.
#'     \item zlim: Numeric vector specifying colour scale limits.
#'     \item cols: Colour palette for water levels.
#'     \item cex.axis, legend.axis: Text scaling factors.
#'     \item legend_line: Positioning for legend text.
#'     \item add_backdrop: Logical, whether to include the backdrop imagery.
#'     \item add_scale_bar: Logical, whether to include a scale bar.
#'   }
#' @param cities Data frame or NA. Optional city labels with columns: lon, lat, name, cex, colour.
#' @param site_grids Character vector or NA. Optional paths to site-specific grid files.
#'
#' @details
#' The function sets up a plot with a colour scale for water levels, overlays the
#'  multidomain stage data, and optionally adds a backdrop image, city labels,
#' and a scale bar. Aspect ratio checks are performed to warn about potential
#' white space padding.
#'
#' @return Saves a PNG image to the specified directory. No object is returned.
#'
#' @examples
#' make_image_2d(
#'   t = 1,
#'   dir_save = "output/images",
#'   backdrop = backdrop,
#'   multidomain_dir = "model/output",
#'   ll = c(150.0, -35.0),
#'   ur = c(151.0, -34.0),
#'   plot_config = list(
#'     width = 800, height = 600, zlim = c(-10, 10),
#'     cols = ga_ocean_cols(n = 1024, alpha = 1.0), cex.axis = 2,
#'     legend.axis = 2, legend_line = 4, add_backdrop = TRUE,
#'     add_scale_bar = FALSE
#'   )
#' )
make_image_2d <- function(t, dir_save, backdrop, multidomain_dir, ll, ur, plot_config, cities = NA, site_grids = NA) {
  filename <- file.path(dir_save, sprintf("stage%03d.png", t))
  print(paste0("Working on ", filename))
  png(filename = filename, width = plot_config$width, height = plot_config$height)

  xlim <- c(ll[1], ur[1])
  ylim <- c(ll[2], ur[2])
  asp <- 1 / cos(mean(ylim) / 180 * pi)

  request_asp <- abs(diff(ylim)) * asp / abs(diff(xlim))
  plot_dims <- par("pin")
  plot_asp <- plot_dims[2] / plot_dims[1]
  if (plot_asp != request_asp) {
    msg <- paste0("Plot graph aspect ratio is ", plot_asp, ". Data aspect ratio is ", request_asp, ".")
    print(msg)
    warning(msg)
  }

  dry_depth <- 1e-3
  remove_cells_with_maxima_below <- 0.05

  # Configure plot (e.g. margins)
  par(plot_config$par)

  # set-up blank plotting area
  image.plot(
    matrix(0, nrow = 2, ncol = 2),
    asp = asp,
    xlim = xlim,
    ylim = ylim,
    zlim = plot_config$zlim,
    col = plot_config$cols,
    nlevel = length(plot_config$cols) + 1,
    legend.args = list(
      cex = plot_config$cex.axis,
      text = "water level (m)",
      side = 4,
      line = plot_config$legend_line
    ),
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i",
    legend.only = FALSE,
    cex.axis = plot_config$cex.axis,
    legend.cex = plot_config$legend.axis
  )

  # add the background imagery
  if (!is.null(backdrop)) {
    plot(
      backdrop,
      xlim = xlim,
      ylim = ylim,
      asp = asp,
      add = TRUE,
      xlab = "",
      ylab = ""
    )
  }

  # add the stage
  multidomain_image(
    multidomain_dir,
    variable = "stage",
    time_index = t,
    xlim = xlim,
    ylim = ylim,
    zlim = plot_config$zlim,
    col = plot_config$cols,
    add = TRUE,
    var_transform_function = NULL,
    NA_if_stage_not_above_elev = TRUE,
    NA_if_max_flux_is_zero = TRUE,
    use_fields = TRUE,
    clip_to_zlim = TRUE,
    buffer_is_priority_domain = TRUE,
    asp = asp,
    fields_axis_args = list(),
    dry_depth = dry_depth,
    nc_files = site_grids
  )

  # city labels
  if (!any(is.na(cities))) {
    text(cities$lon, cities$lat, cities$name, pos = 2, adj = 0, cex = cities$cex, col = cities$colour)
  }

  # scale bar
  if (plot_config$add_scale_bar) {
    d_km <- pointDistance(matrix(ll, ncol = 2), cbind(ur[1], ll[2]), lonlat = TRUE) / 2 / 1000
    d_km <- signif(d_km, 1)
    scalebar(d = d_km, type = "bar", divs = 2, lonlat = TRUE, below = "km", adj = c(0, 0.5))
  }

  dev.off()
}
