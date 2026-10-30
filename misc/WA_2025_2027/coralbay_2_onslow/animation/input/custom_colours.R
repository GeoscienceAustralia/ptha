# GA Colourscheme with ocean in the middle
# adapted from Julie Silec's NSW work in 2021
ga_ocean_cols <- function(n = 1024, alpha = 0.7) {
  cols <- c(
    "#b33b3b", "#c26262", "#ca6c37", "#e38e5c", "#e7af91", "#e9dcae",
    "#e5ecee",
    "#6d9d91", "#3c7e6e", "#98c5cf", "#337ca1", "#005a6f", "#082e41"
  )
  # Interpolate with semi-transparent treatment
  cols <- colorRampPaletteAlpha(
    addalpha(rev(cols), alpha = alpha),
    n = 1024
  )
}

addalpha <- function(colors, alpha = 1.0) {
  r <- col2rgb(colors, alpha = TRUE)
  # Apply alpha
  r[4, ] <- alpha * 255
  r <- r / 255.0
  return(rgb(r[1, ], r[2, ], r[3, ], r[4, ]))
}

# colorRampPaletteAlpha()
# From
# https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
colorRampPaletteAlpha <- function(colors, n = 32, interpolate = "linear") {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate = interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha = TRUE)[4, ]
  # Interpolate
  if (interpolate == "linear") {
    l <- approx(a, n = n)
  } else {
    l <- spline(a, n = n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y / 255.0)
  return(cr)
}
