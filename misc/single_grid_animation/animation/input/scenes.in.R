library(dplyr)
library(tibble)

# Optionally label cities in each scene
cities <- data.frame(
  name = c("Hobart"),
  lon = c(147.33),
  lat = c(-42.88),
  colour = c("black"),
  cex = c(2.5)
)

scenes <- bind_rows(
  tibble(
    title = "Mw 8.4 earthquake",
    domains = list(NULL),
    ll_x = 143.0,
    ll_y = -53.0,
    ur_x = 178.5,
    ur_y = -35.0,
    dbox = list(0.8),
    osm_zoom = 6,
    dir = "se_australia",
    cities = NULL,
    plot_initial_conditions = FALSE
  )
)

scenes <- scenes %>%
  mutate(
    asp = 1 / cos((ll_y + ur_y) / 2 / 180 * pi),
    dx = ur_x - ll_x,
    dy = ur_y - ll_y,
    dydx = dy / dx
  )
