library(dplyr)
library(tibble)


cities <- data.frame(
  name = c("Onslow", "Karratha", "Exmouth", "Barrow\nIsland", "Coral Bay", "Geraldton", "Carnarvon"),
  lon = c(115.12, 116.85, 114.17, 115.46, 113.76, 114.58, 113.65),
  lat = c(-21.63, -20.74, -21.95, -20.79, -23.13, -28.80, -24.89),
  colour = c("black", "black", "black", "black", "black", "black", "black"),
  cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5)
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
