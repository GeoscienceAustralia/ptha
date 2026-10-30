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
    title = "Indonesia Region",
    domains = list(NULL),
    ll_x = 85,
    ll_y = -23,
    ur_x = 128,
    ur_y = 10,
    dbox = list(0.2),
    osm_zoom = 5,
    dir = "indonesia",
    cities = list(NA),
    plot_initial_conditions = TRUE
  ),
  tibble(
    title = "North West Australia",
    domains = list(NULL),
    ll_x = 112,
    ll_y = -23.5,
    ur_x = 117.5,
    ur_y = -19.5,
    dbox = list(0),
    osm_zoom = 7,
    dir = "nw_australia",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Mid West Gascoyne",
    domains = list(NULL),
    ll_x = 107,
    ll_y = -30.5,
    ur_x = 116.2,
    ur_y = -24.0,
    dbox = list(0.5),
    osm_zoom = 8,
    dir = "midwest_gascoyne",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Onslow",
    domains = list(c(68, 69, 73, 74, 117, 118, 119, 120, 121, 122, 123, 124)),
    ll_x = 114.95,
    ll_y = -21.71,
    ur_x = 115.27,
    ur_y = -21.53,
    dbox = list(0.1),
    osm_zoom = 12,
    dir = "onslow",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Cape Range",
    domains = list(c(26, 27, seq(56, 73), seq(105, 116), seq(125, 133))),
    ll_x = 113.36,
    ll_y = -22.36,
    ur_x = 114.66,
    ur_y = -21.58,
    dbox = list(0.1),
    osm_zoom = 9,
    dir = "cape_range",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Coral Bay",
    domains = list(c(53, 147, 148, 149)),
    ll_x = 113.673,
    ll_y = -23.165,
    ur_x = 113.836,
    ur_y = -23.067,
    dbox = list(0.1),
    osm_zoom = 14,
    dir = "coral_bay",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Indian Ocean",
    domains = list(NULL),
    ll_x = 35,
    ll_y = -60,
    ur_x = 155,
    ur_y = 30,
    dbox = list(1),
    osm_zoom = 3,
    dir = "indian_ocean",
    cities = list(NA),
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
