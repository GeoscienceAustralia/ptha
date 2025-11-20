library(dplyr)
library(tibble)


cities <- data.frame(
  name = c("Yeppoon", "Gladstone", "Agnes Waters", "Baffle Creek", "Lady Elliot\nIsland", "Heron Island", "Bundaberg"),
  lon = c(150.74369, 151.25518, 151.90428, 152.02027, 152.715520, 151.913020, 152.3489),
  lat = c(-23.13446, -23.9, -24.21121, -24.53456, -24.113016, -23.441990, -24.8661),
  colour = c("white", "white", "white", "white", "black", "black", "white"),
  cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5)
)

yeppoon_cities <- data.frame(
  name = c("Yeppoon", "Great Keppel\nIsland", "Rosslyn"),
  lon = c(150.79, 150.95262, 150.78),
  lat = c(-23.13065, -23.15, -23.16966),
  colour = c("black", "black", "white"),
  cex = c(3, 3, 3)
)

scenes <- bind_rows(
  tibble(
    title = "Gladstone Harbour",
    domains = list(c(71, 72, 92:105)),
    ll_x = 151.13,
    ll_y = -23.90,
    ur_x = 151.44,
    ur_y = -23.70,
    dbox = list(0.2),
    osm_zoom = 13,
    dir = "gladstone_harbour",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Capricorn Coast",
    domains = list(NULL),
    ll_x = 150.5,
    ll_y = -24.7,
    ur_x = 153,
    ur_y = -23.0,
    dbox = list(0.2),
    osm_zoom = 9,
    dir = "capricornia",
    cities = list(cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Pacific Ocean",
    domains = list(NULL),
    ll_x = 100,
    ll_y = -60,
    ur_x = 300,
    ur_y = 55,
    dbox = list(c(0, 179 - 300, 0, 0)),
    osm_zoom = 2,
    dir = "pacific",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Yeppoon and the\nGreat Keppel Islands",
    domains = list(c(82, 83, 114, 115, 116, 117, 118, 119, 123, 124, 125, 126, 127)),
    ll_x = 150.74,
    ll_y = -23.21,
    ur_x = 151.0,
    ur_y = -23.03,
    dbox = list(0.1),
    osm_zoom = 12,
    dir = "yeppoon_offshore",
    cities = list(yeppoon_cities),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Rosslyn Bay",
    domains = list(c(116, 117)),
    ll_x = 150.7736,
    ll_y = -23.173,
    ur_x = 150.81,
    ur_y = -23.15,
    dbox = list(0.05),
    osm_zoom = 16,
    dir = "rosslyn_bay",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Australia and\nNew Zealand",
    domains = list(NULL),
    ll_x = 110,
    ll_y = -52,
    ur_x = 180,
    ur_y = -8,
    dbox = list(c(30, 0, 0, 0)),
    osm_zoom = 4,
    dir = "sw_pacific",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "SW Pacific Ocean",
    domains = list(NULL),
    ll_x = 135,
    ll_y = -55,
    ur_x = 230,
    ur_y = 10,
    # use negative dbox to go back to anti-meridian
    dbox = list(c(10, 179 - 230, 0, 0)),
    osm_zoom = 4,
    dir = "sw_pacific_full",
    cities = list(NA),
    plot_initial_conditions = TRUE
  ),
  tibble(
    title = "Coral Sea",
    domains = list(NULL),
    ll_x = 140,
    ll_y = -27,
    ur_x = 179,
    ur_y = 0,
    dbox = list(1),
    osm_zoom = 5,
    dir = "coral_sea",
    cities = list(NA),
    plot_initial_conditions = TRUE
  ),
  tibble(
    title = "Boyne Island",
    domains = list(c(89, 90, 91, 94, 95, 96)),
    ll_x = 151.322,
    ll_y = -23.975,
    ur_x = 151.40,
    ur_y = -23.921,
    dbox = list(0.1),
    osm_zoom = 15,
    dir = "boyne_island",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Agnes Waters",
    domains = list(c(121, 122, 66)),
    ll_x = 151.83342,
    ll_y = -24.22,
    ur_x = 151.97554,
    ur_y = -24.13,
    dbox = list(0.1),
    osm_zoom = 14,
    dir = "agnes_waters",
    cities = list(NA),
    plot_initial_conditions = FALSE
  ),
  tibble(
    title = "Lady Elliot Island",
    domains = list(NA),
    ll_x = 152.67016,
    ll_y = -24.14651,
    ur_x = 152.77349,
    ur_y = -24.07051,
    dbox = list(0.001),
    osm_zoom = 11,
    dir = "lady_elliot",
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
