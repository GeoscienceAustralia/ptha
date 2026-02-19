# Script to plot the depth of the two models at a single pixel
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

load_depths <- function(coords, domain_id, sea_level, sea_level_name) {
    path_to_rasters <- Sys.glob(
        paste0(
            "../../../swals/OUTPUTS/v6/ptha18_tidal_check/sea_level_",
            sea_level_name,
            "/random_*/ptha18_random_scenarios_*_HS-full-ambient_sea_level_",
            sea_level,
            "/raster_output_files.tar"
            )
        )
    if (length(path_to_rasters) == 0) {
        stop("No rasters found for sea level ", sea_level)
    }
    # in between ptha18_random_scenarios_ and _HS-full-ambient_sea_level_
    this.scenario_names <- gsub(
        ".*ptha18_random_scenarios_(.*)HS-full-ambient_sea_level_.*", 
        "\\1", 
        path_to_rasters
    )
    vsi_path_to_rasters <- paste0("/vsitar/", path_to_rasters, "/max_stage_domain_", domain_id, ".tif")
    rasters <- rast(vsi_path_to_rasters)
    this.depth <- terra::extract(rasters, coords$value)
    this.sea_level <- rep(sea_level, length(this.scenario_names))
    this.depth <- as.vector(unlist(unname(this.depth)))

    this.df <- data.frame(
        sea_level=this.sea_level,
        scenario=this.scenario_names,
        depth=this.depth
    )
    return(this.df)
}

# Plot the depth of the two models at a single pixel function
plot_depths <- function(depths, coords, residuals=FALSE, offset=2.556) {
    # spread sea levels as columns
    depths_spread <- spread(depths, sea_level, depth)
    col_names <- colnames(depths_spread)

    if (residuals) {
        p <- ggplot(depths_spread, aes(x=depths_spread[, col_names[2]], y=depths_spread[, col_names[3]])) +
            geom_point(aes(y=depths_spread[, col_names[3]] - depths_spread[, col_names[2]]-offset)) +
            geom_abline(intercept = 0, slope = 0) +
            labs(
                x = paste0("Depth (m) at sea level ", col_names[2], " m"),
                y = paste0("Residuals (m). Sea level ", col_names[3], " - Sea level ", col_names[2])) +
            theme_minimal()
    } else {
        # Plot the depth-depth scatter plot
        p <- ggplot(depths_spread, aes(
                x=depths_spread[, col_names[2]],
                y=depths_spread[, col_names[3]] - offset)
            ) +
            geom_point() +
            geom_abline(intercept = 0, slope = 1) +
            labs(
                x = "Max-stage (m) using tidal adjustment technique",
                y = paste0("Max-stage (m) using static tide of ", col_names[3]), " m") +
            theme_minimal()
    }
    return(p)
}


# inputs for the function
setClass(
  "Params",
  slots = list(
    coords = "list",
    domain_id = "character",
    sea_levels = "character"
  )
)
params <- list(
    # new(
    #     "Params",
    #     coords = list(name="lady_elliot_onshore", value=cbind(152.7172, -24.1138)),
    #     domain_id = "109",
    #     sea_levels = c("1.459"="1.459", "vary"="0")
    # ),
    # new(
    #     "Params",
    #     coords = list(name="lady_elliot_offshore", value=cbind(152.7191, -24.0975)),
    #     domain_id = "109",
    #     sea_levels = c("1.459"="1.459", "vary"="0")
    # ),
    # new(
    #     "Params",
    #     coords = list(name="south_trees_channel", value=cbind(151.3200, -23.8500)),
    #     domain_id = "99",
    #     sea_levels = c("2.145"="2.145", "vary"="0")
    # ),
    # new(
    #     "Params",
    #     coords = list(name="south_trees_mangrove", value=cbind(151.2785, -23.8566)),
    #     domain_id = "99",
    #     sea_levels = c("2.145"="2.145", "vary"="0")
    # ),
    new(
        "Params",
        coords = list(name="rosslyn_bay", value=cbind(150.7894, -23.1610)),
        domain_id = "115",
        sea_levels = c("2.556"="2.556", "vary"="0")
    )
    # new(
    #     "Params",
    #     coords = list(name="boyne_river", value=cbind(151.3554,-23.9437)),
    #     domain_id = "90",
    #     sea_levels = c("2.145"="2.145", "vary"="0")
    # )
)

for (param_set in params) {
    coords <- param_set@coords
    domain_id <- param_set@domain_id
    sea_levels <- param_set@sea_levels
    print(coords$name)

    depths <- load_depths(coords, domain_id, sea_levels[1], names(sea_levels)[1])
    depths <- rbind(depths, load_depths(coords, domain_id, sea_levels[2], names(sea_levels)[2]))
    sea_level_names <- names(sea_levels)

    # Plot the depth-depth scatter plot
    plot <- plot_depths(depths, coords)
    dir_save <- "stage-stage_plot"
    if (!dir.exists(dir_save)) {
        dir.create(dir_save)
    }
    save_path <- paste0(dir_save, "/", sea_level_names[1], "_", sea_level_names[2], "_", coords$name, "_", coords$value[1], "_", coords$value[2], ".png")
    ggsave(save_path, plot, width=6, height=6)

    # plot the residuals
    plot <- plot_residuals <- plot_depths(depths, coords, residuals=TRUE)
    save_path <- paste0(dir_save, "/", sea_level_names[1], "_", sea_level_names[2], "_", coords$name, "_", coords$value[1], "_", coords$value[2], "_residuals.png")
    ggsave(save_path, plot, width=6, height=6)
}
