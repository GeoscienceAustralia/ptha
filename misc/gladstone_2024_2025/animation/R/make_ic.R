library(terra)
library(magick)


# make pngs for scenario initial conditions
make_ic_png <- function(k, plot_config, multidomain_dir, scene, dir_ic, ic_files) {
    dir.create(dir_ic, showWarnings = FALSE)

    xlim <- c(scene$ll_x, scene$ur_x)
    ylim <- c(scene$ll_y, scene$ur_y)

    path_ic <- ic_files[[k]]
    tmp_file <- file.path(dir_ic, sprintf("tmp_ic_%03d.png", k))
    print(paste0("Working on ", tmp_file))
    png(filename = tmp_file, width = plot_config$width, height = plot_config$height)

    multidomain_image(
        multidomain_dir,
        variable = "elevation0",
        time_index = 1,
        xlim = xlim,
        ylim = ylim,
        zlim = c(-5000, 1000),
        cols = gray.colors(256),
        add = FALSE,
        var_transform_function = NULL,
        NA_if_stage_not_above_elev = TRUE,
        NA_if_max_flux_is_zero = TRUE,
        use_fields = FALSE,
        clip_to_zlim = TRUE,
        buffer_is_priority_domain = TRUE,
        asp = 1 / cos(mean(ylim) / 180 * pi),
        fields_axis_args = list()
    )

    ic <- rast(path_ic)
    extent <- ext(c(scene$ll_x, scene$ur_x, scene$ll_y, scene$ur_y))
    ic <- extend(ic, extent, fill = 0)

    # Remove initially wet areas that are never really affected by the tsunami
    is_small <- abs(ic) < 0.05
    ic[is_small] <- NA

    # Clip to zlim
    ic <- clamp(ic, lower = plot_config$zlim[1], upper = plot_config$zlim[2])

    # Append the stage
    base::image.plot(
        ic,
        xlim = xlim,
        ylim = ylim,
        zlim = plot_config$zlim,
        col = plot_config$cols,
        add = TRUE,
        useRaster = TRUE,
        cex.axis = 2,
        legend = FALSE
    )

    # clip to plotting region
    par(xpd = FALSE)

    dev.off()
}

combine_ic_frames <- function(tmp_file) {
    # header_file <- file.path(dirname(tmp_file), "header.png")
    # if (!file.exists(header_file)) {
    header <- image_blank(1920, 130, "#00718B") %>%
        image_annotate(
            "Potential scenarios in the SW Pacific",
            font = "C059",
            size = 48,
            color = "white",
            gravity = "center"
        )
    #     image_write(header, header_file)
    # }

    # header <- image_read(header_file)

    side_panel <- image_blank(620, 950, color = "white")
    main <- image_read(tmp_file)

    body <- image_append(c(main, side_panel))
    full <- image_append(c(header, body), stack = TRUE)

    name <- gsub("tmp_", "", basename(tmp_file))
    image_write(full, file.path(dir_ic, name))

    file.remove(tmp_file)
}
