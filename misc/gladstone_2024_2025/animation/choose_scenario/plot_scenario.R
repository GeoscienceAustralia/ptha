#' Quick plot of globbed scenarios' max_stage.
#' Works with tarred files too
library(terra)


md_dirs <- Sys.glob("../swals/OUTPUTS/ptha/sea_level_vary/random_kermadectonga2/ptha18_random_scenarios_kermadectonga2_row_*_Mw_94_HS-full-ambient_sea_level_0")

for (md_dir in md_dirs) {
    tif <- paste("/vsitar", md_dir, "raster_output_files.tar", "max_stage_domain_98.tif", sep = "/")

    print(tif)
    r <- rast(tif)

    fn <- paste0(basename(md_dir), ".png")
    filename <- file.path("plot_max_stage", "domain_98", fn)
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    png(filename = filename)
    plot(r, range = c(0, 1))
    dev.off()
}
