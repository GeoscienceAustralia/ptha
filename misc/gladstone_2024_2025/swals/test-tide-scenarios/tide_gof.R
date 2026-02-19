#' Plot predicted (vary tide) vs predicted (static tide)
#' at several locations
library(terra)

domain_id <- list(
    'south_trees' = 99,
    'rosslyn_bay' = 115,
    'lady_elliot' = 109
)
output_dir <- '../OUTPUTS/tides/'
scenario_names <- c('kermadectonga2', 'large_kermadec', 'solomon2', 'southamerica')
# var <- 'depth_as_max_stage_minus_elevation0'
var <- 'max_speed'
dir_write <- paste0('plot_tide_vary_vs_static/', var)

if (!dir.exists(dir_write)) {
    dir.create(dir_write, recursive=TRUE)
}

stats_log <- file(
    paste0(dir_write, '/stats.log'),
    open='w'
)   

for (scenario_name in scenario_names) {
    for (domain_name in names(domain_id)) {
        print(paste0('Processing ', domain_name, ' ', scenario_name))

        dir_static <- Sys.glob(paste0(output_dir, domain_name, '_', scenario_name, '-full-ambient_sea_level_*/*'))
        dir_vary <- Sys.glob(paste0(output_dir, 'vary_', scenario_name, '-full-ambient_sea_level_*/*'))
        if (length(dir_static) != 1) {
            print(dir_static)
            stop('Expected 1 static directory')
        }
        if (length(dir_vary) != 1) {
            print(dir_vary)
            stop('Expected 1 vary directory')
        }

        tif_static <- paste0(dir_static, '/', var, '_domain_', domain_id[[domain_name]], '.tif')
        tif_vary <- paste0(dir_vary, '/', var, '_domain_', domain_id[[domain_name]], '.tif')
        if (!file.exists(tif_static)) {
            stop(paste0('Expected static tif file at ', tif_static))
        }
        if (!file.exists(tif_vary)) {
            stop(paste0('Expected vary tif file at ', tif_vary))
        }

        static_rast <- rast(tif_static)
        vary_rast <- rast(tif_vary)

        # mask areas where elevation0 is less than -1 in dir_vary/elevation0_domain...
        mask <- rast(paste0(dir_vary, '/elevation0_domain_', domain_id[[domain_name]], '.tif'))
        mask <- mask < -1.0
        static_rast[mask] <- NA
        vary_rast[mask] <- NA

        static <- as.matrix(static_rast)
        vary <- as.matrix(vary_rast)

        pdf(
            paste0(dir_write, '/', domain_name, '_', scenario_name, '_vary_vs_static.pdf'),
            width = 12,
            height = 11
        )

        # get a goodness of fit measure R^2
        gof <- cor(static, vary, use='complete.obs')^2
        eps <- 1e-8
        v_on_s <- vary/(static)
        v_on_s <- v_on_s[!is.na(v_on_s)]
        v_on_s <- v_on_s[!is.infinite(v_on_s)]
        # get the 5th and 95th percentile
        gof2_5 <- quantile(v_on_s, 0.05)
        gof2_95 <- quantile(v_on_s, 0.95)
        rmse <- sqrt(mean((static - vary)^2, na.rm=TRUE))
        median <- median(v_on_s)

        # title
        main <- paste0(gsub('_', ' ', domain_name), ' ', gsub('_', ' ', scenario_name))
        subtitle <- paste0(
            'R^2 = ', round(gof, 3),
            '. RMSE = ', formatC(rmse, format='e', digits=1),
            '. vary/static: 5th percentile = ', round(gof2_5, 3),
            ', 95th percentile = ', round(gof2_95, 3),
            ', median = ', round(median, 3), '.'
        )

        # log stats to log file
        writeLines(main, stats_log)
        writeLines(subtitle, stats_log)

        # set up plot
        max_rast <- max(
            max(static, na.rm=TRUE),
            max(vary, na.rm = TRUE)
        )
        max_val <-  max(max_rast, na.rm=TRUE)
        plot(
            NA,
            NA,
            xlim=c(0, max_val),
            ylim=c(0, max_val),
            type='n',
            xlab='Static tide',
            ylab='Vary tide',
            main=main
        )
        title(sub=subtitle)

        # add x=y line up to max value
        lines(c(0, max_val), c(0, max_val), col='red')

        # plot pixel-wise each value all on same plot. black symbol dot
        points(static, vary, col='black', pch=20)

        dev.off()
    }
}
