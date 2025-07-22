#' Compute the fraction of sea levels that exceed a given threshold (HAT)
#' 
#' Also plot the time series of sea levels and a violin plot of the data
#' See the locations.in.R file for the locations and thresholds
library(dplyr)
library(ggplot2)


source("read_utils.R")

#' Plot the fraction of sea levels that exceed a given threshold for each year
plot_fraction_exceedance <- function(data, threshold, filename){
    data$year <- as.numeric(format(data$date, "%Y"))
    data$value <- as.numeric(data$value)
    print(summary(data))
    fraction <- data %>%
        group_by(year) %>%
        summarise(fraction = sum(value > threshold) / n() * 100)

    ggplot(fraction, aes(x=year, y=fraction)) +
        geom_line() +
        geom_point() +
        labs(title=paste("Fraction of sea levels exceeding modelled sea level"),
             x="Year",
             y="Fraction %") +
        theme_minimal()
    ggsave(filename)
}

#' Same as above but one line for each location
plot_all_fraction_exceedance <- function(data, filename=NA){
    data$year <- as.numeric(format(data$date, "%Y"))
    print(summary(data))

    fraction_tidal_plane <- data %>%
        group_by(year, name) %>%
        summarise(fraction = sum(value > threshold_tidal_plane_LAT) / n() * 100)
    fraction_tidal_plane$threshold_name <- "MSQ HAT"
    fraction_TPXO9 <- data %>%
        group_by(year, name) %>%
        summarise(fraction = sum(value > threshold_TPXO9_LAT) / n() * 100)
    fraction_TPXO9$threshold_name <- "TPXO9 HAT"
    fraction <- rbind(fraction_TPXO9, fraction_tidal_plane)

    ggplot(fraction, aes(x=year, y=fraction, color=name)) +
        geom_line() +
        # more space between the facets
        facet_wrap(~threshold_name) +
        geom_point() +
        labs(title=paste("Fraction of sea levels exceeding modelled sea level"),
             x="Year",
             y="Fraction %") +
        theme_minimal() +
        guides(fill = guide_legend(ncol = 3)) +
        theme(
            legend.position="top",
            legend.title = element_blank(),
            legend.text = element_text(size=14),
            axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            plot.title = element_text(size=24),
            panel.margin = unit(5, "lines"),
            strip.text = element_text(size=18),
        )
    if(!is.na(filename)) {
        ggsave(filename, width=10, height=5)
    }
}

#' Plot a time series of the data with a threshold line, group by name
plot_time_series <- function(data, location, filename){
    # use points not lines cause the data is too dense.
    p <- ggplot(data, aes(x=date, y=value)) +
        geom_point(size=0.001) +
        # plot both threshold lines
        geom_hline(yintercept=location$threshold_tidal_plane_LAT, linetype="dotted", color="red") +
        geom_hline(yintercept=location$threshold_TPXO9_LAT, linetype="dotted", color="blue") +
        labs(title="Time series of sea levels",
             x="Date",
             y="Sea level (m LAT)") +
        theme_minimal() +
        theme(legend.position="none") +
        geom_hline(yintercept=location$AHD_in_LAT, linetype="dashed", color="black")
        # annotate("text", x=3, y=location$AHD_in_LAT, label=paste("AHD"), vjust=-1, color="black")
    ggsave(filename, plot=p, width=10, height=10)
}

#' Violin plot of the data for each year with a threshold line
plot_violin <- function(data, location, datum, filename){

    stopifnot(datum %in% c("AHD", "LAT", "MSL", "MODEL"))
    alt_datum = 0
    if (datum == "AHD") {
        data$value <- data$value - location$AHD_in_LAT
        alt_datum <- location$msl_in_LAT - location$AHD_in_LAT
        y_threshold_tidal_plane <- location$threshold_tidal_plane_LAT - location$AHD_in_LAT
        y_threshold_TPXO9 <- location$threshold_TPXO9_LAT - location$AHD_in_LAT
    } else if (datum == "MSL") {
        data$value <- data$value - location$msl_in_LAT
        y_threshold_tidal_plane <- location$threshold_tidal_plane_LAT - location$msl_in_LAT
        y_threshold_TPXO9 <- location$threshold_TPXO9_LAT - location$msl_in_LAT
        alt_datum <- location$AHD_in_LAT - location$msl_in_LAT
    } else if (datum == "MODEL") {
        # as for AHD case, but TPXO9 values are left in MSL
        data$value <- data$value - location$AHD_in_LAT
        alt_datum <- location$msl_in_LAT - location$AHD_in_LAT
        y_threshold_tidal_plane <- location$threshold_tidal_plane_LAT - location$AHD_in_LAT
        y_threshold_TPXO9 <- location$threshold_TPXO9_LAT - location$msl_in_LAT
    }
    data$year <- as.numeric(format(data$date, "%Y"))

    # Filter to data above the minimum threshold
    threshold_in_dataum <- min(location$threshold_tidal_plane_LAT, location$threshold_TPXO9_LAT, na.rm=TRUE)
    data_top <- data %>% filter(value > threshold_in_dataum)

    if(all(is.null(location$bounds_AHD))) {
        location$bounds_AHD <- c(-Inf, Inf)
    }

    datum_label <- datum
    if (datum_label == "MODEL") datum_label <- "AHD"

    p <- ggplot(data, aes(x=factor(year), y=value)) +
        geom_violin(bounds = location$bounds_AHD) +
        labs(title=paste("Sea level distribution by year at", location$name),
             x="Year",
             y=paste0("Sea level (m ", datum_label, ")")) +
        geom_jitter(data=data_top, aes(x=factor(year), y=value), width=0.2, height=0, color="grey", alpha=0.5, size=0.1) +
        ylim(-3.3, 3.3) +
        # white (not transparent) background
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        # plot mean and median and color them for the legend
        stat_summary(aes(color="Mean"), fun=mean, geom="point", shape=15, size=3, color="black") +
        stat_summary(aes(color="Median"), fun=median, geom="point", shape=17, size=3, color="red") +
        # add a legend
        theme(legend.position="top") +
        scale_color_manual(values=c("Mean"="black", "Median"="red")) +
        # plot both threshold lines
        geom_hline(yintercept=y_threshold_tidal_plane, linetype="dotted", color="red") +
        annotate("text", x=1, y=y_threshold_tidal_plane, label=paste("MSQ HAT"), vjust=-0.8, hjust=0, color="red", size=7) +
        geom_hline(yintercept=y_threshold_TPXO9, linetype="dotted", color="blue") +
        annotate("text", x=1, y=y_threshold_TPXO9, label=paste("TPXO9 HAT"), vjust=+1.5, hjust=0, color="blue", size=7) +
        theme(
            axis.title.x = element_text(size=24),
            axis.title.y = element_text(size=24),
            axis.text.x = element_text(size=20),
            axis.text.y = element_text(size=24),
            plot.title = element_text(size=24),
        )
        if (datum == "MSL") {
            p <- p + geom_hline(yintercept=0, linetype="dashed", color="black") +
            geom_hline(yintercept=alt_datum, linetype="dashed", color="blue") +
            annotate("text", x=1, y=alt_datum, label=paste("AHD"), vjust=-1, hjust=0, color="blue", size=7)
        } else if (datum %in% c("AHD", "MODEL")) {
            # plot MSL line
            p <- p + geom_hline(yintercept=0, linetype="dashed", color="black") +
            geom_hline(yintercept=alt_datum, linetype="dashed", color="blue") +
            annotate("text", x=1, y=alt_datum, label=paste("MSQ MSL"), vjust=-1, hjust=0, color="blue", size=7)
        }
    ggsave(filename, plot=p, width=10, height=8)
}

#' Main function to run the analysis
#' 
#' locations: a list of named lists, each containing the data and metadata files for a location
#' filename_all: the name of the file to save the plot of all locations
main <- function(locations, filename_all="all_locations_fraction_exceedance.png", dir_write=".") {
    data <- data.frame(date=as.POSIXct(character()), value=numeric(), name=character(), thresdhold=numeric())

    dir_write = "results_MODEL"
    dir.create(dir_write)

    for(location in locations){
        # read metadata
        metadata <- read_metadata(location$metadata_file)
        location$AHD_in_LAT <- as.numeric(metadata['Australian Height Datum (AHD) level (metres) relative to station datum',])
        if(is.na(location$AHD_in_LAT) || location$AHD_in_LAT == 0.0) {
            stop("location$AHD_in_LAT is NA or 0.0")
        }

        # read data
        location_data <- read_data(location$files)

        # convert threshold to LAT
        if (is.na(location$msl_in_LAT)) {
            location$msl_in_LAT <- location$AHD_in_LAT
        }
        stopifnot(!is.na(location$msl_in_LAT))
        location$threshold_tidal_plane_LAT <- location$threshold_tidal_plane_MSL + location$msl_in_LAT
        location$threshold_TPXO9_LAT <- location$threshold_TPXO9_MSL + location$msl_in_LAT

        # 1 plot time series for this location
        filename <- paste0(dir_write, "/", location$file_write_name_stem, "_tidal_planes_time_series.png")
        plot_time_series(location_data, location_data$threshold_tidal_plane_MSL, location$AHD_in_LAT, filename=filename)
        plot_time_series(location_data, location_data, filename=filename)

        # 2 plot HAT exceedence for this location
        filename <- paste0(dir_write, "/", location$file_write_name_stem, "_fraction_exceedance.png")
        plot_fraction_exceedance(location_data, location, filename=filename)

        # 3 plot violin plot for this location
        filename <- paste0(dir_write, "/", location$file_write_name_stem, "_violin.png")
        print(location)
        plot_violin(location_data, location, datum="MODEL", filename=filename)

        # add the names, both thresholds and append
        location_data$name <- location$name
        location_data$threshold_tidal_plane_LAT <- location$threshold_tidal_plane_LAT
        location_data$threshold_TPXO9_LAT <- location$threshold_TPXO9_LAT
        data <- rbind(data, location_data)
    }

    # plot all locations
    filename_all <- paste(dir_write, "locations_fraction_exceedance.png", sep="/")
    plot_all_fraction_exceedance(data, filename_all)
}

source("locations.in.R")
main(locations)
