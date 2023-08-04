#
# Basic plot of the pressure data, and the high-frequency component (removing periods > 2h)
#
source('parse_gauge_data.R')
source('global_variables.R')
source('three_panel_plot.R')

#
# MSL Pressure
#
pressure_metadata = read.csv(MSLP_METADATA_TABLE_FILE)

all_pressure_files = pressure_metadata$original_data_file

# Read all the pressure data
all_pressure_data = lapply(all_pressure_files, function(x) try(read_pressure_gauge_any(x)))
names(all_pressure_data) = gsub('.csv', '', basename(all_pressure_files))

for(i in 1:length(all_pressure_data)){
    # Remove any data outside a specified window. This matters for the tide gauge data (for speed and for artefact removal) 
    # -- for the MSLP data here it is just applied for consistency. 
    start_limit = julian(START_TIME_LIMIT_GAUGE_DATA)
    end_limit   = julian(END_TIME_LIMIT_GAUGE_DATA)
    k = which((all_pressure_data[[i]]$juliant < start_limit) | (all_pressure_data[[i]]$juliant > end_limit))
    if(length(k) > 0){
        all_pressure_data[[i]] = all_pressure_data[[i]][-k,]
    }
}


# Plot every station (if possible)
pdf(paste0(OUTPUT_GRAPHICS_DIR, '/highpass_filtered_plot_pressure_gauges.pdf'), width=10, height=8)
for(i in 1:length(all_pressure_data)){
    if(is(all_pressure_data[[i]], 'try-error')) next 

    # Make a plot
    three_panel_plot(all_pressure_data[[i]], names(all_pressure_data)[i], add_highpass_filtered_series=TRUE, 
        filter_threshold=1/(2*3600), site_lon = pressure_metadata$lon[i], site_lat = pressure_metadata$lat[i])

    # Add filtered pressures to the outputs
    data_t = as.numeric(all_pressure_data[[i]]$juliant - all_pressure_data[[i]]$juliant[1])*3600*24
    filt = try(spectral_highpass_filter(data_t, all_pressure_data[[i]]$pressure, interp_dt = 15, cutoff_frequency=1/(2*3600)))
    if(!is(filt, 'try-error')){
        highfreq_pressure = filt$highfreq
        # If the data is NA, ensure the filtered data is NA too
        k = which(is.na(all_pressure_data$pressure))
        if(length(k) > 0) highfreq_pressure[k] = NA
        all_pressure_data[[i]]$resid2h = highfreq_pressure
    }
}
dev.off()

# For convenience later
#saveRDS(all_pressure_data, paste0(OUTPUT_DIR, 'all_pressure_data.RDS'))


#
# Export all gauges to csv
#
dir.create(OUTPUT_MSLP_DIR, showWarnings=FALSE)
for(i in 1:length(all_pressure_data)){
    output_file = pressure_metadata$postprocessed_file[i]
    file_header = c(
        '#',
        paste0('# Name: ', pressure_metadata$station[i]), 
        paste0('# Location (LON/LAT): ', pressure_metadata$lon[i], ', ', pressure_metadata$lat[i]), 
        paste0('# Ground height (m above sea-level, nominal for DES stations): ', pressure_metadata$height_m[i]),
        paste0('# Barometer height (m above sea-level, nominal for DES stations): ', pressure_metadata$height_instrument_m[i]),
        paste0('# Original data (relative to scripts_to_postprocess directory): ', pressure_metadata$original_data_file[i]), 
        '#',
        '# Column variable definitions',
        '#   time: Date and time with UTC timezone ', 
        '#   juliant: Time in days since 1970-01-01 00:00:00 UTC', 
        '#   pressure: mean sea level pressure (hectopascals).',
        '#   quality_flag: A quality flag provided in the original Bureau of Meteorology data, with interpretation listed below.',
        '#       Y - quality controlled and acceptable',
        '#       N - not quality controlled',
        '#       S - quality controlled and considered suspect',
        '#       W - quality controlled and considered wrong',
        '#       I - quality controlled and inconsistent with other known information',
        '#       U - element not observed',
        '#       R - element removed',
        '#       F - no quality control performed',
        '#       NA - (for DES data) no quality control flag was provided in the original data',
        '#   resid2h: High-frequency component of pressure, retaining wave periods of 2 hours or less (hectopascals).',
        '#')
    cat(paste0(file_header,'\n'), file=output_file, sep="", fill=FALSE)

    output_data = data.frame(time=all_pressure_data[[i]]$time, juliant=all_pressure_data[[i]]$juliant,
        pressure=all_pressure_data[[i]]$pressure, quality_flag = all_pressure_data[[i]]$quality_flag, 
        resid2h=all_pressure_data[[i]]$resid2h)

    write.table(output_data, file=output_file, sep=",", row.names=FALSE, quote=FALSE, append=TRUE)
}

