#
# Plotting, cleaning and de-tiding the tide-gauge data
#
source('parse_gauge_data.R')
source('global_variables.R')
source('three_panel_plot.R')
source('detiding.R')

# Read the gauge metadata
gauge_metadata = read.csv(TIDEGAUGE_METADATA_TABLE_FILE)

# Read the gauge time-series
all_gauges = lapply(gauge_metadata$original_data_file, function(x) try(read_tide_gauge_any(x)))
names(all_gauges) = gsub('.csv', '', basename(gauge_metadata$original_data_file), fixed=TRUE)

for(i in 1:length(all_gauges)){
    # For speed, remove any data outside a specified window
    # The start limit is also convenient in removing some artefacts in some BOM port data.
    start_limit = julian(START_TIME_LIMIT_GAUGE_DATA)
    end_limit   = julian(END_TIME_LIMIT_GAUGE_DATA)
    k = which((all_gauges[[i]]$juliant < start_limit) | (all_gauges[[i]]$juliant > end_limit))
    if(length(k) > 0){
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }
}

#
# Make a simple plot with the data, and a 'spectrally detided' time-series, around the time
# of the tsunami
#
pdf(paste0(OUTPUT_GRAPHICS_DIR, '/highpass_filtered_plot_tide_gauges_before_editing_out_artefacts.pdf'), width=10, height=8)
for(i in 1:length(all_gauges)){
    if(is(all_gauges[[i]], 'try-error')) next
    three_panel_plot(all_gauges[[i]], names(all_gauges)[i], add_highpass_filtered_series=TRUE,
        filter_threshold=1/(3*3600), site_lon = gauge_metadata$lon[i], site_lat = gauge_metadata$lat[i])

}
dev.off()

# 
# Cleaning and de-tiding the tide-gauge data.
#
#   Cleaning involves removing some spurious records
#
#   For the detiding, we try 2 approaches at all sites.
#   1) combining TPX072 with spectral removal of low-frequencies. 
#   2) spectral removal of low-frequencies. 
#   Note tidal predictions may not always work well (e.g. local bathymetric
#   effects in some estuarine or shallow coastal sites, transient weather, etc).
#
all_gauges_detided = vector(mode='list', length=length(all_gauges))
names(all_gauges_detided) = names(all_gauges)
for(i in 1:length(all_gauges)){
    print(c('Detiding gauge ', i))

    #
    # CLEANING OF TIDE GAUGES
    # A few gauges have "obviously problematic" entries, which we
    # remove in a case-by-case way. In earlier iterations there were some other
    # gauges that have since been removed (due to double-ups with other sites).
    #
    #if(names(all_gauges_detided)[i] == "BOMQC_thevenard"){
    #    # There is an obviously spurious depression here, not present in the neighbouring IOC site
    #    k = which((all_gauges[[i]]$juliant > 19006.8) & (all_gauges[[i]]$juliant < 19007.2))
    #    all_gauges[[i]]$stage[k] = NA
    #}else if(names(all_gauges_detided)[i] == "BOMQC_esperance"){
    #    # Spurious constant patch here, remove it.
    #    t0 = strptime('2022-01-07 22:41:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    #    t1 = strptime('2022-01-11 23:11:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    #    k = which((all_gauges[[i]]$time > t0) & (all_gauges[[i]]$time < t1))
    #    all_gauges[[i]]$stage[k] = NA
    #}else if(names(all_gauges_detided)[i] == "BOMPorts_ardrossan"){
    if(names(all_gauges_detided)[i] == "BOMPorts_ardrossan"){
        # This one is very noisy, only take a limited part
        t0 = strptime('2022-01-12 01:12:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        t1 = strptime('2022-01-17 08:30:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        k = which((all_gauges[[i]]$time < t0) | (all_gauges[[i]]$time > t1))
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }else if(names(all_gauges_detided)[i] == "Birkdale_2022-01-01T00_00-2022-01-31T23_59_tidedata"){
        # There are some gaps -- take only the interior 
        t0 = strptime('2022-01-08 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        t1 = strptime('2022-01-19 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        k = which((all_gauges[[i]]$time < t0) | (all_gauges[[i]]$time > t1))
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }else if(names(all_gauges_detided)[i] == "Brisbane River (Pinkenba Wharf)_2022-01-01T00_00-2022-01-31T23_59_tidedata"){
        # There are some negative spikes (? zeros ?)
        k = which(all_gauges[[i]]$stage == 0)
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }else if(names(all_gauges_detided)[i] == "ferg"){
        # There are some positive spikes
        k = which(all_gauges[[i]]$stage > 4)
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }else if(names(all_gauges_detided)[i] == "groo"){
        # Noisy after 28th
        t0 = strptime('2022-01-28 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        k = which((all_gauges[[i]]$time > t0))
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }else if(names(all_gauges_detided)[i] == "20220114-20220120_Newcastle Primary - (East) Tide"){
        # Noise begins sometime on 17th
        t0 = strptime('2022-01-17 10:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
        k = which((all_gauges[[i]]$time > t0))
        all_gauges[[i]] = all_gauges[[i]][-k,]
    }

    #
    # Adjustment of gauge coordinate for tidal prediction. A few gauges are at
    # inland locations not recognized as water by TPXO. For those cases, move
    # the gauge coordinate to somewhere recognized as water, and extract the
    # TPXO prediction from there.
    #
    gauge_coord = c(gauge_metadata$lon[i], gauge_metadata$lat[i])
    if(names(all_gauges_detided)[i] %in% c("GB_STHELENS_aest_ahd")){
        gauge_coord = gauge_coord + c(0.05, 0.025)
    }else if(names(all_gauges_detided)[i] %in% c("DE_NEWNORFOLK_aest_ahd")){
        gauge_coord = gauge_coord + c(0.3, -0.1)
    }else if(names(all_gauges_detided)[i] == "BOMPorts_port-pirie"){
        gauge_coord = gauge_coord + c(-0.2, -0.2)
    }

    # Detiding using BOTH tpxo AND a simple spectral approach. Graphically the
    # de-tided series turn out to be pretty much the same (as expected if the
    # tidal-prediction is low frequency).
    all_gauges_detided[[i]] = try(extract_tsunami_from_stage(all_gauges[[i]], 
        gauge_coord, low_frequency_limit = 1/(3*3600)))

}
# Quick plot
pdf(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauges_comparison_multiple_detiding_techniques.pdf'), width=10, height=8)
plot_detided_results(all_gauges_detided)
dev.off()

#
# Export all gauges to csv -- and remove the original 'stage' value for data purchased
# from the NSW Port Authority by Geoscience Australia (in this case, we have permission to
# redistribute the 'detided data'). 
#
dir.create(OUTPUT_TIDE_DIR, showWarnings=FALSE)
for(i in 1:length(all_gauges_detided)){
    output_file = gauge_metadata$postprocessed_file[i]
    file_header = c(
        '#',
        paste0('# Name: ', gauge_metadata$name[i]), 
        paste0('# Location (LON/LAT): ', gauge_metadata$lon[i], ', ', gauge_metadata$lat[i]), 
        paste0('# Original data (relative to scripts_to_postprocess directory): ', gauge_metadata$original_data_file[i]), 
        '#',
        '# Column variable definitions',
        '#   time: Date and time with UTC timezone ', 
        '#   juliant: Time in days since 1970-01-01 00:00:00 UTC', 
        '#   stage: Waterlevel (m) measured on the tide gauge (with a local vertical datum). This is missing for files with PANSW in the name, due to data licencing issues. In some instances we have removed parts of the time-series that were judged to be spurious (but users can consult the original files as needed).',
        '#   resid3h: High-frequency component of stage, retaining wave periods of 3 hours or less (m). This is provided even in the cases for which we cannot provide the stage data (due to licencing issues).',
        '#')
    cat(paste0(file_header,'\n'), file=output_file, sep="", fill=FALSE)

    output_data = data.frame(time=all_gauges_detided[[i]]$time, juliant=all_gauges_detided[[i]]$juliant,
        stage=all_gauges_detided[[i]]$stage, resid3h=all_gauges_detided[[i]]$tsunami_highpass)
    # Remove stage for Port Authority gauges, due to licencing restrictions
    if(grepl('PANSW', gauge_metadata$name[i])) output_data$stage = NA

    write.table(output_data, file=output_file, sep=",", row.names=FALSE, quote=FALSE, append=TRUE)
}
