library(rptha)

##############################################################################
# INPUT DATA 
#
interface_shapefile = '../SOURCE_CONTOURS_2_UNIT_SOURCES/source_shapefiles/alaska.shp'
# sourcename (occurs in input mux_file directory names, and in output files)
sourcename = 'alaska'
# A string which will glob all the mux2 files for the source zone
mux_files_glob = 'alaska_mux2_test/alaska_*/*mux2'
# Subfault width/length (must be the same as used to make the unit sources)
desired_subfault_length = 100
desired_subfault_width = 50
# Value of approx_dx/approx_dy passed to 'discretized_source_summary_statistics'
discrete_source_stats_dx = 4000
# Read a chunk of stations at once. If the chunksize is too large we might run
# out of memory. Must be an integer
station_chunksize = 10000
# Time between mux2 file tide gauge records
mux_timestep = 20 # FIXME: This is contained in the output file
# Minimum and maximum magnitude we should consider
Mmax = 9.2
Mmin = 7.8
#
# END INPUT
##############################################################################

## Get the discrete source geometric information

discrete_source = discretized_source_from_source_contours(interface_shapefile, 
    desired_subfault_length, desired_subfault_width, make_plot=FALSE)

unit_source_statistics = discretized_source_summary_statistics(discrete_source, 
    approx_dx = discrete_source_stats_dx, approx_dy = discrete_source_stats_dx)



## Get all earthquake events
all_eq_events = lapply(as.list(seq(Mmin, Mmax, 0.1)), 
    f<-function(x) get_all_events_of_magnitude_Mw(x, unit_source_statistics))


## Convert to a single table which holds all the events + their subfaults
add_unit_source_indices_to_event_table<-function(event){
    event_statistics = event$event_statistics

    # Make a character vector, where each entry is a string
    # with all unit source indices for that event, separated by '-'
    event_index_string = unlist(lapply(event$event_indices, 
        f<-function(x) paste0(x, sep="-", collapse="")))

    event_statistics = cbind(event_statistics, 
        data.frame(event_index_string = event_index_string))

    return(event_statistics)
}

all_eq_tables = lapply(all_eq_events, add_unit_source_indices_to_event_table)

# big_eq_table + unit_source_statistics hold everything we need about the event geometry
# (but not probability)
big_eq_table = all_eq_tables[[1]]
for(i in 2:length(all_eq_tables)) big_eq_table = rbind(big_eq_table, all_eq_tables[[i]])

nearthquakes = length(big_eq_table[,1])


## Get all the mux files (which each hold the tide-gauge results for a single
## unit source)
all_mux = Sys.glob(mux_files_glob)
# Find the 'subfault number' in the unit_source_statistics table which 
# is associated with each mux file
all_mux_subfault_numbers = match(basename(dirname(all_mux)), 
    paste0(sourcename, '_', unit_source_statistics$downdip_number, '_', 
        unit_source_statistics$alongstrike_number))


# Figure out the number of stations, and store the location information for
# later usage.
tmp = read_mux2_data(all_mux[1])
mux_data_loc = tmp$loc
nstations = length(mux_data_loc[,1])
rm(tmp)

# Make an array which can store the max height and wave period for each
height_period_array = array(NA, dim=c(nstations, nearthquakes, 2)) 

# If the number of stations is < station_chunksize, read all at once
station_chunksize = min(station_chunksize, nstations)

# Make a vector defining the start/end of each chunk
station_starts = seq(0, nstations, by = station_chunksize)
if(max(station_starts) != nstations) station_starts = c(station_starts, nstations)


## Loop over all chunks, read the data and sum the unit sources, then compute 
## max stage and wave period
for(i in 1:(length(station_starts) - 1)){

    print(paste0('Chunk ', i))

    all_data = list()
    # Get the indices of this chunk
    inds = (station_starts[i]+1):station_starts[i+1]

    ## Read all the data
    for(mux in all_mux){
        print(paste0('Reading ', mux))
        all_data[[mux]] = read_mux2_data(mux, inds=inds)
    }

    ## Check that all files are consistent with each other
    for(mux in all_mux[-1]){
        if(any(all_data[[1]]$loc != all_data[[mux]]$loc)){
            stop('mux-data location tables are not identical')
        }

        if(any(dim(all_data[[1]]$wave) != dim(all_data[[mux]]$wave))){
            stop('mux-data wave timeseries do not have the same dimensions')
        }
    }

    ## For each earthquake, add the unit source tsunami in the correct proportions
    for(j in 1:nearthquakes){
        print(paste0('Event ', j, ' of ', nearthquakes))
        slip = big_eq_table$slip[j]
        # The following gives us a vector of the 'subfault_numbers' for each
        # unit source in the event
        unit_sources = as.numeric(unlist(
            strsplit(as.character(big_eq_table$event_index_string[j]), '-')[[1]]
            ))

        wave_field = all_data[[1]]$wave*0
        # Main loop
        for(us in unit_sources){
            wave_field = wave_field + slip*all_data[[which(all_mux_subfault_numbers == us)]]$wave
        }

        # Store max wave height
        wave_max = apply(wave_field, 1, max)
        height_period_array[inds,j,1] = wave_max

        # Store max wave period
        # FIXME: Consider the effect of only working on part of the series
        # around the peak
        period_zc = apply(wave_field, 1, f<-function(x) zero_crossing_period(x, dt=mux_timestep))
        height_period_array[inds, j, 2] = period_zc
    }
}


saveRDS(mux_data_loc, file=paste0(sourcename, '_tide_gauge_loc.RDS'))
saveRDS(height_period_array, file=paste0(sourcename, '_height_period_array.RDS'))
saveRDS(big_eq_table, file=paste0(sourcename,'_earthquake_events_table.RDS'))

#' Make a quick plot of 'z' at tsunami points (indicated by height)
plot_waves_with_unit_sources<-function(x, y, z, zcol='black', z_threshold=0.1, ...){

    # Points and bars
    plot(x,y, asp=1, pch='.', col='black', ...)
    # Only put bars where z>threshold
    k = which(z>z_threshold)
    arrows(x[k], y[k], x[k], y[k]+z[k], lwd=0.2, col=zcol, length=0)

    # Unit sources
    unit_sources = as.numeric(unlist(
        strsplit(as.character(big_eq_table$event_index_string[i]), '-')[[1]]
        ))

    usg = discrete_source$unit_source_grid

    red_trans = rgb(1,0,0)
    for(us in unit_sources){
        dip_index = unit_source_statistics$downdip_number[us]
        strike_index = unit_source_statistics$alongstrike_number[us]
    
        polygon(c(usg[dip_index + 0:1, 1, strike_index], 
                  usg[dip_index + 1:0, 1, strike_index+1]),
                c(usg[dip_index + 0:1, 2, strike_index], 
                  usg[dip_index + 1:0, 2, strike_index+1]),
                col=red_trans, border=red_trans)
    }
    
    # Rupture area
    plot(discrete_source$depth_contours, add=TRUE, col=rgb(0,1,0, alpha=0.5), lwd=0.2)

}

pdf(paste0(sourcename, '_quick_plot.pdf'), width=18, height=8)
for(i in 1:nrow(big_eq_table)){
    par(mfrow=c(1,2))
    plot_waves_with_unit_sources(mux_data_loc$geolong, mux_data_loc$geolat, 
        height_period_array[,i,1]*10, 
        zcol='grey', main=paste0("Mw: ", big_eq_table$Mw[i], ' , wave height'))

    plot_waves_with_unit_sources(mux_data_loc$geolong, mux_data_loc$geolat, 
        sqrt(height_period_array[,i,2])/10 - 1, 
        zcol='grey', main=paste0("Mw: ", big_eq_table$Mw[i], ' , period'))
}
dev.off()
