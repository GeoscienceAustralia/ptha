library(rptha)

##############################################################################
# INPUT DATA 
#

# sourcename (occurs in input mux_file directory names, and in output files, and 
# is the list-entry-name of the discrete_source in the list held in discrete_source_RDS
sourcename = 'alaska'

# Discrete sources RDS filename. This holds a list with the discrete sources, with
# the name of sourcename associated with the discrete source for this source-zone
discrete_source_RDS = '../source_contours_2_unit_sources/all_discretized_sources.RDS'

# A string which will glob all the mux2 files for the source zone
mux_files_glob = 'alaska_mux2_test/alaska_*/*mux2'
#mux_files_glob = '../test/honshu/honshu_mux2_test/*/*mux2'
#mux_files_glob = paste0('../OUTPUTS/', sourcename, '/*/*mux2')

# Value of approx_dx/approx_dy passed to 'discretized_source_summary_statistics'
# Smaller values can be more accurate, because statistics are computed by 
# filling the source with integration points that have this spacing (on the surface)
discrete_source_stats_dx = 4000

# Read a chunk of stations at once. If the chunksize is too large we might run
# out of memory. Must be an integer
station_chunksize = 10000

# Time between mux2 file tide gauge records
mux_timestep = 20 # FIXME: This is contained in the output file

# Minimum and maximum magnitude we should consider
Mmax = 9.6
Mmin = 7.5
dMw = 0.1

# Basename for output_folder. Outputs will go in ./output_folder/sourcename/
# The trailing slash is important
output_folder = paste0('outputs/', sourcename, '/')

#
# END INPUT
###############################################################################

dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)

## Get the discrete source geometric information

discrete_source = readRDS(discrete_source_RDS)[[sourcename]]

## FIXME: The following lines of code should probably be replaced with reads from files
## since the same variables are required when computing event probabilities. Practically,
## the user will probably want to produce the variables at that stage, so they can
## conveniently understand the modelled earthquake probabilities and compare
## with data.
if(TRUE){
    ## Get summary statistics for unit sources (both to save for later use, and to 
    ## pass to 'get_all_earthquake_events'.
    unit_source_statistics = discretized_source_summary_statistics(discrete_source, 
        approx_dx = discrete_source_stats_dx, approx_dy = discrete_source_stats_dx)

    # Get all earthquake events
    big_eq_table = get_all_earthquake_events(discrete_source, unit_source_statistics, 
        Mmin, Mmax, dMw=dMw)
}
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
        # Note the potential for treating multiple events at once if they have the same
        # value of 'us'
        for(us in unit_sources){
            wave_field = wave_field + all_data[[which(all_mux_subfault_numbers == us)]]$wave
        }
        wave_field = slip*wave_field

        # Store max wave height
        height_period_array[inds,j,1] = apply(wave_field, 1, max) 

        # Store max wave period
        # FIXME: Consider the effect of only working on part of the series
        # around the peak
        height_period_array[inds, j, 2] = 
            apply(wave_field, 1, f<-function(x) zero_crossing_period(x, dt=mux_timestep))
    }
}


saveRDS(mux_data_loc, file=paste0(output_folder, sourcename, '_tide_gauge_loc.RDS'))
saveRDS(height_period_array, file=paste0(output_folder, sourcename, '_height_period_array.RDS'))
saveRDS(big_eq_table, file=paste0(output_folder, sourcename,'_earthquake_events_table.RDS'))
saveRDS(unit_source_statistics, file=paste0(output_folder, sourcename,'_unit_source_statistics.RDS'))
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
