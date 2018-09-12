#
# Make the quick_station_stage_exceedance_rates plots for all hazard points
#
library(ncdf4)

# Get station lon/lat
fid = nc_open('tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc', readunlim=FALSE)
all_hazard_point_lonlat = cbind(
    as.numeric(ncvar_get(fid, 'lon')),
    as.numeric(ncvar_get(fid, 'lat')))
nc_close(fid)

# The command-line arguments should tell us which percentiles to run.
# They range from 0-100, and are integers
# e.g. Rscript quick_station_plots_all_sites.R 0 100
command_arguments = as.numeric(commandArgs(trailingOnly=TRUE))
stopifnot(length(command_arguments) == 2)
stopifnot(all(is.finite(command_arguments)))
stopifnot(command_arguments[1] > -1 & command_arguments[2] < 101 & 
    command_arguments[1] < command_arguments[2])
stopifnot(all(round(command_arguments) == command_arguments))


# Set up the parallel cluster
library(parallel)
MC_CORES=16 # Memory limited
cl = makeCluster(MC_CORES)

# Read the code that is used to make the plots
parLapply(cl, as.list(1:MC_CORES), f<-function(x){
    source('quick_station_stage_exceedance_rates.R')
    })

# This function runs a chunk of sites in parallel
run_chunk_of_sites<-function(lower_longitude, upper_longitude){

    lower_longitude = floor(lower_longitude*100)/100   #round(lower_longitude, 2)
    upper_longitude = ceiling(upper_longitude*100)/100 #round(upper_longitude, 2)

    # Make an output directory
    output_dir = paste0('site_results/station_summary_plots_longitudes_', 
        lower_longitude, '_', upper_longitude)
    dir.create(output_dir, recursive = TRUE, showWarnings=FALSE)
    output_dir = paste0(output_dir, '/')

    site_indices = which(
        (all_hazard_point_lonlat[,1] >= lower_longitude) &
        (all_hazard_point_lonlat[,1] <=  upper_longitude) )

    # Get the coordinates we work with, in a list
    sites = lapply(as.list(site_indices), f<-function(x) all_hazard_point_lonlat[x,])

    # All the nodes need to know about the output directory
    clusterExport(cl, c('output_dir'), envir=environment())

    # Do the plots at all "sites"
    parLapply(cl, sites, f<-function(x){
        result = try(quick_source_deagg(round(x[1],4), round(x[2],4), output_dir=output_dir))
        if(class(result) == 'try-error'){
            # It is possible that a pdf device is still open, which might mess up
            # future plotting
            try(dev.off())
        }
        gc()
        return(0)
        })

    # Now, let's zip up the output_dir
    zip_dir_command = paste0('zip -r ', output_dir, '.zip output_dir')
    system(zip_dir_command)
    # In future, delete the directory that held the original files
    # unlink(output_dir, recursive=TRUE)
}

# Run all the plots, in batches where each corresponds to 1% of the
# stations
for(i in seq(command_arguments[1], command_arguments[2]-1, by=1)){

    # Longitude ranges to plot in one batch and save in one folder
    lower_longitude = quantile(all_hazard_point_lonlat[,1], p=i/100)    
    upper_longitude = quantile(all_hazard_point_lonlat[,1], p=(i+1)/100)    

    run_chunk_of_sites(lower_longitude, upper_longitude)
}


stopCluster(cl)

