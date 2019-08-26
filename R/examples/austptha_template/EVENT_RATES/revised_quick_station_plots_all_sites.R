#
# Make the revised_quick_station_stage_exceedance_rates plots for all hazard points
#
# Call this script with 4 command line arguments, e.g.
#
#     Rscript revised_quick_station_plots_all_sites.R 0 100 2 16
#
# The first two command-line arguments should tell us which percentiles to run,
# while the third and fourth control the parallel distribution of work
#
# The first two commandline arguments range from 0-100, and are integers. For
# instance, to run the 'middle 20%' of hazard points (in terms of their sorted
# longitude), do:
#     Rscript revised_quick_station_plots_all_sites.R 40 60 2 16
#
# The third and fourth command line arguments are like 'this_image()' and 'num_images()'
# in Fortran2008 (or in MPI, this would be 'mpi rank + 1' and 'mpi comm size'). In
# practice the final commandline argument is equal to the number of CPUs on the
# machine, and the second-last varies from 1 up to this number. We manually run
# a separate R job for each case (in the job start-up script) It would be
# better to do this using R's paralle routines, but they are causing me memory
# problems, hence this approach

# Parse the commandline arguments
command_arguments = as.numeric(commandArgs(trailingOnly=TRUE))
stopifnot(length(command_arguments) == 4)
stopifnot(all(round(command_arguments) == command_arguments))
stopifnot(all(is.finite(command_arguments)))

longitude_percentiles = command_arguments[1:2]
stopifnot(longitude_percentiles[1] > -1 & longitude_percentiles[2] < 101 & 
    longitude_percentiles[1] < longitude_percentiles[2])

THIS_IMAGE = command_arguments[3]
NUM_IMAGES = command_arguments[4]

# Get station lon/lat
library(ncdf4)
library(parallel)
fid = nc_open('revised1_tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc', readunlim=FALSE)
all_hazard_point_lonlat = cbind(
    as.numeric(ncvar_get(fid, 'lon')),
    as.numeric(ncvar_get(fid, 'lat')))
nc_close(fid)

# This function runs a chunk of sites in parallel
run_chunk_of_sites<-function(lower_longitude, upper_longitude){

    lower_longitude = round(lower_longitude, 2)
    upper_longitude = round(upper_longitude, 2)

    # Make an output directory
    output_dir = paste0('revised1_site_results/revised1_station_summary_plots_longitudes_', 
        lower_longitude, '_', upper_longitude)
    dir.create(output_dir, recursive = TRUE, showWarnings=FALSE)
    output_dir = paste0(output_dir, '/')

    # Find which stations to run on THIS_IMAGE
    site_indices = which(
        (all_hazard_point_lonlat[,1] >= lower_longitude) &
        (all_hazard_point_lonlat[,1] <  upper_longitude) )
    my_inds = splitIndices(length(site_indices), NUM_IMAGES)[[THIS_IMAGE]]
    my_sites = site_indices[my_inds]


    for(site_ind in my_sites){
        # To prevent R's memory from growing and growing until the job is killed, submit
        # each job separately via 'system'
        my_program = paste0(" source('revised_quick_station_stage_exceedance_rates.R');", " quick_source_deagg(", 
            round(all_hazard_point_lonlat[site_ind,1],4), ' , ', 
            round(all_hazard_point_lonlat[site_ind,2],4), ' , ',
            "output_dir = '", output_dir, "')")
        run_command = paste0(' Rscript -e "', my_program, '"')
        job_output = system(run_command, intern=TRUE, ignore.stdout=TRUE, ignore.stderr=TRUE)
    }

}

# Run all the plots, in batches where each corresponds to 1% of the
# stations
for(i in seq(longitude_percentiles[1], longitude_percentiles[2]-1, by=1)){

    # Longitude ranges to plot in one batch and save in one folder
    lower_longitude = quantile(all_hazard_point_lonlat[,1], p=i/100)    
    upper_longitude = quantile(all_hazard_point_lonlat[,1], p=(i+1)/100)    

    run_chunk_of_sites(lower_longitude, upper_longitude)
}

