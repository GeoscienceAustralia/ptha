#
# When computing stochastic slip events by linear summation, it is common that
# the initial call to 'run_make_all_earthquake_tsunami_events.PBS' will not
# complete in the 48 hour time limit [e.g. especially for large source-zones].
#
# To complete the runs, we use run_unfinished_stochastic_slip_tsunami.R, which
# splits the work up onto multiple nodes. To avoid doing unsupported parallel
# netcdf writes on those runs, they save their outputs to the R_images_tmp folder.
# This script reads those outputs and puts them in the main netcdf output file.
#

library(ncdf4)
all_R_images = Sys.glob('R_images_tmp/*.RDS') 

if( length(all_R_images) > 0 ){
   
    for(i in 1:length(all_R_images)){
        # Get the R image
        r_env = new.env()
        load(all_R_images[i], envir=r_env)

        output_nc_file = nc_open(r_env$output_file_name, readunlim=FALSE, write=TRUE)

        # Put each variable, only in the contiguous part of my_events
        ncvar_put(output_nc_file, output_nc_file$var$max_stage, r_env$gauge_event_max_stage, 
                  start=c(r_env$my_events[1], 1), count=c(length(r_env$my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$period, r_env$gauge_event_reference_period,
                  start=c(r_env$my_events[1], 1), count=c(length(r_env$my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$stage_range, r_env$gauge_event_peak_to_trough,
                  start=c(r_env$my_events[1], 1), count=c(length(r_env$my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$arrival_time, r_env$gauge_event_arrival_time,
                  start=c(r_env$my_events[1], 1), count=c(length(r_env$my_events), -1))
        gc()

        ncvar_put(output_nc_file, output_nc_file$var$initial_stage, r_env$gauge_event_initial_stage,
                  start=c(r_env$my_events[1], 1), count=c(length(r_env$my_events), -1))
        gc()

        nc_close(output_nc_file)

        rm(r_env)
        gc()
    }

}

# At this point, check that everything has finished
r_env = new.env()
load(all_R_images[i], envir=r_env)
output_nc_file = nc_open(r_env$output_file_name, readunlim=FALSE, write=TRUE)
# Find 'un-written data' flag, which is -999.999
# Because netcdf will only store to float precision, we just check for values which are
# less than the following [which should all be float(-999.999) ]
nul_r_less_than = r_env$nul_r + 1

# Get the max-stage data, which should no longer have un-written values
max_stage = ncvar_get(output_nc_file, 'max_stage')
nc_close(output_nc_file)

sum_missing = sum(max_stage < nul_r, na.rm=TRUE)
if(sum_missing > 0){
    stop('ERROR: There are still unwritten values in the file')
}else{
    print("PASS")
}
    
