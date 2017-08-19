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

#
# Get local parameters code
#
config_env = new.env()
source('config.R', local=config_env)

command_arguments = commandArgs(trailingOnly=TRUE)

library(ncdf4)

# The code can work with either stochastic slip, or variable uniform slip
if(any(grepl('-stochastic_slip', command_arguments))){
    all_R_images = Sys.glob(paste0(config_env$tmp_RDS_dir, '/*stochastic_slip*.RDS'))
}else if(any(grepl('-variable_uniform_slip', command_arguments))){
    all_R_images = Sys.glob(paste0(config_env$tmp_RDS_dir, '/*variable_uniform_slip*.RDS'))
}else{
    stop('Must provide either -stochastic_slip or -variable_uniform_slip as commandline arguments')
}

if( length(all_R_images) > 0 ){
   
    for(i in 1:length(all_R_images)){
        # Get the R image
        r_env = new.env()
        load(all_R_images[i], envir=r_env)

        output_nc_file = nc_open(r_env$output_file_name, readunlim=FALSE, write=TRUE)

        # Writing by looping over gauges seems faster than writing all gauges at once
        for(j in 1:(dim(r_env$gauge_event_max_stage)[2])){
            # Put each variable, only in the contiguous part of my_events
            ncvar_put(output_nc_file, output_nc_file$var$max_stage, r_env$gauge_event_max_stage[,j], 
                      start=c(r_env$my_events[1], j), count=c(length(r_env$my_events), 1))

            ncvar_put(output_nc_file, output_nc_file$var$period, r_env$gauge_event_reference_period[,j],
                      start=c(r_env$my_events[1], j), count=c(length(r_env$my_events), 1))

            ncvar_put(output_nc_file, output_nc_file$var$stage_range, r_env$gauge_event_peak_to_trough[,j],
                      start=c(r_env$my_events[1], j), count=c(length(r_env$my_events), 1))

            ncvar_put(output_nc_file, output_nc_file$var$arrival_time, r_env$gauge_event_arrival_time[,j],
                      start=c(r_env$my_events[1], j), count=c(length(r_env$my_events), 1))

            ncvar_put(output_nc_file, output_nc_file$var$initial_stage, r_env$gauge_event_initial_stage[,j],
                      start=c(r_env$my_events[1], j), count=c(length(r_env$my_events), 1))
        }

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
# less than the following [they should all be float(-999.999) ]
nul_r_less_than = r_env$nul_r + 1

# Get the max-stage data, which should no longer have un-written values
max_stage = ncvar_get(output_nc_file, 'max_stage')
nc_close(output_nc_file)

sum_missing = sum(max_stage < nul_r_less_than, na.rm=TRUE)
if(sum_missing > 0){
    stop('ERROR: There are still unwritten values in the file')
}else{
    print("PASS")
}
    
