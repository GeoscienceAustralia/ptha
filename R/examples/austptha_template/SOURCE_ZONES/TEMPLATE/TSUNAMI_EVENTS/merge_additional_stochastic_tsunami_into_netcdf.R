#
#
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
