#
# Code to create new netcdf files containing only the "max_stage" variable from files of the form
#    all_SLIPTYPE_slip_earthquake_events_tsunami_SOURCENAME.nc
# where SLIPTYPE and SOURCENAME are replaced with the relevant variable
#
# This can make remote reads of 'max_stage' faster and more reliable
#

library(rptha)

all_source_zones = dirname(Sys.glob('*/TSUNAMI_EVENTS'))

# These files include the tsunami 'max_stage' variable, but they have many other fields, and this causes
# poor performance. 
# I typically only need to read max_stage, so it is good to make a file containing only that
all_tsunami_max_stage_nc_files = c(
    paste0(all_source_zones, '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_', all_source_zones, '.nc'),
    paste0(all_source_zones, '/TSUNAMI_EVENTS/all_variable_uniform_slip_earthquake_events_tsunami_', all_source_zones, '.nc'),
    paste0(all_source_zones, '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_', all_source_zones, '.nc'))

if(!all(file.exists(all_tsunami_max_stage_nc_files))) stop('Problem matching files')

make_limited_copy<-function(nc_file){
    # This will make a file with ONLY the max_stage variable of the nc_file
    job_command = paste0('nccopy -V max_stage ', nc_file, ' ', substring(nc_file, 1, nchar(nc_file)-3), '_MAX_STAGE_ONLY.nc')
    print(job_command)
    system(job_command)
}

library(parallel)
mclapply(all_tsunami_max_stage_nc_files, make_limited_copy, mc.cores=24, mc.preschedule=FALSE)
