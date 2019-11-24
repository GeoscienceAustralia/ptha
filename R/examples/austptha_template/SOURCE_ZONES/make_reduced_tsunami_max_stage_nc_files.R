#
# Code to create new netcdf files containing only the "max_stage" variable from files of the form
#    all_SLIPTYPE_slip_earthquake_events_tsunami_SOURCENAME.nc
# where SLIPTYPE and SOURCENAME are replaced with the relevant variable
#
# This can make remote reads of 'max_stage' faster and more reliable.
#
# The new filenames begin and end with MAX_STAGE_ONLY, i.e.:
#    MAX_STAGE_ONLY_all_SLIPTYPE_slip_earthquake_events_tsunami_SOURCENAME_MAX_STAGE_ONLY.nc
#
# When I first tried this I only put MAX_STAGE_ONLY at the end -- however, that is a bad idea,
# because it will break prior code that relies on matching only the original files with 
#    Sys.glob("TSUNAMI_EVENTS/all_*_earthquake_events_tsunami_*.nc'")
# To work around this, we begin AND end the names with MAX_STAGE_ONLY
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
    new_file = paste0(substring(nc_file, 1, nchar(nc_file)-3), '_MAX_STAGE_ONLY.nc')
    new_file = paste0(dirname(new_file), '/MAX_STAGE_ONLY_', basename(new_file))
    job_command = paste0('nccopy -V max_stage ', nc_file, ' ', new_file)
    print(job_command)
    #system(job_command)
}

library(parallel)
mclapply(all_tsunami_max_stage_nc_files, make_limited_copy, mc.cores=24, mc.preschedule=FALSE)

## Fix for the first iteration of the code
# old_names = Sys.glob('*/TSUNAMI_EVENTS/all_*MAX_STAGE_ONLY.nc')
# new_names = paste0(dirname(old_names), '/MAX_STAGE_ONLY_', basename(old_names))
# file.rename(old_names, new_names)
