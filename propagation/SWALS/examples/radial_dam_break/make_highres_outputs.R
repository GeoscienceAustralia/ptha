#
# After running with a refined mesh, specified with:
#    integer(ip), parameter, dimension(2):: global_nx = [400, 400] * 8 + 1
# this code was used to save some outputs for regression testing
#
source('../../plot.R')
ti = 11
j = 1601

x = get_multidomain(Sys.glob('OUTPUTS/RUN*')[1])

vel = sqrt((x[[1]]$ud[,,ti]**2 + x[[1]]$vd[,,ti]**2)/(x[[1]]$stage[,,ti] - x[[1]]$elev0 + 1.0e-03)**2)

library(ncdf4)
fid = nc_open(Sys.glob(paste0(x[[1]]$output_folder, '/Grid*.nc')))
max_flux = ncvar_get(fid, 'max_flux')
max_speed = ncvar_get(fid, 'max_speed')
arrival_time = ncvar_get(fid, 'arrival_time')
max_stage = ncvar_get(fid, 'max_stage')

output_data = data.frame(
     x=x[[1]]$xs, 
     vel_last_time_step = vel[,j], 
     stage_last_time_step = x[[1]]$stage[,j,ti],
     uh_last_time_step = x[[1]]$ud[,j,ti],
     vh_last_time_step = x[[1]]$vd[,j,ti],
     max_stage = max_stage[,j],
     max_speed = max_speed[,j],
     arrival_time = arrival_time[,j],
     max_flux = max_flux[,j])

options('digits'=5)
write.csv(output_data, file='reference_results.csv', row.names=FALSE)
zip('reference_results.csv.zip', 'reference_results.csv')
