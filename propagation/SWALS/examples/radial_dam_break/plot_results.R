# Compare with high-res outputs. Note this is not an analytical result,
# just a regression test. 

source('../../plot.R')
ti = 11
j = 401

x = get_multidomain(Sys.glob('OUTPUTS/RUN*')[1])

# Speed at last timestep
vel = sqrt((x[[1]]$ud[,,ti]**2 + x[[1]]$vd[,,ti]**2)/(x[[1]]$stage[,,ti] - x[[1]]$elev0 + 1.0e-03)**2)

# Get various flow-maxima variables
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


output_data_reference = read.csv(unzip('reference_results.csv.zip'))


for(var in names(output_data)[2:ncol(output_data)]){

    if(var == 'arrival_time'){
        # Remove 'large-negative-values' that represent areas where the wave has not arrived
        k = which(output_data[[var]] < -1e+06)
        output_data[[var]][k] = NA
        k = which(output_data_reference[[var]] < -1e+06)
        output_data_reference[[var]][k] = NA
    }

    png(paste0('radial_dam_break_reference_solution_', var, '.png'), width=7, height=5, units='in', res=200)
    plot(output_data$x, output_data[[var]], t='p', main=var, xlab='x', ylab=var)
    points(output_data_reference$x, output_data_reference[[var]], t='l', col='red')
    legend('topleft', c('Coarse', 'Fine reference'), col=c('black', 'red'), pch=c(1, NA), lty=c(NA, 'solid'), bty='n')
    if(var == 'vh_last_time_step') title(sub='Should be approximately zero')
    dev.off()

    # Check the errors
    difference = approx(output_data_reference$x, output_data_reference[[var]], xout=output_data$x)$y - output_data[[var]]

    relerr = mean(abs(difference), na.rm=TRUE)/mean(abs(output_data[[var]]), na.rm=TRUE)

    if(var == 'arrival_time'){
        # Relative tolerance of typical error
        tol = 0.1
        if(relerr < tol){
            print('PASS')
        }else{
            print('FAIL')
        }

    }else if(var == 'vh_last_time_step'){

        # This is near-zero -- relative tolerance is not appropriate -- just check it is very small
        if(mean(abs(output_data[[var]]), na.rm=TRUE) < 1.0e-10){
            print('PASS')
        }else{
            print('FAIL')
        }

    }else{
        # Relative tolerance of typical error
        tol = 0.03
        if(relerr < tol){
            print('PASS')
        }else{
            print('FAIL')
        }
    }
}

