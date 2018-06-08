#
# Convert a tsunami stage exceedance rate netcdf file to a shapefile
#
# Combined, this information defines the return period variables we wil store
exceedance_rates = c(1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000, 1/2500, 1/5000, 1/10000)
shear_modulus_type = c('', 'variable_mu_')
event_type = c('uniform', 'stochastic', 'variable_uniform')
value_type = c('', '_upper_ci', '_lower_ci')

input_file = commandArgs(trailingOnly=TRUE)[1]
if(!is.na(input_file)){
    if(!file.exists(input_file)){
        stop(paste0('The following input file does not exist. Check the file name spelling. ', input_file))
    }
}else{
    stop('Must provide input file as command line argument, e.g. \n Rscript "tsunami_stage_exceedance_rate_global.nc" ')
}

library(rptha)

# Read the site coordinates
fid = nc_open(input_file, readunlim=FALSE)

# Prepare a data.frame to hold the outputs

stages = ncvar_get(fid, 'stage')

output = data.frame(lon=ncvar_get(fid, 'lon'), lat=ncvar_get(fid, 'lat'), 
    elev=ncvar_get(fid, 'elev'), gaugeID=ncvar_get(fid, 'gaugeID'))
for(k in 1:length(shear_modulus_type)){
    for(l in 1:length(event_type)){
        for(j in 1:length(value_type)){

            myvar_nc_name = paste0(shear_modulus_type[k], event_type[l], 
                '_slip_rate', value_type[j])

            all_rates = ncvar_get(fid, myvar_nc_name)

            myvar_output_names = paste0(shear_modulus_type[k], event_type[l], 
                '_slip_stage', value_type[j], '_', (1/exceedance_rates))
            print(myvar_output_names)
            local_rate = matrix(NA, nrow=ncol(all_rates), ncol=length(myvar_output_names))
            for(m in 1:nrow(local_rate)){
                rr = all_rates[,m]
                if(all(is.na(rr))) next
                if(all(rr == 0)) next
                local_rate[m,] = approx(rr, stages, xout=exceedance_rates, ties='min')$y
            }   
 
            temp_df = as.data.frame(local_rate)
            names(temp_df) = myvar_output_names
            output = cbind(output, temp_df)
        }
    }
}
write.csv(output, 
    paste0('tsunami_stages_at_fixed_return_periods_from_file_', basename(input_file), '.csv'),
    row.names=FALSE)

nc_close(fid)

#
# Apply Green's law based wave height adjustment
# Assume 1m nearshore
#
output_greens_law = output
green_adjust = (pmax(0, -output$elev))**0.25
for(i in 5:ncol(output_greens_law)){
    output_greens_law[,i] = output_greens_law[,i] * green_adjust
}
write.csv(output_greens_law, 
    paste0('tsunami_stages_greens_law_to_1m_at_fixed_return_periods_from_file_', basename(input_file), '.csv'),
    row.names=FALSE)

