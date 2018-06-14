#
# Convert a tsunami stage exceedance rate netcdf file to a shapefile
#
# Combined, this information defines the return period variables we wil store
exceedance_rates = c(1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000, 1/2500, 1/5000, 1/10000)
#shear_modulus_type = c('', 'variable_mu_')
shear_modulus_type = 'variable_mu_'
#event_type = c('uniform', 'stochastic', 'variable_uniform')
event_type = 'stochastic'
value_type = c('', '_upper_ci', '_lower_ci')

# We cannot put all earthquake-type variables in the shapefile, due to severe
# name-mangling, caused by the format's limit to short attribute names. So only
# return a single type, defined here
model_for_shapefile = 'variable_mu_stochastic_slip_stage'
# To avoid putting the above name in the csv file (which is so specific it will sound
# strange to inexperienced users, who are the target of this csv output), we replace
# column names containing that with the following variable in the csv file
substitute_name_model_for_csv = 'STAGE'

input_file = 'tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc'
output_dir = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/'


#input_file = commandArgs(trailingOnly=TRUE)[1]
#if(!is.na(input_file)){
#    if(!file.exists(input_file)){
#        stop(paste0('The following input file does not exist. Check the file name spelling. ', input_file))
#    }
#}else{
#    stop('Must provide input file as command line argument, e.g. \n Rscript "tsunami_stage_exceedance_rate_global.nc" ')
#}

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
# Force to be numeric (netcdf comes out as array)
for(i in 1:ncol(output)) output[,i] = as.numeric(output[,i])
# Store gaugeID's to 1 decimal places
output$gaugeID = round(output$gaugeID, 1)

# Store hazard info to 3 significant figures
output_signif = output 
output_signif[,5:ncol(output)] = signif(output[,5:ncol(output)], 3)
new_names  = gsub(model_for_shapefile, substitute_name_model_for_csv, names(output_signif))
write.csv(output_signif, 
    paste0(output_dir, 'tsunami_stages_at_fixed_return_periods.csv'),
    row.names=FALSE, col.names=new_names)

nc_close(fid)

#
store_as_shapefile<-function(output, name){
    
    # Because of severe name-mangling, we need to reduce the range of columns provided
    output = output
    names_output = names(output)
    priority_vars = c(1:4, grep(model_for_shapefile, names_output))
    if(length(priority_vars) <= 4){
        stop(paste0('Cannot find any ', model_for_shapefile, ' variables to put in shapefile'))
    }
    output = output[,priority_vars]

    names_output = names(output)
    names_reduced = gsub(paste0(model_for_shapefile, '_lower_ci'), 'STGl', names_output)
    names_reduced = gsub(paste0(model_for_shapefile, '_upper_ci'), 'STGu', names_reduced)
    names_reduced = gsub(model_for_shapefile, 'STG', names_reduced)

    names(output) = names_reduced

    output_sp = SpatialPointsDataFrame(coords=output[,c('lon', 'lat')], 
        data=output, proj4string=CRS("+init=epsg:4326"))

    writeOGR(output_sp, dsn=paste0(output_dir, name), layer=name, 
        driver='ESRI Shapefile', overwrite=TRUE)

}
store_as_shapefile(output_signif, 'tsunami_stages_at_fixed_return_periods')


#
# Apply Green's law based wave height adjustment Assume 1m nearshore. This is
# theoretically invalid, but often used as a 'quick and very very inaccurate'
# estimator of onshore impacts. 
#
adjust_using_greens_law = FALSE
if(adjust_using_greens_law){
    output_greens_law = output
    green_adjust = (pmax(0, -output$elev))**0.25
    for(i in 5:ncol(output_greens_law)){
        output_greens_law[,i] = output_greens_law[,i] * green_adjust
    }
    output_greens_law_signif = output_greens_law
    output_greens_law_signif[,5:ncol(output)] = signif(output_greens_law[,5:ncol(output)], 3)
    write.csv(output_greens_law_signif, 
        paste0(output_dir, 'tsunami_stages_greens_law_to_1m_at_fixed_return_periods.csv'),
        row.names=FALSE)
    store_as_shapefile(output_greens_law_signif, 'tsunami_stages_greens_law_to_1m_at_fixed_return_periods')
}
