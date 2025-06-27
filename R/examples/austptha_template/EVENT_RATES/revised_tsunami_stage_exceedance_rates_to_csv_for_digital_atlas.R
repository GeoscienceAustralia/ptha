#
# Convert a tsunami stage exceedance rate netcdf file to a CSV file for the digital atlas
# The idea is to provide a few variables, to serve as an 'advertisement' for serious access
# via github: https://github.com/GeoscienceAustralia/ptha/tree/master/ptha_access
#

# Combined, this information defines the return period variables we wil store
#exceedance_rates = c(1/100, 1/250, 1/500, 1/1000, 1/2500, 1/5000, 1/10000)
exceedance_rates = c(1/1000)
shear_modulus_type = "" #c('', 'variable_mu_')
#event_type = c('uniform', 'stochastic', 'variable_uniform')
event_type = 'stochastic'
#value_type = c('', '_upper_ci', '_lower_ci', '_median', '_16pc', '_84pc')
value_type = c('', '_84pc')
value_type_name = c('_logictreemean', '_84thpercentile')

input_file = 'revised1_tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc'
output_dir = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/'


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
                '_slip_stage', value_type_name[j], '_1in', (1/exceedance_rates))
            print(myvar_output_names)
            local_stage = matrix(NA, nrow=ncol(all_rates), ncol=length(myvar_output_names))
            for(m in 1:nrow(local_stage)){
                rr = all_rates[,m]
                if(all(is.na(rr))) next
                if(all(rr == 0)) next
                local_stage[m,] = approx(rr, stages, xout=exceedance_rates, ties='min')$y
            }   
 
            temp_df = as.data.frame(local_stage)
            names(temp_df) = myvar_output_names
            output = cbind(output, temp_df)
        }
    }
}
# Force to be numeric (netcdf comes out as array)
for(i in 1:ncol(output)) output[,i] = as.numeric(output[,i])

# Store gaugeID's to 1 decimal places. This doesn't loose information
# but prevents floating point truncation issues from affecting the gaugeID's
output$gaugeID = round(output$gaugeID, 1)

# Store hazard info to 3 significant figures
output_signif = output 
output_signif[,5:ncol(output)] = signif(output[,5:ncol(output)], 3)

# Store location info to 6 significant figures
output_signif[,1:3] = signif(output_signif[,1:3], 6)

#
# Remove points that are not near Australia, or that are too shallow
#
ELEV_UPPER_BOUND = -40
CLIP_POLY = structure(c(160.397151550204, 160.397151550204, 43.2425387928432, 
    62.1668028598893, 89.6509418716872, 92.5196567180365, 107.92743161859, 
    119.703043366911, 122.1090622703, 127.129313251412, 129.651006140541, 
    135.249627050352, 140.420254212925, 141.935583426118, 143.300536457848, 
    144.896837461059, 149.500662093506, 160.651634318832, 166.250255228642, 
    175.828061247905, 164.723358616876, 160.397151550204, -62.4470529104392, 
    -73.0890595985091, -69.2949528662407, -48.7049834045402, -15.5759538886357, 
    -6.32203502944444, -9.0982106872018, -11.6199035763314, -11.6199035763314, 
    -9.42209784727349, -8.91313231001797, -9.12134548434977, -9.15604768007174, 
    -9.64187842017927, -9.15604768007173, -8.25379059130059, -11.7240101634973, 
    -15.9808128387253, -25.5123492636922, -27.0855154697547, -37.0797478376812, 
    -62.4470529104392), dim = c(22L, 2L))

keep = which(output_signif$elev < ELEV_UPPER_BOUND & point.in.polygon(output_signif$lon, output_signif$lat, CLIP_POLY[,1], CLIP_POLY[,2]))

output_signif = output_signif[keep,]

write.csv(output_signif, file='PTHA18_max_stage_1in1000_heterogeneous_slip_constant_rigidity_near_Australia_deeper_than_40m_for_digital_atlas.csv', row.names=FALSE)

