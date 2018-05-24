#
# Combine subuction zone trace information from multiple sources to a single table
#

# Make sure traces_table.csv is the last one -- since we need to apply some hacks to it
sources  = c('bird_traces_table.csv', 'Arakan_traces_table.csv', 'Seramsouth_traces_table.csv', 'traces_table.csv')

#
# Data to help tweak outer-rise parameters from traces_table.csv
#

#
# Original annual_slip_x_area derived for outer-rise sources
#
outer_rise_slip_x_area = list(
    outer_rise_sunda          = 113999950,
    outer_rise_puysegur       = 5640785  ,
    outer_rise_tonga_kermadec = 141119555,
    outer_rise_new_hebrides   = 68660116,
    outer_rise_timor          = 5179095,
    outer_rise_solomons       = 67016122)

# Annual slip_x_area derived for 'corresponding subduction zones' in units of m^3 / year
subduction_slip_x_area = list(
    outer_rise_sunda         = 48458909870,
    outer_rise_puysegur      = 1101709520,
    outer_rise_tonga_kermadec= 42457361492,
    outer_rise_new_hebrides  = 14296631957,
    outer_rise_timor         = 2433991591,
    outer_rise_solomons      = 17372316500)

#
# Ratio of outer-rise moment rate to subduction-zone moment rate
#
outer_rise_multiplier = 0.11

# Nominal shear modulus values
sub_shear_mod = 3e+10
ors_shear_mod = 6e+10

# Nominal dip values
sub_dip = 15/180*pi
ors_dip = 45/180*pi

#
# Read each of the sources csv files into a big array
# When we pass through 'traces_table.csv', apply a correction to the outer-rise events
#
for(i in 1:length(sources)){
    newtab = read.csv(sources[i])

    if(sources[i] == 'traces_table.csv'){
        print('FIXING traces_table.csv')
        # Initially, JG assigned outer-rise source-zones a moment rate of ~
        # 0.4% of the corresponding subduction-zone rate, based on estimates in
        # Sleep (2012). However, comparison with observed earthquake rates
        # indicated that much higher values were required 

        # Hence, here we apply a 'fix' to the data
        tofix = names(outer_rise_slip_x_area)    
        for(j in 1:length(tofix)){

            current_moment_ratio = (outer_rise_slip_x_area[[tofix[j]]]*ors_shear_mod/cos(ors_dip)) / 
                                   (subduction_slip_x_area[[tofix[j]]]*sub_shear_mod/cos(sub_dip))

            # Assume 'relative' rates of velocities are correct, but they need to be rescaled
            indices_to_fix = which(grepl(tofix[j], newtab$name))

            velocity_scale = outer_rise_multiplier / current_moment_ratio
            print(c(tofix[j], current_moment_ratio, outer_rise_multiplier, velocity_scale, length(indices_to_fix)))

            newtab$RL_vel[indices_to_fix]  = newtab$RL_vel[indices_to_fix]  * velocity_scale
            newtab$Div_vel[indices_to_fix] = newtab$Div_vel[indices_to_fix] * velocity_scale
            newtab$Vel_L2R[indices_to_fix] = newtab$Vel_L2R[indices_to_fix] * velocity_scale
        }

    }

    if(i == 1){
        bigtab = newtab
    }else{
        bigtab = rbind(bigtab, newtab)
    }
}

# Flip Azimuth in some regions. This is done because Bird (2003) does not always put the arrows
# in the up-dip direction. We want consistency for plotting, so make this change.
library(rgdal)
library(geosphere)
flip_regions = readOGR('Bird_reverse_vector_regions', layer='Bird_reverse_vector_regions')
pointloc = cbind(0.5*(bigtab$Long1 + bigtab$Long2), 0.5*(bigtab$Lat1 + bigtab$Lat2))
pointloc = SpatialPoints(coords=pointloc, proj4string=CRS(proj4string(flip_regions)))

to_flip = which(!is.na(over(pointloc, flip_regions)))

bigtab_test = bigtab
bigtab_test$Azi_Vel[to_flip] = (bigtab_test$Azi_Vel[to_flip] - 180)%%360

write.csv(bigtab_test, 'sourcezone_traces_table_merged.csv', row.names=FALSE)
