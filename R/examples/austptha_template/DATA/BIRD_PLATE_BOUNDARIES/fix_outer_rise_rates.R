#
# Original annual_slip_x_area derived for outer-rise sources
# These should be derived using the same file that is read into
# 'newtab' below. See the bottom of 'compute_rates_all_sources.R'
# for one example of how to derive these numbers.
#
outer_rise_slip_x_area = list(
    outer_rise_sunda          = 2214785041.78252,
    outer_rise_puysegur       = 59675721.9798216,
    outer_rise_tonga_kermadec = 1789131746.12161,
    outer_rise_new_hebrides   = 642728261.709481,
    outer_rise_timor          = 40031364.5301319,
    outer_rise_solomons       = 826092262.831024)

# Annual slip_x_area derived for 'corresponding subduction zones' in units of m^3 / year
subduction_slip_x_area = list(
    outer_rise_sunda         = 57183575618.9499,
    outer_rise_puysegur      = 1167011875.12453,
    outer_rise_tonga_kermadec= 50494228701.5377,
    outer_rise_new_hebrides  = 16898390874.3854,
    outer_rise_timor         = 974715412.491742,
    outer_rise_solomons      = 21101607868.6022 )


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

## Previously we fixed the outer-rise rates with this file. But after the
## slab2.0 update, we needed to update the 'already updated' file
#newtab = read.csv('sourcezone_traces_table_merged_uncorrected.csv')
newtab = read.csv('sourcezone_traces_table_merged_BEFORE_SLAB2.csv')

print('FIXING outer-rise convergence rates')
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

#
write.csv(newtab, 'sourcezone_traces_table_merged.csv', row.names=FALSE)
