#
# Get strike/dip/rake of scenarios to enable plotting.
# These results are not used for modelling.
#

# Scenarios we use
target_scenarios = read.csv('target_scenarios_data_frame.csv')
target_scenarios_time = strptime(target_scenarios$event_time, format='%Y-%m-%d %H:%M:%S', tz='UTC')

# ISC-GEM earthquake catalogue
gem_cat = read.csv('/media/gareth/Windows7_OS/Users/gareth/Documents/work/My Documents/Documents/Gareth/projects/tsunami/references/GENERATION/EARTHQUAKE_GENERATION/EARTHQUAKE_SIZE_SCALING_AND_BACKGROUND/GEM_EARTHQUAKE_CATALOGUE/isc-gem-up-to-2020/isc-gem/parsed_gem_catalogue.csv')
gem_cat_time = strptime(gem_cat$date, format='%Y-%m-%d %H:%M:%S', tz='UTC')

# GCMT monthly data
gcmt_monthly = read.csv('/media/gareth/Windows7_OS/Users/gareth/Documents/work/My Documents/Documents/Gareth/projects/tsunami/references/GENERATION/EARTHQUAKE_GENERATION/EARTHQUAKE_SIZE_SCALING_AND_BACKGROUND/CMT_CATALOGUE/GCMT_Monthly_2005_2023/GCMT_monthly_data.csv')
gcmt_cat_time = strptime(gcmt_monthly$datetime, format='%Y/%m/%d %H:%M:%S', tz='UTC')

# Find an event in the catalogue with time close to the target_scenario time.
# If more than one event is "close" then choose the one with the largest magnitude.
match_event<-function(target_scenario_time, cat_time, cat_mw, timethresh_secs=180){

    timediff = as.numeric(target_scenario_time - cat_time, units='secs')
    near_time_inds = which(abs(timediff) < timethresh_secs)

    if(length(near_time_inds) > 1){
        # This is sufficient for our examples
        ii = which.max(cat_mw[near_time_inds])
        near_time_inds = near_time_inds[ii]
    }

    # Avoid empty return
    if(length(near_time_inds) == 0) near_time_inds = NA

    return(near_time_inds)
}

# Extract focal mechanism data from GCMT
get_mech_gcmt<-function(ind){
    if(is.na(ind)) return(NA)
    out = data.frame(gcmt_monthly[c('Mw', 'hypo_lon', 'hypo_lat', 'strk1', 'dip1', 'rake1', 'strk2', 'dip2', 'rake2')][ind,])
    names(out) = c('Mw', 'hypo_lon', 'hypo_lat', 'strk1', 'dip1', 'rake1', 'strk2', 'dip2', 'rake2')
    out$catalogue = 'GCMT'
    return(out)
    } 
# Extract focal mechanism data from GEM catalogue, and use names consistent with above.
get_mech_gem<-function(ind){
    if(is.na(ind)) return(NA)
    out = data.frame(gem_cat[c('mw', 'lon', 'lat', 'str1', 'dip1', 'rake1', 'str2', 'dip2', 'rake2')][ind,])
    names(out) = c('Mw', 'hypo_lon', 'hypo_lat', 'strk1', 'dip1', 'rake1', 'strk2', 'dip2', 'rake2')
    out$catalogue = 'ISC-GEM'
    return(out)
    }


# Find nearby events in ISC-GEM
# Some of our events are more recent (2020+) and do not match 
gem_event_id = lapply(target_scenarios_time, function(x) match_event(x, gem_cat_time, gem_cat$mw))
# As above for GCMT -- here we are missing events before 2005
gcmt_event_id = lapply(target_scenarios_time, function(x) match_event(x, gcmt_cat_time, gcmt_cat$Mw))

gem_event_foc = lapply(gem_event_id, get_mech_gem)
gcmt_event_foc = lapply(gcmt_event_id, get_mech_gcmt)

# Define the focal mechanism using GCMT for events 2005+, and ISC-GEM elsewhere
combined_foc = vector(mode='list', length=length(gem_event_foc))
for(i in 1:length(gcmt_event_foc)){
    if(!all(is.na(gcmt_event_foc[[i]]))){
        # In practice this includes events since 2005
        combined_foc[[i]] = gcmt_event_foc[[i]]
    }else{
        # Sumatra 2004, Chile 1960. Of course we could get GCMT for Sumatra
        # 2004, and it has lower Mw due to their point-source approx.
        combined_foc[[i]] = gem_event_foc[[i]]
    }
}

# Write choice of data to file.
target_scenarios_in_catalogue = do.call(rbind, combined_foc)
target_scenarios_in_catalogue = cbind(
    target_scenarios_in_catalogue,
    data.frame(event_name = target_scenarios$event_name))

output_file = 'target_scenarios_earthquake_catalogue_mechanisms_for_plotting_only.csv'
cat(c('#',
      '# This file contains GCMT or ISC-GEM information on the focal mechanism for each of the historical earthquakes.',
      '# The information is ONLY used for plotting - our "similar magnitude range" is sometimes defined differently.',
      '# Longitude and latitude refer to the hypocentre because ISC-GEM does not provide centroid coordinates.',
      '#'),
    file=output_file,
    sep="\n")

write.table(target_scenarios_in_catalogue, output_file, append=TRUE, row.names=FALSE, sep=",")
    
