# Function to "manually" despike data at a few sites where that was judged to be important.
manually_despike_data<-function(gauge_data, gauge_data_file_name){

    if(grepl('puysegur_2009', gauge_data_file_name, ignore.case=TRUE) | # Different naming convention for non-random vs random scenarios
       grepl('puysegur2009', gauge_data_file_name, ignore.case=TRUE) ){
        # For Puysegur 2009, there is noise in the Twofold Bay 1 min record, 
        # especially at the beginning, but also repeated spikes later on.

        # This 'de-spiking' isn't very good.
        # SOLUTION: I removed this gauge from the good_nearshore_data - anyway now we have another Eden gauge.
        obs = gauge_data$TwofoldBay_1min_BOM$event_data$obs
        old_resid = obs$resid

        # Identify spikes as points that differ from the running median
        # by more than 3cm. Delete those points.
        old_resid_filtered = runmed(obs$resid, 3)
        p = which(abs(old_resid_filtered - old_resid) > 0.03) # Spikes
        if(length(p) > 0){
            gauge_data$TwofoldBay_1min_BOM$event_data$obs = 
                gauge_data$TwofoldBay_1min_BOM$event_data$obs[-p,]
        }

    }else if(grepl('solomon_2007', gauge_data_file_name, ignore.case=TRUE) | # Different naming convention for non-random vs random scenarios
             grepl('solomon2007', gauge_data_file_name, ignore.case=TRUE)){
        # For solomon 2007, there looks to be a spike at the end of the
        # Fort Denison tide-gauge record. 

        obs = gauge_data$Sydney_FortDenison_1min_PA$event_data$obs
        fil_resid = runmed(obs$resid, 3)
        # Subjectively, this captures the spikes I want to remove
        p = which(abs(obs$resid - fil_resid) > 0.025)
        if(length(p) > 0) obs = obs[-p,]
        gauge_data$Sydney_FortDenison_1min_PA$event_data$obs = obs

    }else if(grepl('kermadec2021', gauge_data_file_name, ignore.case=TRUE)){
        # This event shows a late-time seiche at many gauges in NSW that will inflate the late time observations.
        # Remove obs near this time
        all_sites = names(gauge_data)
        # This will match the names of gauges that should be edited
        to_edit = c( grep('1min_DPIE', all_sites), grep('1min_BOM', all_sites), grep('1min_PA', all_sites))
        sites_to_edit = all_sites[to_edit]
        for(nm in sites_to_edit){
            obs = gauge_data[[nm]]$event_data$obs
            k = which(obs$juliant < 18692.2)
            gauge_data[[nm]]$event_data$obs = obs[k,]
        }

    }else if(grepl('newhebrides2021', gauge_data_file_name, ignore.case=TRUE)){
        # This event shows late-time waves at a few gauges near Sydney and
        # Hawkesbury that I interpret as weather related, not tsunami.
        # The timings of the spurious waves differ to a degree. Truncate everything
        # after a user-provided time.
        time_truncate = data.frame( 
            site = c('Sydney_FortDenison_1min_PA_b', 'Sydney_MiddleHarbour_1min_DPIE', "Hawkesbury_Patonga_1min_DPIE", "Sydney_Botany_Bay_Pilot_Jetty_1min_PA"),
            trunctime = c(18670.1                  ,         18670.4                 ,          18670.5              ,            18670.25                    ))
        for(i in 1:nrow(time_truncate)){
            site = time_truncate$site[i]
            trunctime = time_truncate$trunctime[i]
            # Remove data at the site when julian time > trunctime
            obs = gauge_data[[site]]$event_data$obs
            k = which(obs$juliant < trunctime)
            gauge_data[[site]]$event_data$obs = obs[k,]
        }
    }
    return(gauge_data)
}


