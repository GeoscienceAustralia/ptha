#
# Get the PTHA18 access codes, needed for functions below
#
get_ptha_script = './get_PTHA_results.R'
ptha18 = new.env()
source(get_ptha_script, local=ptha18, chdir=TRUE)

#
# Get the R session resulting from "compute_rates_all_sources.R", needed for functions below
#

# Get the saved R-session associated with the source-zone event-rate computation from PTHA18. 
compute_rates_session = './compute_rates_all_sources_session.RData'
if(!file.exists(compute_rates_session)){
    # If we are on NCI we might be able to get a copy from here
    compute_rates_session_NCI = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/compute_rates_all_sources_session.RData'
    if(file.exists(compute_rates_session_NCI)){
        file.copy(compute_rates_session_NCI, compute_rates_session)
    }
}
if(!file.exists(compute_rates_session)){
    # If the file wasn't found in the above locations, then download it locally
    compute_rates_session_download = 
        'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/compute_rates_all_sources_session.RData'
    download.file(compute_rates_session_download, compute_rates_session)
}
crs_data = new.env()
load(compute_rates_session, envir=crs_data)


#' For a given source-zone and segment, get the PTHA18 scenario rates and their conditional
#' probabilities [conditional on Mw]. Note in regular ptha_access outputs we provide a mixture
#' of segmented and unsegmented treatments, whereas this function can isolate each segment, as 
#' well as the unsegmented branch
#'
#' @param source_zone name of source-zone in PTHA18
#' @param segment the name of a segment associated with the source-zone, or '' if referring to the unsegmented source.
#'
#' Note that in PTHA18, the 'full' source is often combination of segmented and
#' unsegmented branches -- whereas here we separate them.
get_PTHA18_scenario_conditional_probability_and_rates_on_segment<-function(source_zone, segment=''){

    if(segment != ''){
        source_zone_segment = paste0(source_zone, '_', segment)
    }else{
        source_zone_segment = source_zone
    }

    # Convenient shorthand to refer to the environment with the target source-zones data
    sz = crs_data$source_envs[[source_zone_segment]]

    # Get available Mw values [constant rigidity case]. Should be
    #     7.2, 7.3, ...., 9.6, 9.7, 9.8 
    # but not all of these will have a non-zero rate.
    mws = sort(unique(sz$event_table$Mw))
    if(!all(abs(diff(mws) - 0.1) < 1.0e-06)) stop('discretization of mws are not as expected')
    if(!all(abs(range(mws) - c(7.2, 9.8)) < 1.0e-06)) stop('range of mws not as expected')


    # Get the 'uniform slip' scenario probabilities conditional on Mw.
    # Constant rigidity.
    FAUS_event_rates = sz$event_rates # Logic-tree mean
    FAUS_conditional_probabilities = rep(0, nrow(sz$event_table))
    FAUS_mw = sz$event_table$Mw
    for(i in 1:length(mws)){
        k = which(FAUS_mw == mws[i])
        prob_with_Mw = sum(FAUS_event_rates[k])
        if(prob_with_Mw > 0){
            FAUS_conditional_probabilities[k] = FAUS_event_rates[k]/prob_with_Mw
        }
    }

    # For each 'parent' uniform slip event, we have a set (N=15) of 'child'
    # heterogeneous-slip events (HS) and variable_area_uniform_slip events (VAUS).
    # Their unequal rates will add up to give the associated uniform_slip (FAUS)
    # rates. 
    # If we divide by the associated FAUS rate (so the child-scenario
    # numbers add to one), then those numbers will not vary depending on the
    # segment [but beware NA values caused by division by zero, associated with
    # impossible scenarios]
    get_rates_and_uniform_event_row<-function(slip_type){

        if(!any(slip_type %in% c('stochastic', 'variable_uniform'))){
            stop('slip_type must be either "stochastic" or "variable_uniform"')
        }

        # Read the data from NCI
        nc_file  = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
            'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', slip_type,
            '_slip_earthquake_events_', source_zone, '.nc')
        fid = nc_open(nc_file, readunlim=FALSE)
        rates_full_source = ncvar_get(fid, 'rate_annual')
        uniform_event_row = ncvar_get(fid, 'uniform_event_row')
        mw = ncvar_get(fid, 'Mw')
        nc_close(fid)

        unique_uniform_event_row = sort(unique(uniform_event_row))
        child_conditional_prob = rates_full_source*0
        parent_uniform_scenario_rate = rates_full_source*0
        for(i in 1:length(unique_uniform_event_row)){
            k = which(uniform_event_row == unique_uniform_event_row[i])
            parent_uniform_scenario_rate[k] = sum(rates_full_source[k]) # Deliberately all the same
            if(parent_uniform_scenario_rate[k[1]] > 0){
                child_conditional_prob[k] = rates_full_source[k]/parent_uniform_scenario_rate[k]
            }
        }

        return(data.frame(
            rates=rates_full_source, 
            uniform_event_row=uniform_event_row,
            parent_uniform_scenario_rate=parent_uniform_scenario_rate,
            child_conditional_prob=child_conditional_prob,
            Mw = mw))
    }

    HS_data = get_rates_and_uniform_event_row('stochastic')
    VAUS_data = get_rates_and_uniform_event_row('variable_uniform')

    # Get scenario probabilities on our source/segment combination, conditional
    # on a scenario with the same magnitude having occurred.
    HS_prob_given_Mw = HS_data$child_conditional_prob * 
        FAUS_conditional_probabilities[HS_data$uniform_event_row]
    VAUS_prob_given_Mw = VAUS_data$child_conditional_prob *
        FAUS_conditional_probabilities[VAUS_data$uniform_event_row]
    # Get the logic-tree mean rates applied to the individual scenarios.
    HS_event_rates = HS_data$child_conditional_prob * 
        FAUS_event_rates[HS_data$uniform_event_row]
    VAUS_event_rates = VAUS_data$child_conditional_prob * 
        FAUS_event_rates[VAUS_data$uniform_event_row]

    ## These should either be '1' (or '0' for impossible magnitudes).
    # aggregate(HS_prob_given_Mw, by=list(HS_data$Mw), sum)
    # aggregate(VAUS_prob_given_Mw, by=list(VAUS_data$Mw), sum)

    # Below are some plots used for ad-hoc check of the code [should ONLY work
    # for TOTALLY UNSEGMENTED source-zones]. It will not work where PTHA18 puts
    # partial weight on segmentation.
    if(FALSE){

        # For a totally UNSEGMENTED source-zone [e.g. puysegur2, newguinea2 ],
        # the following should plot on the 1:1 line. I confirmed it does for
        # "puysegur2" and "newguinea2".
        plot(HS_data$parent_uniform_scenario_rate, FAUS_event_rates[HS_data$uniform_event_row])
        grid(); abline(0, 1, col='red')
        # If the source-zone includes partial weight on a segmented
        # interpretation (e.g. as do most of our large source-zones), then it
        # won't plot on the 1:1 line in general, because the "file-values" read
        # above represent a mixture of the different segment interpretations

        # Absolute error -- should be within round-off on totally UNSEGMENTED source-zones
        # I confirmed this holds for puysegur2 and newguinea2
        plot((HS_data$parent_uniform_scenario_rate - FAUS_event_rates[HS_data$uniform_event_row]))

        ## Relative error -- should be consistent with round-off on totally
        ## UNSEGMENTED source-zones, beware near-zero division.
        ## I confirmed this holds for puysegur2 and newguinea2
        plot((HS_data$parent_uniform_scenario_rate - FAUS_event_rates[HS_data$uniform_event_row])/
             HS_data$parent_uniform_scenario_rate)

        # Plots as above, for VAUS
        plot(VAUS_data$parent_uniform_scenario_rate, FAUS_event_rates[VAUS_data$uniform_event_row])
        grid(); abline(0, 1, col='red')

        plot((VAUS_data$parent_uniform_scenario_rate - FAUS_event_rates[VAUS_data$uniform_event_row]))
        plot((VAUS_data$parent_uniform_scenario_rate - FAUS_event_rates[VAUS_data$uniform_event_row])/
             VAUS_data$parent_uniform_scenario_rate)
    }

    output = list(FAUS_prob_given_Mw = FAUS_conditional_probabilities,
                  FAUS_event_rates = FAUS_event_rates,
                  FAUS_mw = FAUS_mw,
                  HS_prob_given_Mw = HS_prob_given_Mw,
                  HS_event_rates = HS_event_rates,
                  HS_mw = HS_data$Mw,
                  HS_uniform_event_row = HS_data$uniform_event_row,
                  VAUS_prob_given_Mw = VAUS_prob_given_Mw,
                  VAUS_event_rates = VAUS_event_rates,
                  VAUS_mw = VAUS_data$Mw,
                  VAUS_uniform_event_row = VAUS_data$uniform_event_row)

    return(output)

}


