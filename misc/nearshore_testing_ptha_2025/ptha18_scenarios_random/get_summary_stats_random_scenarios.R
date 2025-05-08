library(rptha)
# Get R image files storing the simulated scenaros
r_image_files = Sys.glob('set_range_of_mw*/find_scenarios_near_historical_events_R_image_C.RData')

extract_scenario_metadata<-function(r_image_file){
    load(r_image_file)
    return(random_scenarios)    
}

all_metadata = lapply(r_image_files, extract_scenario_metadata)

# Unpack the data for sampled scenarios
get_metadata_of_slip_type<-function(slip_type){
    do.call(rbind, lapply(all_metadata, 
        function(x){
            # Main metadata
            event_data = do.call(rbind, lapply(x, function(y) y[[slip_type]]$events))
            # Row indices
            event_rows = do.call(rbind, lapply(x, function(y) data.frame(desired_event_rows=y[[slip_type]]$desired_event_rows)))
            # Number of times the row was sampled (sampling with replacement)
            event_counts = do.call(rbind, lapply(x, function(y) data.frame(desired_event_rows_count=y[[slip_type]]$desired_event_rows_count)))
            # Name of the 'batch'
            batch_names = names(x)
            batch_lengths = unlist(lapply(x, function(y) nrow(y[[slip_type]]$events)))
            batch_names_rep = rep(batch_names, times=batch_lengths)
            # Name of the historical event 
            historical_event_names = unlist(lapply(batch_names_rep, function(x) strsplit(x, '-')[[1]][1]))

            stopifnot(nrow(event_data) == nrow(event_rows))
            stopifnot(nrow(event_data) == nrow(event_counts))
            stopifnot(nrow(event_data) == length(batch_names_rep))
            output = cbind(event_data, event_rows, event_counts, data.frame(batch_names=batch_names_rep, historical_event_names=historical_event_names))
        }))
}
FAUS_metadata = get_metadata_of_slip_type('uniform_slip') 
VAUS_metadata = get_metadata_of_slip_type('variable_area_uniform_slip')
HS_metadata = get_metadata_of_slip_type('heterogeneous_slip') 
metadata_row_sums = nrow(FAUS_metadata) + nrow(VAUS_metadata) + nrow(HS_metadata)

count_tifs = Sys.glob('set_range_of_mw_and_centroi*/*/*.tif')
stopifnot(length(count_tifs) == metadata_row_sums)

# Ensure we can find a matching tif in the simulations
matching_ptha18_random_scenarios = Sys.glob("/media/gareth/Windows7_OS/Users/gareth/Documents/work/Inundation_tsunami/2019_07_australia_wide_model/clean_model_files/analysis_ptha18_scenarios_2024/REDUCED_OUTPUTS/ptha18_random_like_historic*" )

match_pattern = gsub('.tif', '', basename(count_tifs))
match_count = unlist(lapply(match_pattern, function(x) length(grep(x, matching_ptha18_random_scenarios))))
missing_runs = basename(count_tifs[match_count==0])
if(length(missing_runs) > 0){
    print(missing_runs)
    stop('RUNS ABOVE ARE MISSING!!!')
}


#
# Summary info that is helpful in the paper
#
paper_summary_stats<-function(){

    # Chile 1960
    k = which(FAUS_metadata$sourcename == 'southamerica' & FAUS_metadata$Mw > 9.0)
    print(paste0('Chile 1960 FAUS slip range: ', min(FAUS_metadata$slip[k]), 
        ' ', max(FAUS_metadata$slip[k])))

    #
    # Solomon 2007
    #

    k = which(FAUS_metadata$sourcename == 'solomon2')
    print(paste0('Solomon 2007 FAUS slip range: ', min(FAUS_metadata$slip[k]), 
        ' ', max(FAUS_metadata$slip[k])))

    k = which(HS_metadata$sourcename == 'solomon2')
    HS_slip_max = unlist(lapply(HS_metadata$event_slip_string, 
        function(x) max(as.numeric(strsplit(x, split="_")[[1]]))))
    print(paste0('Solomon 2007 HS slip range: ', 
        min(HS_slip_max[k]), ' ', max(HS_slip_max[k])))

    k = which(VAUS_metadata$sourcename == 'solomon2')
    VAUS_slip_max = unlist(lapply(VAUS_metadata$event_slip_string, 
        function(x) max(as.numeric(strsplit(x, split="_")[[1]]))))
    print(paste0('Solomon 2007 VAUS slip range: ', 
        min(VAUS_slip_max[k]), ' ', max(VAUS_slip_max[k])))


}
