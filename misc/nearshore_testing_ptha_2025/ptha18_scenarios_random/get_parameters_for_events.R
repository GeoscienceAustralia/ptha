# For each event, the parameters that determine how we sample scenarios are
# defined in scripts like 
#     ./*/find_scenarios_near_historic_events.R
#
# Those scripts also do other things (which in hindsight isn't the best design
# choice!).
#
# Here we extract the parameters of interest from the scripts.
#

scripts_to_scrape = Sys.glob('./*/find_scenarios_near_historic_events.R')

all_scripts = lapply(scripts_to_scrape, readLines)

# Find the lines we want
start_ind = lapply(all_scripts, function(x) grep('target_events = list(', x, fixed=TRUE))
end_ind = lapply(all_scripts, function(x) grep('# Get the corresponding events for each', x, fixed=TRUE)-2)

code_fragments = vector(mode='list', length=length(all_scripts))
for(i in 1:length(all_scripts)){
    si = start_ind[[i]]
    ei = end_ind[[i]]
    code_fragments[[i]] = all_scripts[[i]][si:ei]
}

# Evaluate the code
target_scenarios = lapply(code_fragments, function(x) eval(parse(text=x)))

# Edit the result so it contains variables we need
for(i in 1:length(target_scenarios)){
    original_script = scripts_to_scrape[i]

    # Extract "NSCENARIOS" from the script
    nscenarios = as.numeric(gsub('NSCENARIOS = ', "", 
        all_scripts[[i]][grep("NSCENARIOS = ", all_scripts[[i]], fixed=TRUE)]))

    for(j in 1:length(target_scenarios[[i]])){
        # Info on the event
        target_scenarios[[i]][[j]]$original_script = original_script
        target_scenarios[[i]][[j]]$event_name = names(target_scenarios[[i]])[j]
        target_scenarios[[i]][[j]]$NSCENARIOS = nscenarios

        # Flatten the coordinates (for later conversion to data.frame)
        ep1 = target_scenarios[[i]][[j]]$event_point_1
        target_scenarios[[i]][[j]]$event_point_1_lon = ep1[1]
        target_scenarios[[i]][[j]]$event_point_1_lat = ep1[2]
        ep2 = target_scenarios[[i]][[j]]$event_point_2
        target_scenarios[[i]][[j]]$event_point_2_lon = ep2[1]
        target_scenarios[[i]][[j]]$event_point_2_lat = ep2[2]
        # Remove the length-2 coordinates
        target_scenarios[[i]][[j]]$event_point_1 = NULL
        target_scenarios[[i]][[j]]$event_point_2 = NULL
    }
}

# Make a data.frame with the information
target_scenarios_df = target_scenarios
for(i in 1:length(target_scenarios)){
    for(j in 1:length(target_scenarios[[i]])){
        target_scenarios_df[[i]][[j]] = as.data.frame(target_scenarios[[i]][[j]])
    }
    target_scenarios_df[[i]] = do.call(rbind, target_scenarios_df[[i]])
}
target_scenarios_df = do.call(rbind, target_scenarios_df)

# Store it for later
write.csv(target_scenarios_df, 'target_scenarios_data_frame.csv', row.names=FALSE)
