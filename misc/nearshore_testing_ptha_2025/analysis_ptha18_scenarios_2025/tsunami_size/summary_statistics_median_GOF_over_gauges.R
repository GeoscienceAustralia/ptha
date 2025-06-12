#
# For each event and slip type, find the best fit scenarios according to a GOF metric used in "analysis_good_scenarios_plot.R"
# - For each event, rank the slip types
# - Take the average rank over all events. 
#   - This is an indication of the quality of the best scenarios (smaller is better)
#

GOF_results_to_check = c(
    'best_scenarios_good_nearshore_time_varying_hybrid_norm3_median/best_scenarios_GOF.RDS',
    'best_scenarios_good_nearshore_time_varying_hybrid_norm4_median/best_scenarios_GOF.RDS')

for(GOF_results in GOF_results_to_check){
    print('##########################')
    print('#')
    print('##########################')
    print(paste0('GOF_results: ', GOF_results))
    best_scenarios = readRDS(GOF_results)

    unique_events = names(best_scenarios)

    slip_types = names(best_scenarios[[1]])

    ranks = vector(mode='list', length=length(unique_events))
    names(ranks) = unique_events
    for(event in unique_events ){

        slip_scores = unlist(sapply(slip_types, function(st) best_scenarios[[event]][[st]]$x[1]))
        ranks[[event]] = rank(slip_scores)

    }

    combined_ranks = do.call(rbind, ranks)
    print('# Rank for each slip type and event (lower is better)')
    print(combined_ranks)
    print('# Mean rank for each slip type (lower is better)')
    print(colMeans(combined_ranks))

    print('')
}
