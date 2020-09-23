#
# Run this on NCI to find the 'good fitting' scenarios for each
# 'dart buoy test event'
#
# Run in folder ......./AustPTHA_c/SOURCE_ZONES/
#
# The code does not write any files, but at the end it uses 'dput'
# to print the object with the scenarios of interest. This can be
# used in other scripts ( here it was put in best_fitting_HS.R and
# best_fitting_VAUS.R ).
#
#
#


slip_type = 'stochastic' 
#slip_type = 'variable_uniform'

number_best_scenarios = 3

event_stats_rds_files = Sys.glob('*/TSUNAMI_EVENTS/plots/*.Rdata')
# Remove the defunct source-zones {where geometries were revised, often due to SLAB2.0 release}.
to_remove = c(grep("^sunda/", event_stats_rds_files), 
              grep("^puysegur/", event_stats_rds_files),
              grep("^kermadectonga/", event_stats_rds_files),
              grep("^newhebrides/", event_stats_rds_files),
              grep("^solomon/", event_stats_rds_files))
event_stats_rds_files = event_stats_rds_files[-to_remove]


outputs = vector(mode='list', length=length(event_stats_rds_files))

for(i in 1:length(event_stats_rds_files)){

    # Read the RDS file that has the "corresponding family of scenarios" for a
    # given event, and associated DART_buoy comparison statistics.
    event_env = new.env()
    load(event_stats_rds_files[i], envir=event_env)

    if(slip_type == 'stochastic'){
        slip_stats = event_env$stochastic_slip_stats
    }else if(slip_type == 'variable_uniform'){
        slip_stats = event_env$variable_uniform_slip_stats
    }else{
        stop('unrecognized slip type')
    }

    ngauges = length(slip_stats)
    nscenarios = length(slip_stats[[1]])

    # Get the goodness of fit statistic at each gauge
    gf_matrix = matrix(NA, ncol=ngauges, nrow=nscenarios)
    for(g in 1:ngauges){
        for(j in 1:nscenarios){
            gf_matrix[j,g] = slip_stats[[g]][[j]]$model_data_similarity_time
        }
    }

    # Find the top few scenarios
    if(ncol(gf_matrix) > 1){
        # Multi-gauge case -- use the 'median over gauges' goodness of fit statistic
        scenarios_summary_stats = apply(gf_matrix, 1, median)
    }else{
        # Single gauge cause
        scenarios_summary_stats = gf_matrix[,1] 
    }

    # Rank the goodness-of-fit
    #
    # Here ensure the ranks cover all positive integer values, with repetition for ties which
    # will almost surely be the same scenario. For instance "rank" might lead to sorted values
    #     1 1 3 4 5 5 7 8
    # which we transform to
    #     1 1 2 3 4 4 5 6
    #scenarios_ranks = rank(scenarios_summary_stats)
    sss = sort(unique(scenarios_summary_stats))
    scenarios_ranks = sapply(scenarios_summary_stats, f<-function(x) match(x, sss))
    top_few = which(scenarios_ranks < (number_best_scenarios+0.5))

    # Extract the row indies for the top few
    desired_event_rows = sapply(top_few, f<-function(x){
        as.numeric(rownames(slip_stats[[1]][[x]]$events_with_Mw))})

    # Pack the desired information into the output datastructure -- we will dput this
    outputs[[i]] = list(
        stats_file = event_stats_rds_files[i],
        source_zone = basename(dirname(dirname(dirname(event_stats_rds_files[i])))),
        desired_event_rows = desired_event_rows,
        score = scenarios_summary_stats[top_few])

    # Clean up
    rm(event_env, ngauges, nscenarios, gf_matrix, scenarios_summary_stats, scenarios_ranks, slip_stats)
    gc()

}

dput(outputs)
