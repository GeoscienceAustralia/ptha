all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*.Rdata')

uniform_store = list()
stochastic_store = list()
variable_uniform_store = list()

for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    load(all_Rdata[i], envir=event_env)

    ngauges = length(event_env$uniform_slip_stats)

    uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    stochastic_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    variable_uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )

    for(gauge_ind in 1:ngauges){

        # Uniform
        uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))

        uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(uniform_store[[i]]$model[,gauge_ind]))

        # Stochastic
        stochastic_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        stochastic_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        stochastic_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(stochastic_store[[i]]$model[,gauge_ind]))

        # variable_uniform
        variable_uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        variable_uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        variable_uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(variable_uniform_store[[i]]$model[,gauge_ind]))
        
    }

}

#
# Compute the fraction of model events that had stage-range exceeding the data.
# If multiple gauges exist, take the median value over all gauges.
#
meds_uniform          = unlist(lapply(uniform_store         , f<-function(x) median(colMeans(x$model > x$data))))
meds_stochastic       = unlist(lapply(stochastic_store      , f<-function(x) median(colMeans(x$model > x$data))))
meds_variable_uniform = unlist(lapply(variable_uniform_store, f<-function(x) median(colMeans(x$model > x$data))))


meds2_uniform          = (lapply(uniform_store         , f<-function(x) (colMeans(x$model > x$data))))
meds2_stochastic       = (lapply(stochastic_store      , f<-function(x) (colMeans(x$model > x$data))))
meds2_variable_uniform = (lapply(variable_uniform_store, f<-function(x) (colMeans(x$model > x$data))))


save.image('model_data_envelope_summary_statistics.Rdata')

#
# Info on number of unique models [since double-ups can occur, e.g.
# if 2 variable_uniform slip models have the same area and Mw. This
# even happens for single unit-source stochastic models
#
for(i in 1:12){
    print(paste0('################', i))
    print(basename(all_Rdata[i]))
    print(c(uniform_store[[i]]$nunique_models, '/', length(uniform_store[[i]]$model[,1])))
    print(c(stochastic_store[[i]]$nunique_models, '/', length(stochastic_store[[i]]$model[,1])))
    print(c(variable_uniform_store[[i]]$nunique_models, '/', length(variable_uniform_store[[i]]$model[,1])))
}
