library(parallel)


setwd("/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access")
source('get_PTHA_results.R')

find_min_time <- function(source_zone_event_data, point_ID, tol=0.000001, sample_size=NULL) {
    num_cores <- detectCores() - 1  # Use one less than the total number of cores
    total_events <- nrow(source_zone_event_data$events)

    # Determine the sample size
    if (is.null(sample_size) || sample_size > total_events) {
        sample_size <- total_events
    }

    # Randomly sample event IDs
    sampled_event_ids <- sample(1:total_events, sample_size)

    event_times <- mclapply(sampled_event_ids, function(id) {
        model_ts <- get_flow_time_series_at_hazard_point(source_zone_event_data, event_ID=id, hazard_point_ID=point_ID)
        flow_mat <- model_ts$flow[[1]]
        model_ts <- data.frame(time=model_ts$time, max_stage=flow_mat[,,1], uh=flow_mat[,,2], vh=flow_mat[,,3])
        ind <- which(abs(model_ts$max_stage) > tol | 
                     abs(model_ts$uh) > tol | 
                     abs(model_ts$vh) > tol)
                     
        if (length(ind) > 0) {
            return(model_ts$time[ind[1]])
        } else {
            return(Inf)
        }
    }, mc.cores = num_cores)
    
    min_time <- min(unlist(event_times))
    print(min_time / 3600)
    return(min_time)
}


# Selected point IDs on the edge of the nested domains
point_ID <- c(2130.5, 1903.5, 2118.5, 1626.5, 1437.1)
sample_size <- 100
tol=1e-8

print("Alaska")
alaska <- get_source_zone_events_data('alaskaaleutians')
for (i in 1:length(point_ID)) {
    min_time <- find_min_time(alaska, point_ID[i], tol=tol, sample_size=sample_size)
    print(min_time / 3600)
}

print("Japan")
japan <- get_source_zone_events_data('kurilsjapan')
for (i in 1:length(point_ID)) {
    min_time <- find_min_time(japan, point_ID[i], tol=tol, sample_size=sample_size)
    print(min_time / 3600)
}

print("Chile")
chile <- get_source_zone_events_data('southamerica')
for (i in 1:length(point_ID)) {
    min_time <- find_min_time(chile, point_ID[i], tol=tol, sample_size=sample_size)
    print(min_time / 3600)
}
