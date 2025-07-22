# Get the detailed source-zone info
ptha18_detailed = new.env()
source('../../get_detailed_PTHA18_source_zone_info.R', chdir=TRUE, local=ptha18_detailed)

source_zones = names(ptha18_detailed$crs_data$source_envs)


#source_zone = source_zones[1]

get_source_zone_Mw_max_CDF<-function(source_zone){

    # Object with data for the source-zone
    se = ptha18_detailed$crs_data$source_envs[[source_zone]]

    # Quick exit for defunct sources (these have rate=0)
    if(all(se$event_rates == 0) ) return(NULL)

    # Get the rate curves
    all_rate_curves = se$mw_rate_function(NA, return_all_logic_tree_branches=TRUE)

    # Get the Max-max CDF
    var = 'Mw_max'
    vars = sort(unique(all_rate_curves$all_par[[var]]))
    vars_cdf = sapply(vars, f<-function(x){
        sum(all_rate_curves$all_par_prob*(all_rate_curves$all_par[[var]] <= x))} )

    # Invert at the target percentiles
    target_percentiles = c(0, 0.01, 0.05, seq(0.1, 0.9, by=0.1), 0.95, 0.99, 1)
    Mw_max_at_target_percentiles = approx(vars_cdf, vars, xout=target_percentiles, method='constant', f=1, rule=2)

    # Output as a data.frame
    output_df = as.list(signif(Mw_max_at_target_percentiles$y, 4))
    names(output_df) = paste0('p_', as.character(Mw_max_at_target_percentiles$x))
    output_df = as.data.frame(output_df)

    output_df = cbind(data.frame(source_zone = source_zone, is_a_segment=se$is_a_segment), output_df)
    return(output_df)
}

all_results = lapply(source_zones, get_source_zone_Mw_max_CDF)
names(all_results) = source_zones

# Remove 0 length results [source-zones that were not used in PTHA18]
result_length = unlist(lapply(all_results, length))
keep_results = all_results[result_length > 0]

# Convert to data.frame and order by source-zone
output_df = do.call(rbind, keep_results)
order = order(output_df$source_zone)
output_df = output_df[order,]

write.csv(output_df, file='Mw_max_percentiles_from_PTHA18.csv', row.names=FALSE)
