#
# Please don't work in the PTHA directories -- this script should run fine from
# anywhere on Raijin
#
library(rptha)

# This file contains an R session, which has objects containing the rate curves
# for all source-zones. It was produced with the code 'compute_rates_all_sources.R'
# in the same directory
RDS_input_file = '/short/w85/tsunami/MODELS/AustPTHA_c/EVENT_RATES/compute_rates_all_sources_session.RData'
# The above file will automatically be updated if I make any changes to the source-zones.
# If you get to a point where you don't want more updates, just copy the RDS_input_file
# to a safe location, and point to that instead.


# This loads the R session. The memory a few GB, enough that you'll need to get
# an interactive session on NCI.
load(RDS_input_file)

# Type 'ls()' to see all the objects -- it's like typing 'dir()' in ipython

# The source-zone rate curves are stored in a list 'source_envs'
# To see the names of all sources, type 'names(source_envs)'

#
# Example getting all the logic tree curves for a single source zone
#
args = commandArgs(trailingOnly=TRUE)
source_zone = args[1]
print(source_zone)
all_rate_curves_my_source = source_envs[[source_zone]]$mw_rate_function(NA, 
    return_all_logic_tree_branches=TRUE)

# all_rate_curves_my_source is a list, containing:
#     Mw_seq -- vector of Mw values at which the logic-tree Mw-rate cuves are evaluated, e.g. (7.2, 7.21, 7.22, ...., 9.59, 9.6)
Mw_seq = all_rate_curves_my_source$Mw_seq
#     all_par -- data.frame with rate curve parameters. Rows correspond to distinct logic-tree paths
all_par = all_rate_curves_my_source$all_par
#     all_par_prob -- vector with (post-Bayes-update) weights for each rate curve
all_par_prob = all_rate_curves_my_source$all_par_prob
#     all_par_prob_prior -- vector with prior weights for each rate curve (i.e. before considering GCMT).
all_par_prob_prior = all_rate_curves_my_source$all_par_prob_prior
#     a_parameter -- vector with GR 'a' parameter for each row in all_par
a_parameter = all_rate_curves_my_source$a_parameter
#     all_rate_matrix -- matrix, giving the rate of exceedance for each Mw in
#         Mw_seq. Number of rows = nrow(all_par) = length(all_par[,1]). Number of columns = length(Mw_seq)
all_rate_matrix = all_rate_curves_my_source$all_rate_matrix

#
# Example: Calculation of mean rate = sum[ rate_curve_i * weight_of_rate_curve_i]
#
mean_rate_my_source = Mw_seq*0
for(i in 1:nrow(all_par)){
    mean_rate_my_source = mean_rate_my_source + all_par_prob[i]*all_rate_matrix[i,]
}

plot(Mw_seq, mean_rate_my_source, log='y', xlab='Mw', ylab='Exceedance rate (events/year)')


#
# Example: Calculation of the 90th percentile (Bayesian viewpoint) of the rate
# of events above Mw=8.0
#
desired_Mw = 8.0
desired_quantile = 0.9

ind = which.min(abs(Mw_seq-desired_Mw)) # Column index
rate_values = all_rate_matrix[,ind] # Rate values for all logic trees
rate_values_order = order(rate_values) # Index of rates from lowest .... highest 
rate_values_sorted = rate_values[rate_values_order]
prob_sorted = all_par_prob[rate_values_order] # Curve weights, sorted like the rates
prob_sorted_cum = cumsum(prob_sorted) # Cumulative weights -- should sum to 1.0
stopifnot(all.equal(max(prob_sorted_cum), 1.0)) # If this is not true, I've made a mistake

index_nearest_quantile = which.min(abs(prob_sorted_cum-desired_quantile))
the_answer = rate_values_sorted[index_nearest_quantile]


#
# Example: Calculation of new mw_rate_curves, with mw_min taking a lower value
#
#
new_mw_min = 6.0
new_Mw_seq = seq(new_mw_min, 9.8, by=0.1)
new_all_rate_matrix = matrix(NA, ncol=length(new_Mw_seq), nrow=nrow(all_par))
for(i in 1:nrow(new_all_rate_matrix)){
    a = a_parameter[i]
    b = all_par$b[i]
    mw_max = all_par$Mw_max[i]
    dist_type = all_par$Mw_frequency_distribution[i]
    if(dist_type == 'truncated_gutenberg_richter'){
        new_all_rate_matrix[i,] = Mw_exceedance_rate_truncated_gutenberg_richter(
            new_Mw_seq, a = a, b = b, Mw_min = new_mw_min, Mw_max = mw_max)
    }else if(dist_type == 'characteristic_gutenberg_richter'){
        new_all_rate_matrix[i,] = Mw_exceedance_rate_characteristic_gutenberg_richter(
            new_Mw_seq, a = a, b = b, Mw_min = new_mw_min, Mw_max = mw_max)
    }else{
        stop('Error, wrong dist type, BUG')
    }
}

#
# plot(new_Mw_seq, new_all_rate_matrix[10100,])
# points(Mw_seq, all_rate_matrix[10100,], col='red', t='l')
#

#
# Example: Calculation of mean rate = sum[ rate_curve_i * weight_of_rate_curve_i]
#
mean_rate_my_source = new_Mw_seq*0
for(i in 1:nrow(all_par)){
    mean_rate_my_source = mean_rate_my_source + all_par_prob[i]*new_all_rate_matrix[i,]
}

pdfname = paste0('rate_curves_', source_zone, '.pdf')
pdf(pdfname, width=9, height=7)

# Plot the curve
plot(new_Mw_seq, mean_rate_my_source, log='y', t='o',  xlab='Mw', ylab='Exceedance rate (events/year)')

#plot(0:1, col='white')
title(source_zone)
grid(col='orange')
abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), col='orange', lty='dotted')
dev.off()

# Export to a file
output_data = data.frame(Mw = new_Mw_seq, exceedance_rate = mean_rate_my_source)
output_file = paste0( 'output_mw_rate_curve_', source_zone, '.csv')
write.csv(output_data, output_file, row.names=FALSE)

