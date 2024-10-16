#
# Quick code to get exceedance-rate of 2.4m at Hillaries (corresponding to a
# wharf with buildings in the middle of the harbour).
#
# This is hacked out of "extract_max_stage_at_a_point.R"
#

scenario_max_stages = readRDS('rdata/scenario_max_stages_115.7393_-31.8225.RDS') # This is Hillarys

# Get the logic-tree-mean importance-sampling scenario weights and rates
importance_sampling_scenarios_logic_tree_mean = list()
importance_sampling_scenarios_logic_tree_mean$outerrisesunda =
    read.csv('../../sources/hazard/random_outerrisesunda/random_scenarios_outerrisesunda_logic_tree_mean_curve_HS.csv')
importance_sampling_scenarios_logic_tree_mean$sunda2 =
    read.csv('../../sources/hazard/random_sunda2/random_scenarios_sunda2_logic_tree_mean_curve_HS.csv')
STAGE_POINTS_FOR_EXRATE_CURVE = seq(0.01, 5.01, len=1001)
nonlinear_model_MSL = 0.6


# Get PTHA18 functions
ptha18 = new.env()
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R'
file_home = '../../../../../../AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

output_max_stages_IS = STAGE_POINTS_FOR_EXRATE_CURVE + nonlinear_model_MSL

# Get the nonlinear-model max-stage for each importance_sampling_scenario
nonlinear_model_curves = list()
for(nm in names(importance_sampling_scenarios_logic_tree_mean)){
    k = which(scenario_max_stages$nonlinear_source_name == nm)
    mtch = match(importance_sampling_scenarios_logic_tree_mean[[nm]]$inds,
        scenario_max_stages$nonlinear_scenario_row[k])

    # Need one max-stage for every event in the PTHA18 -- but non NA values only for modelled events
    max_stages = rep(NA, max(importance_sampling_scenarios_logic_tree_mean[[nm]]$inds))
    max_stages[importance_sampling_scenarios_logic_tree_mean[[nm]]$inds] = scenario_max_stages$max_stage[k[mtch]]

    # Get the importance-sampling-based exceedance rates, and an approximate variance
    exrates_with_uncertainty = matrix( unlist(lapply(output_max_stages_IS, function(x){
        ptha18$estimate_exrate_uncertainty(
            importance_sampling_scenarios_logic_tree_mean[[nm]],
            max_stages, threshold_stage=x)
        })), ncol=2, byrow=TRUE)

    # Pack the data into a data-frame with useful names
    output_df = data.frame(max_stage = output_max_stages_IS,
        exrate=exrates_with_uncertainty[,1],
        exrate_variance = exrates_with_uncertainty[,2])

    nonlinear_model_curves[[nm]] = output_df
}

exrate_2.4m = approx(nonlinear_model_curves[[1]][,1], nonlinear_model_curves[[1]][,2], xout=2.4)$y + 
              approx(nonlinear_model_curves[[2]][,1], nonlinear_model_curves[[2]][,2], xout=2.4)$y

# > 1/exrate_2.4m
# [1] 1707.657

