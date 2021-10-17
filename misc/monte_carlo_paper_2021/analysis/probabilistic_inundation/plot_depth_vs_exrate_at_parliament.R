#
# Get codes for working with PTHA databases
#
library(rptha)
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
source(file_here, local=ptha18, chdir=TRUE)

# Read the RDS file that was created with "depth_vs_exrate_at_gauge.R"
run_series_name = commandArgs(trailingOnly=TRUE)[1]
rds_name = paste0(run_series_name, '_depth_and_stage_exrate_curve_at_parliament.RDS')

# Get stage-vs-exrate curves at parliament
parliament_data = readRDS(rds_name)
x = parliament_data$stage_exrate_curves
coord = parliament_data$results_df[1,c('lon', 'lat')]
      
# Make arrays that nominally hold stage/elev for ALL SCENARIOS ON THIS
# SOURCE-ZONE IN THE PTHA. Actually they will only be non-NA for randomly
# sampled scenarios, however, the function we use them in
# (estimate_exrate_uncertainty) requires this structure
all_model_results = parliament_data$results_df
random_scenarios = parliament_data$scenarios_databases$unsegmented_HS # Could alternatively use any segment, no change to results.
scenarios_to_results_inds = match(random_scenarios$md_dir, all_model_results$md_dir)
random_scenario_stage = rep(NA, length(max(random_scenarios$inds)))
random_scenario_elev  = rep(NA, length(max(random_scenarios$inds)))
# The random_scenario stage and elevation which follow are only "!is.na()" at
# scenarios we actually simulated -- but that is all that is required.
random_scenario_stage[random_scenarios$inds] = all_model_results$max_stage[scenarios_to_results_inds] 
random_scenario_elev[random_scenarios$inds]  = all_model_results$elev[scenarios_to_results_inds] 
rm(random_scenarios, scenarios_to_results_inds, all_model_results) # Emphasise latter variables were just temporary 
random_scenario_depth = random_scenario_stage - random_scenario_elev

# Store exceedance-rate curves at these points
peak_stages = 10**seq(log10(0.01), log10(10), by=0.005)
peak_depths = 10**seq(log10(0.01), log10(20), by=0.005)


# Get data on "random_scenarios" -- while these all contain the same scenarios, the weight/rate information differs
# depending on whether we are working with the unsegmented, union-of-segments, or 50:50 weight of each. 
scenario_base = ifelse(grepl('_repeated', run_series_name),
    '../../sources/random_repeated/',
    '../../sources/random/')

mean_random_scenarios           = read.csv(paste0(scenario_base, 'random_scenarios_kermadectonga2_logic_tree_mean_curve_HS.csv'))
segments_union_random_scenarios = read.csv(paste0(scenario_base, 'random_scenarios_kermadectonga2_logic_tree_segmented_union_HS.csv'))
unsegmented_random_scenarios    = read.csv(paste0(scenario_base, 'random_scenarios_kermadectonga2_unsegmented_HS.csv'))

# For an existing plot, add depth (or stage) vs exrate curves with 95% MC errors
add_exrate_curves_with_mc_2sd<-function(random_scenarios, curve_type='depth', COL='black', LWD=1){
    
    if(curve_type == 'depth'){

        var_vals = peak_depths
        scenarios_exrate_var = do.call(rbind, 
            lapply(peak_depths, function(x){ 
                ptha18$estimate_exrate_uncertainty(random_scenarios, 
                    random_scenario_depth, threshold=x)}))

    }else if(curve_depth == 'stage'){

        var_vals = peak_stages
        scenarios_exrate_var = do.call(rbind, 
            lapply(peak_stages, function(x){ 
                ptha18$estimate_exrate_uncertainty(random_scenarios, 
                    random_scenario_stage, threshold=x)}))
    }

    points(var_vals, scenarios_exrate_var[,1], t='l', col=COL, lwd=2*LWD)
    for(prcntl in c(0.025, 0.975)){
        two_sigma = qnorm(prcntl)*sqrt(scenarios_exrate_var[,2])
        points(var_vals, scenarios_exrate_var[,1] + two_sigma,
            t='l', lty='dashed', col=COL, lwd=LWD)
    }

    return(invisible(scenarios_exrate_var))
}


# Here we use the basic importance-sampling weights. 
unsegmented_vals = cbind(x$unsegmented_HS$depth, x$unsegmented_HS$exrate_depth)
segmented_vals = cbind(x$unsegmented_HS$depth, 
    (x[["hukurangi_segment_HS"]]$exrate_depth + 
     x[["kermadec_segment_HS" ]]$exrate_depth +
     x[["tonga_segment_HS"    ]]$exrate_depth ))


png(paste0('Tsunami_inundation_depth_vs_frequency_parliament_', rds_name, '.png'),
    width=7.5, height=6, units='in', res=300)
par(mar=c(4.5, 6.5, 3, 2))
par(family='serif')
plot(c(0.1, 8), c(1.0e-04, 1.0e-02),
     xlab='', ylab='', log='y', col='white',
     xlim=c(0.1, 10), ylim=c(1.0e-04, 1/100), axes=FALSE)
mixed_curve = add_exrate_curves_with_mc_2sd(mean_random_scenarios, 
    curve_type='depth', COL='black', LWD=1.5)
unsegmented_curve = add_exrate_curves_with_mc_2sd(unsegmented_random_scenarios, 
    curve_type='depth', COL='red')
segmented_curve = add_exrate_curves_with_mc_2sd(segments_union_random_scenarios, 
    curve_type='depth', COL='green')
axis(side=1, cex.axis=1.3)
axis(side=2, at=10**seq(-5, -2), cex.axis=1.3,
     labels=c('1/100000', '1/10000', '1/1000', '1/100'), 
     las=1)
add_log_axis_ticks(side=2)
abline(v=seq(0,10),  col='darkgrey', lty='dotted')
mtext('Tsunami inundation depth (m)', side=1, line=2.5, cex=1.4)
mtext('Mean exceedances per year', side=2, line=4.5, cex=1.4)

abline(h=1/c(500, 2500, 10000), col='orange')

title(paste0('Onshore depth exceedance-rate \n (', 
   'lon=', round(coord[1],4), ', lat=', round(coord[2],4), ', sea-level = MSL)'), 
    cex.main=1.5)

legend('topright', 
    c('Logic-tree mean (50% unsegmented, 50% union-of-segments)', 
      'Logic-tree mean (100% unsegmented)', 
      'Logic-tree mean (100% union-of-segments)'),
    lty=c('solid', 'solid', 'solid'), 
    lwd=c(2, 2, 2), 
    col=c('black', 'red', 'green'),
    cex=1.2, 
    bg=rgb(1,1,1,alpha=0.5), 
    box.col=rgb(1,1,1,alpha=0.5))

text(0.8, 2.5e-03, 
     'Dashed lines give 95% CI for the \n exceedance-rate in Equation 1 (single Monte-Carlo sample method)', 
     adj=c(0,0), cex=1.2)

dev.off()

# Combined segmented/unsegmented exceedance-rates, with 50% weight on each case
mean_vals = 0.5*(unsegmented_vals + segmented_vals)
target_exrate = uniroot(f<-function(x){(1 - exp(-50*x)) - 0.02}, lower=1/3000, upper=1/2000, tol=1e-12)$root

# Some useful stats
depth_at_target_exrate = approx(mean_vals[,2], mean_vals[,1], xout=target_exrate, ties='min')$y
print(c('depth-at-target-exrate (basic importance-sampling): ', depth_at_target_exrate))

exrate_10cm = approx(mean_vals[,1], mean_vals[,2], xout=0.1)$y
print(c('exrate-of-10cm-inundation (basic importance-sampling): ', exrate_10cm))
print(c('   chance in 50 years    : ', 1 - exp(-50 * exrate_10cm)))


# 2-sigma confidence intervals for the true exceedance-rate, using the single Monte-Carlo sample
# method with the normal approximation.
get_2_sigma_CI_width_at_exrate<-function(curve=mixed_curve, exrate=1/1000){
    # depth uncertainty for fixed exceedance-rate
    c1 = approx(curve[,1] + qnorm(0.975)*sqrt(curve[,2]), peak_depths, xout=exrate, ties='max', rule=2)$y
    c2 = approx(curve[,1] + qnorm(0.500)*sqrt(curve[,2]), peak_depths, xout=exrate, ties='max', rule=2)$y
    c3 = approx(curve[,1] + qnorm(0.025)*sqrt(curve[,2]), peak_depths, xout=exrate, ties='max', rule=2)$y

    # Exceedance-rate uncertainty at the middle-curve depth
    r1 = approx(peak_depths, curve[,1] + qnorm(0.975)*sqrt(curve[,2]), xout=c2, ties='max', rule=2)$y
    r2 = approx(peak_depths, curve[,1] + qnorm(0.500)*sqrt(curve[,2]), xout=c2, ties='max', rule=2)$y
    r3 = approx(peak_depths, curve[,1] + qnorm(0.025)*sqrt(curve[,2]), xout=c2, ties='max', rule=2)$y

    # return(c(c3, c2, c1, r1, r2, r3))
    print(c('exrate_uncertainties (fraction):', r1/r2 - 1, 1 - r3/r2))
    print(c('                     (absolute):', r1-r2 , r3-r2))
    print(c('stage_uncertainties (fraction):', c1/c2 - 1, 1 - c3/c2))
    print(c('                    (absolute):', c1-c2, c3-c2))
}

get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/500)
get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/500)
get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/500)
#> get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/500)
#[1] "exrate_uncertainties (fraction):" "0.0870361669312878"               "0.0870361669312879"              
#[1] "                     (absolute):" "0.000173172352988834"             "-0.000173172352988835"           
#[1] "stage_uncertainties (fraction):" "1.1082754018285"                 "0.906682584069629"              
#[1] "                    (absolute):" "0.163940022371816"               "-0.134119698832324"             
#> get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/500)
#[1] "exrate_uncertainties (fraction):" "0.0928823042779297"               "0.0928823042779297"              
#[1] "                     (absolute):" "0.000160827386718396"             "-0.000160827386718396"           
#[1] "stage_uncertainties (fraction):" "0"                               "0"                              
#[1] "                    (absolute):" "0"                               "0"                              
#> get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/500)
#[1] "exrate_uncertainties (fraction):" "0.0850768894434943"               "0.0850768894434943"              
#[1] "                     (absolute):" "0.000169910839208879"             "-0.000169910839208879"           
#[1] "stage_uncertainties (fraction):" "0.332232999402822"               "0.49371598383845"               
#[1] "                    (absolute):" "0.21503730276619"                "-0.319556918451835"             

get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/1000)
get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/1000)
get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/1000)
#[1] "exrate_uncertainties (fraction):" "0.0888800357867821"               "0.0888800357867818"              
#[1] "                     (absolute):" "8.77043426985337e-05"             "-8.77043426985335e-05"           
#[1] "stage_uncertainties (fraction):" "0.0964359110417132"              "0.0582198438601388"             
#[1] "                    (absolute):" "0.163529473729173"               "-0.0987252603744777"            
#> get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/1000)
#[1] "exrate_uncertainties (fraction):" "0.0971230796959037"               "0.0971230796959036"              
#[1] "                     (absolute):" "9.71230796959037e-05"             "-9.71230796959036e-05"           
#[1] "stage_uncertainties (fraction):" "0.109052800217491"               "0.134203655383933"              
#[1] "                    (absolute):" "0.143077699093755"               "-0.176075719137973"             
#> get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/1000)
#[1] "exrate_uncertainties (fraction):" "0.0865301785998529"               "0.0865301785998528"              
#[1] "                     (absolute):" "8.65301785998528e-05"             "-8.65301785998528e-05"           
#[1] "stage_uncertainties (fraction):" "0.141075455284259"               "0.0849226080968072"             
#[1] "                    (absolute):" "0.306241472996417"               "-0.184346912380055"             


get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/2500)
get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/2500)
get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/2500)

#> get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/2500)
#[1] "exrate_uncertainties (fraction):" "0.104933751266742"                "0.104933751266742"               
#[1] "                     (absolute):" "4.19735005066966e-05"             "-4.19735005066966e-05"           
#[1] "stage_uncertainties (fraction):" "0.0552454927088246"              "0.0665552829045568"             
#[1] "                    (absolute):" "0.217257775972231"               "-0.261734524104207"             
#> get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/2500)
#[1] "exrate_uncertainties (fraction):" "0.10629086769016"                 "0.10629086769016"                
#[1] "                     (absolute):" "4.25163470760642e-05"             "-4.25163470760642e-05"           
#[1] "stage_uncertainties (fraction):" "0.0648682178387623"              "0.0788156182197006"             
#[1] "                    (absolute):" "0.215624838080823"               "-0.26198661660638"              
#> get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/2500)
#[1] "exrate_uncertainties (fraction):" "0.110870714903526"                "0.110870714903526"               
#[1] "                     (absolute):" "4.43482859614103e-05"             "-4.43482859614104e-05"           
#[1] "stage_uncertainties (fraction):" "0.0276820272963525"              "0.0683343129010179"             
#[1] "                    (absolute):" "0.124472088752435"               "-0.307264875118914"          

get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/10000)
get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/10000)
get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/10000)

#> get_2_sigma_CI_width_at_exrate(curve=mixed_curve, exrate=1/10000)
#[1] "exrate_uncertainties (fraction):" "0.15587118454987"                 "0.155871184549869"               
#[1] "                     (absolute):" "1.5587118454987e-05"              "-1.55871184549869e-05"           
#[1] "stage_uncertainties (fraction):" "0.0449676051152661"              "0.0964236326744505"             
#[1] "                    (absolute):" "0.333690811335038"               "-0.71553021639762"              
#> get_2_sigma_CI_width_at_exrate(curve=segmented_curve, exrate=1/10000)
#[1] "exrate_uncertainties (fraction):" "0.154700931364483"                "0.154700931364484"               
#[1] "                     (absolute):" "1.54700931364484e-05"             "-1.54700931364484e-05"           
#[1] "stage_uncertainties (fraction):" "0.0427546222385256"              "0.0612323471831122"             
#[1] "                    (absolute):" "0.257459573466193"               "-0.368728646464158"             
#> get_2_sigma_CI_width_at_exrate(curve=unsegmented_curve, exrate=1/10000)
#[1] "exrate_uncertainties (fraction):" "0.172761539438251"                "0.172761539438251"               
#[1] "                     (absolute):" "1.72761539438251e-05"             "-1.72761539438251e-05"           
#[1] "stage_uncertainties (fraction):" "0.0238681168510235"              "0.0640082058230022"             
#[1] "                    (absolute):" "0.196039326769037"               "-0.525727507350343"             








#
#
#
# Get full probabilistic inundation and uncertainties
#
#
#





# Get the R session resulting from "compute_rates_all_sources.R", needed to 
# work with alternative logic-tree branches
ptha18_source_rate_env = new.env()
source(paste0(dirname(file_here), '/get_detailed_PTHA18_source_zone_info.R'),
       local=ptha18_source_rate_env, chdir=TRUE)

# All logic-tree-branches for the kermadectonga 2 unsegmented, and the segments
all_source_names = c('kermadectonga2', 'kermadectonga2_tonga', 'kermadectonga2_kermadec', 'kermadectonga2_hikurangi')

all_source_samples = list(
    'kermadectonga2' = paste0(scenario_base, 'random_scenarios_kermadectonga2_unsegmented_HS.csv'),
    'kermadectonga2_tonga' = paste0(scenario_base, 'random_scenarios_kermadectonga2_tonga_segment_HS.csv'),
    'kermadectonga2_kermadec' = paste0(scenario_base, 'random_scenarios_kermadectonga2_kermadec_segment_HS.csv'),
    'kermadectonga2_hikurangi' = paste0(scenario_base, 'random_scenarios_kermadectonga2_hukurangi_segment_HS.csv'))
stopifnot(all(all_source_names == names(all_source_samples)))

# Read the random scenarios (note: scenarios are the same in all cases, but the nominal rates differ)
all_samples = list()
for(nm_i in all_source_names){
    all_samples[[nm_i]] = read.csv(all_source_samples[[nm_i]])
}

# Get the depth exceedance rates over all logic-tree branches, for unsegmented and segmented models
all_depth_exrates_all_logic_trees = list()
# For other functions below it is useful to store some components of these results differently.
# The following lists do that.
depth_exceedances_logic_tree_branch = list()
back_computed_mean_curve = list()
logic_tree_posterior_probs = list()
for(nm_i in all_source_names){

    # Make "source_zone" and "segment" to match the ptha18_source_rate_env information

    source_zone = 'kermadectonga2'
    if(source_zone == nm_i){
        # Unsegmented branch
        segment = ''
    }else{
        # Segmented branch -- extract the segment name directly
        segment = strsplit(nm_i, '_')[[1]][2]
    }

    # Compute exceedance rates for all peak_depths, for all logic-tree branches
    all_depth_exrates_all_logic_trees[[nm_i]] = 
        ptha18_source_rate_env$random_scenario_exceedance_rates_all_logic_tree_branches(
            source_zone = source_zone,
            segment = segment,
            random_scenarios = all_samples[[nm_i]],
            all_scenario_stage = random_scenario_depth, 
            threshold_stages=peak_depths)

    #
    # Alternative storage for some results, which is used in code below
    #

    # Exceedance rates for each logic tree branch, and each peak depth
    depth_exceedances_logic_tree_branch[[nm_i]] = 
        all_depth_exrates_all_logic_trees[[nm_i]]$logic_tree_branch_exceedance_rates
    # Average of the depth exceedance-rate curves over all logic tree branches
    # (accounting for posterior probs)
    back_computed_mean_curve[[nm_i]] = 
        all_depth_exrates_all_logic_trees[[nm_i]]$logic_tree_mean_exceedance_rates
    logic_tree_posterior_probs[[nm_i]] = 
        all_depth_exrates_all_logic_trees[[nm_i]]$logic_tree_branch_posterior_prob

}

#
# Checks on the consistency of the "logic-tree-mean" as derived from the branches above, vs the earlier
# results that are indirectly derived.
#

# Test that the mean of the depth-exceedances for the unsegmented "kermadectonga2" source match the previous result
stopifnot(all( abs(back_computed_mean_curve$kermadectonga2 - unsegmented_curve[,1]) <= 1.0e-08*unsegmented_curve[,1]))

# Test that the sum of the depth-exceedences for the segmented sources match the previous result
local_back_computed_mean_curve = rep(0, length(back_computed_mean_curve$kermadectonga2))
for(nm_i in c('kermadectonga2_tonga', 'kermadectonga2_kermadec', 'kermadectonga2_hikurangi')){
    local_back_computed_mean_curve = local_back_computed_mean_curve + back_computed_mean_curve[[nm_i]]
}
stopifnot(all( abs(local_back_computed_mean_curve - segmented_curve[,1]) <= 1.0e-08*segmented_curve[,1]))

#
# We can now compute CDF's for the exceedance-rate of each value of
# peak_depths, on the unsegmented branch, and on each individual segment
#
# To derive the union-of-segments we need to decide how to combine the
# uncertainty in the segments. Co-monotonic? Random? 
#
# To derive the 'unsegmented + union-of-segments', we are just doing a 50:50
# combination of the two distributions.
#
# Below this is all implemented by sampling randomly from the unsegmented and
# segmented models (with sample size from each depending on their respective
# weights), and pooling the samples. 
# 

segment_inds = 2:4 # List entries corresponding to segments
percentile_probs = c(0.025, 0.16, 0.5, 0.84, 0.975)

# Reproducible randomness
set.seed(123)
# Mixed segmented, unsegmented, comonotonic segments
mean_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_depth_exrates_all_logic_trees$kermadectonga2,
    all_depth_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=0.5,
    union_of_segments_wt=0.5,
    segments_copula_type='comonotonic',
    percentile_probs=percentile_probs)

set.seed(123)
# Mixed segmented, unsegmented, independent segments
mean_curve_env_indep = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_depth_exrates_all_logic_trees$kermadectonga2,
    all_depth_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=0.5,
    union_of_segments_wt=0.5,
    segments_copula_type='independent',
    percentile_probs=percentile_probs)

set.seed(123)
# Unsegmented
unseg_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_depth_exrates_all_logic_trees$kermadectonga2,
    all_depth_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=1.0,
    union_of_segments_wt=0.0,
    segments_copula_type='comonotonic',
    percentile_probs=percentile_probs)

set.seed(123)
# Segmented, comonotonic dependence
seg_curve_env = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_depth_exrates_all_logic_trees$kermadectonga2,
    all_depth_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=0.0,
    union_of_segments_wt=1.0,
    segments_copula_type='comonotonic',
    percentile_probs=percentile_probs)

set.seed(123)
# Segmented, independent segments
seg_curve_env_indep = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
    all_depth_exrates_all_logic_trees$kermadectonga2,
    all_depth_exrates_all_logic_trees[segment_inds],
    N=1e+06,
    unsegmented_wt=0.0,
    union_of_segments_wt=1.0,
    segments_copula_type='independent',
    percentile_probs=percentile_probs)


plot_panel<-function(mean_curve_env, plot_title=''){
    # Used to plot the epistemic uncertainties below

    plot(peak_depths, mean_curve_env$mean_exrate, 
        t='l', log='y', ylim=c(1e-04, 1e-02), xlim=c(0.1, 10), lwd=3, las=1, 
        xlab="", ylab="",  cex.lab=1.4, cex.axis=1.4, axes=FALSE)
    axis(side=1, cex.axis=1.5)
    axis(side=2, at=10**seq(-5, -1), cex.axis=1.5,
         labels=c('1/100000', '1/10000', '1/1000', '1/100', '1/10'), 
         las=1)
    add_log_axis_ticks(side=2)
    mtext(side=1, "Tsunami inundation depth (m)", cex=1.7, line=2.5)
    mtext(side=2, "Mean exceedances per year ", cex=1.7, line=4.5)
    title(main=plot_title, cex.main=2.2)
    points(peak_depths, mean_curve_env$percentile_exrate[3,],  t='l', col='brown', lwd=3)
    points(peak_depths, mean_curve_env$percentile_exrate[2,],  t='l', col='blue', lwd=2)
    points(peak_depths, mean_curve_env$percentile_exrate[4,],  t='l', col='blue', lwd=2)
    points(peak_depths, mean_curve_env$percentile_exrate[1,], t='l', col='green', lwd=2)
    points(peak_depths, mean_curve_env$percentile_exrate[5,], t='l', col='green', lwd=2)
    abline(v=seq(0, 20), col='grey', lty='dotted')
    abline(h=1/c(500, 2500, 10000), col='orange')

    legend('right', c('Mean', 'Median', '16/84 %', '2.5/97.5 %'), 
           col=c('black', 'brown', 'blue', 'green'),
           lwd=c(3,3,2,2), pch=c(NA, NA, NA, NA), cex=1.8, bty='n')

    }


#
# Repeat an earlier plot, combined with the segmented/unsegmented results
#

png(paste0('Tsunami_inundation_depth_vs_frequency_parliament_COMBINED_', rds_name, '.png'),
    width=15, height=12, units='in', res=300)
par(mfrow=c(2,2))
par(mar=c(4.5, 6.5, 4, 2))
par(family='serif')
plot(c(0.1, 10), c(1.0e-04, 1.0e-02),
     xlab='', ylab='', log='y', col='white',
     xlim=c(0.1, 10), ylim=c(1.0e-04, 1/100), axes=FALSE)
mixed_curve = add_exrate_curves_with_mc_2sd(mean_random_scenarios, 
    curve_type='depth', COL='black', LWD=1.5)
unsegmented_curve = add_exrate_curves_with_mc_2sd(unsegmented_random_scenarios, 
    curve_type='depth', COL='red')
segmented_curve = add_exrate_curves_with_mc_2sd(segments_union_random_scenarios, 
    curve_type='depth', COL='green')
axis(side=1, cex.axis=1.5)
axis(side=2, at=10**seq(-5, -2), cex.axis=1.5,
     labels=c('1/100000', '1/10000', '1/1000', '1/100'), 
     las=1)
add_log_axis_ticks(side=2)
abline(v=seq(0,10),  col='darkgrey', lty='dotted')
mtext('Tsunami inundation depth (m)', side=1, line=2.5, cex=1.7)
mtext('Mean exceedances per year', side=2, line=4.5, cex=1.7)

abline(h=1/c(500, 2500, 10000), col='orange')

title('Monte-Carlo uncertainity in depth exceedance-rates: \n Example of logic-tree mean curves', 
      cex.main=2.2)

legend('topright', 
    c('50% unsegmented, 50% union-of-segments', 
      '100% unsegmented', 
      '100% union-of-segments'),
    lty=c('solid', 'solid', 'solid'), 
    lwd=c(2, 2, 2), 
    col=c('black', 'red', 'green'),
    cex=2.0, 
    bg=rgb(1,1,1,alpha=0.5), 
    box.col=rgb(1,1,1,alpha=0.5))

text(3.2, 1.0e-03, 
     'Dashed lines show 95% CI for the \nexceedance-rate in Equation 1', 
     adj=c(0,0), cex=2.1)

text(0.5, 1.5e-04, 'A)', cex=3.5)

# Panel 2 -- mixed 50:50 unsegmented:segmented model
plot_panel(mean_curve_env, 
    plot_title='Epistemic uncertainty in depth exceedance-rates: \n 50% Unsegmented, 50% Union-of-segments')
text(0.5, 1.5e-04, 'B)', cex=3.5)
# Illustrate impact of assuming independence between segments
points(peak_depths, mean_curve_env_indep$percentile_exrate[3,],  t='l', col='brown', lwd=1, lty='dashed')
points(peak_depths, mean_curve_env_indep$percentile_exrate[2,],  t='l', col='blue', lwd=1, lty='dashed')
points(peak_depths, mean_curve_env_indep$percentile_exrate[4,],  t='l', col='blue', lwd=1, lty='dashed')
points(peak_depths, mean_curve_env_indep$percentile_exrate[1,], t='l', col='green', lwd=1, lty='dashed')
points(peak_depths, mean_curve_env_indep$percentile_exrate[5,], t='l', col='green', lwd=1, lty='dashed')
text(2, 4.0e-03, adj=c(0,0), 'Solid lines: comonotonic segments \nDashed lines: independent segments', cex=2.1)

# Panel 3 -- 100% unsegmented model
plot_panel(unseg_curve_env, 
    plot_title='Epistemic uncertainty in depth exceedance-rates: \n 100% Unsegmented')
text(0.5, 1.5e-04, 'C)', cex=3.5)

# Panel 4 -- 100% segmented model
plot_panel(seg_curve_env, 
    plot_title='Epistemic uncertainty in depth exceedance-rates: \n 100% Union-of-segments')
text(0.5, 1.5e-04, 'D)', cex=3.5)
# Illustrate impact of assuming independence between segments
points(peak_depths, seg_curve_env_indep$percentile_exrate[3,],  t='l', col='brown', lwd=1, lty='dashed')
points(peak_depths, seg_curve_env_indep$percentile_exrate[2,],  t='l', col='blue', lwd=1, lty='dashed')
points(peak_depths, seg_curve_env_indep$percentile_exrate[4,],  t='l', col='blue', lwd=1, lty='dashed')
points(peak_depths, seg_curve_env_indep$percentile_exrate[1,], t='l', col='green', lwd=1, lty='dashed')
points(peak_depths, seg_curve_env_indep$percentile_exrate[5,], t='l', col='green', lwd=1, lty='dashed')
text(2, 4.0e-03, adj=c(0,0), 'Solid lines: comonotonic segments \nDashed lines: independent segments', cex=2.1)
dev.off()

save.image('R_session_plot_depth_vs_exrate_parliament.RData')

#
# Statistic used in the paper -- how different is the mean curve and the 84th percentile curve?
#
mean_curve_depths = approx(mean_curve_env$mean_exrate,                peak_depths, xout=10**seq(-4, -1, len=100), ties='min')$y
p84_ind = which.min(abs(mean_curve_env$percentile_probs - 0.84))
pc84_curve_depths = approx(mean_curve_env$percentile_exrate[p84_ind,], peak_depths, xout=10**seq(-4, -1, len=100), ties='min')$y
summary(pc84_curve_depths - mean_curve_depths)
#> summary(pc84_curve_depths - mean_curve_depths)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1.146   1.332   1.493   1.525   1.683   2.111      56 



if(FALSE){
    #
    # Some plotting that is nice but I didn't ultimately use.
    #
    plot_single_source_uncertainties<-function(peak_depths, 
        depth_exceedances_logic_tree_branch, curve_weights, back_computed_mean_curve){

        YLIM = c(1.0e-04, 1.0e-02)
        XLIM = c(0.01, 15)

        plot(XLIM, YLIM, xlim=XLIM, ylim=YLIM, col='white', log='y', asp=1, 
            xlab="", ylab="",
            cex.lab=1.4, cex.axis=1.4, las=1, xaxs='i')
        mtext(side=1, 'Inundation Depth (m)', cex=1.4, line=3)
        mtext(side=2, 'Exceedance-rate (events/year)', cex=1.4, line=5)
        
        # Add logic-tree branches, semi-transparent
        greycol = rgb(0.5, 0.5, 0.5, alpha=0.05)
        for(i in 1:ncol(depth_exceedances_logic_tree_branch)){
            points(peak_depths, depth_exceedances_logic_tree_branch[,i], t='l', col=greycol, lwd=0.5)
        }

        points(peak_depths, back_computed_mean_curve, t='l', lwd=2)

        rates_50  = apply(depth_exceedances_logic_tree_branch, 1, function(x) weighted_percentile(x, curve_weights, p=0.500))
        rates_84  = apply(depth_exceedances_logic_tree_branch, 1, function(x) weighted_percentile(x, curve_weights, p=0.840))
        rates_16  = apply(depth_exceedances_logic_tree_branch, 1, function(x) weighted_percentile(x, curve_weights, p=0.160))
        rates_975 = apply(depth_exceedances_logic_tree_branch, 1, function(x) weighted_percentile(x, curve_weights, p=0.975))
        rates_025 = apply(depth_exceedances_logic_tree_branch, 1, function(x) weighted_percentile(x, curve_weights, p=0.025))

        points(peak_depths, rates_975, t='l', col='green', lty='dashed')
        points(peak_depths, rates_025, t='l', col='green', lty='dashed')

        points(peak_depths, rates_84, t='l', col='blue', lty='dashed')
        points(peak_depths, rates_16, t='l', col='blue', lty='dashed')

        points(peak_depths, rates_50, t='l', col='orange')

        output_df = data.frame(
            peak_depths=peak_depths, mean_curve=back_computed_mean_curve, 
            rates_50=rates_50, rates_16=rates_16, rates_84=rates_84,
            rates_975=rates_975, rates_025=rates_025)
        return(invisible(output_df))
    }

    # Plot each individual segment (or unsegmented), and store the curves
    plot_titles = list(
        'kermadectonga2' = 'Unsegmented',
        'kermadectonga2_tonga' = 'Tonga segment',
        'kermadectonga2_kermadec' = 'Kermadec segment',
        'kermadectonga2_hikurangi' = 'Hikurangi segment')
    store_results = vector(mode='list', length=length(plot_titles)) # Used to store summary statistics
    names(store_results) = names(plot_titles)
    png(paste0('Tsunami_inundation_depth_vs_frequency_parliament_uncertainty_', rds_name, '.png'),
        width=12, height=9, units='in', res=300)
    par(mfrow=c(2,2))
    par(mar = c(5,6.5,2,1))
    options(scipen=5)
    for(nm_i in names(depth_exceedances_logic_tree_branch)){
        store_results[[nm_i]] = plot_single_source_uncertainties(
            peak_depths, 
            depth_exceedances_logic_tree_branch[[nm_i]], 
            logic_tree_posterior_probs[[nm_i]], 
            back_computed_mean_curve[[nm_i]])
        title(main=plot_titles[[nm_i]], cex.main=2)
    }
    dev.off()
 

    #
    # Separate plot of unsegmented/segmented models
    #
    unsegmented_model = store_results$kermadectonga2
    segmented_model_comonotonic = store_results$kermadectonga2_tonga
    NN = ncol(segmented_model_comonotonic)
    segmented_model_comonotonic[,2:NN] = segmented_model_comonotonic[,2:NN] + 
        store_results$kermadectonga2_kermadec[,2:NN] + store_results$kermadectonga2_hikurangi[,2:NN]

    par(mfrow=c(1,2))
    par(mar=c(5,6.5,3,2))
    plot(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$mean_curve, 
        t='l', log='y', ylim=c(1e-04, 2e-02), xlim=c(0, 15), lwd=2, las=1, 
        xlab="", ylab="",
        cex.lab=1.4, cex.axis=1.4)
    mtext(side=1, "Inundation depth (m)", cex=1.4, line=2.5)
    mtext(side=2, "Exceedance-rate (events/year)", cex=1.4, line=5)
    title(main='Union-of-segments \n Epistemic Uncertainty', cex.main=1.5)
    points(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$rates_50,  t='l', col='brown', lwd=2)
    points(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$rates_84,  t='l', col='blue')
    points(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$rates_16,  t='l', col='blue')
    points(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$rates_975, t='l', col='green')
    points(segmented_model_comonotonic$peak_depths, segmented_model_comonotonic$rates_025, t='l', col='green')
    abline(v=seq(0, 20), col='grey', lty='dotted')
    abline(h=1/c(500, 2500, 10000), col='orange')
    legend('topright', c('Mean', 'Median', '68% CI', '95% CI'), 
           col=c('black', 'brown', 'blue', 'green'),
           lwd=c(2,2,1,1), pch=c(NA, NA, NA, NA), cex=1.3, bty='n')

    plot(unsegmented_model$peak_depths, unsegmented_model$mean_curve, 
        t='l', log='y', ylim=c(1e-04, 2e-02), xlim=c(0, 15), lwd=2, las=1, 
        xlab="", ylab="",
        cex.lab=1.4, cex.axis=1.4)
    mtext(side=1, "Inundation depth (m)", cex=1.4, line=2.5)
    mtext(side=2, "Exceedance-rate (events/year)", cex=1.4, line=5)
    title(main='Unsegmented \n Epistemic Uncertainty', cex.main=1.5)
    points(unsegmented_model$peak_depths, unsegmented_model$rates_50,  t='l', col='brown', lwd=2)
    points(unsegmented_model$peak_depths, unsegmented_model$rates_84,  t='l', col='blue')
    points(unsegmented_model$peak_depths, unsegmented_model$rates_16,  t='l', col='blue')
    points(unsegmented_model$peak_depths, unsegmented_model$rates_975, t='l', col='green')
    points(unsegmented_model$peak_depths, unsegmented_model$rates_025, t='l', col='green')
    abline(v=seq(0, 20), col='grey', lty='dotted')
    abline(h=1/c(500, 2500, 10000), col='orange')
    legend('topright', c('Mean', 'Median', '68% CI', '95% CI'), 
           col=c('black', 'brown', 'blue', 'green'),
           lwd=c(2,2,1,1), pch=c(NA, NA, NA, NA), cex=1.3, bty='n')

}

