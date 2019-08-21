#
# This script assumes:
#    'revised_station_hazard_curves_PREPROCESSING.R' and 
#    'submit_all_PBS_revised_station_hazard_curves_FINAL.R' and 'revised_station_hazard_curves_FINAL_MERGE.R'
# have already been run.
#
# The script here was originally developed to process results where the
# logic-trees were sampled (3000 branches) to reduce the computational effort
# (i.e.  revised_station_hazard_curves.R) -- so those results had some (small)
# monte-carlo error. But at 99% of stations the mean ARI=500 max-stage changed
# by < 1%. 
#


# These files have the calculations for all source-zones. They were done in
# parallel (chunks of approx. 5 points at a time), and the result is a list of
# length 20185/5 = 4037. Note sometimes 4 or 6 points were done, because the
# function 'splitIndices' does not produce an even partition in this case, even
# though it could!
# all_stage_calc_files = Sys.glob('./preprocessed_source_rate_revised_stage_exrates/stage_exceedance_rate_percentiles_*.RDS')
#
# These files were created using the 'revised_station_hazard_curves_FINAL_MERGE.R'
all_stage_calc_files = Sys.glob('./preprocessed_source_rate_revised_stage_exrates_FULL_MERGED/stage_exceedance_rate_percentiles_*.RDS')
all_stage_calcs = lapply(all_stage_calc_files, f<-function(x) readRDS(x))

# These files have the Mw-rate curves, including the subset of 3000 randomly sampled curves.
# all_sampled_rate_calc_files = Sys.glob('./preprocessed_source_rate_revised_stage_exrates/stage_exceedance_rate_SEG_*.RDS')

# The stage_seq (stages at which rates are stored) is stored in files like this, and is the same everywhere 
one_rate_calc = readRDS('./preprocessed_source_rate_revised_stage_percentiles/preprocessed_rate_info_puysegur2.RDS')
STAGE_SEQ = one_rate_calc$stage_seq

get_mean_exrates<-function(one_stage_calc, slip_type = 'stochastic', mu_type = 'fixed_mu'){
    all_mean_exrates = lapply(one_stage_calc, 
        f<-function(x) x[[slip_type]][[mu_type]]$stage_mean_exrate)
    mean_exrates = do.call(cbind, all_mean_exrates)
    return(mean_exrates)
}

get_percentile_exrates<-function(one_stage_calc, pc, slip_type = 'stochastic', mu_type = 'fixed_mu'){
    row = switch(pc, '0.025' = 1, '0.16' = 2, '0.5' = 3, '0.84' = 4, '0.975' = 5)
    all_pc_exrates = lapply(one_stage_calc, 
        f<-function(x) x[[slip_type]][[mu_type]]$stage_percentile_exrates[row,,])
    pc_exrates = do.call(cbind, all_pc_exrates)
    return(pc_exrates)
}


get_ari500_stages<-function(mean_exrates, stage_seq=STAGE_SEQ, target_rate=1/500){

    # Interpolate the ari500 stage. Note it is possible for all exrates to be
    # 0, which will break the interpolation. So we add a 'cap' with a negative
    # rate + impossibly large stage. Combined with 'ties='min'' and the fact
    # that STAGE_SEQ extends to 20, this cap will never be used, but allows the
    # interpolator to not fail
    ari500_stages = rep(0, ncol(mean_exrates))
    for(i in 1:length(ari500_stages)){
        if(all(is.na(mean_exrates[,i]))){
            ari500_stages[i] = NA
        }else{
            ari500_stages[i] = approx(c(mean_exrates[,i], -1), c(STAGE_SEQ, 1e+06), 
                                      xout=target_rate, rule=2, ties='min')$y
        }
    }
    return(ari500_stages)
}

# Get the sum of the mean exceedance rates for the new calculations
mean_exrates_summed = get_mean_exrates(all_stage_calcs[[1]])
for(i in 2:length(all_stage_calcs)){
    mean_exrates_summed = mean_exrates_summed + get_mean_exrates(all_stage_calcs[[i]])
}
ari500_stages = get_ari500_stages(mean_exrates_summed)

# As above for the 16th and 84th percentiles
# Co-monotonic
pc025_exrates_como = get_percentile_exrates(all_stage_calcs[[1]], '0.025')
pc16_exrates_como  = get_percentile_exrates(all_stage_calcs[[1]], '0.16')
pc50_exrates_como  = get_percentile_exrates(all_stage_calcs[[1]], '0.5')
pc84_exrates_como  = get_percentile_exrates(all_stage_calcs[[1]], '0.84')
pc975_exrates_como = get_percentile_exrates(all_stage_calcs[[1]], '0.975')
for(i in 2:length(all_stage_calcs)){
    pc025_exrates_como = pc025_exrates_como + get_percentile_exrates(all_stage_calcs[[i]], '0.025')
    pc16_exrates_como  = pc16_exrates_como  + get_percentile_exrates(all_stage_calcs[[i]], '0.16')
    pc50_exrates_como  = pc50_exrates_como  + get_percentile_exrates(all_stage_calcs[[i]], '0.5')
    pc84_exrates_como  = pc84_exrates_como  + get_percentile_exrates(all_stage_calcs[[i]], '0.84')
    pc975_exrates_como = pc975_exrates_como + get_percentile_exrates(all_stage_calcs[[i]], '0.975')
}

ari500_stages_025 = get_ari500_stages(pc025_exrates_como)
ari500_stages_16  = get_ari500_stages(pc16_exrates_como)
ari500_stages_50  = get_ari500_stages(pc50_exrates_como)
ari500_stages_84  = get_ari500_stages(pc84_exrates_como)
ari500_stages_975 = get_ari500_stages(pc975_exrates_como)

output = data.frame(
    revised_ari500_stage_lower=ari500_stages_025,
    revised_ari500_stage_16pc=ari500_stages_16,
    revised_ari500_stage_median=ari500_stages_50,
    revised_ari500_stage_mean=ari500_stages,
    revised_ari500_stage_84pc=ari500_stages_84,
    revised_ari500_stage_upper=ari500_stages_975)
write.csv(output, file='revised_ari500_stages.csv', row.names=FALSE)

#
# Original results (PTHA18 GA Record)
#
library(rptha)
#fid = nc_open('../../FIG/spatial_hazard/tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc')
fid = nc_open('./tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc')
# Hazard point data
pt_var = c('lon', 'lat', 'elev', 'gaugeID')
hp = lapply(pt_var, f<-function(x) as.numeric(ncvar_get(fid, x)))
names(hp) = pt_var
hp = as.data.frame(hp)
# Get the equivalent of mean_exrates_summed, with the old results
ssr   = ncvar_get(fid, 'stochastic_slip_rate')
ssr16 = ncvar_get(fid, 'stochastic_slip_rate_16pc')
ssr84 = ncvar_get(fid, 'stochastic_slip_rate_84pc')
ssr50 = ncvar_get(fid, 'stochastic_slip_rate_median')
nc_close(fid)

old_ari500_stages    = get_ari500_stages(ssr)
old_ari500_stages_16 = get_ari500_stages(ssr16)
old_ari500_stages_84 = get_ari500_stages(ssr84)
old_ari500_stages_50 = get_ari500_stages(ssr50)

## The results reported below were obtained using the "logic-tree-sampling" approach for
## computational efficiency (i.e. revised_station_hazard_curves.R), and so have slightly more
## variability than the 'full' calculations.
##     
##     
##     ## The ari500 stage based on the mean exrate should be "the same" in the new
##     ## and old results, except for the monte-carlo sampling error (due to our use
##     ## of sampled logic-trees to reduce computational effort)
##     ## Indeed, the differences are very small -- 99% of values differ by less than 1%
##     summary(ari500_stages/old_ari500_stages)
##     #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     # 0.9516  0.9993  1.0014  1.0011  1.0035  1.0293      15 
##     quantile(ari500_stages/old_ari500_stages, p=c(0.005, 0.995), na.rm=TRUE)
##     #     0.5%     99.5% 
##     #0.9882142 1.0086763 
##     mean(abs(ari500_stages/old_ari500_stages - 1) < 0.01, na.rm=TRUE)
##     #[1] 0.9909767
##     
##     #
##     # Percentiles can differ, because we are using a different method to compute them
##     # The effect seems small, especially at higher percentiles.
##     #
##     summary(ari500_stages_16/old_ari500_stages_16)
##     #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     # 0.5956  1.0451  1.0806  1.0470  1.1029  1.4171      15 
##     summary(ari500_stages_50/old_ari500_stages_50)
##     #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     # 0.9391  1.0223  1.0324  1.0383  1.0493  1.1658      15 
##     summary(ari500_stages_84/old_ari500_stages_84)
##     #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     # 0.8878  0.9935  1.0056  1.0058  1.0176  1.1700      15 
##     
