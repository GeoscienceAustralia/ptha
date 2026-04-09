#
# Create qsub scripts to run jobs that failed
# This will need hacking in general.
#

# After the first fix with alternate nesting, the failed dirs did not show up as plotting errors, so I manually added them here
fix_type = 'low_ts'
failed_dirs_full = c(
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0107982_Mw_94_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0109377_Mw_95_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0109396_Mw_95_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0109584_Mw_95_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0110660_Mw_96_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0110734_Mw_96_HS-full-ambient_sea_level_0",
    "ptha18-kalbarri2coralbay-hazard/random_sunda2/ptha18_random_scenarios_sunda2_row_0110754_Mw_96_HS-full-ambient_sea_level_0")

failed_dirs = unlist(lapply(failed_dirs_full, function(x) strsplit(x, '-full')[[1]][1]))

## Read the run scripts that were previously submitted. They will be used to determine the run commands
run_commands = do.call(c, lapply(Sys.glob('../submitted_kalbarri2coralbay_hazard/ptha*.sh'), readLines))

job_run_commands = vector(mode='list', length=length(failed_dirs))
for(i in 1:length(failed_dirs)){
    k = grep(failed_dirs[i], run_commands, fixed=TRUE)
    stopifnot(length(k) == 2)
    inds = (k[1]-3):(k[1]+4)

    # Models with low timestep
    job_run_commands[[i]] = c(
        '#!/bin/bash',
        '#PBS -P w85',
        '#PBS -q normalsr',
        '#PBS -l walltime=17:00:00', ## These ones need a longer runtime because the timestep is reduced
        '#PBS -lmem=1000GB',
        '#PBS -lncpus=208',
        '#PBS -l wd',
        '#PBS -l storage=gdata/w85+scratch/w85',
        '',
        'source SWALS_ifx_modules_2025_llvm.sh',
        '# Load R as well (just for tarring directories)',
        'module load R/4.3.1',
        paste0('rm OUTPUTS/', failed_dirs_full[i], '/RUN*.tar'),
        paste0('rm OUTPUTS/', failed_dirs_full[i], '/multi*.log'),
        run_commands[inds]
    )

    # An alternate namelist makes the models run with reduced timestep.
    job_run_commands[[i]] = gsub('multidomain_kalbarri2coralbay_B_hazard.nml', 'multidomain_kalbarri2coralbay_B_hazard_lowts.nml', job_run_commands[[i]])

    cat(job_run_commands[[i]], file=paste0('run_failed_job_low_ts_', 1000 + i, '.sh'), sep="\n")

}
