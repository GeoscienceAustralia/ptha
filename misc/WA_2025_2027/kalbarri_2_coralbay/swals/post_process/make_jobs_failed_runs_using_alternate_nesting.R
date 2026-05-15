#
# Create qsub scripts to run jobs that failed
# This will need hacking in general.
#

fix_type = 'alternate_nesting_with_debug'


get_failed_dirs<-function(){
    run_failed_info = readLines('Energy_stats_kalbarri2coralbay_first_runs_with_some_failures.txt')
    k = grep('Try error', run_failed_info)
    run_failed_info = run_failed_info[k]
    failed_dirs = unlist(lapply(run_failed_info, function(x) dirname(strsplit(x, '/OUTPUTS/')[[1]][2])))
    return(failed_dirs)
}
failed_dirs_full = get_failed_dirs()
failed_dirs = unlist(lapply(failed_dirs_full, function(x) strsplit(x, '-full')[[1]][1]))

## Read the run scripts that were previously submitted. They will be used to determine the run commands
run_commands = do.call(c, lapply(Sys.glob('../submitted_kalbarri2coralbay_hazard/ptha*.sh'), readLines))

job_run_commands = vector(mode='list', length=length(failed_dirs))
for(i in 1:length(failed_dirs)){
    k = grep(failed_dirs[i], run_commands, fixed=TRUE)
    stopifnot(length(k) == 2)
    inds = (k[1]-3):(k[1]+4)

    # Compile the model to use the alternate nesting scheme, and run in 'debug' mode so that it will report the locations of any failures.
    # This is only a bit slower than a regular run, and if it doesn't work, at least we get more info on the error.
    job_run_commands[[i]] = c(
        '#!/bin/bash',
        '#PBS -P w85',
        '#PBS -q normalsr',
        '#PBS -l walltime=8:00:00',
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

    job_run_commands[[i]] = gsub('model_build', 'model_both_debug_and_old_nesting', job_run_commands[[i]])

    cat(job_run_commands[[i]], file=paste0('run_failed_job_debug_and_old_nesting_', 1000 + i, '.sh'), sep="\n")

}
