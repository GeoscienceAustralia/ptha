get_failed_dirs<-function(){
    run_failed_info = readLines('check_log_files_plots/ptha_runs_December_2025.txt')
    k = grep('Try error', run_failed_info)
    run_failed_info = run_failed_info[k]
    failed_dirs = unlist(lapply(run_failed_info, function(x) dirname(strsplit(x, '/OUTPUTS/')[[1]][2])))
    return(failed_dirs)
}
failed_dirs_full = get_failed_dirs()
failed_dirs = unlist(lapply(failed_dirs_full, function(x) strsplit(x, '-full')[[1]][1]))

run_commands = do.call(c, lapply(Sys.glob('../run_scripts_haza*/ptha*.sh'), readLines))
#print(run_commands)

job_run_commands = vector(mode='list', length=length(failed_dirs))
for(i in 1:length(failed_dirs)){
    k = grep(failed_dirs[i], run_commands, fixed=TRUE)
    #print(c(failed_dirs[i], k))
    stopifnot(length(k) == 2)
    inds = (k[1]-3):(k[1]+4)
    #print('')
    #print(run_commands[inds])
    #print('')
    job_run_commands[[i]] = c(
        '#!/bin/bash',
        '#PBS -P w85',
        '#PBS -q normalsr',
        '#PBS -l walltime=24:00:00',
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
    job_run_commands[[i]] = gsub('multidomain_kalbarri2onslow_20251218_hazard.nml', 'multidomain_kalbarri2onslow_20251218_hazard_lowts.nml', job_run_commands[[i]])

    cat(job_run_commands[[i]], file=paste0('run_failed_job_', 1000 + i, '.sh'), sep="\n")
}
