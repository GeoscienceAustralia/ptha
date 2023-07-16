ALL_INVERSIONS_FOLDERS = Sys.glob('../OUTPUTS/*/RUN*')

#for(i in 1:length(all_INVERSIONS_folders)){
parallel_job<-function(i){
    ALL_INVERSIONS_FOLDERS = Sys.glob('../OUTPUTS/*/RUN*')
    run_folder = ALL_INVERSIONS_FOLDERS[i]
    # Get the text denoting the event, which comes just after INVERSIONS_
    event_flag = strsplit(basename(dirname(run_folder)), '_')[[1]][2]

    # Map between 
    event2script = list('Chile1960'   = 'plot_southamerica1960.R',
                        'Chile2010'   = 'plot_southamerica2010.R', 
                        'Chile2015'   = 'plot_southamerica2015.R', 
                        'andaman2004' = 'plot_sumatra2004.R', 
                        'Tohoku2011'  = 'plot_tohoku2011.R')

    job_command = paste0('Rscript ', event2script[[event_flag]], ' ', run_folder)
    print(job_command)
    system(job_command)
}

library(parallel)
mclapply(1:length(ALL_INVERSIONS_FOLDERS), parallel_job, mc.cores=24, mc.preschedule=FALSE)
