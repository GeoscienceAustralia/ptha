source('../config.R') # Get logical variable NEARSHORE_TESTING_PAPER_2020_ONLY
  
ALL_INVERSIONS_FOLDERS = Sys.glob('../OUTPUTS/*/RUN*')

#for(i in 1:length(all_INVERSIONS_folders)){
parallel_job<-function(i){
    ALL_INVERSIONS_FOLDERS = Sys.glob('../OUTPUTS/*/RUN*')
    run_folder = ALL_INVERSIONS_FOLDERS[i]
    # Get the text denoting the event, which comes just after INVERSIONS_
    event_flag = strsplit(basename(dirname(run_folder)), '_')[[1]][1]


    if(!NEARSHORE_TESTING_PAPER_2020_ONLY){
        # Include events where we've run PTHA18 scenarios, but not yet inversions
        event2script = list('Chile1960'   = 'plot_southamerica1960.R',
                            'Chile2010'   = 'plot_southamerica2010.R',
                            'Chile2015'   = 'plot_southamerica2015.R',
                            'Chile2014'   = 'plot_southamerica2014.R',
                            'Sumatra2004' = 'plot_sumatra2004.R',
                            'Tohoku2011'  = 'plot_tohoku2011.R',
                            'Puysegur2009' = 'plot_puysegur2009.R',
                            'Solomon2007' = 'plot_solomon2007.R')
    }else{
        # Include only the published inversions, for the 2020 nearshore testing paper
        event2script = list('Chile1960'   = 'plot_southamerica1960.R',
                            'Chile2010'   = 'plot_southamerica2010.R',
                            'Chile2015'   = 'plot_southamerica2015.R',
                            'Sumatra2004' = 'plot_sumatra2004.R',
                            'Tohoku2011'  = 'plot_tohoku2011.R')
    }

    job_command = paste0('Rscript ', event2script[[event_flag]], ' ', run_folder)
    print(job_command)
    system(job_command)
}

library(parallel)
mclapply(1:length(ALL_INVERSIONS_FOLDERS), parallel_job, mc.cores=24, mc.preschedule=FALSE)

