PBS_text = "#!/bin/bash
#PBS -P w85
#PBS -q normal 
#PBS -l walltime=30:00:00
#PBS -lmem=32GB
#PBS -lncpus=16
#PBS -l wd

source ~/R_351_NCI_modules.sh

Rscript revised_station_hazard_curves_FINAL.R STARTINDEX ENDINDEX
"

make_PBS_text<-function(STARTINDEX, ENDINDEX){

    new_text = gsub('STARTINDEX', as.character(STARTINDEX), PBS_text, fixed=TRUE)
    new_text = gsub('ENDINDEX', as.character(ENDINDEX), new_text, fixed=TRUE)

    return(new_text)
}

# Make a whole bunch of jobs, each having a subset of the points
number_of_jobs = 50
number_of_points = 20185
library(parallel)
point_ranges = splitIndices(number_of_points, number_of_jobs)

for(i in 1:length(point_ranges)){
    lower = min(point_ranges[[i]])
    upper = max(point_ranges[[i]])
    chunk_PBS_text = make_PBS_text(lower, upper)
    output_file = paste0('run_revised_station_hazard_curves_FINAL_', lower, '_', upper, '.PBS')
    cat(chunk_PBS_text, file=output_file)
    system(paste0("qsub ", output_file))
}
