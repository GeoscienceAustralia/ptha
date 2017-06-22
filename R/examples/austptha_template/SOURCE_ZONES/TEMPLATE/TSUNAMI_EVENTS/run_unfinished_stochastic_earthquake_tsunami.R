#
# qsub additional stochastic_tsunami runs if they did not all finish 
#
# This is run after run_make_all_tsunami_events.sh, in case that did not finish
#

library(ncdf4)

source_zone = basename(dirname(getwd()))

# Open (potentially incomplete) stochastic slip earthquake events file
nc_file = paste0('all_stochastic_slip_earthquake_events_tsunami_', source_zone, '.nc')
fid = nc_open(nc_file, readunlim=FALSE)

# Number of events
nevents = fid$var$event_index_string$varsize[2]

# How many batches did we originally split into -- note the '4500' must match the 
# run_make_all_tsunami_events.sh script
nbatch = floor(nevents/4500 + 1)

batch_inds = parallel::splitIndices(nevents, nbatch)

# Find those batches that we missed [i.e. the job was killed before they finished ]
missed = c()
for(i in 1:length(batch_inds)){
    # Get any stage series in the i'th batch -- here just take the first one in
    # batch_inds
    stg = ncvar_get(fid, 'max_stage', start=c(batch_inds[[i]][1],1), count=c(1,1))
    # Recall the missing data value was -999.999
    if(stg < -999| is.na(stg)) missed = c(missed, i)
}


#
# Template PBS script text -- note REPLACEWITHMID in final line, which will be auto-replaced
# with 'missed' integers
#
pbs_text = "#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=32GB
#PBS -lncpus=16
#PBS -l wd

# Source key R modules -- not that you will need the right packages installed
# as well (see comments in the script that is sourced)
# NOTE THIS IS ONLY FOR NCI, COMMENT OUT OTHERWISE
source R_modules.sh

nevents=$( ncdump -h all_stochastic_slip_earthquake_events_*.nc | grep 'table_rows = UNLI' | awk '{print $6}' | tr '(' ' ' )
nsplit=$( expr $nevents / 4500 + 1 )

Rscript make_all_earthquake_tsunami.R --stochastic_slip --subset REPLACEWITHMYID $nsplit --save_as_RDS
"

# Make scripts for each unfinished job, and run them
if(length(missed) > 0){
    for(ii in missed){
        pbs_local = gsub('REPLACEWITHMYID', ii, pbs_text)
        sub_script_name = paste0('temp_qsub_stochastic_extra_', ii, '.PBS')
        cat(pbs_local, file=sub_script_name)
        qsub_command = paste0('qsub ', sub_script_name)
        system(qsub_command)
    }
}


