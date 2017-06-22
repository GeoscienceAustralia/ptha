#
# Script to 'qsub' additional stochastic_tsunami runs if they did not all finish
# when run from 'run_make_all_tsunami_events.PBS'. The latter script executes
# a number of jobs to do make the stochastic slip tsunami, in serial, and for
# large source-zones they will not finish in the 48 hours available.
#
# This script will probably only work on NCI [or systems with similar setup +
# using PBS for job submission].
#
# This script should be run after 'run_make_all_tsunami_events.sh'. If the
# latter script did finish, then this script will not do anything.
#
# For each 'batch' of events we tried to treat earlier, we check 'max_stage' in
# the netcdf file to see if the batch has missing data. If it does, we rerun
# the batch (on its own PBS job), and save as an RDS file [to avoid the
# possibility of parallel writes if running more than one job, which can cause
# netcdf to fail].
#
# Once all the above jobs have finished, we'll need to run another script to
# insert their data into the netcdf file
#
###############################################################################

#
# INPUTS: THESE MUST BE CONSISTENT WITH OTHER SCRIPTS FOR EVERYTHING TO WORK
# DO NOT CHANGE WITHOUT UNDERSTANDING THOSE LINKAGES [commented below]
#

# This number was used to determine the number of batches in
# run_make_all_tsunami_events.PBS (computation of nsplit)
batch_chunk_size = 4500 

# This number is greater than the number used as a missing data value in
# 'make_all_earthquake_tsunami.R'
# See nul_r = -999.999 in that script.
# Note -- to identify missing data, we check that "max_stage < missing_data_less_than".
# We do not try an exact check - because the netcdf data is stored as float,
# and will not be able to exactly represent -999.999. 
# -999 is well outside the range of max_stage vlaues, so it's no problem.
missing_data_less_than = -999.0 

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

#
# END INPUTS [do not change them without understanding cross-script linkages!]
#

#############################################################################

#
# MAIN CODE
#

library(ncdf4)

source_zone = basename(dirname(getwd()))

# Open (potentially incomplete) stochastic slip earthquake tsunami events file
nc_file = paste0('all_stochastic_slip_earthquake_events_tsunami_', source_zone, '.nc')
fid = nc_open(nc_file, readunlim=FALSE)

# Number of events
nevents = fid$var$event_index_string$varsize[2]

# How many batches did we originally split into? -- note the 'batch_chunk_size'
# must match the run_make_all_tsunami_events.sh script
nbatch = floor(nevents/batch_chunk_size + 1)

batch_inds = parallel::splitIndices(nevents, nbatch) # Same split as in make_all_earthquake_tsunami.R

# Find those batches that we missed [i.e. the job was killed before they finished ]
missed = c()
for(i in 1:length(batch_inds)){
    # Get any stage series in the i'th batch -- here just check the first entry
    # in the first row
    stg = ncvar_get(fid, 'max_stage', start=c(batch_inds[[i]][1],1), count=c(1,1))
    # Recall the missing data value was -999.999
    if(stg < missing_data_less_than | is.na(stg)) missed = c(missed, i)
}


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


