#
# Code to run */TSUNAMI_UNIT_SOURCE/check_runs_complete.R on a user-defined set of source-zones.
#

#source('TEMPLATE/TSUNAMI_UNIT_SOURCE/check_runs_complete.R')

source_zones = c('izumariana', 'kermadectonga', 'kurilsjapan', 'newhebrides',
    'puysegur', 'solomon', 'sunda', 'southamerica')

# DEBUG
#source_zones = source_zones[-1]

mydir = getwd()
for(sz in source_zones){
    setwd(paste0(sz, '/TSUNAMI_UNIT_SOURCE'))
    source('check_runs_complete.R')
    print(paste0('Checking ', sz, ' .....'))
    try(check_models_have_been_run(sz, verbose=TRUE))
    try(check_model_gauge_integrity(sz, verbose=TRUE))
    setwd(mydir)
}
