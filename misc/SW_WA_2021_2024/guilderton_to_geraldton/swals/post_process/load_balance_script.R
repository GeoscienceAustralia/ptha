#
# Run this from inside a multidomain directory to compute a load_balance file,
# and report on expected time-savings. If the results are good you can them manually
# copy that file to 
#

mydir = getwd()
if(!grepl('RUN_', basename(mydir))){
    stop('This script must be run from inside a multidomain directory')
}

source('/g/data/w85/tsunami/CODE/gadi/ptha_mm/propagation/SWALS/plot.R')
# Assume we want every processor to have part of domain 1. This makes sense
# if domain 1 is a global domain, and the other domains may utilise static_before_time.
x = make_load_balance_partition('.', domain_index_groups=list(1, 2:9999))

