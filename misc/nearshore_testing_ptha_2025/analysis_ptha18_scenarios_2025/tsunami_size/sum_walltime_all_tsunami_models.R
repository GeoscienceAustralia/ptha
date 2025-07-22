#
# Sum the model times in all runs
#

all_runs = Sys.glob('../REDUCED_OUTPUTS/*/RUN*/*.log.zip')

readfun<-function(file_log_zip){
    x = readLines(unzip(file_log_zip))
    time_line = x[length(x)-1]
    numeric_time = as.numeric(gsub("Total WALLCLOCK time: ", "", time_line))
    return(numeric_time)
}
try_readfun<-function(x) try(readfun(x))
library(parallel)
all_times = mclapply(all_runs, readfun, mc.cores=16)

Ncores = 384
compute_time_core_hours = sum(unlist(all_times))/3600 * Ncores
#> compute_time_core_hours/1e+03
#[1] 1512.913

