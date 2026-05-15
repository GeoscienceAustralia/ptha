#
# Run from inside the multidomain directory of interest, e.g.:
#   Rscript ../../../post_process/report_domain_run_times.R
# By default it reports the top 10 domains - to report the top 
# 25 (or whatever) domains, provide the number:
#   Rscript ../../../post_process/report_domain_run_times.R 25
#

file_ptha_plot = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(file_ptha_plot)

all_logs = Sys.glob('multi*.log')

domain_run_times = lapply(all_logs, function(x) get_domain_wallclock_times_in_log(x))

all_domain_run_times = do.call(rbind, domain_run_times)

print('Stem and leaf plot of domain run times')
stem(all_domain_run_times$time)

top_few = order(all_domain_run_times$time, decreasing=TRUE)
top_N = 10
# Optionally set top_N from command line
ca = commandArgs(trailingOnly=TRUE)
if(length(ca) > 0) top_N = as.numeric(ca[1])

if(length(top_few) > top_N) top_few = top_few[1:top_N]
print('Here are the slowest domains:')
print(all_domain_run_times[top_few,])
