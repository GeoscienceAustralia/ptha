#
# Run from inside the multidomain directory of interest, e.g.
#   Rscript ../../../post_process/report_domain_run_times.R
#

file_ptha_plot = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/propagation/SWALS/plot.R'
source(file_ptha_plot)

all_logs = Sys.glob('multi*.log')

domain_run_times = lapply(all_logs, function(x) get_domain_wallclock_times_in_log(x))

all_domain_run_times = do.call(rbind, domain_run_times)

print('Stem and leaf plot of domain run times')
stem(all_domain_run_times$time)

top_few = order(all_domain_run_times$time, decreasing=TRUE)
if(length(top_few) > 10) top_few = top_few[1:10]
print('Here are the slowest domains:')
print(all_domain_run_times[top_few,])
