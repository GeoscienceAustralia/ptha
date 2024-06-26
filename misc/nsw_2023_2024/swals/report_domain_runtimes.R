#
# Run from inside the multidomain directory of interest, e.g.
#   Rscript ../../../report_domain_run_times.R
#

file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

all_logs = Sys.glob('multi*.log')

domain_run_times = lapply(all_logs, function(x) get_domain_wallclock_times_in_log(x, wallclock_time_line_spacing=17))

all_domain_run_times = do.call(rbind, domain_run_times)

print('Stem and leaf plot of domain run times')
stem(all_domain_run_times$time)

top_few = order(all_domain_run_times$time, decreasing=TRUE)
if(length(top_few) > 10) top_few = top_few[1:10]
print('Here are the slowest domains:')
print(all_domain_run_times[top_few,])
