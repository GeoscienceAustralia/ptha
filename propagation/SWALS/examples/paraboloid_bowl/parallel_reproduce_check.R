# Run this script after running 
# ./parallel_reproduce.sh
# To compare results computed with coarrays and openmp.

source('../../plot.R')

md_dirs = Sys.glob('OUTPUTS/RUN*')

# Get VH at a late time, in OMP/MPI runs
x = merge_domains_nc_grids(multidomain_dir=md_dirs[1], domain_index=1, desired_var='vh', desired_time_index=300)
xp = merge_domains_nc_grids(multidomain_dir=md_dirs[2], domain_index=1, desired_var='vh', desired_time_index=300)

# Should be exactly the same
run_range = range(x$vh - xp$vh)
if(all(run_range == 0)){
    print('PASS')
}else{
    print('FAIL')
}

