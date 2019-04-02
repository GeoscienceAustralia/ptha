# Run all the validation tests
# 
# Before running, ensure OMP_NUM_THREADS is NOT set. If it is set, I can't use OMP 
# in the system calls inside R.
#

basedir = getwd()
test_files = c(Sys.glob('../../examples/*/run_model.sh'), Sys.glob('../../examples/*/*/run_model.sh'))

test_example<-function(test_file){

    t0 = Sys.time()

    # Move to the test_file directory, and run the job
    setwd(dirname(test_file))
    print('')
    print(paste0('Testing ', test_file))
    run_results = system('bash ./run_model.sh', intern=TRUE) 
    # Filter out system messages (e.g. often when R closes a graphics device)
    k = grep('[1]', run_results, fixed=TRUE)
    cat(run_results[k], sep="\n")
    setwd(basedir)

    # Record the time taken to run the job 
    t1 = Sys.time()
    time_taken = t1 - t0
    print(paste0('Time taken (build, run, plot) : ', format(time_taken)))

    return(list(result=run_results, elapsed_time=time_taken, file=test_file))
}

test_results = lapply(test_files, test_example)

test_results_filename = paste0('test_results_', format(Sys.time(), '%Y_%m_%d_%H_%M_%S'), '.RDS')
saveRDS(test_results, test_results_filename)
