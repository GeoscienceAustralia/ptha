# Run all the validation tests
# 
# Before running, ensure OMP_NUM_THREADS is NOT set. If it is set, I can't use OMP 
# in the system calls inside R.
#
library(parallel)
basedir = getwd()
test_files = c(Sys.glob('../../examples/*/run_model.sh'), Sys.glob('../../examples/*/*/run_model.sh'))
#test_files = c(Sys.glob('../../examples/dispersive/*/run_model.sh'))

test_example<-function(test_file){

    t0 = Sys.time()

    # Move to the test_file directory, and run the job
    setwd(dirname(test_file))
    print('')
    print(paste0('Testing ', test_file))
    run_results = system(paste0('export OMP_NUM_THREADS=', detectCores(), ' ; bash ./run_model.sh'), intern=TRUE) 
    # Filter out system messages (e.g. often when R closes a graphics device)
    k = grep('[1]', run_results, fixed=TRUE)
    cat(run_results[k], sep="\n")
    setwd(basedir)

    # Record the time taken to run the job 
    t1 = Sys.time()
    time_taken = t1 - t0
    print(paste0('Time taken (build, run, plot) : ', format(time_taken)))

    return(list(result=run_results, elapsed_time=time_taken, test_file = test_file))
}

test_results = lapply(test_files, test_example)

test_results_filename = paste0('test_results_', format(Sys.time(), '%Y_%m_%d_%H_%M_%S'), '.RDS')
saveRDS(test_results, test_results_filename)

#
# If a 'reference test results' file exists, then use it to check for speed regressions (only)
#
reference_test_file = 'test_results_REFERENCE.RDS'
if(file.exists(reference_test_file)){

    test_referen = readRDS(reference_test_file)

    # Match model names     
    test_referen_model_names = unlist(lapply(test_referen, f<-function(x) x$test_file))
    test_results_model_names = unlist(lapply(test_results, f<-function(x) x$test_file))

    index_link = match(test_results_model_names, test_referen_model_names)

    # Get speed for corresponding tests
    test_speed_results = unlist(lapply(test_results, f<-function(x) as.double(x$elapsed_time, units='secs')))
    test_speed_referen = sapply(index_link, f<-function(x) as.double(test_referen[[x]]$elapsed_time, units='secs'))

    # Speed regressions: Allow 15% on any one test, and 10% overall?
    if(any(test_speed_results > 1.15 * test_speed_referen) | sum(test_speed_results) > 1.1*sum(test_speed_referen)){
        print('WARNING: Speed regressions')
    }else{
        print('PASS no speed regressions')
    }

    print('Total time:')
    print(sum(test_speed_results))
    print('Reference time:')
    print(sum(test_speed_referen))

}
