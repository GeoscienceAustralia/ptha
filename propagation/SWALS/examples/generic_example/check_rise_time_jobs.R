source('../../plot.R')


# Assume there are 2 jobs only here
model_runs = Sys.glob('OUTPUTS/test_tohoku/RUN*') 
if(length(model_runs) != 2) stop('Did not find exactly 2 model runs')

# The first job has a zero rise-time, the second has a rise time of 1000 seconds
# Theoretically (for linear equations), the stage-time-series of the second job
# should equal a time-average of the first (over the rise_time window)
# For nonlinear equations this will be close to true as well, but perhaps less exact.
# Allow a higher error threshold for the nonlinear equations,
solver_type = commandArgs(trailingOnly=TRUE)
if(length(solver_type) == 0){
    print("FAIL (you didn't pass the solver_type to check_rise_time_jobs.R)")
    stop('Deliberate halt')
}
if(solver_type == 'linear'){
    tol = 1e-04
    delay_forcing_time = 0.0
}else if(solver_type == 'rk2'){
    tol = 1e-03
    delay_forcing_time = 30.0
}else{
    stop('FAIL -- unknown solver type')
}

x = lapply(model_runs, get_all_recent_results)

approx_f = approxfun(x[[1]]$gauges$time, x[[1]]$gauges$time_var$stage[10,], rule=2)
rise_time = 1000
smoother = sapply(x[[2]]$gauges$time, f<-function(t) mean(approx_f(t + seq(-rise_time, 0) - delay_forcing_time)) )

#plot(x[[1]]$gauges$time, x[[1]]$gauges$time_var$stage[10,],t='l')
#points(x[[2]]$gauges$time, x[[2]]$gauges$time_var$stage[10,],t='l', col='red')
#points(x[[2]]$gauges$time, smoother, t='l', col='green')



# Difference between model with rise time, vs theoretical results
err = mean(abs(x[[2]]$gauges$time_var$stage[10,] - smoother))
if(err < tol){
    print('PASS')
}else{
    print('FAIL ')
}
