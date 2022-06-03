#
# Code to convergence test model runs
#
# I assume OUTPUTS contains only model runs made with:
#   source run_convergence_test.sh 
# 
#

# Number of models that were run
NM = 4

source('../../plot.R')

all_md = Sys.glob('OUTPUTS/RUN*/RUN*')

xs = lapply(as.list(all_md), get_all_recent_results)

x_values = lapply(xs, f<-function(x) x$xs)
y_values = lapply(xs, f<-function(x) x$ys)

# Find the lowest-res domain
lowres_ind = which.max(unlist(lapply(xs, f<-function(x) x$dx[1])))

# Interpolate results at x-values that are on the lowest-res domain, and not in
# the periodic-halo
desired_x = x_values[[lowres_ind]]
desired_x = desired_x[which(desired_x > 0 & desired_x < 1)]
desired_y = y_values[[lowres_ind]]
desired_y = desired_y[which(desired_y > 0 & desired_y < 1)]

desired_xy = expand.grid(desired_x, desired_y)

# Interpolation function
get_values_at_desired_xy<-function(domain, var='stage', time=0.05){
    suppressPackageStartupMessages(library(fields))
  
    desired_time = which.min(abs(domain$time - time)) 

    vals = interp.surface(list(x=domain$xs, y=domain$ys, z=domain[[var]][,,desired_time]), loc=as.matrix(desired_xy))

    out = matrix(vals, ncol=length(desired_x), byrow=TRUE)

    return(out)
}

get_order_accuracy<-function(var, time=0.05){
    all_stages = lapply(xs, f<-function(x) get_values_at_desired_xy(x, var=var, time=time))
    errs = unlist(lapply(all_stages[1:(NM-1)], f<-function(x) max(abs(x - all_stages[[NM]]))))
    dx = unlist(lapply(xs[1:(NM-1)], f<-function(x) x$dx[1]))

    # Order of accuracy ~ log2(coarser_error/finer_error)
    order = c(NA, log(errs[-(NM-1)]/errs[-1], 2))

    if( all(order[2:(NM-1)] > 1.5) & (mean(order[2:(NM-1)]) > 1.8) ){
        print('PASS')
    }else{
        print('FAIL: Convergence less rapid than previously')
        print(order[2:(NM-1)])
    }
    
    cbind(dx, errs, order) 
}

# Second order accuracy with rk2
tmp = get_order_accuracy('stage', time=0.05)
tmp = get_order_accuracy('ud', time=0.05)
tmp = get_order_accuracy('vd', time=0.05)

# Clearly 2nd order accuracy
tmp = get_order_accuracy('stage', time=0.056)
tmp = get_order_accuracy('ud', time=0.056)
tmp = get_order_accuracy('vd', time=0.056)
