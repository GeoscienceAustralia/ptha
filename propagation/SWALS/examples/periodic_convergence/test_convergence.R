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

    if(all(order[2:(NM-1)] > 1.5)){
        print('PASS')
    }else{
        print('FAIL: Convergence less rapid than previously')
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


## > # Second order accuracy with rk2
## > get_order_accuracy('stage', time=0.05)
##               dx         errs    order
## [1,] 0.009900990 0.0308245329       NA
## [2,] 0.004975124 0.0091937866 1.745348
## [3,] 0.002493766 0.0021909298 2.069116
## [4,] 0.001248439 0.0004498888 2.283903
## > get_order_accuracy('ud', time=0.05)
##               dx        errs    order
## [1,] 0.009900990 0.272771857       NA
## [2,] 0.004975124 0.081466958 1.743408
## [3,] 0.002493766 0.020086947 2.019957
## [4,] 0.001248439 0.004090322 2.295972
## > get_order_accuracy('vd', time=0.05)
##               dx         errs    order
## [1,] 0.009900990 0.0610027632       NA
## [2,] 0.004975124 0.0180286382 1.758584
## [3,] 0.002493766 0.0043889836 2.038332
## [4,] 0.001248439 0.0008708502 2.333390
## > 
## > # Clearly 2nd order accuracy
## > get_order_accuracy('stage', time=0.055)
##               dx         errs    order
## [1,] 0.009900990 0.0916853064       NA
## [2,] 0.004975124 0.0146657571 2.644239
## [3,] 0.002493766 0.0038584719 1.926350
## [4,] 0.001248439 0.0007738359 2.317930
## > get_order_accuracy('ud', time=0.055)
##               dx        errs    order
## [1,] 0.009900990 0.838361174       NA
## [2,] 0.004975124 0.133864304 2.646801
## [3,] 0.002493766 0.035777026 1.903666
## [4,] 0.001248439 0.007132469 2.326560
## > get_order_accuracy('vd', time=0.055)
##               dx        errs    order
## [1,] 0.009900990 0.144194847       NA
## [2,] 0.004975124 0.024678463 2.546695
## [3,] 0.002493766 0.005916193 2.060512
## [4,] 0.001248439 0.001311882 2.173031
## > 
