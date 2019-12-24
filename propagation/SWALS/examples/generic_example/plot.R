#
# Read the darts and make a plot
#
source('../../plot.R')

model_type = commandArgs(trailingOnly=TRUE)[1]

# Model data
x = get_all_recent_results(rev(Sys.glob('OUTPUTS/test_tohoku/RUN_*'))[1])

# Dart data into a list
dart_data_files = Sys.glob('dart_test_data/[0-9]*.csv')
dart_data = lapply(dart_data_files, f<-function(x) read.csv(x))
names(dart_data) = gsub('.csv', '', basename(dart_data_files))

# Plot each dart as a panel
png(paste0('test_tohoku_DARTs_', model_type, '.png'), width=4, height=8, units='in', res=200)
par(mfrow=c(length(dart_data), 1))
par(mar=c(2,2,2,2))
gof_stats = rep(NA, length(dart_data))
model_var_stat = rep(NA, length(dart_data)) # Use this to check for changes in results
for(i in 1:length(dart_data)){
    dart_id = as.numeric(names(dart_data)[i])

    gauge_ind = match(dart_id, x$gauges$gaugeID)
    if(length(gauge_ind) != 1) stop('Could not fine exactly one gauge matching the dart')

    plot(dart_data[[i]][,1:2], t='l')
    points(x$gauges$time, x$gauges$time_var$stage[gauge_ind, ], t='l', col='red')
    title(dart_id, line=-1, cex=2, adj=1)

    model_fun = approxfun(x$gauges$time, x$gauges$time_var$stage[gauge_ind,])

    model_at_obs_times = model_fun(dart_data[[i]][,1])

    # Use this only to check for changes in results. 
    model_var_stat[i] = var(model_at_obs_times)

    # Make some kind of crude GOF statistic, which might detect catestrophic failure
    # Skip the first 1000 seconds for comparison, to filter the seismic waves
    zero_filter = (dart_data[[i]][,1] > 1000)
    model_range = range(model_at_obs_times * zero_filter)
    data_range = range(dart_data[[i]][,2] * zero_filter)
    range_ratio_log = log10(diff(model_range)/diff(data_range))
    gof_stats[[i]] = range_ratio_log
}
dev.off()

options(digits=12)

# Compare model_var_stat to the values from previous runs, to flag changes
expected_model_var_stats = list(
    "linear" =  c(0.0117662541101088, 0.0170184727742347, 0.00203174633168809, 
                  0.00198896568619196, 0.0319943378663546, 0.00788146806215795),

    "linear_with_nonlinear_friction" = c(0.0109163588083476, 0.0160886427199833, 
                                         0.00195355853476261, 0.00192765017953425, 
                                         0.0296076938837027, 0.00730800462575785), 

    "almost_linear" = c(0.0121756007780842, 0.0166737089216131, 0.0019850065701512, 
                        0.00204578160004785, 0.030819169692472, 0.00768198621715397), 

    "almost_linear_with_nonlinear_friction" = c(0.011672406201134, 0.0160287341397014, 
                                                0.00195195196755947, 0.001975799999565, 
                                                0.0297337459556888, 0.00737455656792283),

    ## Older results with froude_limit = 3.0                                
    #"leapfrog_nonlinear" = c(0.01018940190279, 0.01538559938733, 0.00186367844663, 
    #                         0.00189781202883, 0.02824797105766,  0.00700442470616),

    "leapfrog_nonlinear" = c(0.01019171450267, 0.01538531135840, 0.00186367837222, 
                             0.00189782061909, 0.02824770303451, 0.00700402286858),

    "linear_no_coriolis" = c(0.0114685360428089, 0.0171835712601261, 0.00194199504293318, 
                             0.00201682360561679, 0.0323221991569785, 0.00772623170987776),

    "linear_with_nonlinear_friction_no_coriolis" = c(0.0106521254715517, 0.0162701293387362, 
                                                     0.00186321670672357, 0.00195998848015663, 
                                                     0.0297601680229403, 0.0071694923930676), 

    "almost_linear_no_coriolis" = c(0.011748617982088, 0.0168087915480528, 0.00188845720196856, 
                                    0.00207088111467101, 0.0312334402938191, 0.00751337437756913), 

    "almost_linear_with_nonlinear_friction_nocoriolis" = c(0.0113546147788896, 0.0162295791353805, 
                                                           0.00186283720042138, 0.00200637061290488, 
                                                           0.0297637680102435, 0.00724945950376884), 

    #"leapfrog_nonlinear_nocoriolis" = c(0.00986683022037, 0.01547212705389, 0.00176974300331, 
    #                                    0.00192589322547, 0.02837528543449, 0.00685083032816)
    "leapfrog_nonlinear_nocoriolis" = c(0.00986497511932, 0.01547296598008, 0.00176974292866, 
                                        0.00192589106903, 0.02836895758817, 0.00685149956726)

    )

# Basic regression test. Note that compiler options (and probably compiler choice) can lead to slight
# differences (e.g. compiling with/without -march=native on gfortran)
model_var_stat_change = abs(model_var_stat/expected_model_var_stats[[model_type]] - 1)
if(all( model_var_stat_change < 1.0e-03)){
    print('PASS')
}else{
    print(c('FAIL (change in results) ', model_type, model_var_stat_change))
    #options(digits=12)
    #print(model_var_stat)
}

# Check on the error vs DART. Note this is a PTHA scenario, so will not be as close as a good inversion.
# From examining the results, these thresholds seem reasonable.
# They are ad-hoc choices, but have a chance of catching significant introduced errors.
err_tol = c( 0.15 , 0.1, 0.1, 0.2, 0.32, 0.153)
if(all(abs(unlist(gof_stats)) < err_tol)){
    print('PASS')
}else{
    print(c('FAIL', unlist(gof_stats)))
}
