source('../../../plot.R')
# Since we call SWALS from inside R, let's pass the openmp commands.
omp_run_command = Sys.getenv('OMP_RUN_COMMAND')

#
# First standard test case
#
# h/d ~ 0.0185
data_caseA = list(
    d = 0.298,
    g = 9.8,
    h_on_d = 0.0185,
    plot_ylim = c(-0.03, 0.07),
    plot_xlim = c(-5, 20),
    DEPTH_EPS = 1.0e-03,
    snapshot_files = Sys.glob('../test_repository/BP04-JosephZ-Single_wave_on_simple_beach/profs/Case0_0185*'),
    png_filename = 'Model-vs-data_0.0185.png',
    target_max_runup_dimensionless = 0.076,
    # Percentage error in max runup
    allowed_runup_error = 0.10
)
#
# Second standard test case
#
data_caseB = list(
    d = 0.1562,
    g = 9.8,
    h_on_d = 0.3,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.55,
    snapshot_files = Sys.glob('../test_repository/BP04-JosephZ-Single_wave_on_simple_beach/profs/Case0_3*'),
    png_filename = 'Model-vs-data_0.3.png',
    # Percentage error in max runup
    allowed_runup_error = 0.10
    )

#
# Other cases I choose to run (for which we don't have snapshot data)
# The h/d, d, and target_max_runup_dimensionless are based on the Lab_runup.txt data
#

# h/d ~ 0.1
data_caseC = list(
    d = 0.0981,
    g = 9.8,
    h_on_d = 0.097,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.274,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.10
    )

# h/d ~ 0.04 (reported to be near experimental breaking)
data_caseD = list(
    d = 0.2855,
    g = 9.8,
    h_on_d = 0.04,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.156,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.10
    )

# h/d ~ 0.2
data_caseE = list(
    d = 0.3535,
    g = 9.8,
    h_on_d = 0.193,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.426,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.12
    )

# h/d ~ 0.2 again, but slightly larger
data_caseF = list(
    d = 0.0989,
    g = 9.8,
    h_on_d = 0.236,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.467,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.12
    )

# h/d ~ 0.6
data_caseG = list(
    d = 0.1454,
    g = 9.8,
    h_on_d = 0.607,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.805,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.12
    )

# h/d ~ 0.4
data_caseH = list(
    d = 0.2349,
    g = 9.8,
    h_on_d = 0.394,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.641,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.12
    )

# h/d ~ 0.16
data_caseI = list(
    d = 0.2144,
    g = 9.8,
    h_on_d = 0.16,
    plot_ylim = c(-0.1, 0.5),
    plot_xlim = c(-10, 20),
    DEPTH_EPS = 1.0e-03,
    target_max_runup_dimensionless = 0.384,
    snapshot_files = NA,
    png_filename = NA,
    # Percentage error in max runup
    allowed_runup_error = 0.12
    )

# Run all cases
MODEL_CASES = list(data_caseA, data_caseB, data_caseC, data_caseD, data_caseE, data_caseF, data_caseG, data_caseH, data_caseI)

ncase = length(MODEL_CASES)
model_runup_test_store = matrix(NA, ncol=5, nrow=ncase)
colnames(model_runup_test_store) = c('d', 'h_on_d', 'R_on_d', 'target_r_on_d', 'error_fraction')

for(NJOB in 1:ncase){

    MODEL_CASE = MODEL_CASES[[NJOB]]

    # Get the config data
    attach(MODEL_CASE, warn.conflicts=FALSE)

    # Run the model
    #system(paste0(omp_run_command, ' ./BP4_testcases "rk2" ', d, ' ', h_on_d))
    #system(paste0('./BP4_testcases "leapfrog_nonlinear" ', d, ' ', h_on_d))
    system(paste0(omp_run_command, ' ./BP4_testcases "midpoint" ', d, ' ', h_on_d))

    # Get time-series at 2 locations
    # Avoid noisy messages about IO
    sink('tmp_sink_file')
    md = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])
    x = md[[1]] #get_all_recent_results()
    sink()

    if(!is.na(png_filename)){

        # Times should be scaled by sqrt(g/d)
        snapshot_file_times_dimensionless = unlist(lapply(strsplit(snapshot_files, 't='), f<-function(x) as.numeric(x[2])))
        snapshot_file_times = snapshot_file_times_dimensionless * sqrt(d/g)

        # Relevant model indices
        nearest_time_index = sapply(snapshot_file_times, f<-function(y) which.min(abs(y - x$gauges$time)))

        # Plot time snapshots
        png(png_filename, width=10, height=5, units='in', res=200)
        par(mfrow=c(2,3))
        for(i in 1:length(snapshot_files)){
            model_ind = nearest_time_index[i]
            stg = x$gauges$time_var$stage[,model_ind]
            elv = x$gauges$static_var$elevation0
            k = which(stg > elv + DEPTH_EPS)
            plot(x$gauges$lon[k]/d, stg[k]/d, t='o', ylim=plot_ylim, cex=0.1, pch=19,
                 xlim = plot_xlim, xlab='x/d', ylab='Stage / d', cex.axis=1.3, cex.lab=1.3)
            points(x$gauges$lon/d, elv/d, t='l', col='grey')
            dat = read.table(snapshot_files[i])
            points(dat, col='red', t='o', pch=19, cex=0.1)
            title(main = paste0('t/T = ', snapshot_file_times_dimensionless[i]), cex.main=1.5)
            grid(col='brown')
            points(x$gauges$lon/d, stg/d, t='l')
        }
        dev.off()
    }

    #
    # Check max runup
    #
    runup_model = apply(x$gauges$time_var$stage, 2, f<-function(y){
        elev = x$gauges$static_var$elevation0
        ind = min(which(y > elev + DEPTH_EPS))
        return(y[ind])
             })
    max_runup_model_dimensionless = max(runup_model)/d
    
    # Allow for some error.
    #print(c(paste('TEST h/d=', h_on_d), 'Observed Runup', 'Target runup'))
    error_fraction = abs(max_runup_model_dimensionless - target_max_runup_dimensionless)/target_max_runup_dimensionless 
    if(error_fraction < allowed_runup_error){
        #print(c('PASS', max_runup_model_dimensionless, target_max_runup_dimensionless))
        print('PASS')
    }else{
        #print(c('FAIL', max_runup_model_dimensionless, target_max_runup_dimensionless))
        print('FAIL')
    }

    # Store the stats
    model_runup_test_store[NJOB,] = c(d, h_on_d, max_runup_model_dimensionless, 
                                      target_max_runup_dimensionless, error_fraction)
}

#
# Scaling plot
#
lab_data = read.table('../test_repository/BP04-JosephZ-Single_wave_on_simple_beach/Lab_runup.txt')
png('Runup_scaling_plot.png', width=8, height=6, units='in', res=300)
plot(lab_data[,1:2], xlab='h/d', ylab='R/d', main='NTHMP Benchmark problem 4: Runup scaling', 
     cex.main=1.7, cex.lab=1.5, cex.axis=1.5)
grid()
points(model_runup_test_store[,2], model_runup_test_store[,3], col='red', pch=19)
legend('topleft', c('Data', 'Model'), col=c('black', 'red'), pch=c(1, 19), cex=1.2, bty='n')
dev.off()
write.csv(model_runup_test_store, 'runup_model_test_results.csv', row.names=FALSE)

