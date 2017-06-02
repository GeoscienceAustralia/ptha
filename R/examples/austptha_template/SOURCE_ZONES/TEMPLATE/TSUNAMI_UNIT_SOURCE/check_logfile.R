#' Function to read the model logfile, and summarise output
#'
#' It can plot the max stage envelope and mass balance over time,
#' and flag warnings if the peak stage is too extreme
#'
check_logfile<-function(logfile_name, plot=FALSE, extreme_stage_warning=Inf){
   
    log_data = readLines(logfile_name) 

    log_time_ind = grep('Time:', log_data)

    times = as.numeric(log_data[log_time_ind+1])
    
    log_stage_ind = grep('Stage:', log_data)
    stage_max = as.numeric(log_data[log_stage_ind+1])
    stage_min = as.numeric(log_data[log_stage_ind+2])

    mass_bal_ind = grep('Mass balance:', log_data)
    mass_bal = as.numeric(gsub('Mass balance:', '', log_data[mass_bal_ind]))


    stage_global_max = max(stage_max)
    stage_global_min = min(stage_min)

    if(plot){
        ylims = c(stage_global_min, stage_global_max)

        kk = min(length(times), length(stage_max), length(stage_min))
       
        par(mfrow=c(2,1)) 
        plot(times[1:kk], stage_max[1:kk], t='l', 
            main=paste0('Stage envelope \n ', logfile_name), 
            ylim=ylims, xlab='Time (s)', ylab = 'Stage (m)')
        points(times[1:kk], stage_min[1:kk], t='l')
        grid()

        kk = min(length(times), length(mass_bal))
        plot(times[1:kk], mass_bal[1:kk], t='l', xlab='Time (s)', ylab='Volume (m^3)')

        mass_balance_mean = mean(mass_bal)
        mass_balance_range = diff(range(mass_bal))

        title(main=paste0('Volume balance: (relative range = ', 
            signif(mass_balance_range/mass_balance_mean, 3), ')'))
        grid()
    }

    if(abs(stage_global_max) > extreme_stage_warning){
        print(paste0('Extreme max stage warning: ', logfile_name))
    }

    if(abs(stage_global_min) > extreme_stage_warning){
        print(paste0('Extreme min stage warning: ', logfile_name))
    }

}

### Main
source_name = basename(dirname(getwd()))
all_logs = Sys.glob(
    paste0('/g/data/w85/tsunami/AustPTHA/version1/unit_sources/', 
        source_name, '/unit_source_tsunami/RUN_*/*/log*'))

# Try ordering them along-strike/down-dip
us_name = basename(dirname(dirname(all_logs)))
dip_ind = sapply(us_name, f<-function(x) strsplit(x, '_')[[1]][[4]])
strike_ind = sapply(us_name, f<-function(x) strsplit(x, '_')[[1]][[5]])
dip_strike = cbind(as.numeric(dip_ind), as.numeric(strike_ind))
o1 = order(dip_strike[,2], dip_strike[,1])
all_logs = all_logs[o1] # Now in correct order

pdf('tsunami_log_check.pdf', width=10, height=8)
for(logi in all_logs){
    try(check_logfile(logi, plot=TRUE, extreme_stage_warning=5))
}
dev.off()


