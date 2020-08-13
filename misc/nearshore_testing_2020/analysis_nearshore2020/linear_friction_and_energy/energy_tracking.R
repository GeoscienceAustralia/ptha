#
# Have a look at the global-scale energy decay in various models.
#
ONE_HOUR = 3600

## These log-files are the same as those used above (for which we have data)
## There is no reason to restrict the result in this way here
#manning_log_files = paste0(dirname(files_manning0.035), "/multidomain_log_image_00000000000000000001.log")

## This is 'ALL' the manning_0.035 runs with high-res areas
manning_log_files = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_manning*-highres_NSW/*/*.log'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_manning*-highres_australia/*/*.log'))
# Drop HIGHRES scenarios
manning_log_files = manning_log_files[-grep('HIGHRES', manning_log_files) ]
# Drop 'Chile 1960' scenarios that had highres_australia has high-resolution, because they are repeats of the ones with highres_NSW
# (done out of interest for WA colleagues).
manning_log_files = manning_log_files[-which(grepl('Chile1960', manning_log_files) & grepl('highres_australia', manning_log_files)) ]

# This is 'ALL' the manning_0.035 runs, even those without high-res areas
# This shows there is little change to the energy budget induced by using high-res areas; as expected
# because they cover a small area of the domain
#manning_log_files = c(
#    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_manning*/*/*.log'))


run_info = basename(dirname(dirname(manning_log_files)))

get_time_and_global_energy_from_log<-function(log_file){
    x = readLines(log_file)
    # Find log-lines just above the reported total energy
    k = grep('Global energy-total', x)
    energy_on_rho = as.numeric(x[k+1])
    # In these log files, the desired time (secs) happens to be in index "k-33"
    # (detail of the sequence of print-outs, likely to change in future)
    if(all(substring(x[k-34], 1, 5) == 'Time:')){
        time = as.numeric(x[k-33])
    }else{
        stop('Did not find "Time:" in index (k-34) -- has the log file format changed?')
    }

    return(data.frame(time=time, energy_on_rho=energy_on_rho))
}

manning_energy = lapply(manning_log_files, get_time_and_global_energy_from_log)

# Look at the "end - start" energy relative to the initial energy
manning_energy_change = unlist(lapply(manning_energy, 
    f<-function(x) (x[length(x[,2]),2] - x[1,2])/x[1,2]))
#> manning_energy_change
# [1] -0.8420079 -0.7907710 -0.7901383 -0.8019070 -0.6263134 -0.6567077
# [7] -0.8452437 -0.7775568 -0.8666060 -0.9271099 -0.9587252 -0.9543358
#> range(manning_energy_change)
#[1] -0.9587252 -0.6263134


energy_decay_constant<-function(time, energy, fit_start_time = -1){
    # Suppose energy ~ energy_time_zero * exp(CONSTANT * time)
    # Here we estimate CONSTANT. We can optionally fit only to later times
    # by setting e.g. fit_start_time = 3600*12  (start fit at 12 hours)
    k = which(time > fit_start_time)
    energy = energy[k]
    time = time[k]
    time = time - time[1]

    energy_norm = energy/energy[1]
    fit = lm(log(energy_norm) ~ time + 0)
    return(fit)
}

energy0 = unlist(lapply(manning_energy, f<-function(x) x$energy_on_rho[1]))

# Decay-time from log-linear fit to energy-vs-time
decay_time_hours = unlist(lapply(manning_energy, f<-function(x){
    fit = energy_decay_constant(x$time, x$energy_on_rho)
    return(1/coef(fit) * 1/ONE_HOUR)
             }), use.names=FALSE)
decay_constant = 1/(decay_time_hours * ONE_HOUR)

# Decay-time from log-linear fit to energy-vs-time, AFTER 36 hours of simulation
decay_constant_skip36hr = unlist(lapply(manning_energy, 
    f<-function(x) coef(energy_decay_constant(x$time, x$energy_on_rho, fit_start_time=36*ONE_HOUR))))
decay_time_hours_skip36hr = 1/decay_constant_skip36hr * 1/ONE_HOUR

# Another estimation technique -- line passing through energy at t0 and tlast
decay_constant_firstlast = unlist(lapply(manning_energy, f<-function(x){
    e0 = x$energy_on_rho[1]
    l = length(x$energy_on_rho)
    decay_coef = ( log(x$energy_on_rho[l]) - log(e0) )/(x$time[l] - x$time[1])
    }))
decay_time_firstlast = 1/decay_constant_firstlast * 1/ONE_HOUR

# 
png('Energy_decay_time_fit_to_global_manning_energy_dissipation.png', width=15, height=6, units='in', res=300)
par(mfrow=c(1,3))
plot(energy0, -decay_time_hours, 
     main='Energy e-folding time, based on a log-linear fit to \n the modelled dissipation with manning=0.035', 
     log='xy', ylim=c(10, 100))
text(energy0, -decay_time_hours, substring(run_info, 1, 12))
grid()

plot(energy0, -1/decay_constant_skip36hr * 1/ONE_HOUR, 
     main='Energy e-folding time: as per left but ignoring first 36h', 
     log='xy', ylim=c(10, 100))
text(energy0, -1/decay_constant_skip36hr * 1/ONE_HOUR, substring(run_info, 1, 12))
grid()

plot(energy0, -1/decay_constant_firstlast * 1/ONE_HOUR, 
     main='Energy e-folding time: match first/last energies in 60h', 
     log='xy', ylim=c(10, 100))
text(energy0, -1/decay_constant_firstlast * 1/ONE_HOUR, substring(run_info, 1, 12))
grid()
dev.off()

pdf('Energy_decay_manning_models.pdf', width=8, height=6)
for(i in 1:length(manning_energy)){
    plot(manning_energy[[i]]$time, manning_energy[[i]]$energy, t='l', log='y')
    points(manning_energy[[i]]$time, energy0[i] * exp(decay_constant[i]*manning_energy[[i]]$time), t='l', col='red')
    grid(col='orange')
    title(substring(run_info[i], 1, 30))
    title(paste0('Time decay from ', c('start :', '36 hours:', '(First/last)'), 
                 round(-1/c(decay_constant[i], decay_constant_skip36hr[i], decay_constant_firstlast[i]) * 1/ONE_HOUR)), 
          line=0)
    points(manning_energy[[i]]$time, energy0[i] * exp(manning_energy[[i]]$time*(-1e-05)), t='l', col='brown', lty='dotted')
    points(manning_energy[[i]]$time, energy0[i] * exp(manning_energy[[i]]$time*(decay_constant_skip36hr[i])), t='l', col='green')
    points(manning_energy[[i]]$time, energy0[i] * exp(manning_energy[[i]]$time*(decay_constant_firstlast[i])), t='l', col='purple')
}
dev.off()

#
# Summarise energy decay for 'linear-friction-offshore' models. 
#
linear_friction_log_files = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_linear_friction*-highres_NSW/*/*.log'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_linear_friction*-highres_australia/*/*.log'))
linear_friction_energy = lapply(linear_friction_log_files, get_time_and_global_energy_from_log)
# Look at the "end - start" energy relative to the initial energy
linear_friction_energy_change = unlist(lapply(linear_friction_energy, 
    f<-function(x) (x[length(x[,2]),2] - x[1,2])/x[1,2]))
#> linear_friction_energy_change
# [1] -0.8852562 -0.8847494 -0.8855192 -0.8858640 -0.8856396 -0.8856412
# [7] -0.8848984 -0.8846983 -0.8844612 -0.8875214 -0.8871057 -0.8857192

# For this model we know the expected energy decay over 2.5 days -- it should
# closely match the values above
theoretical_energy_change = exp(2.5 * 24 * 3600 * (-1e-05))  - 1
#> theoretical_energy_change
#[1] -0.8846749
summary(theoretical_energy_change - linear_friction_energy_change )
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0002137  0.0001863  0.0009045  0.0009146  0.0010805  0.0028466 



##
## SUMMARISE ENERGY CONSERVATION FOR THE "FRICTIONELSS-OFFSHORE" MODELS 
##
frictionless_log_files = c(
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_no_friction*-highres_NSW/*/*.log'),
    Sys.glob('../gauge_RDS_files/OUTPUTS/*full-linear_with_no_friction*-highres_australia/*/*.log'))

frictionless_energy = lapply(frictionless_log_files, get_time_and_global_energy_from_log)
# Look at the "end - start" energy relative to the initial energy
frictionless_energy_change = unlist(lapply(frictionless_energy, 
                                           f<-function(x) (x[length(x[,2]),2] - x[1,2])/x[1,2]))
summary(frictionless_energy_change)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.027245 -0.009810 -0.004792 -0.007952 -0.002363  0.002539 

#
# Check "start - max-energy" and "start - min-energy" to catch any 'drift' not identified above
#
frictionless_energy_change_positive = unlist(lapply(frictionless_energy, 
                                           f<-function(x) (max(x[,2]) - x[1,2])/x[1,2]))
summary(frictionless_energy_change_positive)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0005300 0.0007745 0.0008320 0.0010504 0.0010452 0.0025391 

frictionless_energy_change_negative = unlist(lapply(frictionless_energy, 
                                           f<-function(x) (min(x[,2]) - x[1,2])/x[1,2]))
summary(frictionless_energy_change_negative)
#     Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0272448 -0.0098104 -0.0048077 -0.0084020 -0.0027305 -0.0001631 



